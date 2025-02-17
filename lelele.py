BACKEND_FLATN: str = 'flatn'
BACKEND_FPYLLL: str = 'fpylll'
BACKEND_DEFAULT: str = BACKEND_FPYLLL

from hashlib import sha256
from typing import Iterator, Generator

def _wrap_lin(ctx: 'LeLeLe', val: 'LinearCombination | Variable | int') -> 'LinearCombination':
    assert isinstance(ctx, LeLeLe)

    # linear combination
    if isinstance(val, LinearCombination):
        return val

    # variable
    if isinstance(val, Variable):
        return val.lin()

    # constant
    try:
        return LinearCombination(ctx=ctx, combine={ctx.one(): int(val)})
    except ValueError:
        raise ValueError('failed to covert %r to linear combination' % val)

def _zero_matrix(n: int, m: int) -> list[list[int]]:
    M = []
    for _ in range(n):
        M.append([0] * m)
    return M

def _identity(n: int) -> list[list[int]]:
    '''
    Returns an identity matrix of size n x n
    '''
    M = _zero_matrix(n, n)
    for i in range(n):
        M[i][i] = 1
    return M

def _concat(A: list[list[int]], B: list[list[int]]) -> list[list[int]]:
    assert len(A) == len(B), 'Matrices must have the same number of rows'
    return [A[i] + B[i] for i in range(len(A))]

class LeLeLe:
    def __init__(self):
        self.vars: list['Variable'] = []
        self.vone: 'Variable | None' = None
        self.constraints: list[tuple['LinearCombination', int]] = []

    def one(self) -> 'Variable':
        '''
        Returns a "variable" which should the constant value one.
        '''
        if self.vone is None:
            self.vone = self.bit(name='1')
        return self.vone

    def bit(self, name: str | None = None) -> 'Variable':
        return self.var(name, 'bit').short(norm=1)

    def byte(self, name: str | None = None) -> 'Variable':
        return self.var(name, 'byte').short(norm=0xff)

    def word(self, width: int, name: str | None = None) -> 'Variable':
        assert width > 0
        return self.var(name, 'word').short(norm=(1 << width)-1)

    def var(self, name: str | None = None, prefix: str = 'var') -> 'Variable':
        idx = len(self.vars)
        var = Variable(
            self,
            name if name else '%s%d' % (prefix, idx),
            idx
        )
        self.vars.append(var)
        return var

    def add_constraint(self, lin: 'LinearCombination', norm: int) -> None:
        assert isinstance(lin, LinearCombination)
        self.constraints.append((lin, int(norm)))

    def system(self) -> list[list[int]]:
        '''
        Returns the Matrix representing the LLL system.

        This does not require fpylll.
        '''

        # calculathe largest norm
        max_norm = max([norm for _, norm in self.constraints], default=1)

        rows = len(self.vars)
        cols = len(self.constraints)

        M = _zero_matrix(rows, cols)

        for (i, (lin, norm)) in enumerate(self.constraints):
            rescale = max_norm // norm # rescale factor
            for (var, scl) in lin.combine.items():
                M[var.index][i] = scl * rescale

        return M

    def solve_raw(self, backend: str | None = None, solver_args: dict | None = None) -> list[list[int]]:
        '''
        Solves the system and returns the raw reduced matrix.

        This method requires either fpylll or flatn to be installed.

        Args:
            backend (str, optional):
                The backend to use for matrix reduction. Can be 'fpylll', 'flatn' or None.
                If None, will use flatn if available, otherwise fpylll.
            solver_args (dict, optional): Additional arguments to pass to the chosen solver backend.

        Returns:
            list[list[int]]: The reduced matrix

        Raises:
            ImportError: If neither fpylll nor flatn is installed
            ValueError: If an invalid backend is specified
        '''

        # generate the matrix
        M = self.system()

        # Default empty args if none provided
        if solver_args is None:
            solver_args = {}

        # reduce the matrix with the chosen backend
        if backend == 'fpylll':
            from fpylll import IntegerMatrix, LLL
            R = IntegerMatrix.from_matrix(M)
            LLL.reduction(R, **solver_args)
            return [list(row) for row in R]
        elif backend == BACKEND_FLATN:
            import flatn
            return flatn.reduce(M, **solver_args)
        elif backend is None:
            return self.solve_raw(backend=BACKEND_DEFAULT, solver_args=solver_args)
        else:
            raise ValueError('Invalid backend')

    def solve(self, backend: str | None = None, solver_args: dict | None = None) -> 'Solution':
        '''
        Solves the system and assigns solution values to all variables.

        This method solves the linear system using matrix reduction and assigns
        the resulting values to each variable in the system.

        Args:
            backend (str, optional): The backend solver to use ('fpylll' or 'flatn').
                Defaults to 'flatn' if available, otherwise 'fpylll'.
            solver_args (dict, optional): Additional arguments to pass to the solver backend.

        Returns:
            Solution: A Solution object containing the variable assignments and reduced matrix.

        Raises:
            ImportError: If neither fpylll nor flatn solver backend is installed.
            ValueError: If an invalid backend is specified.
        '''
        return Solution(
            self.solve_raw(
                backend=backend,
                solver_args=solver_args
            ),
            self.vone,
            self.vars,
            self.constraints
        )

    def __repr__(self) -> str:
        cons = []
        for (lin, norm) in self.constraints:
            cons.append('    0x%x >= |%s|' % (norm,lin))
        return 'LeLeLe(\n' + '\n'.join(cons) + '\n)'

class Solution:
    def __init__(
        self,
        R: list[list[int]],
        vone: 'Variable | None',
        vars: list['Variable'],
        cons: list[tuple['LinearCombination', int]],
    ):
        self.R = [list(row) for row in R]
        self.vone = vone
        self.vars = list(vars)
        self.cons = list(cons)
        self.rels = [rel for (rel, _norm) in cons]

        assert len(self.R) == len(self.vars)
        assert len(self.R[0]) == len(self.rels)

        # extract solution from matrix
        self.all_sols = self.solve()
        self.assign_vars: dict['Variable', int] = {}
        self.assign_rels: dict['LinearCombination', int] = {}
        self.next_solution()

    def reduced(self) -> list[list[int]]:
        return self.R

    def __str__(self) -> str:
        longest_var = max(len(str(var)) for var in self.assign_vars)
        lines = []
        lines.append("Solution(")
        for var, value in self.assign_vars.items():
            if var is None:
                continue
            padding = " " * (longest_var - len(str(var)))
            lines.append(f"    {padding}{var} = {hex(value)}")
        lines.append(")")
        return "\n".join(lines)

    def __call__(self, obj: 'Variable | LinearCombination') -> int:
        if isinstance(obj, Variable):
            return self.assign_vars[obj]

        if isinstance(obj, LinearCombination):
            # try to lookup directly
            if obj in self.assign_rels:
                return self.assign_rels[obj]

            # compute from variables
            result = 0
            for (var, scl) in obj:
                try:
                    result += scl * self.assign_vars[var]
                except KeyError:
                    raise ValueError(f"Variable {var} not assigned")
            return result

        raise TypeError('Expected Variable or LinearCombination')

    def next_solution(self) -> None:
        assign_vars, assign_rels = next(self.all_sols)
        self.assign_vars = assign_vars
        self.assign_rels = assign_rels

    def solve(self) -> Generator[tuple[dict['Variable', int], dict['LinearCombination', int]], None, None]:
        for sol in self.R:
            # map: var -> int (assignment)
            assign_rels: dict['LinearCombination', int] = {}
            assign_vars: dict['Variable', int] = {}

            # assign linear combinations
            for (rel, val) in zip(self.rels, sol):
                assign_rels[rel] = val

            # use a basic peeling algorithm to recover variables:
            # map: linear relation -> unknown terms
            unknown = {rel: len(rel) for rel in self.rels}

            # map: linear relation -> value of linear combine of unknown terms
            values = {rel: val for rel, val in zip(self.rels, sol)}

            # map: var -> LinearCombination's with var in it
            appear = {var: [] for var in self.vars}
            for rel in self.rels:
                for var in rel.vars():
                    appear[var].append(rel)

            while True:
                # find singleton
                for (rel, n) in unknown.items():
                    if n == 1: break
                else:
                    break

                # find the single unknown variable
                val = values[rel]
                for (var, scl) in rel:
                    if var not in assign_vars:
                        break
                else:
                    assert False, "unreachable"

                # var * scl = val
                print(rel)
                print(val)
                print(scl)
                assert val % scl == 0
                assigned = val // scl

                # remove from other combinations
                if var in assign_vars:
                    assert assign_vars[var] == assigned
                    continue

                for orel in appear[var]:
                    values[orel] -= assigned * orel.coeff(var)
                    assert unknown[orel] > 0
                    unknown[orel] -= 1

                # removed from every relation
                appear[var] = []

                # add to assignments
                assign_vars[var] = assigned

            yield (assign_vars, assign_rels)

class LinearCombination:
    '''
    An immutable linear combination of variables.
    '''
    def __init__(self, ctx: 'LeLeLe', combine: dict['Variable', int]) -> None:
        assert isinstance(combine, dict), f"combine must be a dictionary, got {type(combine)}"
        self.ctx = ctx
        self.combine = combine
        self.solutions: list[int] | None = None

        # allows quick comparison:
        # without O(n) complexity for every equality check
        self.identifier = sha256(str(self).encode()).hexdigest()

    def vars(self) -> set['Variable']:
        return set(self.combine.keys())

    def coeff(self, var: 'Variable') -> int:
        return self.combine.get(var, 0)

    def __len__(self) -> int:
        return len(self.combine)

    def __iter__(self) -> Iterator[tuple['Variable', int]]:
        return iter(sorted(self.combine.items()))

    def __eq__(self, other: object) -> bool:
        # check the type
        if not isinstance(other, LinearCombination):
            return False

        # check the context
        if self.ctx != other.ctx:
            return False

        # check the linear combination
        return self.identifier == other.identifier

    def __hash__(self) -> int:
        return hash((id(self.ctx), self.identifier))

    def __repr__(self) -> str:
        lin = ['%r * %s' % (v, hex(s)) if s != 1 else '%r' % v for (v, s) in self]
        return ' + '.join(lin)

    def __mod__(self, other) -> 'LinearCombination':
        try:
            n = int(other)
        except ValueError:
            raise ValueError('Modulo of linear combination only defined for integers')
        combine = dict(self.combine)
        combine[self.ctx.var()] = n
        return LinearCombination(ctx=self.ctx, combine=combine)

    def __neg__(self) -> 'LinearCombination':
        return LinearCombination(
            ctx=self.ctx,
            combine={s: -v for (s, v) in self.combine.items()}
        )

    def __sub__(self, other) -> 'LinearCombination':
        return self + (-other)

    def __rsub__(self, other) -> 'LinearCombination':
        return self.__sub__(other)

    def __add__(self, other) -> 'LinearCombination':
        if other == 0: return self # this is convenient

        # convert to linear combination
        other = _wrap_lin(self.ctx, other)
        assert self.ctx == other.ctx, 'linear combinations belong to different systems'

        # combine terms
        combine = dict(self.combine)
        for (var, scl) in other.combine.items():
            try:
                combine[var] += scl
            except KeyError:
                combine[var] = scl

        return LinearCombination(ctx=self.ctx, combine=combine)

    def __radd__(self, other) -> 'LinearCombination':
        return self.__add__(other)

    def __mul__(self, other) -> 'LinearCombination':
        try:
            n = int(other)
        except ValueError:
            raise ValueError('Can only mul. linear combination by integer, not %r' % other)
        return LinearCombination(
            ctx=self.ctx,
            combine={v: s * n for (v, s) in self.combine.items()}
        )

    def __rmul__(self, other) -> 'LinearCombination':
        return self.__mul__(other)

    def short(self, norm: int = 1) -> 'LinearCombination':
        '''
        Constrain the linear combination to have small norm.

        The "norm" parameter should be understood
        as the "max-value" of the variable/linear combination.
        '''
        if norm <= 0:
            raise ValueError('Norm must be positive')
        self.ctx.add_constraint(self, norm)
        return self

    def __getitem__(self, n: int) -> int:
        '''
        Return the n'th solution
        '''

        # check if constraint and solution set
        if self.solutions is not None:
            return self.solutions[n]

        # otherwise compute from variable assignments
        return sum([s * v[n] for (v, s) in self.combine.items()])

    def __int__(self) -> int:
        return int(self[0])

    def __index__(self) -> int:
        return int(self)

class Variable:
    def __init__(self, ctx: 'LeLeLe', name: str, index: int):
        self.ctx = ctx
        self.name = name
        self.index = index
        self.solutions: list[int] | None = None # not solved

    def lin(self) -> LinearCombination:
        return LinearCombination(ctx=self.ctx, combine={self: 1})

    def key(self) -> tuple[str, int, int]:
        '''
        Returns the ordering key for consistent comparisons.
        '''
        return (self.name, self.index, id(self.ctx))

    def __lt__(self, other: 'Variable') -> bool:
        if not isinstance(other, Variable):
            return NotImplemented
        return self.key() < other.key()

    def __le__(self, other: 'Variable') -> bool:
        if not isinstance(other, Variable):
            return NotImplemented
        return self.key() <= other.key()

    def __gt__(self, other: 'Variable') -> bool:
        if not isinstance(other, Variable):
            return NotImplemented
        return self.key() > other.key()

    def __ge__(self, other: 'Variable') -> bool:
        if not isinstance(other, Variable):
            return NotImplemented
        return self.key() >= other.key()

    def __neg__(self) -> LinearCombination:
        return - self.lin()

    def __sub__(self, other) -> LinearCombination:
        return self.lin() - _wrap_lin(self.ctx, other)

    def __rsub__(self, other) -> LinearCombination:
        return self.lin() - _wrap_lin(self.ctx, other)

    def __add__(self, other) -> LinearCombination:
        return self.lin() + _wrap_lin(self.ctx, other)

    def __radd__(self, other) -> LinearCombination:
        return self.__add__(other)

    def __mul__(self, other) -> LinearCombination:
        return self.lin() * int(other)

    def __rmul__(self, other) -> LinearCombination:
        return self.__mul__(other)

    def __repr__(self) -> str:
        return self.name

    def __index__(self) -> int:
        return int(self)

    def __mod__(self, other) -> LinearCombination:
        return self.lin() % int(other)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Variable):
            return False
        return self.key() == other.key()

    def __hash__(self) -> int:
        return hash(self.key())

    def short(self, norm: int = 1) -> 'Variable':
        '''
        Constrain the variable to have small norm.
        '''
        self.lin().short(norm)
        return self
