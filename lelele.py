from __future__ import annotations

BACKEND_FLATN: str = 'flatn'
BACKEND_FPYLLL: str = 'fpylll'
BACKEND_DEFAULT: str = BACKEND_FPYLLL

from math import gcd
from typing import Iterator, Generator
from functools import reduce, total_ordering

def _wrap_lin(ctx: LeLeLe, val: LinearCombination | Variable | int) -> LinearCombination:
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
        raise ValueError(f'failed to convert {val!r} to linear combination')

def _zero_matrix(n: int, m: int) -> list[list[int]]:
    return [[0] * m for _ in range(n)]

def _reduce_row_gcd(row: list[int]) -> None:
    g = reduce(gcd, row, 0)
    if g > 1:
        for i in range(len(row)):
            row[i] //= g

class LeLeLe:
    def __init__(self):
        self.vars: list[Variable] = []
        self.vone: Variable | None = None
        self.constraints: list[tuple[LinearCombination, int]] = []

    def one(self) -> Variable:
        '''
        Returns a "variable" which should be the constant value one.
        '''
        if self.vone is None:
            self.vone = self.bit(name='1')
        return self.vone

    def bit(self, name: str | None = None) -> Variable:
        return self.var(name, 'bit').short(norm=1)

    def byte(self, name: str | None = None) -> Variable:
        return self.var(name, 'byte').short(norm=0xff)

    def word(self, width: int, name: str | None = None) -> Variable:
        assert width > 0
        return self.var(name, 'word').short(norm=(1 << width)-1)

    def var(self, name: str | None = None, prefix: str = 'var') -> Variable:
        idx = len(self.vars)
        var = Variable(
            self,
            name if name else f'{prefix}{idx}',
            idx
        )
        self.vars.append(var)
        return var

    def add_constraint(self, lin: LinearCombination, norm: int) -> None:
        assert isinstance(lin, LinearCombination)
        self.constraints.append((lin, int(norm)))

    def system(self) -> list[list[int]]:
        '''
        Returns the Matrix representing the LLL system.

        This does not require fpylll.
        '''

        # calculate the largest norm
        max_norm = max([norm for _, norm in self.constraints], default=1)

        rows = len(self.vars)
        cols = len(self.constraints)

        M = _zero_matrix(rows, cols)

        for i, (lin, norm) in enumerate(self.constraints):
            rescale = max_norm // norm
            for var, scl in lin.combine.items():
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
        if backend == BACKEND_FPYLLL:
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

    def solve(self, backend: str | None = None, solver_args: dict | None = None) -> Solutions:
        '''
        Solves the system and returns an iterable of solutions.

        This method solves the linear system using matrix reduction and returns
        a Solutions object that iterates over non-degenerate variable assignments.

        Args:
            backend (str, optional): The backend solver to use ('fpylll' or 'flatn').
                Defaults to 'flatn' if available, otherwise 'fpylll'.
            solver_args (dict, optional): Additional arguments to pass to the solver backend.

        Returns:
            Solutions: An iterable of Solution snapshots. Supports direct calls
                via sol(var) which delegates to the first non-degenerate solution.

        Raises:
            ImportError: If neither fpylll nor flatn solver backend is installed.
            ValueError: If an invalid backend is specified.
        '''
        return Solutions(
            self.solve_raw(
                backend=backend,
                solver_args=solver_args
            ),
            self.vone,
            self.vars,
            self.constraints,
            self.system(),
        )

    def __repr__(self) -> str:
        cons = []
        for (lin, norm) in self.constraints:
            cons.append(f'    {norm:#x} >= |{lin}|')
        return 'LeLeLe(\n' + '\n'.join(cons) + '\n)'

class Solution:
    '''
    Immutable snapshot of a single solution's variable assignments.
    '''
    def __init__(
        self,
        assign_vars: dict[Variable, int],
        assign_rels: dict[LinearCombination, int],
    ):
        self.assign_vars = assign_vars
        self.assign_rels = assign_rels

    def __str__(self) -> str:
        longest_var = max((len(str(var)) for var in self.assign_vars), default=0)
        lines = []
        lines.append("Solution(")
        for var, value in self.assign_vars.items():
            if var is None:
                continue
            padding = " " * (longest_var - len(str(var)))
            lines.append(f"    {padding}{var} = {hex(value)}")
        lines.append(")")
        return "\n".join(lines)

    def __call__(self, obj: Variable | LinearCombination) -> int:
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


class Solutions:
    '''
    Iterable over solutions from lattice reduction. Skips degenerate rows.
    '''
    def __init__(
        self,
        R: list[list[int]],
        vone: Variable | None,
        vars: list[Variable],
        cons: list[tuple[LinearCombination, int]],
        M: list[list[int]],
    ):
        self.R = [list(row) for row in R]
        self.vone = vone
        self.vars = list(vars)
        self.cons = list(cons)
        self.rels = [rel for (rel, _norm) in cons]
        self.M = M
        self._first: Solution | None = None

        assert len(self.R) == len(self.vars)
        assert len(self.R[0]) == len(self.rels)

    def reduced(self) -> list[list[int]]:
        '''Returns the raw reduced matrix.'''
        return self.R

    @staticmethod
    def _solve_row(M: list[list[int]], sol: list[int], num_vars: int, num_cons: int):
        '''
        Solve M^T * x = sol for x using Gauss-Jordan elimination over Z.

        M is num_vars x num_cons, so M^T is num_cons x num_vars.
        We build the augmented matrix [M^T | sol] and reduce.

        Returns a list mapping variable-index -> value,
        or None if the row is inconsistent / doesn't determine all variables.
        '''

        # build augmented matrix [M^T | sol]
        nrows = num_cons
        ncols = num_vars + 1  # +1 for augmented column
        A = []
        for c in range(num_cons):
            row = [M[v][c] for v in range(num_vars)]
            row.append(sol[c])
            A.append(row)

        # Gauss-Jordan elimination over Z
        pivot_row = [None] * num_vars  # maps column -> row with pivot
        pr = 0  # next pivot row
        for col in range(num_vars):
            # find row with smallest nonzero absolute value in this column (from pr onward)
            best = None
            for r in range(pr, nrows):
                if A[r][col] != 0:
                    av = abs(A[r][col])
                    if best is None or av < best[0]:
                        best = (av, r)
            if best is None:
                continue  # no pivot in this column

            _, br = best
            # swap to pivot position
            A[pr], A[br] = A[br], A[pr]

            # make pivot positive
            if A[pr][col] < 0:
                for j in range(ncols):
                    A[pr][j] = -A[pr][j]

            # eliminate all other rows (both above and below)
            for r in range(nrows):
                if r == pr or A[r][col] == 0:
                    continue
                # eliminate: row[r] = row[r] * pivot - row[pr] * A[r][col]
                factor = A[r][col]
                pivot_val = A[pr][col]
                for j in range(ncols):
                    A[r][j] = A[r][j] * pivot_val - A[pr][j] * factor
                _reduce_row_gcd(A[r])

            pivot_row[col] = pr
            pr += 1

        # extract variable values
        result = [None] * num_vars
        for col in range(num_vars):
            r = pivot_row[col]
            if r is None:
                return None  # underdetermined
            pv = A[r][col]
            rhs = A[r][num_vars]
            if pv == 0:
                return None
            if rhs % pv != 0:
                return None  # not an integer solution
            result[col] = rhs // pv

        return result

    def __iter__(self) -> Generator[Solution, None, None]:
        num_vars = len(self.vars)
        num_cons = len(self.rels)

        for sol in self.R:
            # assign linear combinations
            assign_rels: dict[LinearCombination, int] = {}
            for (rel, val) in zip(self.rels, sol):
                assign_rels[rel] = val

            # solve M^T * x = sol via Gaussian elimination
            result = self._solve_row(self.M, sol, num_vars, num_cons)
            if result is None:
                # could not determine all variables from this row
                continue

            # skip degenerate: vone=0 means all constants vanish
            if self.vone is not None and result[self.vone.index] == 0:
                continue

            # normalize sign: LLL may return -v instead of v;
            # if vone is negative, flip the entire solution
            if self.vone is not None and result[self.vone.index] < 0:
                result = [-v for v in result]
                assign_rels = {rel: -val for rel, val in assign_rels.items()}

            assign_vars: dict[Variable, int] = {}
            for var, val in zip(self.vars, result):
                assign_vars[var] = val

            yield Solution(assign_vars, assign_rels)

    def _get_first(self) -> Solution:
        '''Returns the first non-degenerate solution, caching it.'''
        if self._first is None:
            try:
                self._first = next(iter(self))
            except StopIteration:
                raise ValueError('No non-degenerate solution found')
        return self._first

    def __call__(self, obj: Variable | LinearCombination) -> int:
        '''Convenience: delegates to the first non-degenerate solution.'''
        return self._get_first()(obj)

    def __str__(self) -> str:
        '''Convenience: delegates to the first non-degenerate solution.'''
        return str(self._get_first())

class LinearCombination:
    '''
    An immutable linear combination of variables.
    '''
    def __init__(self, ctx: LeLeLe, combine: dict[Variable, int]) -> None:
        assert isinstance(combine, dict), f"combine must be a dictionary, got {type(combine)}"
        self.ctx = ctx
        self.combine = combine
        self._key = frozenset(combine.items())

    def vars(self) -> set[Variable]:
        return set(self.combine.keys())

    def coeff(self, var: Variable) -> int:
        return self.combine.get(var, 0)

    def __len__(self) -> int:
        return len(self.combine)

    def __iter__(self) -> Iterator[tuple[Variable, int]]:
        return iter(sorted(self.combine.items()))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, LinearCombination):
            return False
        if self.ctx is not other.ctx:
            return False
        return self._key == other._key

    def __hash__(self) -> int:
        return hash((id(self.ctx), self._key))

    def __repr__(self) -> str:
        lin = [f'{v!r} * {hex(s)}' if s != 1 else f'{v!r}' for (v, s) in self]
        return ' + '.join(lin)

    def __mod__(self, other) -> LinearCombination:
        try:
            n = int(other)
        except ValueError:
            raise ValueError('Modulo of linear combination only defined for integers')
        combine = dict(self.combine)
        combine[self.ctx.var()] = n
        return LinearCombination(ctx=self.ctx, combine=combine)

    def __neg__(self) -> LinearCombination:
        return LinearCombination(
            ctx=self.ctx,
            combine={s: -v for (s, v) in self.combine.items()}
        )

    def __sub__(self, other) -> LinearCombination:
        return self + (-other)

    def __rsub__(self, other) -> LinearCombination:
        return _wrap_lin(self.ctx, other) - self

    def __add__(self, other) -> LinearCombination:
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

    def __radd__(self, other) -> LinearCombination:
        return self.__add__(other)

    def __mul__(self, other) -> LinearCombination:
        try:
            n = int(other)
        except ValueError:
            raise ValueError(f'Can only mul. linear combination by integer, not {other!r}')
        return LinearCombination(
            ctx=self.ctx,
            combine={v: s * n for (v, s) in self.combine.items()}
        )

    def __rmul__(self, other) -> LinearCombination:
        return self.__mul__(other)

    def short(self, norm: int = 1) -> LinearCombination:
        '''
        Constrain the linear combination to have small norm.

        The "norm" parameter should be understood
        as the "max-value" of the variable/linear combination.
        '''
        if norm <= 0:
            raise ValueError('Norm must be positive')
        self.ctx.add_constraint(self, norm)
        return self

@total_ordering
class Variable:
    def __init__(self, ctx: LeLeLe, name: str, index: int):
        self.ctx = ctx
        self.name = name
        self.index = index

    def lin(self) -> LinearCombination:
        return LinearCombination(ctx=self.ctx, combine={self: 1})

    def key(self) -> tuple[str, int, int]:
        '''
        Returns the ordering key for consistent comparisons.
        '''
        return (self.name, self.index, id(self.ctx))

    def __lt__(self, other: Variable) -> bool:
        if not isinstance(other, Variable):
            return NotImplemented
        return self.key() < other.key()

    def __neg__(self) -> LinearCombination:
        return - self.lin()

    def __sub__(self, other) -> LinearCombination:
        return self.lin() - _wrap_lin(self.ctx, other)

    def __rsub__(self, other) -> LinearCombination:
        return _wrap_lin(self.ctx, other) - self.lin()

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

    def __mod__(self, other) -> LinearCombination:
        return self.lin() % int(other)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Variable):
            return False
        return self.key() == other.key()

    def __hash__(self) -> int:
        return hash(self.key())

    def short(self, norm: int = 1) -> Variable:
        '''
        Constrain the variable to have small norm.
        '''
        self.lin().short(norm)
        return self
