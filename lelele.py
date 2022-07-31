SOLVE = True

try:
    '''
    If fpylll is installed LeLeLe can also solve the system.
    '''
    from fpylll import *
except ImportError:
    SOLVE = False

def _wrap_lin(ctx, val):
    assert isinstance(ctx, LeLeLe)

    # linear combination
    if isinstance(val, LinearCombination):
        return val

    # variable
    if isinstance(val, Variable):
        return val.lin()

    # constant
    try:
        c = int(val)
        return LinearCombination(ctx=ctx, combine=[(c, ctx.one())])
    except ValueError:
        raise ValueError('failed to covert %r to linear combination' % val)

def _zero_matrix(n, m):
    M = []
    for _ in range(n):
        M.append([0] * m)
    return M

class LeLeLe:
    def __init__(self):
        self.vars = []
        self.vone = None
        self.constraints = []

    def one(self):
        '''
        Returns a "variable" which should the constant value one.
        '''
        if self.vone is None:
            self.vone = self.bit(name='1')
        return self.vone

    def bit(self, name=None):
        return self.var(name, 'bit').short(norm=1)

    def byte(self, name=None):
        return self.var(name, 'byte').short(norm=0xff)

    def word(self, width, name=None):
        assert width > 0
        return self.var(name, 'word').short(norm=(1 << width)-1)

    def var(self, name=None, prefix='var'):
        var = Variable()
        var.ctx = self
        if name is not None:
            var.name = name
        else:
            var.name = '%s%d' % (prefix, len(self.vars))
        var.index = len(self.vars)
        self.vars.append(var)
        return var

    def add_constraint(self, lin, norm):
        assert isinstance(lin, LinearCombination)
        self.constraints.append((lin, int(norm)))

    def system(self):
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
            for (scale, var) in lin.combine:
                M[var.index][i] = scale * rescale

        return M

    def solve(self):
        '''
        Solves the system and assigns the solution to the variables
        '''

        if not SOLVE:
            raise ImportError('You must install fpylll to use the solve method')

        M = self.system()

        # transformation (used to derieve assignment of variables)
        U = IntegerMatrix.identity(len(M))

        # reduced basis (used to derieve assignment of constraints)
        R = IntegerMatrix.from_matrix(M)

        # run LLL and save the result (for debugging)
        LLL.reduction(R, U)

        self.U = U
        self.R = R
        self.M = IntegerMatrix.from_matrix(M)

        # reset solutions for each variable
        for var in self.vars:
            var.solutions = []

        # reset solutions for each constraint
        for (con, _) in self.constraints:
            con.solutions = []

        # assign values to variables/constraints
        for row, assign in zip(R, U):
            if all([r == 0 for r in row]):
                # trivial solution
                continue

            # sanity check
            assert len(assign) == len(self.vars)
            assert len(row) == len(self.constraints)

            # add solutions to variables
            for var, val in zip(self.vars, assign):
                var.solutions.append(val)

            # add solutions to constraints
            for (con, _), val in zip(self.constraints, row):
                con.solutions.append(val)

        if self.vone and int(self.vone) != 1:
            raise ValueError('.one() constant not assigned 1 in solution (assigned 0x%x); check the norms assigned with .short' % int(self.vone))

        return (self.R, self.U, self.M)

    def __repr__(self):
        cons = []
        for (lin, norm) in self.constraints:
            cons.append('    0x%x = |%s|' % (norm,lin))
        return 'LLLSystem(\n' + '\n'.join(cons) + '\n)'

class LinearCombination:
    def __init__(self, ctx, combine):
        self.ctx = ctx
        self.combine = combine # linear combination
        self.solutions = None

    def __eq__(self, other):
        if not isinstance(other, LinearCombination):
            return False
        if self.ctx != other.ctx:
            return False
        if self.combine != other.combine:
            return False
        return True

    def __repr__(self):
        lin = ['%r * %s' % (v, hex(s)) if s != 1 else '%r' % v for (s, v) in self.combine]
        return ' + '.join(lin)

    def __mod__(self, other):
        try:
            n = int(other)
        except ValueError:
            raise ValueError('Modulo of linear combination only defined for integers')

        return LinearCombination(
            ctx=self.ctx,
            combine=self.combine + [(n, self.ctx.var())]
        )

    def __neg__(self):
        return LinearCombination(
            ctx=self.ctx,
            combine=[(-s, v) for (s, v) in self.combine]
        )

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return self.__sub__(other)

    def __add__(self, other):
        if other == 0: return self # this is convenient

        # combine terms
        # i.e. [(s1, var)] + [(s2, var)] would become [(s1 + s2, var)]
        other = _wrap_lin(self.ctx, other)

        assert self.ctx == other.ctx, 'linear combinations belong to different systems'

        le_vars = {v: s for (s, v) in self.combine}

        for (s, v) in other.combine:
            try:
                le_vars[v] += s
            except KeyError:
                le_vars[v] = s

        return LinearCombination(
            ctx=self.ctx,
            combine= [(s,v) for (v,s) in le_vars.items()]
        )

    def __radd__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        try:
            n = int(other)
        except ValueError:
            raise ValueError('Can only mul. linear combination by integer, not %r' % other)

        return LinearCombination(
            ctx=self.ctx,
            combine=[(s * n, v) for (s, v) in self.combine],
        )

    def __rmul__(self, other):
        return self.__mul__(other)

    def short(self, norm=1):
        '''
        Constrain the linear combination to have small norm.

        The "norm" parameter should be understood
        as the "max-value" of the variable/linear combination.
        '''
        assert norm > 0
        self.ctx.add_constraint(self, norm)
        return self

    def __getitem__(self, n):
        '''
        Return the n'th solution
        '''

        # check if constraint and solution set
        if self.solutions is not None:
            return self.solutions[n]

        # otherwise compute from variable assignments
        return sum([s * v[n] for (s, v) in self.combine])

    def __int__(self):
        return int(self[0])

    def __index__(self):
        return int(self)

class Variable:
    def __init__(self):
        self.solutions = None # not solved

    def lin(self):
        return LinearCombination(ctx=self.ctx, combine=[(0x1, self)])

    def __neg__(self):
        return - self.lin()

    def __sub__(self, other):
        return self.lin() - _wrap_lin(self.ctx, other)

    def __rsub__(self, other):
        return self.lin() - _wrap_lin(self.ctx, other)

    def __add__(self, other):
        return self.lin() + _wrap_lin(other)

    def __radd__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        return self.lin() * int(other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __repr__(self):
        return self.name

    def __index__(self):
        return int(self)

    def __mod__(self, other):
        return self.lin() % int(other)

    def short(self,norm=1):
        '''
        Constrain the variable to have small norm.
        '''
        self.lin().short(norm)
        return self

    def __getitem__(self, n):
        '''
        Return the n'th solution
        '''
        if self.solutions is None:
            raise ValueError('Must solve the system first')
        return self.solutions[n]

    def __int__(self):
        return int(self[0])
