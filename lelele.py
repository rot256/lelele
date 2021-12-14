SOLVE = True

try:
    '''
    If fpylll is installed LeLeLe can also solve the system.
    '''
    from fpylll import *
except ImportError:
    SOLVE = False

def _wrap_lin(var):
    if isinstance(var, LinearCombination):
        return var
    if isinstance(var, Variable):
        return LinearCombination(ctx=var.ctx, combine=[(0x1, var)])

    return int(var)

def _zero_matrix(n, m):
    M = []
    for _ in range(n):
        M.append([0] * m)
    return M

class LeLeLe:
    def __init__(self):
        self.vars = []
        self.constraints = []

    def var(self, name=None):
        var = Variable()
        var.ctx = self
        if name is not None:
            var.name = name
        else:
            var.name = 'v%d' % len(self.vars)
        var.index = len(self.vars)
        self.vars.append(var)
        return var

    def add_constrain(self, lin):
        self.constraints.append(lin)

    def system(self):
        '''
        Returns the Matrix representing the LLL system.

        This does not require fpylll.
        '''

        rows = len(self.vars)
        cols = len(self.constraints)

        M = _zero_matrix(rows, cols)

        for (i, con) in enumerate(self.constraints):
            for (scale, var) in con.combine:
                M[var.index][i] = scale

        return M

    def solve(self, norm=None):
        '''
        Solves the system and assigns the solution to the variables
        '''

        if not SOLVE:
            raise ImportError('You must install fpylll to use the solve method')

        M = self.system()

        # transformation (used to derieve assignment of variables)
        U = IntegerMatrix.identity(len(M))

        # reduced basis (used to derieve assignment of constrains)
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
        for con in self.constraints:
            con.solutions = []

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
            for con, val in zip(self.constraints, row):
                con.solutions.append(val)

        return (self.R, self.U, self.M)

    def __repr__(self):
        cons = []
        for con in self.constraints:
            cons.append('    [%s]' % con)
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
        lin = ['%s * %r' % (hex(s), v) if s != 1 else '%r' % v for (s, v) in self.combine]
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
        assert isinstance(other, LinearCombination), 'can only add linear combination to linear combination, not %s' % type(other)
        assert self.ctx == other.ctx, 'linear combinations belong to different systems'
        return LinearCombination(
            ctx=self.ctx,
            combine=self.combine + other.combine,
        )

    def __radd__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        try:
            n = int(other)
        except ValueError:
            raise ValueError('Can only mul. linear combination by integer')

        return LinearCombination(
            ctx=self.ctx,
            combine=[(s * n, v) for (s, v) in self.combine],
        )

    def __rmul__(self, other):
        return self.__mul__(other)

    def short(self):
        '''
        Constrain the linear combination to have small norm.
        '''
        self.ctx.add_constrain(self)
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
        return self[0]

class Variable:
    def __init__(self):
        self.solutions = None # not solved

    def __sub__(self, other):
        return _wrap_lin(self) - _wrap_lin(other)

    def __neg__(self):
        return - _wrap_lin(self)

    def __rsub__(self, other):
        return _wrap_lin(self) - _wrap_lin(other)

    def __add__(self, other):
        return _wrap_lin(self) + _wrap_lin(other)

    def __radd__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        return _wrap_lin(self) * _wrap_lin(other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __repr__(self):
        return self.name

    def short(self):
        '''
        Constrain the variable to have small norm.
        '''
        _wrap_lin(self).short()
        return self

    def __getitem__(self, n):
        '''
        Return the n'th solution
        '''
        if self.solutions is None:
            raise ValueError('Must solve the system first')
        return self.solutions[n]

    def __int__(self):
        return self[0]
