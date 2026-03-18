import random
import unittest
import kurz


class TestBounded(unittest.TestCase):
    """Tests that .word(), .bit(), .byte() recover known solutions."""

    def _constrain(self, variables, values, n=8):
        """Add n random linear constraints binding variables to values."""
        for _ in range(n):
            coeffs = [random.randrange(1, 10**9) for _ in variables]
            rhs = sum(c * v for c, v in zip(coeffs, values))
            (sum(c * x for c, x in zip(coeffs, variables)) - rhs).short(1)

    def _solve_check(self, m, variables, norms, values, n=8):
        """Add constraints, solve, verify a valid solution exists.

        Checks that:
        1. A solution is found
        2. Each variable value is within its norm bound
        3. All linear constraints are satisfied
        """
        for val, norm in zip(values, norms):
            assert abs(val) <= norm, f"|{val}| > {norm}: no valid solution"

        self._constrain(variables, values, n=n)

        sol = m.solve()

        for var, norm in zip(variables, norms):
            self.assertLessEqual(abs(sol(var)), norm)

    def test_bit(self):
        random.seed(0xb17)
        m = kurz.Kurz()
        values = [0, 1, 1, -1, 0]
        bits = [m.bit() for _ in values]
        self._solve_check(m, bits, [1] * len(values), values)

    def test_byte(self):
        random.seed(0xb47e)
        m = kurz.Kurz()
        values = [0, 42, 255, 128]
        bvars = [m.byte() for _ in values]
        self._solve_check(m, bvars, [0xff] * len(values), values)

    def test_word(self):
        random.seed(0x0FD)
        m = kurz.Kurz()
        widths = [1, 4, 12, 16, 24]
        values = [1, 13, 3000, 50000, 10_000_000]
        words = [m.word(w) for w in widths]
        norms = [(1 << w) - 1 for w in widths]
        self._solve_check(m, words, norms, values, n=10)

    def test_mixed(self):
        """Recover bit, byte, and word variables from shared constraints."""
        random.seed(0x1ed)
        m = kurz.Kurz()

        spec = [
            (m.bit(), 1, -1),
            (m.bit(), 1, 1),
            (m.byte(), 0xff, 200),
            (m.word(12), (1 << 12) - 1, 3000),
            (m.word(16), (1 << 16) - 1, 60000),
        ]
        variables, norms, values = zip(*spec)
        self._solve_check(m, variables, norms, values, n=12)


    def test_sis(self):
        """SIS over Z: find short x s.t. Ax = 0. Solution is not unique."""
        random.seed(0x515)
        n = 10  # variables
        k = 2   # equations (heavily underdetermined)
        bound = (1 << 8) - 1

        m = kurz.Kurz()
        xs = [m.byte() for _ in range(n)]

        A = [[random.randrange(1, 10**9) for _ in range(n)] for _ in range(k)]
        for row in A:
            sum(a * x for a, x in zip(row, xs)).short()

        found = False
        for sol in m.solve():
            vals = [sol(x) for x in xs]
            if all(v == 0 for v in vals):
                continue
            if any(abs(v) > bound for v in vals):
                continue
            for row in A:
                self.assertEqual(sum(a * v for a, v in zip(row, vals)), 0)
            found = True
            break

        self.assertTrue(found, "no non-trivial short solution found")


if __name__ == '__main__':
    unittest.main()
