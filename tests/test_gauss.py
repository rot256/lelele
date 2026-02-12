"""
Test that Gaussian elimination can recover variables
that only appear in multi-variable constraints.
"""

import unittest
import lelele as le


class TestGaussianElimination(unittest.TestCase):

    def test_two_unknowns(self):
        m = le.LeLeLe()
        x = m.var()
        y = m.var()
        x_true, y_true = 42, 17

        for (a, b) in [(3, 7), (5, 11), (13, 2), (1, 1)]:
            rhs = a * x_true + b * y_true
            (a * x + b * y - rhs).short(norm=1)

        sol = m.solve()
        self.assertEqual(sol(x), x_true)
        self.assertEqual(sol(y), y_true)

    def test_three_unknowns(self):
        m = le.LeLeLe()
        x = m.var()
        y = m.var()
        z = m.var()
        x_true, y_true, z_true = 10, 20, 30

        for (a, b, c) in [(1, 2, 3), (4, 5, 6), (7, 8, 10), (2, 3, 1)]:
            rhs = a * x_true + b * y_true + c * z_true
            (a * x + b * y + c * z - rhs).short(norm=1)

        sol = m.solve()
        self.assertEqual(sol(x), x_true)
        self.assertEqual(sol(y), y_true)
        self.assertEqual(sol(z), z_true)


if __name__ == '__main__':
    unittest.main()
