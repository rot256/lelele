import random
import unittest
import lelele as le


class TestTiny(unittest.TestCase):

    def test_tiny(self):
        v0 = 0
        v1 = 1

        m = le.LeLeLe()

        v = m.var()
        b0 = m.bit()
        b1 = m.bit()

        for _ in range(5):
            s0 = random.randrange(0, 10000000000)
            s1 = random.randrange(0, 10000000000)
            rs = s0 * v0 + s1 * v1
            (s0 * b0 + s1 * b1 - rs).short(1)

        (1000 * v - 1000 * b0).short(1)

        sol = m.solve()

        self.assertEqual(sol(b0), v0)
        self.assertEqual(sol(b1), v1)
        self.assertEqual(sol(v), sol(b0))


if __name__ == '__main__':
    unittest.main()
