"""
Seccon 2021 "oooooo" challenge.
"""

import random
import unittest

from Crypto.Util.number import long_to_bytes, bytes_to_long, getPrime

from lelele import LeLeLe


class TestOooooo(unittest.TestCase):

    def test_oooooo_seccon_2021(self):
        fails = 10

        for _ in range(128):
            message = b""
            for _ in range(128):
                message += b"o" if random.getrandbits(1) == 1 else b"O"

            M = getPrime(len(message) * 5)
            S = bytes_to_long(message) % M

            LEN = 128

            le = LeLeLe()

            v1 = ord('o')
            v0 = ord('O')

            D = v1 - v0

            b = [le.bit() for _ in range(LEN)]

            v0s = bytes_to_long(bytes([v0] * LEN))
            con = v0s - S

            d = sum([b[i] * 2 ** (i * 8) for i in range(LEN)]) * D + con

            (d % M).short()

            try:
                sol = le.solve()
            except ValueError:
                continue

            s = [chr(v1) if sol(bi) else chr(v0) for bi in b][::-1]
            s = ''.join(s)
            if s == message.decode('utf-8'):
                return

            self.assertGreater(fails, 0, 'too many failed attempts')
            fails -= 1

        self.fail('test failed: no successful solve in 128 iterations')


if __name__ == '__main__':
    unittest.main()
