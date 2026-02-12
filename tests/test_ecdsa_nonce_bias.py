"""
ECDSA nonce bias attack (Hidden Number Problem).

Ported from https://github.com/bitlogik/lattice-attack

When a few bits of each ECDSA nonce k are leaked, the private key
can be recovered via lattice reduction. Given ECDSA signatures where:

    s_i = k_i^{-1} * (h_i + r_i * d)  (mod n)

and the lower `l` bits of each k_i are known (call them kp_i), we have:

    k_i = 2^l * u_i + kp_i

where u_i is the unknown upper part, bounded by n / 2^l.
Rearranging gives the Hidden Number Problem:

    u_i = B^{-1} * (r_i * s_i^{-1} * d + s_i^{-1} * h_i - kp_i)  (mod n)

which is solvable by LLL lattice reduction.
"""

import random
import hashlib
import unittest

from lelele import LeLeLe

# secp256k1 parameters
P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
GX = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
GY = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8


def _modinv(a, m):
    return pow(a, -1, m)


def _ec_add(p1, p2):
    if p1 is None:
        return p2
    if p2 is None:
        return p1
    x1, y1 = p1
    x2, y2 = p2
    if x1 == x2 and y1 == y2:
        lam = (3 * x1 * x1) * _modinv(2 * y1, P) % P
    elif x1 == x2:
        return None
    else:
        lam = (y2 - y1) * _modinv(x2 - x1, P) % P
    x3 = (lam * lam - x1 - x2) % P
    y3 = (lam * (x1 - x3) - y1) % P
    return (x3, y3)


def _ec_mul(k, point):
    result = None
    addend = point
    while k:
        if k & 1:
            result = _ec_add(result, addend)
        addend = _ec_add(addend, addend)
        k >>= 1
    return result


G = (GX, GY)


def _ecdsa_sign(privkey, msg_hash, nonce):
    R = _ec_mul(nonce, G)
    r = R[0] % N
    s = _modinv(nonce, N) * (msg_hash + r * privkey) % N
    return r, s


class TestEcdsaNonceBias(unittest.TestCase):

    def _recover_privkey(self, sigs, leaked_bits, pubkey):
        """Build HNP lattice and recover private key using LeLeLe.

        Args:
            sigs: list of (r, s, msg_hash, leaked_nonce_lsbs) tuples
            leaked_bits: number of known LSBs per nonce
            pubkey: target public key for verification

        Returns:
            Recovered private key, or None on failure.
        """
        le = LeLeLe()
        d = le.var(name='privkey')

        B = 1 << leaked_bits
        Binv = _modinv(B, N)

        for r_i, s_i, h_i, kp_i in sigs:
            sinv = _modinv(s_i, N)
            a_i = Binv * r_i * sinv % N
            b_i = Binv * (sinv * h_i - kp_i) % N
            ((a_i * d + b_i) % N).short(N >> leaked_bits)

        d.short(N)

        for sol in le.solve():
            candidate = sol(d) % N
            if candidate > 0 and _ec_mul(candidate, G) == pubkey:
                return candidate
        return None

    def _gen_sigs(self, privkey, num_sigs, leaked_bits):
        """Generate ECDSA signatures with leaked nonce LSBs."""
        B = 1 << leaked_bits
        sigs = []
        for _ in range(num_sigs):
            msg = random.randbytes(32)
            h = int(hashlib.sha256(msg).hexdigest(), 16)
            k = random.randrange(1, N)
            r, s = _ecdsa_sign(privkey, h, k)
            sigs.append((r, s, h, k % B))
        return sigs

    def test_ecdsa_nonce_bias_8bit(self):
        """Recover secp256k1 private key from 50 sigs with 8-bit nonce leak."""
        random.seed(1337)
        privkey = random.randrange(1, N)
        pubkey = _ec_mul(privkey, G)
        sigs = self._gen_sigs(privkey, num_sigs=50, leaked_bits=8)
        recovered = self._recover_privkey(sigs, leaked_bits=8, pubkey=pubkey)
        self.assertEqual(recovered, privkey)

    def test_ecdsa_nonce_bias_7bit(self):
        """Recover secp256k1 private key from 60 sigs with 7-bit nonce leak."""
        random.seed(42)
        privkey = random.randrange(1, N)
        pubkey = _ec_mul(privkey, G)
        sigs = self._gen_sigs(privkey, num_sigs=60, leaked_bits=7)
        recovered = self._recover_privkey(sigs, leaked_bits=7, pubkey=pubkey)
        self.assertEqual(recovered, privkey)


if __name__ == '__main__':
    unittest.main()
