"""
Seccon 2021 "oooooo" challenge.
"""

import random

from Crypto.Util.number import bytes_to_long, getPrime

from lelele import LeLeLe

random.seed(42)

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

# print the system for debugging
print(le)

# iterate over solutions until we find one with valid bit values
for sol in le.solve():
    vals = set(sol(bi) for bi in b)
    if not (vals <= {0, 1} or vals <= {0, -1}):
        continue
    s = [chr(v1) if sol(bi) else chr(v0) for bi in b][::-1]
    s = ''.join(s)
    if s == message.decode('utf-8'):
        print("recovered:", s)
        break
