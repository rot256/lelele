'''
A very basic example: the oooooo challenge at Seccon 2021
'''

## Challenge ###

from Crypto.Util.number import long_to_bytes, bytes_to_long, getPrime

import random

message = b""
for _ in range(128):
    message += b"o" if random.getrandbits(1) == 1 else b"O"

M = getPrime(len(message) * 5)
S = bytes_to_long(message) % M

# We are given M and S (and must recover "message")

## Solve Script ###

from lelele import *

LEN = 128

le = LeLeLe()

v1 = ord('o')
v0 = ord('O')

D = v1 - v0

b = [le.var().short() for _ in range(LEN)]

v0s = bytes_to_long(bytes([v0] * LEN))
con = v0s - S # the constant

C = le.var().short()

d = sum([b[i] * 2**(i * 8) for i in range(LEN)]) * D
v = (C * con + d) % M

v.short()

le.solve()

s = [chr(v1) if int(bi) else chr(v0) for bi in b][::-1]

s = ''.join(s)

assert s == message.decode('utf-8')
