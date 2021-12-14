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

l = LeLeLe()

v1 = ord('O')
v2 = ord('o')

D = v1 ^ v2
C = v1 & v2

com = bytes([C] * LEN)
com_M = bytes_to_long(com)

b = [l.var().short() for _ in range(LEN)]

cons = com_M - S

con_var = l.var().short()

d = sum([b[i] * 2**(i * 8) for i in range(LEN)]) * D
v = (d + con_var * cons) % M

v.short()

l.solve()

s = [chr(C + D) if int(bi) else chr(C) for bi in b][::-1]

s = ''.join(s)

assert s == message.decode('utf-8')
