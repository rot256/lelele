'''
A very basic example: the oooooo challenge at Seccon 2021
'''

## Challenge ###

from Crypto.Util.number import long_to_bytes, bytes_to_long, getPrime

import random

for _ in range(100):

    message = b""
    for _ in range(128):
        message += b"o" if random.getrandbits(1) == 1 else b"O"

    M = getPrime(len(message) * 5)
    print(M)
    S = bytes_to_long(message) % M

    # We are given M and S (and must recover "message")

    ## Solve Script ###

    from lelele import *

    LEN = 128

    le = LeLeLe()

    v1 = ord('o')
    v0 = ord('O')

    D = v1 - v0

    b = [le.bit() for _ in range(LEN)]

    v0s = bytes_to_long(bytes([v0] * LEN))
    con = v0s - S # the constant

    d = sum([b[i] * 2**(i * 8) for i in range(LEN)]) * D + con

    (d % M).short() # affine conbination is short mod M

    try:
        print(le)
        le.solve() # unfortunately only about 1/2 of runs of this challenge are solvable
    except ValueError:
        continue

    s = [chr(v1) if int(bi) else chr(v0) for bi in b][::-1]

    s = ''.join(s)

    assert s == message.decode('utf-8')
    break

else:
    assert False, 'test failed'
