import random

LEN = 128

message = b""
for _ in range(LEN):
    message += b"o" if random.getrandbits(1) == 1 else b"O"

M = 2811830987405905292751966156324695448083623705675212665673615457818396286443991422909672526682947555018391066902157822920734594475696540429928336527302857491989936012335568308434637664589366489
S = int.from_bytes(message, 'big') % M

from lelele import *

le = LeLeLe()

v1 = ord('o')
v0 = ord('O')

D = v1 - v0

b = [le.var().short() for _ in range(LEN)]

v0s = int.from_bytes(bytes([v0] * LEN), 'big')
con = v0s - S # the constant

C = le.var().short()

d = sum([b[i] * 2^(i * 8) for i in range(LEN)]) * D
v = (C * con + d) % M

v.short()

le.solve()

s = [chr(v1) if int(bi) else chr(v0) for bi in b][::-1]

s = ''.join(s)

assert s == message.decode('utf-8')
