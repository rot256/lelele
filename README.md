# LeLeLe

`LeLeLe` is a very simple library (<300 lines) to help you more easily implement lattice attacks, the library is inspired by `Z3Py` (python interface for Z3).
Manually constructing lattices for LLL attacks is usually a messy process of debugging list comprehensions,
`LeLeLe` solves this by allowing you to simply require that a linear combination of variables is `.short()` and then `.solve()` for concrete values,
the solution is assigned to the variables and can be retrieved by using `int(var)`.
`LeLeLe` turns a hard to understand/debug mess like (example from [H1@ Google 2021 Writeup](https://rot256.dev/post/h1/)):

```sage
cols = (L // B) * 2 + 1
M = []

# short mod n, so first column should contain a vector (n, 0, ..., 0)
M.append([n] + (cols - 1) * [0])

# require that |v_i| are short and add ti[i] * v to the short linear combination
# using a vector (ti[i], 0, ..., 0, 1, 0, ..., 0)
for i, v in enumerate(ti[1:]):
    M.append([v] + [0] * i + [1] + [0] * (cols - i - 2))

# add the final u term which should occure at most once
# to do this add (u*inv, 0, ..., 0, 2^8)
M.append([int(u * inv)] + [0] * (cols - 2) + [K])

# print the matrix for debugging
M = Matrix(M)
print(M)

# run LLL
row = M.LLL()[0]

# print solution
row[0] = -row[0]
print(row)
```

Into a more readable:

```python3
from lelele import *

le = LeLeLe()

q = le.var()
V = [le.var().short() for _ in range(len(ti))] # short variables

# define short linear combination mod n
w = sum([t*v for (v, t) in zip(V, ti)]) + inv * u * q
w %= n
w.short()

# q should be taken at most once: require that q * <<large number>> is small
(q * 0x100).short()

# prints a description of the system
print(le)

# find a solution
le.solve()

# print values assigned in solution
print(-int(w), [int(v) for v in V])
```

## Installation

Rather than installing the library using `pip`,
`LeLeLe` is simply intended to be copy-pasted to the same folder as your exploit.
This is because I do not want to commit to a stable API, so if your exploits should keep working
please copy-paste whenever you need it (essentially vendoring it).

## Requirements

It is recommended to install `fpylll`, such that `LeLeLe` can also be used to solve the system and automatically assign the solution to all the free variables.
`LeLeLe` does not require SageMath.

Without `fpylll`, `LeLeLe` can still be used to construct the lattices using `.system()` and you can then apply LLL to the resulting lattice using another tool:

```python3
from lelele import *

le = LeLeLe()

q = le.var()
V = [le.var().short() for _ in range(len(ti))] # short variables

# define short linear combination mod n
w = sum([t*v for (v, t) in zip(V, ti)]) + inv * u * q
w %= n
w.short()

# q should be taken at most once: require that q * <<large number>> is small
(q * 0x100).short()

# export lattice, a list of lists of ints: [[int]]
M = le.system()
```
