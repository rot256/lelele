# LeLeLe

`LeLeLe` is a super simple library to help you implement lattice attacks inspired by `Z3Py` (python interface for Z3):
manually constructing lattices for LLL attacks is usually a messy process of debugging list comprehensions,
`LeLeLe` solves this by allowing you to simply require that a linear combination of variables are `.is_short()` and then `.solve()` for concrete values.
Turning a hard to understand/debug mess like (example from [H1@ Google 2021 Writeup](https://rot256.dev/post/h1/), using Sage):

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

# run LLL
M = Matrix(M)
row = M.LLL()[0]
row[0] = -row[0]

print(row)
```

Into:

```python3
q = le.var()
V = [le.short() for _ in range(len(ti))] # short variables (sugar for .is_short on a var)

# define short linear combination mod n
w = sum([t*v for (v, t) in zip(V, ti)]) + inv * u * q
w %= n
w.short()

# q should be taken at most once: require that q * <<large number>> is small
(q * 0x100).short()

# find a solution
le.solve()

# print values assigned in solution
print(-int(w), [int(v) for v in V])
```

## Requirements

It is recommended to install `fpylll`, such that `LeLeLe` can also be used to solve the system and automatically assign the solution to all the free variables.
`LeLeLe` does not require SageMath.

Without `fpylll`, `LeLeLe` can still be used to construct the lattices using `.system()` and then apply LLL to the resulting lattice using another tool e.g.

```
q = le.var()
V = [le.short() for _ in range(len(ti))] # short variables (sugar for .is_short on a var)

# define short linear combination mod n
w = sum([t*v for (v, t) in zip(V, ti)]) + inv * u * q
w %= n
w.short()

# q should be taken at most once: require that q * <<large number>> is small
(q * 0x100).short()

# export lattice, a list of lists of ints: [[int]]
M = le.system()
```

