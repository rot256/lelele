# LeLeLe

`LeLeLe` is a very simple library (~500 lines) to help you more easily implement lattice attacks, the library is inspired by `Z3Py` (python interface for Z3).
Manually constructing lattices for LLL attacks is usually a messy process of debugging list comprehensions,
`LeLeLe` solves this by allowing you to simply require that a linear combination of variables is `.short()` and then `.solve()` for concrete values,
the solution is returned as an object and can be queried using `sol(var)`.
`LeLeLe` turns a hard to understand/debug mess like:

```python3
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
from lelele import LeLeLe

le = LeLeLe()

V = [le.byte() for _ in range(len(ti))]

# define short linear combination mod n
w = sum([t * v for (v, t) in zip(V, ti)]) + inv * u
(w % n).short()

# prints a description of the system for debugging
print(le)

# find a solution
sol = le.solve()

# print values assigned in solution
print(-sol(w), [sol(v) for v in V])
```

Where `print(le)` shows the system of constraints:

```
LeLeLe(
    0xff >= |byte0|
    0xff >= |byte1|
    ...
    0xff >= |byte30|
    0x1 >= |1|
    0x1 >= |1 * <constant> + byte0 * 0x100000000 + ... + var32 * <modulus>|
)
```

## Installation

Rather than installing the library using `pip`,
`LeLeLe` is simply intended to be copy-pasted to the same folder as your exploit.
This is because I do not want to commit to a stable API, so if your exploits should keep working
please copy-paste whenever you need it (essentially vendoring it).

## Requirements

To automatically solve systems, `LeLeLe` needs a lattice reduction backend.
Two backends are supported:

- [`fpylll`](https://github.com/fplll/fpylll) (default)
- [`flatn`](https://github.com/rot256/flatn)

Install one of them (e.g. `pip install fpylll`).
You can select the backend explicitly:

```python3
sol = le.solve(backend='fpylll')  # or backend='flatn'
```

Without either backend, `LeLeLe` can still be used to construct the lattice matrix using `.system()` and you can then apply LLL to the resulting lattice using another tool:

```python3
M = le.system()
```

## Iterating over solutions

`le.solve()` returns a `Solutions` object. You can query the first solution directly:

```python3
sol = le.solve()
print(sol(v))  # value of variable v in the first solution
```

Or iterate over all non-degenerate solutions from the reduced lattice:

```python3
for sol in le.solve():
    print(sol(v))
```

Each `Solution` can be queried for variable values via `sol(var)` or for linear combination values via `sol(expr)`.

## Examples

See the [`examples/`](examples/) directory for full worked examples:

- [`h1_google_2021.py`](examples/h1_google_2021.py) -- H1@Google 2021 CTF challenge
- [`h1_google_2021_sage.sage`](examples/h1_google_2021_sage.sage) -- calling LeLeLe from Sage.
- [`oooooo_seccon_2021.py`](examples/oooooo_seccon_2021.py) -- Seccon 2021 "oooooo" challenge
