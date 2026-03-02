# Kurz

`Kurz` ("short" in German) is a very simple library (~500 lines) to help you more easily implement lattice attacks, the library is inspired by `Z3Py` (python interface for Z3).
Manually constructing lattices for LLL attacks is usually a messy process of debugging list comprehensions,
`Kurz` solves this by allowing you to simply require that a linear combination of variables is `.short()` and then `.solve()` for concrete values,
the solution is returned as an object and can be queried using `sol(var)`.
`Kurz` turns a hard to understand/debug mess like:

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
from kurz import Kurz

le = Kurz()

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
Kurz(
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
`Kurz` is simply intended to be copy-pasted to the same folder as your exploit.
This is because I do not want to commit to a stable API, so if your exploits should keep working
please copy-paste whenever you need it (essentially vendoring it).

## Requirements

To automatically solve systems, `Kurz` needs a lattice reduction backend.
Two backends are supported:

- [`flatn`](https://pypi.org/project/flatn/) (default)
- [`fpylll`](https://github.com/fplll/fpylll)

Install one of them (e.g. `pip install flatn`).
You can select the backend explicitly:

```python3
sol = le.solve(backend='fpylll')  # or backend='flatn'
```

Without either backend, `Kurz` can still be used to construct the lattice matrix using `.system()` and you can then apply LLL to the resulting lattice using another tool:

```python3
M = le.system()
```

## Tuning Reduction Quality

Both backends accept additional parameters which are forwarded directly to the underlying solver. If a solve fails to find the expected solution, try increasing reduction quality before adding more constraints.

### Flatn

- **`delta`** -- The Lovász condition parameter (0.25 to 1.0). Higher values produce shorter vectors but take longer. Kurz defaults to 0.99; flatter's own default is ~0.62 (rhf 1.0219).
- **`rhf`** -- Root Hermite factor. The shortest vector in the reduced basis satisfies `|b_1| ~ rhf^n * det(L)^{1/n}`. Lower values produce shorter vectors but take longer.
- **`alpha`** -- Flatter's internal slope parameter: `alpha = 2 * log2(rhf)`. Lower values produce shorter vectors but take longer.

```python3
sol = le.solve(rhf=1.01)
```

### Fpylll

- **`delta`** -- The Lovász condition parameter (0.25 to 1.0). Higher values produce shorter vectors but take longer. Defaults to 0.99.
- **`eta`** -- LLL parameter (0 to sqrt(delta)). Defaults to 0.51.

```python3
sol = le.solve(backend='fpylll', delta=0.9999)
```

## Iterating Over Solutions

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
- [`h1_google_2021_sage.sage`](examples/h1_google_2021_sage.sage) -- calling Kurz from Sage.
- [`oooooo_seccon_2021.py`](examples/oooooo_seccon_2021.py) -- Seccon 2021 "oooooo" challenge
