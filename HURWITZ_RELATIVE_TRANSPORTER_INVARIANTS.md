# Relative-transporter invariants on the seven-point Hurwitz scheme

Let `K0=Q(sqrt(-23))`, and let `H_in` be the reduced zero-dimensional inner
Hurwitz scheme of normalized covers with branch-cycle classes
`(2A,23A,23B)`.  Certified continuation of the seven reconstructed maps
identifies its coordinate algebra as

```text
A_H = K0 x L,       [L:K0] = 6,
```

where `L/K0` is a separable field extension.  Thus
`H_in = Spec(A_H)` is finite étale over `K0`, with two connected components:
the degree-one `K0`-point and `Spec(L)`.

For clarity, a function from the seven geometric points to a finite set `S`
is said to *descend over K0* when it is Galois invariant.  Equivalently, it
comes from a `K0`-morphism from `H_in` to the constant finite étale scheme
associated with `S`.  Such a function must be constant on the six geometric
points of the connected scheme `Spec(L)`.

At the prime above `23`, the normalization of the local Hurwitz scheme has
closed-point residue degrees `(1,2,2)`: the point on the degree-one component has degree `1`,
while the two local pieces of the sextic component both have degree `2`.
Consequently the augmentation has the local arithmetic description

```text
epsilon(Theta)(P) = residue_degree(P)-1 mod 2.
```

This equality uses the already identified Hurwitz algebra.  Comparing it
intrinsically with the pointed stable reductions is still the missing
connector; see `HURWITZ_LOCAL_23.md`.

## Definitions of the finite-group invariants

For a geometric Nielsen triple `(x,y,z)`, put `Y=<y>` and `Z=<z>`.  Define

```text
nu(Y,Z) = #{(a,b) in (F_23^times)^2 : order(y^a*z^b)=23}.
```

For `c` equal to `y` or `z` and `d` in `F_23^times`, put

```text
a_c(d) = order(x*x^(c^d)),
kappa_m = sum_{c in {y,z}} #{d : a_c(d)=m},       m in {2,4}.
```

Changing a generator of `Y` or `Z` only reindexes the exponents.  Summing
over both endpoints also makes `(kappa_2,kappa_4)` invariant under endpoint
exchange.  These are cyclic-conjugacy order counts; no additional notion of
a “packet” is required.

Write `N(Y)` and `N(Z)` for the normalizers in `M23`, use
`x^g=g^(-1)*x*g`, and let

```text
Omega = x^N(Y) intersection x^N(Z).
```

The stabilizer of `x` in either normalizer is trivial because its order
divides both `253` and `2688`.  Hence, for every `w` in `Omega`, there are
unique transporter elements `n_Y(w)` and `n_Z(w)` such that

```text
x^n_Y(w) = w = x^n_Z(w).
```

Their relative transporter

```text
c_w = n_Y(w)*n_Z(w)^(-1)
```

centralizes `x`.  Let `A_x` be the eight two-element orbits of `<x>` on its
moved letters.  Modulo `<x>`, the element `c_w` induces a permutation
`cbar_w` of `A_x`; the acting quotient is

```text
C_M23(x)/<x> = 2^3:PSL(3,2) = AGL(3,2).
```

Define

```text
mu = #{cbar_w : w in Omega},
Theta = sum_{w in Omega} [cbar_w]
        in F_2[C_M23(x)/<x>].
```

The standard group-algebra augmentation satisfies

```text
epsilon(Theta) = |Omega| mod 2.
```

The Hamming weight `wt(Theta)` is the number of nonzero coefficients in the
natural group basis.  Endpoint exchange sends each basis element to its
inverse, so it preserves `mu`, `wt(Theta)`, and `epsilon(Theta)`.

## Exact descent test

The six embeddings of `L`, in the order used by the Arb certificate, have
Nielsen IDs

```text
(7,4,1,5,3,2).
```

The exact finite-group values are:

| ID | `nu` | `mu` | `(kappa_2,kappa_4)` | `wt(Theta)` | `epsilon(Theta)` |
|---:|---:|---:|---:|---:|---:|
| 1 | 54 | 16 | `(0,16)` | 15 | 1 |
| 2 | 54 | 16 | `(0,16)` | 15 | 1 |
| 3 | 28 | 31 | `(0,8)` | 25 | 1 |
| 4 | 46 | 28 | `(2,6)` | 27 | 1 |
| 5 | 46 | 28 | `(2,6)` | 27 | 1 |
| 6 | 32 | 17 | `(0,12)` | 16 | 0 |
| 7 | 42 | 14 | `(4,12)` | 11 | 1 |

Thus `nu`, `mu`, and `(kappa_2,kappa_4)` do not descend to maps from
`Spec(L)` to constant finite schemes: their values vary on its six geometric
points.  The varying Hamming weights also show that `Theta` cannot descend
as a section of a lisse `F_2`-sheaf whose transition maps permute the natural
group-algebra basis.  This last conclusion is intentionally limited to
basis-compatible descent; arbitrary linear changes of basis need not
preserve Hamming weight.

The augmentation is different:

```text
epsilon(Theta) = |Omega| mod 2 = (1,1,1,1,1,0,1)
```

in deterministic ID order.  It is `0` on the degree-one factor and `1` on
all six embeddings of `L`, so it does descend.

## Boolean component idempotents

Let

```text
e_pub = (1,0),       e_6 = (0,1)
```

be the primitive idempotents of `A_H=K0 x L`.  On all seven geometric points,
the following Boolean functions agree:

```text
[nu=32]
  = [mu=17]
  = [(kappa_2,kappa_4)=(0,12)]
  = [epsilon(Theta)=0]
  = e_pub.
```

Equivalently,

```text
epsilon(Theta) = e_6 = 1-e_pub.
```

The four indicator functions therefore descend to morphisms
`H_in -> {0,1}` and correspond to `e_pub`.  The augmentation corresponds to
`e_6`.  These idempotents separate the two connected components and hence
give the finest clopen partition of `H_in`.

There is an important logical limitation.  This conclusion uses the already
certified `1+6` decomposition.  It does not prove that decomposition or
explain the rational point without the exact maps and their branch cycles.
An independent explanation would have to construct `epsilon(Theta)` in a
Galois-compatible characteristic-23 model and prove directly that its
Boolean morphism cuts out `Spec(L)`, while discarding the richer assignments
that fail to descend.

Run the exact checks with

```text
make verify-hurwitz-relative-transporter
```

The GAP certificate recomputes every finite-group value.  The Sage
certificate verifies irreducibility of the sextic factor, imports the exact
Arb class-to-embedding crosswalk, and proves the component-idempotent
statement.
