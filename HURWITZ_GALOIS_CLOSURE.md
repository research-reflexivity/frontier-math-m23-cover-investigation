# Galois closure of the sextic Hurwitz component

Put `K0=Q(sqrt(-23))`.  The connected degree-six component of the inner
Hurwitz scheme has function field `L/K0`.  Its degree-12 absolute defining
polynomial is

```text
f(a) = a^12 - 6*a^11 + 20*a^10 - 32*a^9 + 44*a^8 - 22*a^7
       + 6*a^6 - 22*a^5 + 44*a^4 - 32*a^3 + 20*a^2 - 6*a + 1.
```

This polynomial is reciprocal.  For `y=a+a^(-1)`, exact elimination gives

```text
g(y) = y^6 - 6*y^5 + 14*y^4 - 2*y^3 - 27*y^2 + 44*y - 44,
f(a) = a^6*g(a+a^(-1)).
```

Let `E=Q[y]/(g)`.  The exact relative-to-absolute field map reconstructed
from `data/hurwitz_algebra_candidate.json` places `sqrt(-23)` in `Q(a)`.
The element `y+sqrt(-23)` has degree 12, so

```text
L = E*K0.
```

The maximal-order discriminants are

```text
disc(E) = 2^4*11*23^4,
disc(L/Q) = 2^8*11^2*23^8,
disc(L/K0) = (2^4*11*23).
```

## The Galois group

The following rational primes split in `K0` and do not divide the
discriminant of `g`.  The factor degrees of `g` modulo them are

| prime | factor degrees | Frobenius cycle type |
|---:|:---|:---|
| `3` | `6` | `6` |
| `139` | `1+5` | `1+5` |
| `2671` | `1+1+1+1+2` | `1+1+1+1+2` |

Irreducibility modulo `3` proves that `g` is irreducible.  An exhaustive
census of the sixteen transitive subgroups of `S6` shows that only `6T16`,
the full symmetric group `S6`, contains all three displayed cycle types.
Sage's exact Galois-group routine and PARI/GP's independent `polgalois`
calculation both return the same group.

Let `N` be the splitting field of `g`.  Since `Gal(N/Q)=S6`, its unique
quadratic subfield is its discriminant field.  The square class of
`disc(g)=2^22*11*23^4` is `11`, so that subfield is `Q(sqrt(11))`, not
`K0`.  Therefore `N` and `K0` are linearly disjoint.  It follows that

```text
the Galois closure of L/K0 is N*K0,
Gal(N*K0/K0) = S6,
Gal(N*K0/Q) = S6 x C2.
```

In the certified sextic-embedding order, the Nielsen IDs are

```text
(7,4,1,5,3,2).
```

Thus the relative Galois group acts as the full symmetric group on the six
maps on the degree-six component.  In particular, the action has no nontrivial block system
and `L/K0` has no proper intermediate fields.  Under the class labeling, the
chosen complex-conjugation comparison has permutation `(1,2)(4,5)`, fixing
IDs `3` and `7` on the sextic component (and ID `6` on the degree-one
component).  Complex conjugation does not belong to the relative group over
`K0`; it also applies the nontrivial automorphism of `K0/Q`.

At the prime above `23`, the relative decomposition group is the Klein four
subgroup with orbit sizes `2+4`, and the two primes of `L/K0` have
`(e,f)=(1,2),(2,2)`.  The exact pointed reductions and their `E8` versus
`A2+A6` singularities are recorded in `HURWITZ_LOCAL_23.md`.

Run the open-source certificates with

```text
make verify-hurwitz-galois-closure
```

The optional independent Magma calculation is

```text
magma -b verification/verify_hurwitz_galois_closure.m
```
