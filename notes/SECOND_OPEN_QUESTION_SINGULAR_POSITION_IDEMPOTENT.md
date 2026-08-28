# The singular-position idempotent at 23

## Outcome

The component-cutting half of the second open question can be carried out
directly on the characteristic-23 canonical models, without using the exact
branch-cycle crosswalk.

Let `P` be one of the three closed points of the normalized local Hurwitz
scheme, let `C_P` be its direct reduced canonical curve, and let

```text
beta_P : C_P -> P1
```

be the reduced degree-23 map.  The normalization is `P1`, the marked zero
and pole give `b=0`, `c=infinity`, and there is a unique rational function
`t` such that

```text
beta_P=t^23.
```

There is no scalar ambiguity: in characteristic 23 the group scheme of
23rd roots of unity has only the identity geometric point, and Frobenius is
injective on the normalization's function field.

Define the distinguished singular point `tau_P` to be the unique
singularity in the `E8` cases and the uniquely characterized `A6`
singularity in the `A2+A6` case.  Put

```text
u(P)=t(tau_P).
```

The exact models give the following table.

| closed Hurwitz point | distinguished singularity | `u(P)` | minimal polynomial over `F_23` |
|:--|:--|:--|:--|
| degree-one component | `E8` | `16` | `u-16` |
| unramified degree-2 piece | `E8` | `sqrt(-1)` | `u^2+1` |
| ramified degree-4 piece | `A6` | primitive cube root | `u^2+u+1` |

For the last row the `A2` singularity has parameter `t=1`; the `A6`
singularity has parameter `r`, where `r^2+r+1=0`.  Frobenius exchanges
`r` and `r^2`.  In the unramified `E8` row Frobenius exchanges the two
roots of `u^2+1`.

Consequently the reduced five-geometric-point special fibre is recovered
from the pointed canonical models by the intrinsic resolvent

```text
R_23(u)=(u-16)*(u^2+1)*(u^2+u+1).
```

The three factors have degrees `1,2,2`, exactly the residue degrees of the
three closed points.  They are pairwise coprime.  Thus `u` is a primitive
coordinate for the reduced special fibre, not merely a statistic with the
right degree pattern.

## The Boolean component morphism

In

```text
B=F_23[u]/(R_23)
```

the polynomial

```text
e_sing(u)=2*u^4+2*u^3+4*u^2+2*u+3
```

is idempotent and has values

```text
0 mod (u-16),
1 mod (u^2+1),
1 mod (u^2+u+1).
```

Indeed, if

```text
A(u)=(u^2+1)*(u^2+u+1),
```

then `A(16)=11` and `11^(-1)=21 mod 23`, so

```text
e_sing=1-21*A=1+2*A.
```

Idempotents lift uniquely through nilpotent thickenings and through the
complete local normalization.  Therefore `e_sing` also separates the
degree-one local factor from both local factors of the sextic component,
including the multiplicity-two ramified factor.  Under the already
reconstructed algebra

```text
A_H=K0 x L,
```

its unique lift is the sextic component idempotent `(0,1)`.  This constructs
the Boolean morphism cutting out the degree-six component from the
characteristic-23 models themselves; no Nielsen labels or exact continued
branch cycles enter this calculation.

Equivalently, on the three closed points,

```text
e_sing(P)=f(P)-1 mod 2
```

for residue degrees `f(P)=1,2,2`.

## What remains open

This does **not** yet identify the singular-position idempotent with the
finite-group formula

```text
epsilon(Theta)=|x^N(<y>) intersection x^N(<z>)| mod 2
```

without branch cycles.  The existing exact branch-cycle certificate proves
that both functions equal the sextic idempotent, but that is the comparison
the open question asks to replace by geometry.

The remaining theorem is therefore narrower than the version currently
stated in the paper:

> **Singular-position/normalizer-trace comparison.**  Construct, on the
> pointed wild stable `M23`-map, a finite-branch incidence or vanishing-cycle
> class whose augmentation specializes to `e_sing(u)`.

Resolving the ADE source curve alone cannot prove this theorem.  The
degree-one and unramified quadratic models are both geometrically of type
`E8`, so the resolution graph does not distinguish them.  The distinction
is the arithmetic position of the singularity together with the lifted
Mathieu pointing and gluing.

## Stable-map and patching audit

There are two relevant compactification formalisms, and they should not be
conflated.

1. Abramovich--Oort define a proper mixed-characteristic complete Hurwitz
   stack as the closure of the characteristic-zero Hurwitz locus in a stack
   of stable maps.  In wild characteristic its objects can be nonfinite and
   inseparable.  Properness supplies a carrier for specialization, but not
   the missing internal `M23` pointing.
2. Wewers's special `G`-maps retain the finite `G`-action, the special
   deformation datum, the pointed tails, and the rigid patching data.  This
   is the appropriate formalism for the desired comparison.

For any special `G`-map realizing the certified deformation-datum candidate

```text
a^11=z/(z-2),
omega=4*a^6*da/(a^22-1),
```

the two non-wild conductors are `7` and `15`.  Wewers's patching theorem
then gives

```text
|P(fbar)|=(23-1)*7*15=2310,
Galois orbit length=(23-1)*lcm(7,15)=2310.
```

Thus the unquotiented patching data form one tame orbit; there is no
canonical patching point.  Numerically,

```text
2310=15*154=15*(22*7),
```

and `15` is also the number of relative pointing colors

```text
C_M23(x) backslash M23 / (23:11).
```

This equality is suggestive but is not a construction of a quotient map.
Producing the incidence map from patching data to those fifteen colors, and
then proving that its mod-two pushforward is `e_sing`, is the remaining
geometric task.

## Reproduction

Run

```text
DOT_SAGE=/private/tmp/m23-cover-investigation-sage \
  sage -python verification/verify_hurwitz_pointed_23.py

DOT_SAGE=/private/tmp/m23-cover-investigation-sage \
  sage -python notes/certify_p23_special_deformation_datum.py
```

The first command derives every singularity, the intrinsic normalization
parameter, all four displayed minimal polynomials, the resolvent, and the
Boolean idempotent from the exact canonical models.  The second certifies
the special deformation datum, conductors, and the patching-data counts.
