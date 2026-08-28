# The local Hurwitz scheme and pointed reductions at 23

Put `K0=Q(sqrt(-23))`, let `q=(sqrt(-23))`, and let

```text
g(y)=y^6-6*y^5+14*y^4-2*y^3-27*y^2+44*y-44.
```

The sextic Hurwitz component is `L=E*K0`, where `E=Q[y]/(g)`.  This note
records its exact local structure at `q` and the direct reductions of the
three closed Hurwitz points of the normalized local model.

## Local decomposition group

Modulo `23`,

```text
g(y)=(y-5)^2*(y+1)^4.
```

Over `Q_23`, the two factors have degrees `2` and `4`, discriminant
valuations `1` and `3`, and discriminant-unit residues `9` and `7`.
Therefore their discriminant square classes are respectively `23` and
`-23`.  The quartic factor is totally and tamely ramified.  Its local
Galois closure over `Q_23` is dihedral of order `8`: tame inertia is `C4`, and Frobenius
acts by inversion because `23=-1 mod 4`.  Its discriminant field is
`Q_23(sqrt(-23))`.

After base change to `K0_q`, the two primes of `L` have

```text
(e,f)=(1,2) and (2,2).
```

The decomposition group in the `S6` closure becomes the even subgroup
`V4` of the local `D4`.  In its natural action on the six points of the
degree-six component it has one orbit of size `2` and one regular orbit of size `4`.
Inertia has cycle type `1^2 2^2`; either Frobenius lift has cycle type
`2^3`.

A labeling of these local orbits by the complex Nielsen IDs is not
canonical.  It requires a choice of a prime of the `S6` closure above
`q`; conjugating that prime conjugates the `V4` inside `S6`.

The first residual equations are also explicit.  If `s=sqrt(-23)`, then

```text
U=(y-5)/s       gives U^2=15 mod q,
V=(y+1)^2/s     gives V^2=5 mod q.
```

Both define `F_23^2`, since `15/5=3=7^2 mod 23`.  Thus the normalization
of the full local Hurwitz scheme has three closed points with relative
local data

```text
degree-one component: (e,f)=(1,1),
degree-2 orbit:  (e,f)=(1,2),
degree-4 orbit:  (e,f)=(2,2).
```

## Direct pointed reductions

The exact canonical quadrics, Petri cubics, marked points, and degree-23
map sections are integral at all three local points.  Their direct special
fibres are geometrically integral `(2,3)` complete intersections on smooth
quadrics.  The marked zero `b` and pole `c` remain distinct smooth points.

The singularities are:

| local Hurwitz point | residue field | singularities | Milnor numbers |
|:--|:--|:--|:--|
| degree-one point (class `6`) | `F_23` | `E8` | `8` |
| unramified degree `2` | `F_23^2` | `E8` | `8` |
| ramified degree `4` | `F_23^2` | `A2 + A6` | `2+6` |

The Tjurina numbers equal the displayed Milnor numbers.  In each row the
total delta invariant is `4`, the arithmetic genus of a `(2,3)` complete
intersection.  Hence the geometric normalization is `P1`.

Let `N,D` be the reduced quintic sections.  Their common base scheme has
Hilbert polynomial `7`, it avoids every singular point, and

```text
ord_b(N)=23,  ord_b(D)=0,
ord_c(D)=23,  ord_c(N)=0.
```

Consequently, on the normalization,

```text
div(N/D)=23*b-23*c.
```

After choosing a coordinate `t` with `b=0` and `c=infinity`, and rescaling
it over the perfect residue field, every direct reduced pointed
map is

```text
beta=t^23.
```

The coordinate is unique once this equation is imposed.  It gives an
intrinsic singular-position coordinate on the reduced local Hurwitz fibre.
Choose the unique `E8` singularity in the first two rows and the uniquely
characterized `A6` singularity in the third, and evaluate `t` there.  The
exact values have minimal polynomials

```text
degree-one E8:          u-16,
unramified-degree-2 E8: u^2+1,
ramified-degree-4 A6:   u^2+u+1.
```

(The `A2` point in the last row has `t=1`.)  Hence the reduced five-point
special fibre is recovered from the pointed models by

```text
R23(u)=(u-16)*(u^2+1)*(u^2+u+1).
```

In `F_23[u]/(R23)`, the polynomial

```text
e_sing=2*u^4+2*u^3+4*u^2+2*u+3
```

is idempotent with factorwise values `(0,1,1)`.  It lifts uniquely through
the ramified nilpotent thickening and the complete local normalization, and
its lift is the sextic component idempotent.  This constructs the Boolean
component morphism directly from the characteristic-23 models, without a
complex Nielsen labeling or the exact branch-cycle crosswalk.  See
`notes/SECOND_OPEN_QUESTION_SINGULAR_POSITION_IDEMPOTENT.md`.

This is a useful warning: direct reduction and normalization erase the
separable Mathieu tail.  The positive genus is concentrated in the ADE
singularities.  Recovering the stable `M23`-map requires resolving their
mixed-characteristic deformations and retaining the gluing and pointing
data.  The computation therefore narrows the missing argument but does not
prove the universal pointed parity lemma.

## The binary augmentation locally

On the three closed points of the normalized local Hurwitz scheme, let
`f(P)` be the residue degree.  The singular-position idempotent gives

```text
e_sing(P) = f(P)-1 mod 2,
```

because the two sides have values `(0,1,1)` on residue degrees `(1,2,2)`.
Under the already certified Hurwitz-algebra and branch-cycle identification,
`epsilon(Theta)=e_sing`.  What is still missing is an intrinsic geometric
proof of that equality: the stable Mathieu incidence cycle has not yet been
constructed from the singular position and pointed gluing.

Run the exact certificates with

```text
make verify-hurwitz-local-23
```

The local number-field and decomposition-group calculation is in
`verification/verify_hurwitz_local_23.py`.  The integral canonical models,
singularity types, base divisor, and pointed Frobenius reductions are
checked by `verification/verify_hurwitz_pointed_23.py`.
