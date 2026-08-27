# Certified branch cycles of the seven exact maps

The reconstructed degree-23 maps now have certified branch-cycle
permutations around `0`, `1`, and `infinity`.  This closes the earlier gap
between an exact passport and an identified `M23` Nielsen class.

## Exact degree-23 eliminant

On the sextic Hurwitz component, put the canonical curve in the affine chart
`X0=1` and eliminate `X3=z` from its quadric and Petri cubic.  If the linear
subresultant is

```text
L1(x,y) z + L0(x,y),
```

then `z=-L0/L1`.  The exact quintic numerator and normalized denominator are
linear in `z`.  Writing either section as `S0+z*S1`, its plane form is

```text
S0*L1 - S1*L0.
```

This avoids generic fraction-field simplification and gives a plane sextic
and two plane forms of degree seven in about six seconds on the reference
machine.

For a target value `t`, eliminate `y` from the plane sextic and

```text
Nplane - t*Dplane.
```

The raw resultant `R_t(x)` has degree 42.  The exact gcd of `R_2` and `R_3`
has degree 19 and is independent of `t`; it is the plane-projection
base-point contribution.  Division leaves a degree-23 polynomial.  Exact
resultants at `t=2,...,8` interpolate

```text
P(t,x) in L[t,x],   deg_t(P)=6,   deg_x(P)=23.
```

The unused value `t=9` agrees exactly.  The specializations at `0`, `1`, and
the leading target coefficient have derivative-gcd degrees `22`, `8`, and
`22`; the degree-eight gcd at `1` is squarefree.  Thus the eliminant itself
has exact passport

```text
(23), (2^8 1^7), (23).
```

The exact coefficients are in
`data/hurwitz_monodromy_eliminant_candidate.json`.  Recompute one raw
resultant with

```text
make hurwitz-monodromy-resultant HURWITZ_TARGET=2
```

and assemble parallel results for targets `2,...,9` with
`scripts/assemble_hurwitz_monodromy_eliminant.py`.

## Uniform Arb continuation tubes

The continuation uses the common base target

```text
t_* = 1/2 + 2 i.
```

The loops about `0` and `1` are counterclockwise rational diamonds; the
infinity loop is a clockwise rational diamond enclosing both finite branch
values.  For each straight segment `[a,b]`, set `m=(a+b)/2`.  Arb first
isolates all 23 roots at `m`.  For a root center `c`, a disk of radius `R`,
and a fixed approximate inverse slope `Y`, it evaluates

```text
eta = sup |Y P(t,c)|,
q   = sup |1 - Y P_x(t,x)|
```

for every `t` in the rectangular Arb enclosure of `[a,b]` and every `x` in
a box containing the disk.  A segment is accepted only if

```text
q < 1,                  eta + q*R < R,
```

all 23 disks are pairwise disjoint, and the endpoint root balls lie inside
their unique disks.  The fixed-point map `x -> x-Y*P(t,x)` is therefore a
contraction on each disk uniformly along the whole segment.  If any check
fails, the segment is bisected and retried.  Nearest-neighbor matching is
never used as a proof step.

The six embeddings of the degree-12 field use the exact power-basis
coefficients embedded into Arb balls.  The map on the degree-one component,
which has Nielsen class `6`, is checked independently from the integral optimal equation after the exact target
change

```text
beta = (sqrt(-23)-T)/(sqrt(-23)+T).
```

The completed run contains 150,145 accepted uniform tubes and reaches maximum
dyadic depth 18.  The detailed permutations, precisions, root scalings, and
per-loop tube counts are recorded in
`verification/hurwitz_branch_cycle_summary.json`.

## Nielsen identification

The exact continuation gives cycle shapes

```text
g_0 : 23,
g_1 : 1^7 2^8,
g_infinity : 23,
g_0 g_1 g_infinity = 1
```

for every map.  On the sextic component the Nielsen convention is
`x=g_1`, `y=g_0`.  The independently normalized optimal equation for the
degree-one map
orders the two order-23 endpoints oppositely, so there the convention is
`x=g_1`, `y=g_0^(-1)`.

After conjugating `y` to the fixed 23-cycle, the centralizer translations
match `x` with Nielsen representatives `1,...,7`, exactly once each.  Every
pair generates a permutation group of order `10,200,960`, the fixed natural
copy of `M23`.

Run the quick exact and recorded checks with

```text
make verify-hurwitz-monodromy-eliminant
make verify-hurwitz-branch-cycles
```

To repeat a full continuation, run

```text
make certify-hurwitz-branch-cycles HURWITZ_CLASS=k
make certify-degree-one-branch-cycles
```

The recorded class-`4` run sets `HURWITZ_PRECISION=384`; all other runs use
the default 256-bit precision.

Because the inner Nielsen class has exactly seven elements, the seven exact
maps exhaust it.  Consequently `Spec(K0 x L)` is identified with the reduced
finite inner Hurwitz scheme, and `K0 x L` is its coordinate algebra, without
requiring an a priori LLL height or separation bound.  The LLL reconstruction
remains the discovery route for the displayed coefficients, but is no longer
the logical identification step.
