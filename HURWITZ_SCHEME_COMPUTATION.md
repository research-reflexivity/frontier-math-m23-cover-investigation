# Computing the seven Hurwitz points

## Status

The analytic part of the computation works for all seven inner Nielsen
classes of type `(2A,23A,23B)`.  It produces:

- the four-dimensional weight-two space for each genus-four quotient;
- Taylor expansions at the totally ramified `23A` point in the disk
  coordinate used by Huang--Jackson--Lee--Poonen--Pries--Zhang (HJLPPZ);
- the unique canonical quadric and a Petri cubic in a fixed coefficient
  normalization.

The normalization-invariant ratio `J=q_13/A^2` has now been reconstructed as
an exact rank-seven finite étale algebra over `Q(sqrt(-23))`.  It factors into
a degree-one component carrying Nielsen class `6` and a degree-six field.  Every coefficient
of the normalized canonical quadric and Petri cubic has also been
reconstructed in that algebra.  Exact arithmetic checks verify the field and
factorization, all coefficient plateaus, the coordinate normalization, and
smoothness of the canonical complete intersections on both components.

The two totally ramified marked points and the degree-23 functions have now
also been recovered exactly.  On the sextic component they form one map over
the absolute degree-12 field, whose six embeddings over the chosen embedding
of `Q(sqrt(-23))` are the six missing maps.  Exact saturation proves that the
ratio has degree 23.  The third branch value has an exact CRT/LLL
reconstruction and 24 withheld-prime checks.  The algebraic cover certificate
is now complete as well: an exact characteristic-zero resultant calculation
proves the `2^8 1^7` third fibre.  A second exact elimination removes a
degree-19 target-independent base factor and gives a degree-23 polynomial in
the sheet coordinate.  Arb root continuation around all three branch values
then matches the seven exact maps with the seven Nielsen classes and certifies
`M23` monodromy for every embedding.  Thus the known Nielsen count identifies
the reconstructed algebra unconditionally; an a priori LLL height theorem is
no longer needed for that identification.

The implementation is `scripts/compute_hurwitz_covers.py`.  It uses NumPy and
the explicit seven permutation triples already certified in this repository.

The high-precision path is now implemented in
`scripts/hurwitz_high_precision.py` and `scripts/certify_hurwitz_acb.py`.
It is matrix-free and uses Sage's FLINT/Acb ball fields for the finite
transcendental equations.  Algebra reconstruction is implemented in
`scripts/recognize_hurwitz_algebra.py` and
`scripts/reconstruct_hurwitz_models.py`.  Marked-point and map reconstruction
are implemented in `scripts/reconstruct_hurwitz_marked_points.py` and
`scripts/reconstruct_hurwitz_degree23_maps.py`; exact internal checks are in
`verification/verify_hurwitz_algebra_candidate.py` and
`verification/verify_hurwitz_canonical_models_candidate.py`, together with
the corresponding marked-point, branch-value, and map verifiers.

## Normalization

Use the triangle signature `(2,23,23)`, with `sigma_a=x` and `sigma_b=y`.
The basic triangle in the upper half-plane has

```text
z_a = i,
z_b = mu*i,
mu = 14.482029549242206... .
```

Powers of `delta_b` give all 23 right cosets.  A single Taylor expansion at
`z_a` would be badly conditioned, so the code uses the federalist KMSV
algorithm and expands simultaneously at the 23 points `alpha_i*z_a`.  Each
patch only needs to cover one triangle, whose disk radius is
`0.8708179703675935...`.  The default Cauchy circle has radius
`0.8958179703675935...`; coefficients are balanced as `c_n=b_n*rho^n`.

To compare the curves, the forms are transferred to the unique point above
the `23A` branch.  If `w_b` is the disk coordinate centered at `z_b`, then

```text
q = exp(i*(pi-pi/23))*w_b.
```

This sends the third triangle vertex to the positive real axis and the
order-two vertex to argument `-pi/23`, exactly as in the HJLPPZ disk model.
Echelonization gives `f_i=q^i+O(q^4)` for `i=0,1,2,3`.

The quadric is normalized by the coefficient of `x_0*x_2`; its leading shape
is

```text
x_0*x_2 - x_1^2 + A*(x_0*x_3-x_1*x_2) + ... .
```

The raw coefficient `A` contains the analytic scale of the local parameter
and is not an arithmetic coordinate.  Put `y_i=A^i x_i`.  After multiplying
the quadric by `A^2`, its coefficients transform as
`q_ij -> q_ij A^(2-i-j)`; in particular `J=q_13/A^2` is scale-free.

The Petri cubic is selected modulo the four multiples `x_i Q` by setting the
coefficients of `x_0^2*x_2`, `x_0*x_1*x_2`, `x_0*x_2^2`, and
`x_0*x_2*x_3` to zero.  After the same coordinate change, it is normalized by
`C_222=1`.  These are algebraic coordinate conditions applied identically at
every embedding; unlike an orthogonal-complement convention, they are
Galois-compatible.

## Multi-centre Acb path

For high precision, the 23-chart atlas is replaced by two 23-point orbits of
centres on the long `b--c` edge, with rational Klein parameters `27/500` and
`677/1000`, plus the point over `b`.  On a deterministic order-160 mesh the
largest sampled distance to these 47 centres is `0.45821`; the actual Cauchy
routes used below have target radii below `0.460`.  The separate Arb mesh-cell
certificate proves the all-points bound `0.470869 < 0.471`: Klein mesh cells
are geodesically convex, and the pseudohyperbolic triangle inequality extends
the vertex bounds to every point of both half-triangles.

The resulting Cauchy operator is applied with Horner recurrences and FFTs;
no dense matrix is formed.  For the ball computation, the routewise Horner
recurrences are compiled against FLINT and the Fourier projections use
FLINT's Acb DFT.  (FLINT's product-tree multipoint evaluator was tested but
discarded here because interval widening made the results indeterminate.)
Block fixed-space iteration recovers all four
weight-two forms and avoids the repeated-eigenvalue failure of a scalar
Krylov start.  The code records every triangle-reduction word.  The Acb pass
checks combinatorially that each word lies in the subgroup, then rebuilds
`mu`, the triangle matrices, atlas centres, Möbius maps, and Fourier residual
from the exact signature and rational parameters.

For class `6`, at `N=120`, `Q=320`, two projected-Neumann mixed-precision
corrections put the four 384-bit Acb residual norms between `5.81e-43` and
`1.08e-42`; the canonical-quadric jet residual is `1.39e-40`.  At `N=240`,
`Q=640`, two corrections put the basis residuals below `4.3e-46` and the
quadric residual below `2.0e-44`.  At `N=480`, `Q=1280`, three corrections at
1024 bits put the basis residuals below `8.0e-62` and the quadric residual
below `5.3e-60`; four and five corrections improve the basis residuals to
`1.5e-77` and `3.5e-93`, respectively.  The corresponding quadric residuals
are below `6.3e-76` and `9.3e-92`.

For the six classes on the degree-six component, an independent `N=480`, `Q=1280`,
five-round batch put every basis residual below `8.5e-93`.  A second
`N=360`, `Q=960`, six-round batch at 1024 bits put them below `1.2e-108`.
Its scale-free coefficients agree with the `N=480` batch to roughly 90--92
decimal digits.  The extra `N=360` digits were used for reconstruction; the
`N=480` values were retained as the independent truncation check.

The degree-480 figures above concern output modes only through the polynomial
cutoff.  Evaluating all `Q=1280` output modes exposes a residual of order
`2e-67` in modes `481,...,1279`; this is the honest degree-480 truncation, not
Acb interval noise.  Extending the solved cutoff to `N=700`, while retaining
`Q=1280`, reduces the full all-mode residual below `1e-90` in every class.

The uniform comparison is now recorded in `HURWITZ_TAIL_BOUND.md`.
An Arb mesh-cell proof gives the global radius `0.471`; an a posteriori
IEEE-754/Acb left inverse for modes `0,...,60` has certified Schur margin at
least `6.60e-4`; a Cauchy bootstrap at outer radius `0.99` bounds the exact
branch-normalized forms; and the all-mode residual then certifies at least 75
decimal digits in the first 20 branch-normalized Taylor rows.  The word
`candidate` remains appropriate for algebraic recognition: the tail theorem
certifies the analytic digits, but does not by itself supply an a priori
height/separation theorem identifying every LLL reconstruction with the exact
automorphic invariant.  This no longer blocks the cover theorem: the exact
maps reconstructed from those coefficients independently exhaust the seven
Nielsen classes by the branch-cycle certificate below.

## Calibration

Before running M23, the triangle and Cauchy formulas were tested on KMSV's
degree-five `(5,3,3)` example.  At weight six the operator has exactly three
machine-zero singular values, as required.  This test also caught a
transcription mistake: the formula for `lambda` has denominator
`sin(pi/a)*sin(pi/b)`, without an extra factor of two.  The corrected formula
reproduces KMSV's printed matrix.

For Nielsen class `6`, the class carried by the degree-one component, the
`N=120`, `Q=299` computation
gives

```text
fourth Hejhal singular value = 3.5891e-12
fifth Hejhal singular value  = 6.9170e-3
genus-four gap               = 1.9272e9
quadric jet residual         = 6.95e-9
Petri cubic jet residual     = 1.07e-9
```

At `N=160`, `Q=391`, the quadric and cubic residuals improve to `6.49e-12`
and `8.10e-13`.  The two forbidden leading quadric coefficients are below
`9e-13`, while the coefficients of `x_0*x_2` and `x_1^2` differ from `1` and
`-1` by less than `2e-12`.

## Seven numerical covers

At `N=120`, `Q=299`, the coefficient `A` and quadric residual are:

| ID | `A` | residual |
|---:|:----|---------:|
| 1 | `-2.336837053973 + 0.397912239301 i` | `2.14e-9` |
| 2 | `-0.942589965611 + 0.770679518900 i` | `2.46e-9` |
| 3 | `-0.353948423979 - 0.537014438309 i` | `2.48e-9` |
| 4 | ` 3.500731275443 + 2.571146643331 i` | `7.33e-9` |
| 5 | `-0.259589590908 - 0.604110631773 i` | `3.46e-9` |
| 6 | `-0.631824680770 - 0.044658240056 i` | `6.95e-9` |
| 7 | `-1.258837671244 - 0.715792234714 i` | `4.94e-9` |

Their tentative `A`-polynomial is

```text
X^7
+ ( 2.282896111041 -  1.838162856681 i) X^6
+ (-5.208704266269 - 12.817116112758 i) X^5
+ (-19.331630738856 - 39.691473673913 i) X^4
+ (-17.407879385380 - 62.790502768436 i) X^3
+ ( 2.251072858036 - 54.809820237043 i) X^2
+ (10.922872862830 - 22.306305406471 i) X
+ ( 3.953150636870 -  2.831655733531 i).
```

This is only a floating-point fingerprint of the raw, scale-dependent
coefficient `A`; it is not a coordinate on the Hurwitz scheme.

## Reconstructed coordinate algebra and canonical models

For the normalization-invariant coordinate

```text
J = q_13/A^2,
```

class `6` on the degree-one component reconstructs as

```text
J_6 = (148227 - 34830*sqrt(-23))/142129.
```

The monic polynomial of the seven numerical values of `J` reconstructs over
`K=Q(sqrt(-23))` and factors exactly as this linear factor times an
irreducible sextic.  The full coefficients are recorded in
`data/hurwitz_algebra_candidate.json`.  The absolute degree-12 field of the
sextic component has the small defining polynomial

```text
x^12 - 6*x^11 + 20*x^10 - 32*x^9 + 44*x^8 - 22*x^7
     + 6*x^6 - 22*x^5 + 44*x^4 - 32*x^3 + 20*x^2 - 6*x + 1.
```

This reciprocal polynomial is the compositum of `Q(sqrt(-23))` with the
degree-six trace field defined by

```text
y^6 - 6*y^5 + 14*y^4 - 2*y^3 - 27*y^2 + 44*y - 44.
```

Split-prime factorization types `6`, `1+5`, and `1+1+1+1+2` force the
Galois group of its closure to be `S6`.  Its discriminant field is
`Q(sqrt(11))`, so it is linearly disjoint from `Q(sqrt(-23))`.  Consequently
the relative Galois group on the six maps on the degree-six component is the full `S6`, in
the Nielsen order `(7,4,1,5,3,2)`.  See `HURWITZ_GALOIS_CLOSURE.md`.

At the prime above `23`, the decomposition group in this `S6` is a Klein
four group with orbits `2+4`; the two primes of `L/K0` have
`(e,f)=(1,2),(2,2)`.  Exact reduction of the three local pointed models
gives singularity configurations `E8`, `E8`, and `A2+A6`.  Their
normalizations are rational and all three reduced maps are `t -> t^23`.
The distinguished singular positions have minimal polynomials `u-16`,
`u^2+1`, and `u^2+u+1`, giving an intrinsic residual component idempotent.
See `HURWITZ_LOCAL_23.md` and run `make verify-hurwitz-local-23`.

The JSON filename retains `hurwitz_algebra` as a stable data identifier.  In
the manuscript and this note, `K0 x L` is called the coordinate algebra of
the finite Hurwitz scheme.

All ten quadric coefficients and all twenty cubic coefficients in the fixed
normalization then reconstruct in an integral basis of this field.  Most coordinates are
found by ordinary rational plateaus.  `C_233` requires a higher continued-
fraction bound, while `Q_33` and `C_333` require simultaneous rational
reconstruction of all twelve integral-basis coordinates by LLL:

| coefficient | common denominator | maximum embedding error |
|:--|--:|--:|
| `Q_33` | `79959600090942838037256354600089312590325395259801284` | `1.4e-106` |
| `C_233` | `183891926563909175401913594618178627363643916` | `6.4e-107` |
| `C_333` | `4373159464762185111178624854250694634577724027103703308808` | `7.5e-107` |

For example,

```text
den(Q_33) = 2^2 * 71^4 * 211^4 * 793709999^4,
```

which supplies a useful arithmetic check on the simultaneous reconstruction.
Every reconstructed coefficient agrees with the independent `N=480` model;
the worst discrepancy is `8.9e-90`, at the precision boundary of that batch.

The exact candidate equations are stored in
`data/hurwitz_canonical_models_candidate.json`.  Run

```text
make verify-hurwitz-candidate
make verify-hurwitz-models-candidate
```

to check the degree-seven factorization, the degree-12 coefficient field, all
coefficient plateaus, the fixed normalization, and exact projective smoothness of the
degree-one and degree-six-component quadric--cubic intersections.  These are exact
checks *of the reconstructed candidate*.  The map calculation below removes
the earlier absence of exact degree-23 functions; the automorphic tail issue
remains logically separate from the exact algebraic map and fibre checks.

## Exact degree-23 maps

Write `b=[1:0:0:0]` for one order-23 point and let `c` be the exactly
reconstructed opposite order-23 point.  In the canonical coordinate ring,

```text
dim H^0(5K-23b) = dim H^0(5K-23c) = 4.
```

The degree-five quotient has dimension 27.  Vanishing of the first 23 local
coefficients at either marked point has rank 23, leaving the two
four-dimensional kernels `B` and `C`.  Multiplication is then performed in
the degree-ten quotient, of dimension 57.  Solving

```text
B_i (M C)_j - B_j (M C)_i = 0    for all i<j
```

gives a 16-column system of rank 15.  Its one-dimensional kernel contains an
invertible matrix `M`.  With the deterministic row-zero choice,

```text
u = N/D = B_0/(M C)_0
```

has order 23 at `b` and order `-23` at `c`; moreover `N(c)` and `D(b)` are
nonzero.  Exact projective saturation of the four sections in `B`, away from
`b`, is the unit ideal, and the corresponding saturation of `MC`, away from
`c`, is also the unit ideal.  Thus no residual zero or pole is hidden in the
construction and

```text
div(u) = 23 b - 23 c.
```

This proves that `u` has degree 23 on every geometric fibre of the
reconstructed finite scheme.
The sparse exact quintics are serialized in
`data/hurwitz_degree23_maps_candidate.json`; the sextic component represents
all six missing conjugate maps.

The raw third branch value `lambda=u(a)` is large in the chosen power basis,
so coordinatewise continued fractions are ineffective.  Instead, at each
completely split prime the unique nonzero linear discriminant factor of
multiplicity eight is computed in all twelve residue embeddings and
interpolated back to the power basis.  A 13-dimensional LLL calculation on
71 such primes (a 404-decimal-digit CRT modulus) yields a 1129-bit common-
denominator vector.  The next lattice vectors have at least 1246 bits.  The
selected vector agrees with the independent 1024-bit Arb computation to
`1.7e-101` in integral-basis coordinates and satisfies 24 additional split
primes that were not used by CRT or LLL.  Define

```text
beta = N/(lambda D).
```

Then the branch values are `0,1,infinity`.  The exact value and all
CRT and holdout records are in
`data/hurwitz_degree23_branch_candidate.json`.  Run

```text
make verify-hurwitz-branch-candidate
make verify-hurwitz-maps-candidate
make reconstruct-hurwitz-maps
```

for the quick arithmetic checks and the slower exact canonical-ring
reconstruction.  The generated Magma file
`verification/verify_hurwitz_degree23_branch.m` independently checks the
degree-12 field, primitivity of `lambda`, and all withheld primes.

For the exact characteristic-zero third fibre, project the canonical curve
birationally to a plane sextic.  After cancellation, the numerator and
denominator have plane degree seven.  Eliminating the remaining plane
coordinate gives degree-42 resultants.  At a generic value the resultant and
its derivative have gcd degree 6.  At `lambda` the gcd degree is 14; dividing
by the common degree-6 factor leaves a squarefree degree-8 polynomial.  A
repeated elimination root could still come from two distinct plane points
sharing their projected coordinate.  The Magma certificate rules this out at
the completely split prime `863153`.  In residue embeddings
`2,3,4,5,9,12`, every factor of the residual degree-8 polynomial gives a
degree-1 gcd in the remaining plane coordinate, while all resultant and gcd
degrees are preserved.  Any characteristic-zero common divisor of degree at
least two would retain that degree under such a good reduction.  Consequently
the branch fibre has eight distinct double points and seven simple points:

```text
passport(beta) = (23), (2^8 1^7), (23).
```

The long exact certificate is
`verification/verify_hurwitz_degree23_third_fiber.sage`; its recorded output
is checked by `make verify-all`, while a from-scratch rerun is opt-in:

```text
make verify-hurwitz-third-fiber-exact
```

Independently, `verification/verify_hurwitz_degree23_geometry.m` reduces the
exact data modulo the completely split prime `863153` and checks all twelve
residue embeddings.  Magma obtains constant Hilbert polynomial 15 for the
critical branch fibre in every embedding, accounting for the seven common
base points and eight residual ramification points; six residue embeddings
also give the collision-free plane-projection check above.

## Reproduction

A fast genus-four check is:

```text
DOT_SAGE=/private/tmp/m23-cover-investigation-sage \
  sage -python scripts/compute_hurwitz_covers.py \
  --class-id 6 --terms 30 --samples 92
```

The complete current numerical computation is:

```text
DOT_SAGE=/private/tmp/m23-cover-investigation-sage \
  sage -python scripts/compute_hurwitz_covers.py \
  --all-classes --canonical --terms 120 --samples 299 \
  --output /private/tmp/m23_hurwitz_n120.json
```

For a near-double-precision model of one class, use `--terms 160 --samples
391`.  Beyond this, use Arb/Acb or high-precision Magma with the computed
double-precision four-plane as the initial subspace.

The matrix-free and two-precision Acb checks are:

```text
make hurwitz-multicentre
make hurwitz-acb
```

The mixed-precision canonical-quadric run, using the optional compiled FLINT
bridge when a C compiler and `pkg-config` metadata for FLINT are available,
is:

```text
make hurwitz-acb-model
```

A larger mixed-precision pilot is:

```text
DOT_SAGE=/private/tmp/m23-cover-investigation-sage \
  sage -python scripts/certify_hurwitz_acb.py \
  --class-id 6 --terms 120 --samples 320 \
  --precision-low 192 --precision-high 384 \
  --refine-rounds 1 --neumann-iterations 160
```

## Exact completion and remaining arithmetic questions

The normalization-invariant ratio, rank-seven coordinate algebra, canonical models, marked
points, exact degree-23 map sections, exact third branch value, complete
passport, and all 21 branch-cycle continuations have been checked.  The
permutations match Nielsen IDs `1,...,7` and every pair generates `M23`.
Because the inner Nielsen class has exactly seven elements, the exact maps
exhaust it and identify `Spec(K0 x L)` with the reduced finite inner Hurwitz
scheme.  See `HURWITZ_BRANCH_CYCLES.md`.

What remains is arithmetic interpretation rather than cover identification:

1. construct the group-algebra augmentation intrinsically and in a
   Galois-compatible way from characteristic-23 models, and identify its
   Boolean morphism with the degree-six component without first using the
   exact branch cycles;
2. relate that intrinsic construction to the rational Nielsen point on the
   degree-one component (class `6`) without using its known equation.

The `1+6` connected-component decomposition already resolves the finite
descent test.  The raw values of `nu`, `mu`, the cyclic-conjugacy order counts,
and the group-algebra element `Theta` vary on the degree-six orbit and
therefore do not define the proposed morphisms to constant finite schemes.
Their Boolean singleton predicates all correspond to the primitive
idempotent of the degree-one factor, while the augmentation corresponds to
the complementary degree-six idempotent.  See
`HURWITZ_RELATIVE_TRANSPORTER_INVARIANTS.md`.
