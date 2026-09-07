# Arithmetic and geometry of an M23 Hurwitz scheme

This repository accompanies the paper *Arithmetic and geometry of an
M23 Hurwitz scheme*.
Huang, Jackson, Lee, Poonen, Pries, and Zhang constructed the regular
M23 cover and identified its Galois-fixed Nielsen class. This paper
studies the full seven-cover family and gives a geometric
characterization of the distinguished member.

At the two totally ramified points of each canonical genus-four curve,
take the osculating planes. Their intersection line meets the
distinguished curve in two reduced points and misses each of the
other six curves. The condition uses the curve and the unordered pair
of marked points, so it is preserved by global Galois conjugation.
Its uniqueness therefore gives a model-based proof of fixedness.

The central results are:

- the exact finite Hurwitz algebra K0 x L, where
  K0 = Q(sqrt(-23)) and [L:K0] = 6, with certified maps representing
  all seven Nielsen classes and relative Galois closure S6;
- a relative osculating-incidence scheme supported over the distinguished
  point, checked exactly over both coefficient fields;
- the equivalent relation [C-D] = 3[A-B], with a unique additional
  pair of points on that curve and no such pair on the other six;
- an intrinsic degree-four osculating pencil with geometric monodromy
  S4, distinct from the earlier pencil |K-A-B|;
- the pointed reductions at 23 and the exact comparison of their
  component idempotent with the global incidence and finite-group
  augmentation;
- a harmonic special M23-map built from two explicit characteristic-23
  tail covers, with the marked canonical model and required jets of
  its lift derived from tame descent, effective ramification, and a
  common identification of the sheets;
- independent exact reconstruction and Hensel uniqueness proving
  incidence on that lift, with its identification as the known
  cover checked by an exact isomorphism of maps.

See [the osculating calculation](HURWITZ_OSCULATING.md) for the
short geometric argument, exact inputs, and reproduction commands.
The local construction supplies existence without reading the seven
generic models. The exact seven-model calculation supplies global
uniqueness of the incidence locus. These are different steps: a unique
local lift alone does not establish global rationality.
For the constructed special map, the local model and jets are now
derived, not additional premises. The characteristic-23 tail groups,
finite jet equations, and exact model remain computational inputs.
This is not an equation-free construction from permutations alone,
nor a claim that harmonic position determines every such deformation.

See [the local reconstruction and audit](HARMONIC_RECONSTRUCTION.md).
`make verify-harmonic-reconstruction` reruns the exact model, map,
descent, local comparison and M23 identification. The longer
`make audit-harmonic-hensel` audits all 78 residual identity tests.
`make audit-harmonic-construction` runs the additional special-map,
canonical-lattice, jet, and common-sheet design checks. It does not
mechanically verify the written geometric arguments or rerun Magma.
The independent Magma check is prepared but **not yet run**; these
opt-in targets are separate from `verify-all`.

The local component criterion uses Frobenius and unramified extensions.
On the normalized integral Hurwitz model, the unramified quadratic
extension splits exactly on the degree-six component. Its complement
is the Frobenius-fixed component of the reduced closed fibre.
Comparing with the finite-group augmentation still uses certified
branch-cycle matching.

The relative osculating incidence is a characteristic-zero construction.
It is **not** the earlier proposed Fano--affine specialization or quadratic
orientation construction. That separate comparison remains unproved;
see [the repair record](notes/CORRECTNESS_REPAIRS_2026_09_05.md).
The regular realization and its original descent remain due to
Huang et al.; the present incidence supplies an additional,
computationally established route to the known fixedness.

Two supporting calculations describe the distinguished cover itself: an
explicit minimal-degree `(23,4)` equation and the Fano-plane and affine-cube
arithmetic in the fibre above the rational branch point `T=0`.  These
structures give finite incidence identities, but the branch fibre alone
does not distinguish the fixed Hurwitz point.

A uniform family of `M23` number-field specializations is a separate
application of the minimal equation.  Huang et al. already obtain abstractly
an infinite mutually independent family; the result here gives an explicit
all-members progression and controlled ramification within it.

Here “rational branch point” means that `T=0` is rational on the base; the
points above it need not be rational.  This is a characteristic-zero branch
fibre, not the special fibre of a model over a local base.

Huang et al. prove that a degree-four function exists and is minimal, but do
not compute it.  This paper makes that function and its equation explicit,
derives the branch-fibre description, and compares the two plane coordinates
under the change of generator `W = J0(T,V)/J1(T,V)`.  The fibre of the
minimal-degree equation at `T=0` factors as
an irreducible septic times the square of an irreducible octic.  Evaluating
`W` on the septic and octic factors of the original fibre generates the two
corresponding factors of the optimal fibre.  Thus the Fano point field and
the affine eight-point field are preserved by the change of generator.

Huang et al. use the nonsquare discriminant of the canonical quadric to prove
minimality.  This repository refines their argument by identifying the exact
ruling field: the two rulings are governed by
`Q(sqrt(4873)) = Q(sqrt(11*443))`, not by the branch-orientation field
`Q(sqrt(-23))`.  Consequently the quadratic étale scheme of rulings does not
by itself explain the Galois-fixed Nielsen point.

Building on Huang et al.'s numerical computation of the seven complex
covers, the reconstruction here computes normalized canonical quadrics and
Petri cubics for all seven Hurwitz points, including the six points on the
connected degree-six component.  Exact reconstruction and branch-cycle
certificates identify the resulting covers.  Run `make hurwitz-pilot` for the quick genus-four calibration or
`make hurwitz-numerical` for all seven low-precision models.  A
47-chart matrix-free atlas lowers the observed reduced disk radius from
`0.87082` to less than `0.459`; `make hurwitz-multicentre` repeats class `6`,
the class singled out by the displayed permutation in the paper, in that
atlas.  `make hurwitz-acb` rebuilds its finite
geometry from exact signature and rational-atlas data and checks the residual
in two Acb precisions.  `make hurwitz-acb-model` builds FLINT's Acb FFT bridge
and repeats the canonical quadric on the degree-one component with two mixed-precision
refinement rounds.

The normalization-invariant ratio `J=q_13/A^2` reconstructs the rank-seven
coordinate algebra of the finite inner Hurwitz scheme over
`Q(sqrt(-23))`.  It is the product of a degree-one component and a connected
degree-six component.  The degree-one component carries Nielsen class `6`,
called the distinguished class in the paper because it contains the displayed
permutation.  In this coordinate algebra, all 30
coefficients of the normalized canonical quadric and Petri cubic in the
fixed coefficient normalization, and the two
totally ramified marked points, have stable exact reconstructions.  Exact
canonical-ring linear algebra then recovers the degree-23 ratio: the spaces
`H^0(5K-23b)` and `H^0(5K-23c)` both have dimension four, and their unique
rank-15 multiplier gives sections `N,D` with divisor ratio `23b-23c`.  The
absence of residual base points is checked by exact saturation.

The remaining branch value on the sextic component is reconstructed exactly
from 71 completely split primes and a 13-dimensional LLL calculation, then
passes 24 split primes that were withheld from the reconstruction.  Thus the
six conjugate maps and their normalization `beta=N/(lambda D)` are explicit
in `data/hurwitz_degree23_maps_candidate.json`.  An exact characteristic-zero
plane-sextic calculation proves the third-fibre pattern `2^8 1^7`: the generic
and branch resultants have gcd degrees 6 and 14, and the additional degree-8
factor is squarefree.  An independent Magma calculation obtains critical-fibre
length 15 (seven common base points plus eight ramification points) in all
twelve residue embeddings at a completely split prime, with six
degree-preserving embeddings also excluding projected-point collisions.  The
uniform automorphic comparison is also certified: an Arb mesh-cell proof
covers both half-triangles by the 47 charts with radius below `0.471`, a
finite-`Q` Schur certificate controls all omitted modes, and `N=700`,
`Q=1280` all-mode Acb residuals certify at least 75 decimal digits in the
first 20 branch-normalized Taylor rows.  The word `candidate` now refers to
the reconstruction history of the displayed coefficients, not to an open
cover-identification problem.  An exact plane eliminant has degree 23 after
removing its degree-19 base-point factor.  Arb continuation around `0`, `1`,
and `infinity` certifies 150,145 uniform interval-Newton tubes.  The resulting
triples match Nielsen classes 1 through 7 exactly and every pair generates
`M23`.  Since the inner Nielsen class has seven elements, this identifies
`Spec(K0 x L)` with the finite inner Hurwitz scheme without an a priori LLL
height bound.  See `HURWITZ_SCHEME_COMPUTATION.md`,
`HURWITZ_TAIL_BOUND.md`, and `HURWITZ_BRANCH_CYCLES.md` for the normalization
and certificates.  The `1+6` connected-component decomposition also settles
the descent of the finite relative-transporter invariants.  Their raw values
do not define morphisms from the degree-six component to constant finite
schemes, but four Boolean predicates all correspond to the idempotent of the
degree-one component; the group-algebra augmentation corresponds to the
complementary idempotent.  See
`HURWITZ_RELATIVE_TRANSPORTER_INVARIANTS.md`.

The sextic component has full symmetric Galois closure.  More precisely, its
absolute degree-12 field is the compositum of `Q(sqrt(-23))` with the trace
sextic

```text
y^6 - 6*y^5 + 14*y^4 - 2*y^3 - 27*y^2 + 44*y - 44,
```

whose Galois group is `S6`.  The discriminant field of that closure is
`Q(sqrt(11))`, so the relative Galois group over `Q(sqrt(-23))` is also
`S6`, acting naturally on the six maps on the degree-six component in Nielsen order
`(7,4,1,5,3,2)`.  See `HURWITZ_GALOIS_CLOSURE.md`.

At the unique prime of `Q(sqrt(-23))` above `23`, the local decomposition
group on the sextic component is `V4`, with orbit sizes `2+4` and relative
prime data `(e,f)=(1,2),(2,2)`.  Together with the point on the degree-one
component, the
normalized local integral model of the Hurwitz scheme has residue degrees
`1+2+2`.  Reduction of the exact pointed maps gives an `E8` source for the degree-one
point and the unramified-degree-two point, and an `A2+A6` source for the
ramified-degree-four point.  Each normalization is `P1`, and every reduced
pointed map is the Frobenius map `t -> t^23`.  In the stored-ratio normalization, the distinguished singular
positions have resolvent `(u-16)*(u^2+1)*(u^2+u+1)` and give the
component idempotent. The degree-one stored map has third branch
value congruent to 7. Scaling that value to one changes its singular
coordinate from 16 to `16/7 = -1`, while the other two coordinates
stay unchanged. The uniform branch-normalized resolvent is therefore
`(t+1)*(t^2+1)*(t^2+t+1)`: the distinguished point is the harmonic
position `t=-1` relative to `0, infinity, 1`.  Exact branch-cycle
identification shows that this idempotent and the relative-transporter
augmentation both equal the sextic-component idempotent.  The accompanying
older connector notes record a proposed relative orientation construction,
now explicitly marked unproved.  The revised manuscript instead proves the
Frobenius and unramified-splitting interpretation, while retaining the
finite incidence identities and the obstruction from the vanishing ordinary
coefficient Bockstein.  See [the local calculations](HURWITZ_LOCAL_23.md)
and [the repair record](notes/CORRECTNESS_REPAIRS_2026_09_05.md).

## Reproduce

Requirements for the complete open-source certificate suite are SageMath,
Singular, GAP, PARI/GP, Python 3, and a LaTeX installation with `latexmk`.
Magma provides optional independent repetitions of the main identity,
irreducibility, canonical-quadric, and unpointed `23`-inertia obstruction
checks, as well as the twelve-embedding Hurwitz critical-fibre calculation.
Transcript metadata from the independent public Magma runs is recorded in
`verification/magma_verification_summary.json` and
`verification/gauss_prolongation_magma_summary.json`, and
`verification/canonical_quadric_magma_summary.json`.  The Hurwitz branch and
geometry runs are recorded in the corresponding
`verification/hurwitz_degree23_*_magma_summary.json` files; the Galois and
local runs are recorded in `verification/hurwitz_galois_closure_magma_summary.json`
and `verification/hurwitz_local_23_magma_summary.json`.  Five independent
connector runs are recorded in the corresponding Fano--affine, pinched-tag,
wild-orientation, pointed-Bockstein, and logarithmic-quadratic-line Magma
summary files in `verification/`.

```text
make verify-all
make paper
# focused new exact checks (also included in verify-all)
make verify-hurwitz-osculating
make verify-hurwitz-degree-one-normalization
# optional long exact third-fibre recomputation
make verify-hurwitz-third-fiber-exact
# optional licensed rerun
make verify-magma
```

`make verify-all` also checks the exact Hurwitz map/passport records and that
the recorded Magma runs still match the current deterministically generated
certificates by SHA-256 and byte count.  It also reruns the Arb all-points
atlas cover and checks the compact seven-class tail record.  The expensive
per-class recomputations are available as
`make hurwitz-tail-stability HURWITZ_CLASS=k` and
`make hurwitz-tail-model HURWITZ_CLASS=k`.  The recorded exact branch cycles
are checked by `make verify-hurwitz-branch-cycles`; a full continuation rerun
is available as `make certify-hurwitz-branch-cycles HURWITZ_CLASS=k`.
The recorded class-`4` run uses
`HURWITZ_CLASS=4 HURWITZ_PRECISION=384`; the other classes use the default
256-bit precision.  Run `make certify-degree-one-branch-cycles` for class `6`.
The exact descent test for the relative-transporter invariants is
`make verify-hurwitz-relative-transporter`.
The finite group identities, local arithmetic, divided normalization
telescope, and quadratic-orientation identities relevant to the
characteristic-23 connector are checked by
`make verify-hurwitz-connector`.
The focused exact reconstruction of the formal annulus at the ramified `A6`
node is available separately as `make verify-hurwitz-connector-a6`.
The Galois-closure calculation is checked independently in SageMath and
PARI/GP by `make verify-hurwitz-galois-closure`; its recorded Magma run is
hash-bound into `make verify-all`.

The manuscript is `paper/main.tex`; a successful build produces the final
artifact `output/pdf/m23-cover-investigation.pdf`.

## Private and public histories

Development uses a private canonical repository and a separate public
repository with an independent Git history.  Public releases are generated
from the private working tree using only the paths in `PUBLIC_FILES.txt`:

```text
make export-public EXPORT_DIR=/absolute/path/to/public-repository
```

See `PUBLISHING.md` for the rationale and release checklist.
