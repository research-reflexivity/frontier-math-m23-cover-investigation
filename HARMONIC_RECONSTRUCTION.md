# Exact local reconstruction: result and audit

Status: 7 September 2026. The manuscript now includes the
assembled harmonic construction in `paper/harmonic-construction.tex`,
with its source and supporting checks included in the public export.
The local uniqueness theorem still states its hypotheses
separately; the new preceding construction derives them for the explicit
harmonic special M23-map. The historical investigation entries below
describe the earlier manuscript state, not additional current gaps.

## Construction now included in the manuscript

The argument is organized around distinct mathematical tasks:

1. The two explicit characteristic-23 tails, their recorded function-field
   Galois groups, and the logarithmic datum construct a special M23-map.
2. The direct special-map symmetry count and centerless descent give a
   locally defined marked lift. The marked-distance lemma determines the
   good-reduction field and its tangent character.
3. Canonical differentials construct an integral lattice. The effective
   degree-eight ramification divisor, located by the valuative criterion,
   determines the primitive logarithmic differential and closed equations.
4. Actual node coordinates identify the same wild generator on both tails.
   The two cyclic Witt designs then exclude the negative first jet.
5. Necessary higher jet equations place the lift in the local theorem's
   neighbourhood. Hensel uniqueness and independent exact existence prove
   incidence to all orders. The global seven-cover test supplies the
   separate uniqueness needed for global fixedness.

The manuscript explains why support alone, abstract tail groups, finite
vanishing digits, or local rationality would not suffice. It also retains
the actual third branch value in the all-orders equations; normalizing
the value at the coordinate origin is equivalent only modulo pi^5.

The open-source construction checks can be reproduced with

```sh
make audit-harmonic-construction HARMONIC_WORKERS=8
```

This target does not certify the written geometric arguments. It checks
the existing tail-Magma input hashes, not new Magma executions. A separate
calculator rerun on 7 September 2026 successfully executed the generic-model
input and both tail scripts, as recorded below. The separate Fano–affine
relative specialization question remains open.

### Checks rerun for the 7 September rewrite

The complete `audit-harmonic-construction HARMONIC_WORKERS=8` target
passed, including the two-point first-jet ideal, the parameter elimination
giving the coefficient 3, both Witt designs, and the branch-value regression.
The complete `audit-harmonic-hensel HARMONIC_WORKERS=8` target passed:
the two residual systems have determinants 16 and 11 and the same root.
The deliberately off-locus residual-multiplier test still detects its
higher-order defect, as the manuscript now explains.

The independent exact model checker was rerun with `--ramification
--descent --frame notes/harmonic_intrinsic_frame_residues.json
--checkpoint notes/harmonic_local_recognition_checkpoint.json`, without
`--output`. It passed smoothness, marked jets through pi^33, exact
base divisors, all six cross-products, the reduced eight-point critical
algebra and full three-point passport, and explicit rational descent.
The existing certificate was not overwritten. These tests do not amount
to an independent expert review of the geometric proof.
The separate map-identification checker also passed again, including
the exact marked isomorphism, the identity of rational functions up to
target scaling, and transfer of the recorded M23 branch-cycle certificate.

## What has been established

The specified harmonic deformation and marked jets, together with the
three-point condition and effective critical cluster, determine one
formal lift. An exact global model constructed from that lift has two
reduced osculating-incidence points. Necessity and uniqueness of the
local equations identify the model with the lift to all orders.

This is an exact reconstruction argument, not extrapolation from a
long string of vanishing 23-adic digits. It is also not a new regular
M23 realization: a separate exact comparison identifies the map with
the distinguished cover constructed by Huang, Jackson, Lee, Poonen,
Pries, and Zhang.

There are three deliberately separate inputs:

1. **Local necessity and uniqueness.** The eleven fixed analytic
   equations have affine reduction with Jacobian determinant 16 in
   F23 and one residual root. The 78 degree-at-most-two identity
   evaluations have been rerun successfully. The off-locus test still
   detects a nonzero higher-order defect in the *constant* residual
   incidence multiplier. That constant matrix must not be advertised
   as an exact syzygy.
2. **Independent exact existence.** The small model and map are
   verified without reading the older characteristic-zero data tables.
   All four smoothness charts, both quintic systems and exact base
   divisors, all six global cross-products, and the eight-point reduced
   critical algebra pass. Exact descent to Q passes. The marked model
   agrees with the local canonical slice modulo pi^33, which places it
   in the stipulated neighbourhood; uniqueness then supplies the
   all-orders conclusion.
3. **Identification and M23 monodromy.** An invertible projective
   change identifies the marked curve with the earlier degree-one
   model. An exact cross-multiplication identity identifies their
   rational functions up to target scaling. The existing certified
   branch cycles therefore apply. This comparison intentionally reads
   the old data. It does not replay Arb continuation, and the passport
   alone would not determine M23.

The local argument needs necessity and uniqueness plus exact existence.
It does **not** depend on the separate proposed Picard-functor argument
that the selected eleven equations generate every map equation.

## Reproduction

From the repository root, with SageMath available:

```sh
make verify-harmonic-reconstruction
make audit-harmonic-hensel HARMONIC_WORKERS=8
```

Both component commands of the first target and the full second target
were run successfully on 6 September 2026. The exact model run took
about 102 seconds in this audit; timings vary by machine.
These targets are opt-in, not added to `verify-all`.

The audit target also checks every first-order equation/mark direction,
the universal canonical-ring Groebner basis, and the independent
finite-algebra implementation of the residual equations. Those checks
all passed, including its separate 78-point audit and comparison on
the zero locus. Its Jacobian determinant is 11 rather than 16 because
it uses different coefficient equations; the unique root is the same.

`make verify-hurwitz-osculating` was also rerun successfully: the
degree-one component has common divisor degree two, the sextic
component has degree zero, and both results survive an invertible
projective change of coordinates. The residual divisor relation and
the S4 monodromy test for the degree-four pencil passed as well.

Principal outputs:

- `notes/harmonic_independent_exact_model_certificate.json`: exact
  equations, quintic section witnesses, map, critical value and
  polynomial, and explicit rational descent. Its scope string does
  not claim M23 monodromy by itself.
- `verification/harmonic_model_identification_summary.json`: exact
  source transformation, target scalar, input hashes, and successful
  transfer of the existing M23 certificate.

The original derivation, full residual Jacobian and normalizations are
in [the fixed-system note](notes/HARMONIC_FIXED_HENSEL_SYSTEM_AND_INCIDENCE_IDEAL_2026_09_06.md).
Its historical cautions about exact incidence are superseded by
[the exact reconstruction note](notes/HARMONIC_EXACT_INCIDENCE_FROM_LOCAL_RECONSTRUCTION_2026_09_06.md)
and the manuscript's conditional local theorem. They are retained as
a record of what finite precision did and did not establish.

## Independent Magma check: successfully executed on 7 September 2026

```sh
make prepare-harmonic-magma
make verify-harmonic-magma
```

The first command mechanically emits a standalone Magma input file
from the exact section witnesses; it is not a verification. The second
actually requires a working Magma installation. The generated input
is about 40 KB and checks the witnesses anew: smoothness, marked jets,
basepoint saturation, an independent multiplier kernel, all six
cross-products, and a radical eight-point critical algebra with one
critical value. It does not audit the analytic Hensel argument.

**Execution status: PASS, Magma 2.29-10.** The restored
[Sydney calculator](https://magma.maths.usyd.edu.au/calc/) executed the
unchanged 39,831-byte input successfully on 7 September 2026 in 4.009 seconds.
It checked smoothness, the two-point incidence, marked jets and base divisors,
the independent multiplier kernel and cross-products, the degree-23 function,
and the eight distinct equal critical values giving the full passport.

The [complete rerun report](MAGMA_RERUN_2026_09_07.md) records successful runs
for all 16 public inputs, including both tail function-field Galois groups.
Exact requests and raw responses are retained under
`verification/magma_runs/2026-09-07/`. The first conductor-seven tail attempt
timed out; the unchanged assertions passed after resetting to the previously
recorded random seed. The older records are preserved, not relabeled as new runs.

`make verify-magma-calculator-records` checks the saved evidence locally.
This independent software check is no longer outstanding. It does not
verify the analytic Hensel argument, the written geometric construction,
explicit Q-descent or the pi^33 comparison; those have their separate
SageMath checks and mathematical arguments. No new Arb continuation was run.

## Historical record: what the 6 September manuscript audit left open

This section records the hypotheses left open at the manuscript-audit
stage. The subsequent private proof assembly linked below now gives
a proposed derivation for the constructed harmonic example. It has
not been promoted to the manuscript.

The missing branch-cycle-only implication is now precise. Starting
with an abstract special M23-map and its compatible patching data,
one must construct the degree-23 quotient and its canonical model
so as to obtain **the particular neighbourhood in the theorem**.
That includes:

- the contraction to `q0=HW-U^2-UV-V^2`, `c0=UW^2-V^3` and the
  marked coordinates `z=1,-1`;
- the coefficients C1, C2, C3 and the marked jets, not only the
  harmonic cross-ratio on the reduced normalization;
- the primitive divided logarithmic differential and its effective
  eight-point critical cluster;
- a common labeling of the sheets on overlapping charts, if the
  proposed two-Witt-design argument is used to select C1.

I checked the relevant definitions and lifting construction in
[Wewers, Sections 4.1–4.2, especially Theorem 4.5](https://arxiv.org/html/math/0205026).
They retain identified special fibres and compatible patching data.
They do not identify the resulting quotient with our chosen canonical
slice. Counting lifts or identifying their local fields is not a
substitute for that comparison. The new manuscript theorem does not
invoke the unaudited special-map counts as a way to omit its hypotheses.

The separate automorphism-counting notes also flag Selander's correction
to Wewers's general automorphism description. That full correction
was not available in this audit; we have not represented it as checked.
It is not used in the local uniqueness/exact-existence proof above.

Once the canonical comparison is supplied, the local theorem gives
incidence. The existing global 2-versus-0 calculation then gives its
unique support, and functoriality under conjugation and exchange of
marks gives global fixedness. This last reasoning is unchanged and
does not infer global rationality merely from a local linear factor.

## Historical private investigation record (6 September)

### Further internal scrutiny: valuations and scalar normalization

[The adversarial audit](notes/HARMONIC_ASSEMBLY_ADVERSARIAL_AUDIT_2026_09_06.md)
records a further check of the assembled geometric argument. No fatal
contradiction was found, but this is not independent expert review.
The assembly now explicitly justifies parameter independence of the
marked distance and centers the ramification points on the canonical
model by the valuative criterion, without assuming a contraction map.

The audit also proves that the computational normalization b(0)=1
differs from normalizing the third branch value only from order pi^5
onwards. The active all-orders equations already eliminate that branch
value correctly. A new finite-residue regression test reran the global
function equations modulo pi^9 and found a nonzero discrepancy
u-1=-34*pi^5 modulo pi^9. It also checked stability under perturbing
the two unavailable ramification-factor digits. This makes explicit
why using b-1 is valid in the earlier modulo-pi^5 calculations and
invalid as an all-orders shortcut.

Only private audit notes and a diagnostic checker changed in this
continuation; the manuscript, PDFs, public materials and Magma execution
status are unchanged.

### Subsequent private investigation: assembled geometric construction

[From the harmonic special map to the marked local neighbourhood](notes/HARMONIC_CONSTRUCTION_TO_INCIDENCE_ASSEMBLY_2026_09_06.md)
assembles the construction, local descent, canonical model and jet
selection as one computer-assisted proposition. Its new marked-distance
lemma derives the degree-15 good-reduction extension and the tangent
character from the separation of the two rational marked points.
Consequently that inertia character, the canonical equations, the
effective ramification cluster and the common sheet labels are
conclusions of the written construction, not extra premises within it.

The scope is an existence construction for this particular harmonic
special map. The explicit characteristic-23 tail equations, their
recorded function-field M23 computations, the special-map lifting
theorem and the finite jet/Hensel calculations remain visible inputs.
This is neither a permutations-only algorithm nor a claim that every
tuple of the same passport has harmonic reduction. Global fixedness
still uses the independently checked unique incidence support, not
merely the local field of definition.

Freshly rerun in this continuation: the direct special-map algebra
checker and GAP normalizer checks; the canonical-descent coefficient
and normalization checks; the full first-jet ideal and next two jet
steps; the two cyclic Witt designs and their tail alignment; the full
`audit-harmonic-hensel` target; and the independent exact model checker
with ramification, rational descent and frame comparison through pi^33.
All passed. The exact model checker was run without `--output`, so
it did not overwrite the existing certificate. At the time of this historical
continuation, only stored tail-Magma input hashes were checked and the new
generic-map input had not been run. The successful 7 September calculator
executions are recorded in the current-status section above.

The assembled geometric proof remains a private candidate for independent
review. Its hypotheses are now explicit mathematical/computational
inputs rather than an unexplained assumption of the desired local
model. The current manuscript theorem remains conditional. No manuscript,
PDF, public repository, commit or deployment was changed in this
continuation.

### Subsequent private investigation: the node coefficient comparison

[Leading coefficients at a node and the harmonic sheet comparison](notes/NODE_LEADING_COEFFICIENTS_AND_HARMONIC_SHEET_TRANSPORT_2026_09_06.md)
proves an all-orders coefficient lemma in `S[[a,y]]/(ay-alpha)` and
applies it to the actual regular function `sigma(a)-a`. It also gives
the canonical node-coordinate comparison, including the absence of
horizontal zeros and poles and the exact choice of smoothing roots.
The resulting common sheet labels are i on the conductor-7 tail and
3i on the conductor-15 tail. Thus, once the preceding marked canonical
construction has been established, this comparison is an implication,
not a further independent annular-compatibility hypothesis.

The independent checker `notes/certify_node_leading_coefficients.py`
passed its 117,678 monomial tests, universal valuation identity,
compatible-root identities, formal-slice weight check and exact cyclic
order-23 test model. The enumeration and cyclic example do not certify
the M23 application: the all-orders lemma is proved in the note, and
the upstream special-map construction, descent and group certificates
remain separate inputs to audit together. The paper's hypotheses and
both PDFs remain unchanged in this continuation.

### Subsequent private investigation: tame canonical descent

The new working note
[Tame canonical descent and harmonic jets](notes/TAME_CANONICAL_DESCENT_AND_HARMONIC_JETS_2026_09_06.md)
constructs a rational (3,5)-cuspidal canonical model from a smooth
genus-four tail with the specified tame order-15 action. It includes
a proposed marked-annulus comparison for the constructed harmonic
special G-map and connects it to the existing first- and higher-jet
equations. This reduces the former unspecified-model issue to explicit
descent, node-coordinate, and horizontal-divisor arguments.

The new algebraic checker and fresh first-jet, next-order, and Witt-design
calculations passed. The full geometric application remains a private
proof candidate for a further audit; the new checker does not certify
that application. The manuscript's local hypotheses have not been
removed, and no PDF, public repository, commit, or deployment was
changed during this subsequent investigation.

### Earlier manuscript revision

The private manuscript, its PDF, README, verification targets and
public-export allowlist have been updated. The public working tree,
remote repository and visual app have not been changed in this audit.
The allowlist now contains the files needed for a later deliberate
sync; it does not publish them automatically.
The 46-page PDF builds with resolved references and no overfull boxes;
the revised pages were rendered and visually checked. Its output copy
is byte-identical. No commit has been made during this audit.
