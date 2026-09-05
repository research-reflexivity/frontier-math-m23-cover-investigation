# Correctness repairs and a proved local replacement

Date: 5 September 2026.

This revision concerns the private cover-investigation manuscript. It does
not assess arXiv's moderation decision, claim a proof of the unresolved
relative comparison, or update the public repository or companion apps.

## Mathematical status

The exact algebra K0 x L, its relative S6 closure, the seven reconstructed
maps, and the computational comparison of the component idempotents remain.
The claimed independent relative quadratic-determinant proof has been
withdrawn from the manuscript. The older normalization/gluing and
orientation notes are retained as explicitly marked research history.

### The proved replacement

Theorem 3.13 proves a local component criterion for any finite normal flat
algebra A over the ring of integers R of a nonarchimedean local field.
Write the ramification and residue degrees of its factors as (e_i,f_i).

- The scheme-theoretic fixed locus of q-power Frobenius on the reduced
  closed fibre consists of precisely the factors with f_i=1. It lifts
  uniquely to an open-and-closed subscheme of rank sum(e_i : f_i=1).
- The unramified extension of degree m splits after pullback to a factor
  exactly when m divides f_i.

The proof uses finite-field tensor products and henselian lifting of
idempotents and finite etale algebras. These general facts are standard,
not a claimed new theory. For residue degrees restricted to 1 and 2, the
quadratic split locus is the complement of the Frobenius-fixed locus.

In the M23 case, the pairs are (1,1),(1,2),(2,2). Therefore the complement
of the unique Frobenius-fixed component has generic rank six. The
unramified quadratic extension is R[z]/(z^2+1); its pullback to the actual
normalized integral Hurwitz algebra is split precisely on this complement.
Corollary 3.14 identifies its idempotent with the singular-position
idempotent through the explicit identity

    -(u^23-u)^22 = 2*u^4+2*u^3+4*u^2+2*u+3 mod R23(u).

This is a geometric interpretation on a finite integral moduli scheme,
not a construction of relative incidence correspondences on a family of
curves. Its agreement with epsilon(Theta) still uses the exact branch-cycle
matching. It is not the independent second proof previously claimed.

The proof distinguishes the fixed-locus quotient from the invariant
subalgebra. It also records dependence on the fixed residue field and
gives X^3-2 at 5 as a counterexample to automatic local-to-global descent.
The global decomposition K0 x L is an essential input to the application.

### Repairs to existing arguments

1. Proposition 3.4: the selected quintic canonical sections have a common
   effective divisor E of degree seven, not no common zero. The complete
   four-section systems have fixed divisors 23b and 23c; their
   cross-multiplication identities then prove div(N/D)=23b-23c. The
   equations and the degree-23 maps do not change.
2. The degree-23 cover is distinguished from its M23 Galois closure.
   The two wild branch points are conjugate over Q and individually
   rational over Q(sqrt(-23)).
3. The ramification theorem now lists every constituent prime's (e,f).
   It proves persistence of the quadratic ramification at 23 by the odd
   ramification indices of L0/E8 there. It does not infer persistence
   at 2 or 3 from the octic relative discriminant alone, or infer an
   inertia order from a census restricted to that order.
4. The abstract, introduction, discussion, README, and scope statements
   no longer claim a relative determinant construction or a second
   specialization proof. Finite incidence, the actual finite pinching
   stalk, and elementary quadratic identities are retained with their scope.
5. The certificate section separates recomputation, stored-result checks,
   interval continuation, and noncomputational proof obligations.

### What remains missing

A proof independent of branch-cycle matching must construct a common
relative model and cohomological correspondences extending the finite
incidences, while retaining the specified integral multiplicities
mathcal C and T_D. It must establish the generic and special descriptions
of its invariant by actual maps and filtrations. Ordinary duality,
the finite weight identities, and a telescope of integer differences do
not supply those data.

The old half-rank determinant assignment also lacks the stated monoidal
coherence: swapping two rank-two blocks has ordinary determinant +1,
while the proposed odd graded lines have Koszul swap -1. This obstructs
that assignment with the usual determinant tensor identifications, not
every conceivable quadratic refinement in a different category.

## Verification and reproduction

The new target make verify-hurwitz-frobenius-selector recomputes the local
Hurwitz decomposition and pointed canonical reductions, then checks the
Frobenius equalizer, the polynomial idempotent, the quadratic splitting
in each residue field, and the local-to-global counterexample.

The general theorem is proved in the manuscript; the finite-field script
does not certify all possible local algebras or the missing geometry.
The existing branch-cycle summary target checks stored certified outputs,
not a fresh run of every Arb interval inequality.

The full degree-23 producer can be rerun without overwriting recorded data:

    DOT_SAGE=/private/tmp/m23-cover-investigation-sage sage -python \
      scripts/reconstruct_hurwitz_degree23_maps.py \
      --component all --check-basepoints --output /private/tmp/m23-maps-review.json

No computational equations or certificate input hashes were changed by
the mathematical repairs. The new checker is additional. Public export
lists include it and this note. The publication update exports these repairs
to the public repository and replaces the visual map's withdrawn second-proof
claims with the proved local criterion and the explicit open relative step.

The paper is retitled *An M23 Hurwitz scheme: exact arithmetic and reduction
at 23*. Its abstract and introduction now lead with the exact reconstruction
and the established component comparison. The local interpretation is
separated from the unresolved relative construction; the supporting
projection, branch-fibre arithmetic and specialization results remain in
the main text.

### Reruns completed in this revision

- Local decomposition, all three pointed reductions, and the new
  Frobenius/quadratic-splitting checks: PASS.
- Full degree-23 reconstruction on both components, with residual
  base-locus saturations: PASS. The regenerated JSON agrees with the
  stored equations, with only the additional successful saturation fields.
- S6 closure over K0 in Sage and PARI: PASS.
- Stored branch-cycle consistency, Nielsen matching, and generated M23
  groups: PASS. This was not a fresh Arb continuation run.
- GAP transporter computation and arithmetic descent verifier: PASS.
- PARI constituent-field ramification and octic relative discriminant:
  PASS. The latter's output refers to E'_8/E8; the proof explains
  separately what persists after base change to L0.
- PDF rebuilt; no LaTeX warnings, undefined references, or overfull boxes
  in the final log. The changed mathematical and certificate pages were
  rendered and inspected.

No full verify-all run, fresh Sydney Magma run, or proof-assistant
verification is claimed.
