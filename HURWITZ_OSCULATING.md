# Osculating planes on the seven exact covers

This is the global geometric characterization proved in Section 2 of
the paper. The exact computation on the seven already certified models
establishes existence and global uniqueness. Section 5 additionally
constructs a lift with this incidence from an explicit characteristic-23
special map, without reading the seven generic models. A construction
intrinsic to the branch datum alone, and a proof of global uniqueness
without the seven-model exclusion, remain open.

## Statement and reason for Galois invariance

Let A and B be the two totally ramified points of the degree-23 map.
The canonical vanishing sequence at each is (0,1,2,3). Let h_A and h_B
be their order-three osculating sections, and put

    E_osc = gcd(div(h_A), div(h_B)).

On the distinguished cover E_osc consists of two reduced points.
On each of the other six curves it is empty.
The construction uses canonical jets and the unordered pair {A,B};
it commutes with isomorphisms and global Galois conjugation.
Its one-point support in the Hurwitz scheme is therefore Galois-fixed.
It does not define that support by assuming rationality.

The relative incidence is the common zero scheme of the two
canonical jet-kernel sections on the universal quotient curve.
It is finite over the finite Hurwitz scheme, with geometric fibre
lengths 2,0,0,0,0,0,0.
This characteristic-zero construction is not the separate proposed
relative Fano--affine specialization.

## Equivalent Jacobian criterion

The class xi=[A-B] has exact order 23. For the distinguished curve
there is a unique pair of additional points C,D satisfying

    [C-D] = 3 xi.

There is no such pair over the algebraic closure on any of the other
six curves. Indeed, removing E_osc from the two osculating divisors
gives 3A+D and 3B+C. Conversely, such a degree-four linear equivalence
has at least two sections; Riemann--Roch supplies a complementary
effective divisor of degree two in the canonical system. Uniqueness
of the osculating sections forces that divisor to be E_osc.

The resulting degree-four pencil has geometric monodromy S4.
It is different from the earlier pencil |K-A-B| used to construct
the minimal bidegree-(23,4) equation.

## Reproduce

From the repository root:

    make verify-hurwitz-osculating
    make verify-hurwitz-degree-one-normalization

Both targets are included in make verify-all.
The first executes the maintained certificate
verification/verify_hurwitz_osculating.py twice:

- exact local canonical jets, homogeneous line gcds on both coefficient
  fields, residual-point checks, and the quartic monodromy calculation;
- an independent repetition of the incidence test after an invertible
  projective change of coordinates.

The old notes/certify_global_osculating_separator.py path is a
compatibility wrapper, not a second implementation.

The test does not use numerical approximations or reduction at 23.
Over the degree-one coefficient field the restricted quadric divides
the cubic and has nonzero discriminant. Over the sextic field the two
restrictions are coprime, excluding all six conjugates at once.
The quartic function-field model is checked by a linear subresultant
and exact substitution. Its discriminant's squarefree decomposition
has degree/multiplicity pairs (10,1) and (5,2). Together with the
ramification index three and degree four, this proves geometric S4.

## Input fingerprints

The historic filenames retain the word candidate. Their separate
exact cover and branch-cycle certificates establish their identification
with all seven covers; the new incidence verifier does not replace them.

| Input under data/ | SHA-256 |
| --- | --- |
| hurwitz_algebra_candidate.json | 10052c4ccace3cbeae0af7a307f01b6fc38c0b07404e99dd457c698fba07538a |
| hurwitz_canonical_models_candidate.json | f1594e9e85f3b19927938be1fa391afdb8f4ccc5a3db3ea41a2b2ab426a8e2c4 |
| hurwitz_marked_points_candidate.json | 129fb6c72262ddd14b18dd034d94f662cdf36de7e0eb432a725dcfc1283f0185 |
| hurwitz_degree23_maps_candidate.json | 86cdbd314d01a9123aed0b6890fa4316b4303782b1d69563e020b6ede62a6c2a |

The last file is needed only for the branch-normalization check.

## Relation with reduction at 23

The degree-one stored canonical map has an unnormalized third branch
value lambda. The exact critical algebra has dimension eight and
satisfies N=lambda*D. The unit lambda has residue 7.
Consequently its stored singular coordinate 16 becomes 16/7=-1
when the third branch value is scaled to one. The other two maps
are already normalized. The uniform singular-position polynomial is

    (t+1)(t^2+1)(t^2+t+1).

The harmonic position -1 and the global osculating incidence select
the same component. Theorem 1.4 now constructs a special map with this
marked position and proves incidence on its lift. It uses explicit tail
covers, marked descent, the effective ramification divisor, same-sheet
compatibility, and exact local uniqueness; harmonic position alone is
not sufficient. See [the construction and audit](HARMONIC_RECONSTRUCTION.md).

## Verification record, 6 September 2026

The two new targets passed, including the printed compact incidence
equations and the projective-coordinate repetition.
The existing canonical-model, marked-point, branch-cycle-summary,
Galois-closure, local-23, pointed-reduction, and Frobenius checks also passed.
The branch-cycle-summary check validates stored permutations and
their provenance; the 150,145 Arb continuation tubes were not replayed.
The full verify-all suite was not rerun for this revision.

These new exact computations use SageMath. The existing Magma records
do not independently certify the new osculating theorem or branch scaling.
