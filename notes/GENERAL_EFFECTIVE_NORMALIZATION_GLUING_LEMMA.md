# The effective normalization--gluing lemma

> **Status correction, 2026-09-05.** This is an archived proposed argument,
> not a proved normalization--gluing theorem. The telescope of divided
> integer differences is valid, but it does not construct a relative
> quadratic determinant or its asserted generic and special filtrations.
> The half-rank grading also needs coherence data not supplied here.
> The manuscript no longer invokes this lemma. See
> CORRECTNESS_REPAIRS_2026_09_05.md for the proved replacement and limits.

## Purpose

This note isolates the algebraic statement used by the logarithmic
quadratic-orientation construction.  The statement is deliberately made in
the exact category of **based integral coefficient systems on normalized
graph components**.  It is not a theorem saying that an ordinary nearby-cycle
complex remembers a chosen effective lift.

The distinction matters.  At the finite Mathieu node the two branch
coefficients are `1078` and `112`.  Modules of those ranks do not become
isomorphic after reduction modulo two.  What is true, and what the lemma
uses, is that their coefficient difference is divisible by two.

## 1. Split normalization data

Let `R` be a ring in which `2` is invertible.  The integral normalization
coefficients below live first in `Z[W]`; tensoring their based lattices with
`R` gives the modules on which the quadratic lines are formed.  Let `Gamma`
be a finite tree with a specified set of marked vertices, and fix
two of them, `v_-` and `v_+`.  The marked set is required to contain
every vertex of valence different from two.  At every unmarked vertex the
two restriction maps identify the same basis and the decoration agrees
under that identification.  At every vertex
`v` let `W_v` be a finite based set.  For every edge `e=vw`, choose a finite
based set `W_e` and based restriction maps from the two incident vertex
sets.  In the applications these are the restrictions of
deck-transformation graph labels to the two branches of a split node.

Write

```text
B_v=R[W_v],                   B_e=R[W_e].
```

The associated normalization complex is

```text
N_Gamma(B)=[ direct_sum_v B_v  -->  direct_sum_e B_e ],        (1)
```

in degrees zero and one.  On the `e=vw` component the differential is the
difference of the two branch restrictions.  This is the cellular form of the
normalization sequence

```text
0 -> B_X -> nu_* B_Xtilde -> direct_sum_e i_(e,*) B_e -> 0.    (2)
```

For a split nodal curve, (2) follows directly from the completed local ring

```text
k[[u,v]]/(uv) -> k[[u]] x k[[v]],
```

whose quotient on locally constant coefficients is branch difference.
For a Ferrand pinching it is the module sequence associated with the fibre
product of rings.  In particular, (1) is a finite based perfect complex.

An **effective decoration** is an integral coefficient vector

```text
z_v=sum_(w in W_v) n_(v,w)[w],       n_(v,w)>=0,               (3)
```

on every normalized component.  The coefficients, not merely their
reductions, are part of the object.  Transport the two restrictions of (3)
to the common edge basis and put

```text
delta_e(z)=z_(e,+)-z_(e,-) in Z[W_e].                         (4)
```

We call the decoration **two-divisible at the conductor** when every
coefficient of (4) is even.  Then

```text
kappa_e(z)=delta_e(z)/2 mod 2 in F_2[W_e]                     (5)
```

is its divided conductor class.

For an integral multiplicity vector `m=sum_w n_w[w]`, write
`R[m]=sum_w n_w R[w]` in the Grothendieck group of finite based free
`R`-modules; a negative coefficient contributes the corresponding dual
determinant factor.  Applying such vectors to the normalization filtration
defines a virtual based perfect object.  This is what “decorated complex”
means below.  It is formed before reduction modulo two and need not be a
subcomplex of the undecorated constant-coefficient complex.  We assume that
all unpaired virtual ranks are even whenever `hdet` is used.

## 2. Statement of the lemma

> **Effective normalization--gluing lemma.**  For a two-divisible effective
> decoration on a split nodal marked tree, the following hold.
>
> 1. The local classes (5) are the image of the divided normalized
>    coefficient under the differential of (1).  Hence they glue to one
>    conductor class; they are not independently chosen node invariants.
> 2. Along the unique path from `v_-` to `v_+`, after transporting all edge
>    bases to one base point,
>
>    ```text
>    sum_e kappa_e(z)=(z_(v_+)-z_(v_-))/2 mod 2.              (6)
>    ```
>
>    Every unmarked valence-two transport component on this path occurs
>    twice and cancels.  Marked branching and endpoint terms are retained.
>    Before reduction, the same identity gives the canonical line
>    transport
>
>        hdet(R[z_(v_+)]-R[z_(v_-)])
>          = tensor_(e in P) hdet(R[delta_e(z)]).             (6a)
>
>    The factors are ordered along the oriented path.
> 3. Let `hdet` denote the half-graded determinant of an even-rank based
>    virtual module.  The Knudsen--Mumford determinant isomorphism for the
>    vertex filtration of (1) gives a canonical quadratic-line
>    factorization
>
>    ```text
>    hdet(N_Gamma(B),z)
>       = marked lines tensor (tensor_e hdet(delta_e(z)))
>         tensor H_Gamma,                                      (7)
>    ```
>
>    where the marked factor is the tensor product of the marked vertex
>    lines and `H_Gamma` is the tensor product of the internal based factors
>    occurring once in each of two adjacent normalization triangles.  It
>    has its canonical evaluation trivialization.  Thus only the marked
>    vertex terms and the nonzero divided conductor terms remain.
> 4. The construction is functorial under simultaneous based relabelling,
>    flat base change preserving the split normalization diagram, reversal
>    of a branch, gluing of marked trees, and subdivision of an
>    edge.  Reversal dualizes the relevant virtual line; subdivision inserts
>    the based acyclic complex `[B --1--> B]`.

Here `hdet(N_Gamma(B),z)` means the half-graded determinant of the virtual
based object obtained by applying the integral multiplicity vectors to the
normalization filtration.  It is an object of the quadratic Picard
groupoid, not an ordinary determinant of the reduction of (1).

## 3. Proof

The assertion is local at a node before it is global.  At a split node the
normalization sequence with an arbitrary coefficient module `M` is

```text
0 -> M -> M x M -> M -> 0,
     a |-> (a,a),        (b_+,b_-) |-> b_+-b_-.                (8)
```

Apply (8) coefficient by coefficient to (3).  If (4) is divisible by two,
the reduction of (3) modulo two lies in the kernel of branch difference and
therefore descends across the node.  Write
`z_(e,+)-z_(e,-)=2h_e`.  Division by two in the kernel of
`Z/4 -> F_2` is unique, and `h_e mod 2` is (5): the cochain measuring
the failure of the two branch lifts to agree modulo four.  Thus (5) is the
first obstruction to descending the chosen integral lift modulo four.
This proves the first claim.
It also proves that (5) is unchanged when the two branches are exchanged:
the integral difference changes sign, which has no effect over `F_2`.

Orient the path from `v_-` to `v_+` and write it as

```text
v_0 --e_1-- v_1 -- ... --e_r-- v_r.
```

Use the based edge transports to express every coefficient in the basis at
`v_0`.  Summing (4) gives

```text
(z_1-z_0)+(z_2-z_1)+...+(z_r-z_(r-1))=z_r-z_0.       (9)
```

Division by two and reduction modulo two give (6).  This is the promised
global gluing: the cancellation is an identity of coefficient vectors, not
only an equality of their total parities.  Before reduction, (9) is an
identity of virtual based modules.  Applying hdet in path order proves
(6a), including its quadratic norm and parity.

For the line statement, first cut the tree at all marked vertices.  Each
remaining piece is a chain, and every internal vertex of that chain is
unmarked with identical based restrictions.  Filter each chain by
successively adjoining its vertices and edges.  Each step is a genuine
short exact sequence of based complexes.  Determinant additivity therefore
identifies the determinant of the whole filtered complex with the ordered
tensor product of the determinants of its graded pieces.  These filtrations
are based split: every unmarked transition is a bijection of distinguished
bases.  Clifford volume is multiplicative on such orthogonal based sums;
the doubled-overlap term is the polarization line.  At an unmarked internal
vertex the same based summand occurs in two consecutive
normalization triangles, once covariantly and once contravariantly.  The
two determinant lines are dual, their quadratic norms are inverse, and
evaluation gives their canonical trivialization.  The unpaired graded
pieces are precisely the marked vertex pieces and the virtual branch
differences (4).  Applying the Clifford-volume refinement to the same exact
sequences proves (7).

No coherence choice is hidden here.  Associativity is the associativity of
the determinant functor on a finite filtration.  Gluing two marked trees
amounts to concatenating their filtrations.  If the common marked endpoint
is forgotten after gluing and becomes an unmarked valence-two vertex, its
two determinant factors are dual and evaluate to one; if the mark is
retained, so is its factor.  Subdividing an edge
adds `[B --1--> B]`, which is based contractible and has the canonical unit
determinant.  Flat base change preserves (8), the finite free modules, and
the filtration.  A based relabelling is an isometry of every permutation
module.  These observations prove all the stated functorialities.

A branch point or leaf is marked by definition and its line is retained in
the factorization.  Only valence-two transport vertices with the stated
based-isomorphism property are cancelled.

## 4. Relative form over a trait

Let `S=Spec(R)` now be a strictly henselian trait, let `A` be a finite
torsion-free graph-incidence `R`-algebra with finite normalization `A^nu`,
and let `z_eta` be an effective idempotent of `A^nu tensor_R K`.  Assume
that its conductor quotients carry the based branch identifications of
Section 1, that their coefficient differences are two-divisible, and that
the unpaired virtual ranks are even.  Since `z_eta` satisfies `X^2-X=0`, it
belongs to the integral closure `A^nu`; uniqueness follows from the
injection into the generic algebra.  Its image is a direct summand of the
finite torsion-free `R`-module `A^nu`, hence finite projective.

Suppose these normalized graph components and their conductor branches
come from a finite relative incidence scheme obtained from a pointed
semistable cover.  Write `c` for the conductor.  The Ferrand square

```text
A = A^nu fiber_product_(A^nu/c) A/c
```

and its module sequence supply the two restriction maps to each common
based conductor module.  Their difference defines the relative version of
(1).  When that difference vanishes its kernel is the descended module;
in general its cone records the descent defect.  Because the normalization
terms are finite flat and the conductor terms are finite over the discrete
valuation ring, hence have projective dimension at most one.  The resulting
complex is perfect.  Derived formation of its generic and special fibres commutes
with every flat base change preserving this conductor diagram.

The same construction applies componentwise to a finite ordered family of
effective idempotents.  Keeping this ordered family is essential: only
after taking the conductor filtration do we use the line-valued
polarization law to combine its factors.  In particular, two complementary
packet idempotents may have the same generic symmetric difference while
retaining different, and geometrically meaningful, conductor terms.

Consequently (7) is an isomorphism between two filtrations of one relative
quadratic line.  It is not a scalar comparison between unrelated lines over
the fraction and residue fields.  If local annular identifications carry
all graph labels through one common group isomorphism, they are morphisms of
this relative normalization diagram and the construction descends from a
splitting trait.

This is exactly the amount of geometry supplied by Wewers's patching datum:
one global datum is `(lambda;phi_j)`, with the same `lambda` in every local
pair.  A `G`-equivariant `phi_j` carries the graph of `g` to the graph of
`lambda(g)`.  Thus a change of patching datum simultaneously relabels every
term of (1), and (7) is invariant.

## 5. The Mathieu specialization

For the canonical effective augmented packet, quotient pushforward gives
the two branch coefficients

```text
1078 and 112.
```

The lemma does **not** call these modules isomorphic modulo two.  It applies
the last map in (8) to their coefficient vector and obtains

```text
(1078-112)/2=483=1 mod 2.                            (10)
```

The ordinary divided coefficient on the normalization is the fixed heptad.
The intrinsic degree-three Fano incidence carries it to the seven nonfixed
ramification nodes, so those seven coordinates cancel the corresponding
seven entries of the eight-node conductor packet.  The finite graded term
in (7) is therefore the single fixed-node line.

On every internal wild annulus, the common `lambda` identifies the two
based restrictions, so (4) is zero.  At the marked returned end the two
normalizer traces are equal for a square return and disjoint in adjacent
Mathieu copies for a nonsquare return.  Their symmetric-difference ranks
are `0` and `506`, and their half-graded parities are `0` and `253=1 mod 2`.
Thus the only unpaired terms are the finite unit and the returned trace
line, exactly as asserted in the paper.

## 6. Scope

The lemma proves the global gluing and all coherence statements for the
**effective coefficient normalization complex**.  It does not identify
that complex with ordinary nearby cycles, and it does not claim that the
coefficient boundary (10) is a Chow-theoretic refined intersection.  Those
identifications are neither used nor needed when the connector is defined
in the effective coefficient category.

The generic term is still geometric: on the etale locus of a marked
`G`-cover the normalized self-fibre product is the disjoint union of the
deck-transformation graphs.  Componentwise multiplication of their
idempotents is the common-label pairing.  Additive quotient pushforward
commutes with restriction to the two branches of (8), so the numbers in
(10) are the special branch values of the same distinguished integral graph
coefficient, not values imported from a separate calculation.

## References

- Daniel Ferrand, *Conducteur, descente et pincement*, Bull. Soc. Math.
  France 131 (2003), 553--585, especially Section 2.2.
- Finn Knudsen and David Mumford, *The projectivity of the moduli space of
  stable curves. I: Preliminaries on det and Div*, Math. Scand. 39 (1976),
  19--55.
- Stefan Wewers, *Three point covers with bad reduction*, J. Amer. Math.
  Soc. 16 (2003), 991--1032, Section 4.2.
