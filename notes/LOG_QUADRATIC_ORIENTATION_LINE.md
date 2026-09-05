# The effective logarithmic quadratic-orientation line

> **Status correction, 2026-09-05.** The construction claimed below remains
> unproved. The finite weight identities do not supply a relative
> geometric object, a symmetric monoidal half-graded determinant, or an
> identification of its special parity with the singular-position
> idempotent. In particular, the two alleged filtrations are not established.
> Retained as research history, not as a proof or a certificate.
> See CORRECTNESS_REPAIRS_2026_09_05.md.

## Status

This note constructs the missing orientation object in the category in which
it can actually exist: pointed logarithmic graph-normalization complexes
decorated by their **effective integral graph cycles**.  The construction
does not factor through the ordinary derived category with `F_2`
coefficients.  In particular, it does not turn the zero coefficient
Bockstein into a nonzero one.

The object is a Clifford-volume, or half-graded determinant, line.  Its
specialization is obtained from the normalization--conductor filtration
before reducing the graph cycles modulo two.  This gives a functorial
factorization

```text
generic quadratic-orientation line
    = finite split-node anomaly line
      tensor returned quadratic-orientation line
      tensor canonically cancelling internal factors.
```

For the Mathieu data the three parities are respectively

```text
1+epsilon(Theta),       1,       q(n),
```

so the factorization transports the required quadratic refinement.

## 1. The half-graded determinant

Let `R` be a ring in which `2` is invertible and let `E` be a finite free
permutation module with its standard symmetric form: the distinguished
finite basis is orthonormal.  Write

```text
vol(E)=det(E).
```

If `rank(E)=2m`, give this line the quadratic norm

```text
qvol_E(omega)=(-1)^m det(<,>_E)(omega,omega).        (1)
```

Equivalently, place `E` in its Clifford algebra.  For an orthonormal basis
`e_1,...,e_2m`, the volume element satisfies

```text
(e_1 ... e_2m)^2=(-1)^m.                            (2)
```

Changing the ordering changes the volume element by a sign and therefore
does not change the line or its norm.  A bijection of finite bases is an
isometry and induces the corresponding map of volume lines.  Thus (1) is
functorial without choosing an ordering.

Denote the resulting quadratic line, together with the parity `m mod 2`, by

```text
hdet(E)=(det(E),qvol_E,m mod 2).
```

For even-rank `E` and `F`, the definition extends to their virtual
difference:

```text
hdet(E-F)=hdet(E) tensor hdet(F)^(-1),
parity hdet(E-F)=(rank(E)-rank(F))/2 mod 2.          (3)
```

It also extends by alternating tensor product to a perfect based complex of
even Euler rank.  Exact triangles and finite filtrations give the usual
determinant isomorphisms.  Formula (3) refines them by remembering the Euler
rank modulo four rather than only modulo two.

This is best regarded as a quadratic Picard groupoid, not as the ordinary
Picard groupoid of superlines.  Its commutativity constraint is the one
induced by determinant and Clifford multiplication; the extra sign in a
quadratic sum is recorded by the polarization line below.  Calling it an
ordinary square root of the determinant would suppress exactly the datum
we need.

## 2. The quadratic line of two effective cycles

Let `A` and `B` be two effective `0/1` cycles in a common finite tag set
`W`, of the same parity.  On the split graph algebra they are commuting
idempotents `e_A,e_B`.  The element

```text
s(A,B)=e_A+e_B-2e_A e_B                              (4)
```

is again an idempotent.  Its image is the permutation module on the
symmetric difference `A triangle B`.  Its rank is even, so define

```text
QOr_W(A,B)=hdet(im s(A,B)).                          (5)
```

The parity of (5) is exactly

```text
deg QOr_W(A,B)=|A triangle B|/2 mod 2=Q(A+B).        (6)
```

This formulation matters.  The two idempotents in (4) are retained before
reduction; `A+B` is not first replaced by an anonymous vector in a binary
permutation module.

The polarization is also line-valued.  If `S` and `R` are even subsets,
then there is a canonical based decomposition

```text
Z[S] direct_sum Z[R]
 =Z[S triangle R] direct_sum Z[S intersection R]^(direct_sum 2).
```

Applying `hdet` gives

```text
hdet(S) tensor hdet(R)
 =hdet(S triangle R) tensor Pol(S,R),                (7)
```

where

```text
deg Pol(S,R)=|S intersection R| mod 2.               (8)
```

Taking parities in (7) is precisely

```text
Q(S+R)=Q(S)+Q(R)+<S,R>.
```

The associativity of the line isomorphisms follows from the corresponding
identity for multiplicity functions on `W`; it is not an additional
choice.  This is the quadratic, rather than linear, composition law that an
ordinary determinant line misses.

## 3. The decorated logarithmic normalization complex

Let `S` be a strictly henselian trait and let a pointed semistable
`G`-cover over `S` be equipped with the divisorial log structure from the
special fibre and the horizontal finite branch divisor.  Normalize the
union of the relative deck graphs before pushing it to the natural
degree-`23` quotient.  After the harmless splitting base change, its graph
components form a finite based local system on the marked component tree
joining the finite end and the two wild ends.

The two wild endpoint traces give two effective idempotents on the generic
graph normalization.  Put

```text
P_Dy=C-T_Dy,                  P_Dz=C-T_Dz.
```

These are effective because both traces lie in `C`, and their symmetric
difference is the symmetric difference of the two traces.  Retain the
ordered pair `(P_Dy,P_Dz)` through normalization.  At the horizontal finite
modification this means keeping the full effective class `C` and the
effective trace `T_D` separately.  Do **not** replace their complement by
an arbitrary congruent integral lift.

The two relevant fibre presentations of the resulting ordered decoration
are

```text
QOr(P_Dy,P_Dz)=QOr(T_Dy,T_Dz)                  generically,
hdet(P_D) tensor QOr(T_D,T_D^n)                specially,
```

up to the polarization line.  The latter has degree zero for the Mathieu
packets.  This presentation is obtained only after taking the conductor
filtration: cancelling the two copies of `C` on the generic fibre must not
erase their provenance at the finite conductor.

It is allowed to have a nonzero conductor boundary; that boundary is part
of the decorated object.  The construction and its coherence are proved
abstractly in `GENERAL_EFFECTIVE_NORMALIZATION_GLUING_LEMMA.md`.

For a node `e`, restrict the two effective graph cycles separately to the
two branches of the normalization.  Let `B_(e,+)` and `B_(e,-)` be the
resulting finite based incidence modules after quotient pushforward.  Their
coefficient difference is divisible by two.  This does **not** say that
modules of the two (possibly unequal) ranks become isomorphic modulo two.
Write

```text
D_e=B_(e,+)-B_(e,-)                                 (9)
```

has even virtual rank.  Define the local anomaly line

```text
A_e(C,T_D)=hdet(D_e).                                (10)
```

Let the relative normalization--conductor complex of the decorated graph
closure over the trait be denoted by `mathscr N`.  Its special-fibre
cellular presentation is

```text
N^bullet=[direct_sum_vertices B_v
          --> direct_sum_edges B_e],                (11)
```

with branch-difference differential.  We use it relative to all marked
ends.  The ordered effective idempotents make (11) a **decorated** perfect
complex: its terms remember which summands came from `C` and which came
from `T_D`.  The hyperbolic completion

```text
H(mathscr N)=mathscr N direct_sum RHom(mathscr N,R)[1] (12)
```

is canonically `1`-shifted self-dual, and `mathscr N` is its
distinguished effective Lagrangian.  The Pfaffian-type line of (12) is, by
definition, the
half-graded determinant of this Lagrangian.  Replacing `N` by its
unrestricted mod-two image forgets the Lagrangian and destroys the
construction.

The vertex filtration of (11), or equivalently the
normalization--conductor filtration of the relative graph closure, has
three kinds of graded terms:

1. the symmetric-difference module of the two generic traces;
2. the virtual branch modules `D_e` at the logarithmic nodes; and
3. the symmetric-difference module of the two returned traces at the wild
   marked end.

Let `O_S(C,T_D)` be the half-graded determinant line of this relative
decorated complex over the trait.  The determinant isomorphism for a
filtered perfect complex gives the two fibre descriptions

```text
O_S(C,T_D)|_eta = QOr(T_Dy,T_Dz) tensor U_eta,
O_S(C,T_D)|_s   = QOr(T_D,T_D^n)
                  tensor (tensor_e A_e(C,T_D))
                  tensor U_s.                       (13)
```

Here `U_eta` and `U_s` collect the normalization summands that occur in
opposite cohomological degrees.  The relative normalization contraction
identifies their determinant lines and quadratic norms.  Thus (13) is a
specialization of one line over `S`, rather than a purported canonical
linear map between vector spaces over the different residue fields of the
trait.  After cancelling the paired `U` factors it gives the orientation
factorization written in the Status section.  All terms come from one
decorated normalization complex.  No counit is inserted into a
nearby-cycle pairing, and no tagged dot product is replaced by a bucketwise
product.

### Why internal terms cancel

Cut the marked component tree at all marked vertices and apply the
determinant isomorphism to the one-node charts in each resulting chain.  At
an unmarked valence-two internal vertex the
same based module occurs once in degree zero and once in degree one.  The
two determinant factors are dual and cancel canonically.  Gluing two chain
pieces at a marked vertex is the relative Mayer--Vietoris triangle, so the
cancellation is compatible with composition.

A semistable subdivision inserts the acyclic based complex

```text
[B --identity--> B].                                (14)
```

Its determinant and its quadratic norm have their canonical unit
trivialization.  Thus (13) is invariant under log blowups and does not
depend on a chosen semistable subdivision.

## 4. Functoriality

The construction has the required functorialities.

- **Simultaneous relabelling.**  A bijection of graph tags is an isometry
  of every permutation module and induces (13).  Wewers's common `lambda`
  changes all endpoint and node labels by one such relabelling.
- **Base change.**  Finite projective pullback, determinant, the standard
  permutation form, normalization triangles, their perfect conductor
  resolutions, and their finite filtrations commute with saturated
  log-etale base change.  Galois descent is through the same permutation
  isometries.
- **Branch exchange.**  Exchanging the two branches sends `D_e` to `-D_e`
  and the anomaly line to its dual.  Its parity is unchanged because
  `-kappa=kappa` modulo two.  Reversing the oriented interval dualizes the
  entire isomorphism (13).
- **Composition.**  Gluing pointed marked trees uses relative
  Mayer--Vietoris.  Determinant additivity composes the isomorphisms; the
  doubled-overlap term is exactly the polarization line in (7), so there is
  no hidden associativity sign.
- **Subdivision.**  Formula (14) gives canonical invariance.

These properties show that `O_S` is a functorial relative quadratic line
whose two fibre descriptions give logarithmic specialization of the
decorated orientation.  The specialization is an isomorphism of relative
quadratic-line objects; it does not choose a scalar trivialization of their
unavoidable `+-1` orientation torsor.  It is not a functor of the
underlying ordinary nearby-cycle complex alone.

## 5. Evaluation for the Mathieu cycles

### Generic end

Both traces have size `253`.  Hence

```text
rank(T_Dy triangle T_Dz)
  =2*(253-|T_Dy intersection T_Dz|),
```

and therefore

```text
deg QOr(T_Dy,T_Dz)
  =253-|T_Dy intersection T_Dz|
  =1+epsilon(Theta) mod 2.                           (15)
```

### Returned end

The returned traces are equal for a square return and are disjoint in the
adjacent Mathieu copy for a nonsquare return.  Their symmetric-difference
ranks are `0` and `506`.  Consequently

```text
deg QOr(T_D,T_D^n)=0,253 mod 2=q(n).                 (16)
```

### Finite horizontal node

For the ordered effective pair `(C,T_D)`, the exact branch incidence ranks
are

```text
rank B_+=1078,             rank B_-=112.
```

Thus (10) has parity

```text
deg A_fin(C,T_D)
 =(1078-112)/2=483=1 mod 2.                         (17)
```

The conductor filtration identifies this line with the orientation line
of the effective packet `P_D=C-T_D`: its rank is `3542`, and

```text
3542/2=1771=1 mod 2.
```

At every other internal node, Wewers's common `lambda` identifies the
two based restrictions, so `D_e=0` as a virtual based module.  The
affine-normalizer mismatch is retained at the returned marked end and is
already the line in (16).  Therefore taking parities in (13) gives

```text
1+epsilon(Theta)=1+q(n)=1+e_sing.                   (18)
```

Cancelling the common finite unit proves

```text
epsilon(Theta)=e_sing.
```

Equation (18) is not an interpolation from the seven values.  It is the
parity of the functorial line isomorphism (13), evaluated using the exact
generic trace, finite conductor, and returned trace modules.

This cancellation must not be justified merely from the ordinary parity of
the wild different.  For example, the `A2` chart has different exponent
`38`, whose half is odd.  That different belongs to the ordinary
graph--diagonal trace calculation, whereas the quadratic-orientation
anomaly is the virtual based restriction `D_e`.  Common-`lambda`
patching makes the latter zero and places the nontrivial return in (16).

## 6. Why the effective lifts cannot be discarded

The natural sheet operators satisfy

```text
rho(C)=-I mod 4,       rho(T_D)=I mod 4,
rho(P_D)=2I mod 4.
```

The congruent virtual lift `P_D-2C` maps to zero modulo four.  Under the
half-graded determinant, however, changing the effective lift by `-2C`
changes the parity by

```text
rank(C)=3795=1 mod 2.                               (19)
```

Thus the two lifts define opposite quadratic orientations.  Formula (19)
is not a defect of the construction; it is the proof that the construction
does not factor through the object that ordinary nearby cycles retain.
The ordered effective cycles `(C,T_D)` are part of its domain.

## 7. Scope

The construction proves the orientation transport in the exact category of
effective graph-normalization complexes attached to the pointed log stable
map.  It does not assert that this line is the determinant of ordinary
nearby cycles, and it does not identify the divided class with a coefficient
Bockstein.  The self-dual object is the hyperbolic completion (12); its
extra datum is the effective Lagrangian `N`.  Forgetting that Lagrangian
recovers the earlier vanishing result.

The only geometric inputs are the relative deck-graph closures, their
normalization--conductor complexes, and the common group transport in the
pointed patching datum.  The exact finite certificates evaluate the ranks
of the resulting based modules; they do not substitute for the
specialization isomorphism.

The general gluing proof, including the integral coefficient category,
Ferrand conductor sequence, telescoping identity, and all coherence maps,
is in `GENERAL_EFFECTIVE_NORMALIZATION_GLUING_LEMMA.md`.

## Reproduction

Run

```text
python3 notes/audit_log_quadratic_orientation_line.py
magma -b notes/audit_log_quadratic_orientation_line.m
```

The script checks the half-rank parities, polarization law, change of
effective lift, branch exchange, and subdivision invariance used above.

## References

- Finn Knudsen and David Mumford, *The projectivity of the moduli space of
  stable curves. I: Preliminaries on det and Div*, Math. Scand. 39 (1976),
  19--55.
- Chikara Nakayama, *Nearby cycles for log smooth families*, Compositio
  Math. 112 (1998), 45--75.
- Stefan Wewers, *Three point covers with bad reduction*, J. Amer. Math.
  Soc. 16 (2003), 991--1032.
