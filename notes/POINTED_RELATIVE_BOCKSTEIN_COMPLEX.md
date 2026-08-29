# The pointed quadratic connector and its divided gluing defect

## Status

This note records the corrected construction of the **local divided gluing
defect** and two obstructions to interpreting it as an ordinary Bockstein.
The functorial relative orientation line built from this datum is constructed
in `LOG_QUADRATIC_ORIENTATION_LINE.md`.

The main advance is that the class lives in the relative normalization
complex of the path between the finite and wild markings.  It is not an
ordinary coefficient Bockstein in the natural 23-sheet representation or
in the horizontal logarithmic nearby-cycle stalk, and it is not obtained by
replacing the tagged nearby-cycle pairing with a bucketwise product.  The
surviving invariant is instead the half-weight quadratic refinement attached
to the two distinguished effective graph lifts.

## 1. The unrestricted mapping cone kills the class

Let

```text
rho: Z/4[M23] -> End((Z/4)^23)
```

be the natural permutation representation.  Its adjoint is taken with
respect to the group-basis and matrix-entry pairings.  Double transitivity
gives

```text
|M23|/23      =443520=0 mod 4,
|M23|/(23*22)= 20160=0 mod 4.
```

Consequently `rho*rho^*=0 mod 4`, and the tempting three-term complex

```text
End(V)^* -> Z/4[M23] -> End(V)
```

is indeed self-dual.

Write `C` for the full `2A` graph sum, `T_D` for a normalizer trace, and
`P_D=C-T_D` for the effective complementary packet.  Their integral sheet
operators are

```text
rho(C)  =120 J+1035 I = -I mod 4,
rho(T_D)=  8 J+  69 I =  I mod 4,
rho(P_D)=112 J+ 966 I = 2I mod 4.
```

Thus the distinguished effective lift appears to have Bockstein `I`.
However, the congruent lift

```text
P_D-2C
```

has sheet operator

```text
-128 J-1104 I=0 mod 4.
```

It follows that the apparent class is zero in the cohomology of the
unrestricted mapping cone.  This is not merely a technical annoyance: an
ordinary derived object forgets which of the two congruent integral lifts
is the graph-effective one.

The same calculation explains the positive datum that has to be retained.
The two canonical effective correspondences `C` and `T_D` are separately
fixed.  They are two integral lifts of the same mod-two sheet operator, and
their relative comparison is

```text
rho(T_D)^(-1)*rho(C)=-I=1+2I mod 4.
```

Its first divided layer is canonically `I` once the two graph lifts, rather
than only their mod-two difference, are part of the object.

## 2. Horizontal log nearby cycles do not restore an ordinary Bockstein

Nakayama's log-purity theorem explicitly permits horizontal log structure.
For a log-smooth family over a henselian trait and

```text
L_n=R j_* Z/n
```

on the generic fibre, his tame nearby-cycle formula is

```text
R^q Psi_t(L_n)
  =R^0 Psi_t(L_n) tensor
    wedge^q(M_rel^gp tensor Z/n(-1)).
```

At a smooth point of one horizontal divisor the relative characteristic
group has rank one.  After the splitting base change used for the Mathieu
tags, the local stalk is therefore one free coefficient copy in degrees
zero and one and zero in higher degrees:

```text
R^0 Psi_t(L_n)=Z/n,
R^1 Psi_t(L_n)=Z/n(-1),
R^q Psi_t(L_n)=0 for q>1.
```

For the coefficient sequence

```text
0 -> F_2 -> Z/4 -> F_2 -> 0,
```

reduction from the `Z/4` version of the displayed stalk to the `F_2`
version is surjective
in both nonzero degrees.  Hence every coefficient connecting homomorphism
on this stalk is zero.  A split finite tag set only takes direct sums, so
the conclusion is unchanged.  Naturality of nearby cycles with respect to
the coefficient triangle then shows that a zero Bockstein in the
graph-to-sheet cone cannot become the unit class after specialization.

This rules out the strongest hoped-for shortcut.  The number `483 mod 2`
is not the coefficient Bockstein of a hidden logarithmic nearby-cycle
generator.  It is the divided failure of a **specified effective mod-four
graph lift** to glue.  Forgetting that lift, as the ordinary derived
category does, forgets the class.  Nakayama's formula remains useful for locating
the horizontal degree-one line, but a quadratic or cycle-level refinement
is required to put the Mathieu defect on it.

## 3. The correct split-node square

Let `Q` be a point of a pinched special fibre and let `E_Q` be its
normalization tags.  With `Lambda_n=Z/n`, the specialization unit is

```text
delta: Lambda_n[Q] -> Lambda_n[E_Q],
       [Q] |-> sum_(e in E_Q) [e].
```

The target is the nearby-cycle stalk.  It is essential not to replace it by
the rank-one source.

At a split node `E_Q={+,-}`.  The canonical integral packet acts on the two
normalized branches with coefficients

```text
a_+=1078,       a_-=112.
```

Modulo two both actions are zero, so they define a commutative square with
the zero action on `Lambda_2[Q]`.  The canonical mod-four lift does not
commute with `delta`; its defect is

```text
(a_+,a_-)=(2,0) mod 4.
```

After division by two, its class in

```text
Lambda_2[E_Q]/delta(Lambda_2[Q])
```

is represented by `(1,0)`.  It is the nonzero vanishing-cycle generator,
because the quotient map is branch difference.  Equivalently,

```text
kappa=(a_+-a_-)/2=483=1 mod 2.                 (1)
```

This is a genuine obstruction to lifting the **commutative specialization
square with its distinguished graph maps**.  It does not assert that an
arbitrary diagonal selector acts on the quotient stalk.  In particular,
the two-label counterexample with selectors `(1,0)` and `(0,1)` is excluded:
those selectors do not even preserve the specialization line modulo two.

The local duality is completely explicit.  Put

```text
V=F_2{e_+,e_-},       L=F_2(e_++e_-).
```

For the tagged dot product on `V`, the line `L` is isotropic and

```text
L=L^perp,             V/L=L^vee.
```

The class `(1,0)` obtained above generates `V/L` and pairs to one with the
specialization line.  Thus

```text
0 -> L -> V -> L^vee -> 0
```

is already the required self-dual two-level stalk; no identification of
the tagged and bucketwise pairings is involved.  The divided gluing defect
is the off-diagonal class of its distinguished mod-four lift.  Calling this
class a coefficient Bockstein would be misleading: it depends on retaining
that effective lift.

There is an equivalent quadratic description.  On the even-weight subspace
of a tagged binary permutation module, put

```text
Q(v)=weight(v)/2 mod 2.
```

For even-weight vectors one has the exact polarization identity

```text
Q(v+w)=Q(v)+Q(w)+<v,w>_tagged.                 (2)
```

Every Mathieu packet has even multiplicity in every ordered sheet-pair
bucket.  On the diagonal and transposition branches its quadratic values
are

```text
Q_+=1078/2=539=1 mod 2,
Q_-= 112/2= 56=0 mod 2.
```

Thus (1) is the branch difference of the canonical quadratic refinement of
the original tagged pairing.  This is stronger than a numerical match: it
explains why division by two and the tagged dot product occur in the same
two-level object.

Globally over the ordered sheet pairs, each graph label occurs for 23 source
sheets and `|P_D|=3542`, so

```text
Q(P_D)=23*(3542/2)=1 mod 2.
```

Consequently polarization of `Q` on two packets is their generic tagged
intersection parity.  At the returned end the full-class/trace conductor
has size `98*253`; its half-weight is odd for a square return and zero for a
nonsquare return.  These are exactly the values `1+q(n)`.  Section 7 records
how the effective quadratic determinant transports this same refinement,
rather than merely a function with the same two values, along the pointed
stable map.

## 4. Why the two markings make the class survive

Let `Gamma` be the path in the dual graph from the marked finite component
to the marked wild component, and let

```text
partial Gamma={v_fin,v_wild}.
```

For a path with `r` edges, the absolute cochain map

```text
C^0(Gamma,F_2) -> C^1(Gamma,F_2)
```

is surjective.  Thus an unpointed normalization complex again kills every
edge class.  Relative to the two ends, however,

```text
H^1(Gamma,partial Gamma;F_2)=F_2.              (3)
```

The sum of the edge coordinates detects its generator.  If `h` is a
divided normalization coefficient, telescoping gives

```text
sum_edges dh = h(v_fin)+h(v_wild).             (4)
```

Every internal value occurs twice.  Formula (4) is invariant under
subdivision and is precisely the endpoint formula already proved in the
tree-normalization lemma.  The new point is conceptual: it identifies the
target of that formula as top **relative** cohomology.  Allowing homotopies
at either marked end would destroy the class, exactly as in Section 1.

The relative cochain complex has the expected duality.  For the oriented
interval, cellular Poincare--Lefschetz duality pairs

```text
C^bullet(Gamma,partial Gamma;Lambda)
```

with the absolute chain complex in complementary degree.  Over `F_2` no
orientation sign remains.  Hyperbolic completion therefore gives an
explicit self-dual two-level finite complex.  With graph-label coefficients,
the same construction applies after transporting all labels by the common
group identification in the pointed patching datum.

## 5. Mathieu endpoint values

At the finite end, (1), the fixed-heptad residual, and the intrinsic Fano
incidence bridge give the fixed-node value

```text
h(v_fin)=1.
```

At the wild end, the pointed affine return gives

```text
h(v_wild)=q(n),
```

with values `0,1,1` on the rational `E8`, unramified `E8`, and ramified
`A2+A6` points.  Hence the pointed relative class has value

```text
1+q(n)=1+e_sing.                                  (5)
```

In the raw lifted-trace convention the normalized returned fixed-locus
term has value `1+q(n)` and the ordinary tame node has value `1`; their sum
is `q(n)`.  Equations (4)--(5) are the same endpoint calculation in the two
conventions.

## 6. Explicit candidate packet maps

The preceding calculation also removes the ambiguity in the schematic
notation `[A_graph -> B_quotient]`.  On the graph-normalized generic side,
use the actual inclusion of label modules

```text
K_graph(D_y)=[F_2[T_Dy] -> F_2[C]],
```

where `C=2A(M23)` and `T_Dy=x^D_y`.  Its degree-zero quotient is the
effective packet `P_Dy=C-T_Dy`.  The second endpoint trace `T_Dz` defines
the dual functional by the group-basis dot product.  Its Euler pairing with
the displayed two levels is

```text
|C intersection T_Dz|+|T_Dy intersection T_Dz|
  =1+epsilon(Theta).                                (6)
```

Thus the two packet maps are no longer hypothetical: the first is the cone
of the trace inclusion into the full class, and the second is evaluation on
the returned trace.  Adding the dual cone gives a hyperbolic self-dual
complex.  Equivalently, the self-duality is the polarization (2) of the
half-weight form.

After quotient pushforward at the finite branch, the two levels have the
same mod-two sheet operator but their canonical integral lifts are `-I` and
`+I` modulo four.  The divided differential is therefore the local class
(1).  Pairing its transported relative class with the returned endpoint
trace gives the already certified conductor identity

```text
b_x(P_Dy,T_Dz^n)+b_x(T_Dy,T_Dz^n)
  =b_x(C,T_Dz^n)
  =1+q(n).                                          (7)
```

Equations (6)--(7) give the requested generic and special traces from the
same two separately retained graph levels.  They also explain why working
only with the binary packet `P_D` was insufficient: it discarded the
integral comparison between `C` and `T_D` that defines the divided
differential.

There is a particularly small reformulation of the target.  For two
253-element trace subsets of a common graph-label set, define their relative
half-distance by

```text
Q(T_1,T_2)=|T_1 symmetric_difference T_2|/2 mod 2.
```

Then

```text
Q(T_Dy,T_Dz)=253-|T_Dy intersection T_Dz|
             =1+epsilon(Theta).                     (8)
```

At the returned end the two traces are equal for a square return and are
disjoint in the adjacent Mathieu copy for a nonsquare return.  Hence

```text
Q(T_D,T_D^n)=q(n).                                  (9)
```

The connector is therefore equivalent to the conservation law

```text
Q_generic = 1+Q_returned.                          (10)
```

The unit in (10) is exactly the split-node quadratic defect (1).  This turns
the comparison into a precise orientation statement: transport the
half-graded determinant of the odd trace subspace along the pointed
graph-component normalization and show that its only odd local change is the
finite `2A` node.  The construction in Section 7 does so; the even wild
differents proved in the local parity filter rule out an additional wild
contribution.

## 7. Resolution by the effective quadratic determinant

The comparison is carried by a half-graded Clifford-volume line, not by an
ordinary coefficient Bockstein.  For an even-rank based permutation module
`E`, its half determinant has underlying line `det(E)` and orientation
parity

```text
rank(E)/2 mod 2.
```

Applied to the symmetric-difference module of two effective trace
idempotents, this parity is exactly `Q_W`.  Applied to the virtual branch
module at the finite split node, it is

```text
(1078-112)/2=483=1 mod 2.
```

The pointed graph-normalization complex is hyperbolically completed while
the effective complex is retained as its distinguished Lagrangian.  The
determinant isomorphism for its normalization--conductor filtration then
gives one relative quadratic line with fibre parities

```text
generic:  Q_W(T_Dy+T_Dz)=1+epsilon(Theta),
special:  Q_W(P_D)+Q_W(T_D+T_D^n)=1+q(n).
```

Internal vertex factors occur in dual cohomological degrees, and a
semistable subdivision inserts an acyclic identity complex.  This proves
composition and subdivision invariance.  Simultaneous relabelling and
saturated log-etale base change act by isometries of the based permutation
modules.  The complete construction and its scope are given in
`LOG_QUADRATIC_ORIENTATION_LINE.md`.

The dependence on the effective pair is essential: replacing `P_D` by the
congruent lift `P_D-2C` reverses the half orientation because `|C|=3795` is
odd.  Ordinary nearby-cycle duality forgets this Lagrangian datum and
therefore continues to have zero coefficient Bockstein.

## Reproduction

Run

```text
python3 notes/audit_pointed_relative_bockstein.py
```

The existing group-theoretic endpoint audits remain

```text
sage -gap -A -q notes/audit_universal_branch_graph_packet_vanishing.g
sage -gap -A -q notes/audit_returned_half_graph_conductor.g
sage -gap -A -q notes/audit_fano_affine_incidence_bridge.g
```

## References

- Teruyoshi Yoshida, *On the action of algebraic correspondences on weight
  spectral sequences*, arXiv:1109.2208.
- Chikara Nakayama, *Nearby cycles for log smooth families*, Compositio
  Math. 112 (1998), 45--75.
- Qing Lu and Weizhe Zheng, *Duality and nearby cycles over general bases*,
  Duke Math. J. 168 (2019), 3135--3213, arXiv:1712.10216.
- Qing Lu and Weizhe Zheng, *Categorical traces and a relative
  Lefschetz--Verdier formula*, Forum Math. Sigma 10 (2022), e10,
  arXiv:2005.08522.
