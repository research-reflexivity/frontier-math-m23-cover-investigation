# Re-derivation of the relative cohomological connector

## Status

This note restarts the proposed characteristic-`23` connector from the
definitions of a cohomological correspondence and of nearby cycles.  Its
conclusion is negative but precise:

* the lifted equality incidence has a genuine generic
  Lefschetz--Verdier realization;
* the nearby stalk of the finite Ferrand pinch is the stated permutation
  module;
* the claimed passage from the tagged dot product to the untagged bucket
  product is **not** the nearby-cycle pairing of those correspondences; and
* the mod-four node calculation is a valid secondary Bockstein calculation,
  but it has not yet been realized as the specialization of the generic
  packet pairing.

Consequently the current pinched-tag lemma does not prove the connector.

## 1. A completely named local Lu--Zheng datum

Let `S=Spec(R)` be a strictly henselian trait, with generic point `eta` and
closed point `s`.  Let `p:E->Q` be a finite surjection and let

```text
A_p=R^E fiber_product_(k^E) k^Q,
P=Spec(A_p).
```

For every `e in E`, projection to the `e`-th factor gives a normalization
section

```text
s_e:S -> P.
```

Fix binary endpoint vectors `alpha,beta in Lambda^E`, where
`Lambda=F_2`.  The following are actual data in the notation of
Lu--Zheng:

```text
X=S,                         L=Lambda_eta,
Y=P,                         M=Lambda_(P_eta),
C_alpha=disjoint union_(alpha_e=1) S,
D_beta =disjoint union_(beta_e =1) S.
```

The correspondence `C_alpha` goes from `X` to `Y`: its left map is the
structure map and its right map on the `e`-component is `s_e`.  The
correspondence `D_beta` goes from `Y` to `X` with the two maps reversed.
On the generic fibre every `s_e` is an open-and-closed point of `P_eta`, so
the cohomological maps

```text
u_alpha:c_left^* L -> c_right^! M,
v_beta :d_left^* M -> d_right^! L
```

are the identity on the selected components and zero elsewhere.

The generic pairing support is

```text
F_eta=(C_alpha)_eta fibre_product_(P_eta) (D_beta)_eta.
```

It contains the label `e` exactly when `alpha_e=beta_e=1`.  Therefore

```text
integral <u_alpha,v_beta>_eta
  =sum_e alpha_e beta_e.                              (1)
```

This is the tagged dot product.  This construction is the finite local
model of the lifted equality incidence used in the paper.

## 2. The correct nearby-cycle calculation

For `q in Q`, put `E_q=p^(-1)(q)`.  The strict local geometric generic
fibre of `P` at `q` is the finite discrete set `E_q`.  Hence

```text
(R^0 Psi_P M)_q=Lambda[E_q],
(R^i Psi_P M)_q=0 for i>0.                           (2)
```

After the splitting base change the inertia action on this basis is
trivial.  Without that base change it is the evident permutation action.

Write `V_q=Lambda[E_q]`.  The specialized cohomological maps are

```text
psi(u_alpha):Lambda -> V_q,   1 |-> sum alpha_e[e],
psi(v_beta) :V_q -> Lambda,   [e] |-> beta_e.
```

Their composite and their Lefschetz--Verdier pairing are therefore

```text
psi(v_beta) psi(u_alpha)=sum_e alpha_e beta_e.        (3)
```

Equations (1) and (3) are exactly Corollary 3.10 of Lu--Zheng in this
finite example.

There are also canonical specialization and trace maps

```text
delta_q:Lambda -> V_q,       1 |-> sum_e[e],
sigma_q:V_q -> Lambda,       sum a_e[e] |-> sum a_e.
```

But replacing the maps in (3) by compositions with `delta_q` and
`sigma_q` produces the different scalar

```text
sigma_q(alpha) sigma_q(beta)
  =(sum_e alpha_e)(sum_e beta_e).                    (4)
```

This is the bucket product.  It is not the pairing of
`psi(u_alpha),psi(v_beta)`.

The two-branch test makes the distinction unavoidable.  For

```text
E={1,2}, Q={q}, alpha=(1,0), beta=(0,1),
```

the generic and nearby pairings in (1)--(3) are zero, while (4) is one.
If (4) were the nearby pairing, Lu--Zheng's theorem would fail in this
example.

## 3. Why the conductor triangle does not change the answer

The constant-sheaf normalization--conductor triangle on `P` is

```text
Lambda_P -> nu_*Lambda_(S^E) + i_*Lambda_Q
         -> i_*p_*Lambda_E -> +1.                   (5)
```

On restriction to the generic fibre the two conductor-supported terms
vanish and (5) becomes the identity of `Lambda_E`.  Applying nearby cycles
to that generic triangle therefore gives the identity of `V_q`; it does
not create the map `V_q->Lambda` in (4).

Equivalently, proper pushforward along `p` sends the coefficient object
`Lambda_E` to `p_*Lambda_E`, whose stalk is `V_q`.  It does **not** replace
that object by the constant rank-one sheaf `Lambda_Q`.  Applying the
counit

```text
p_*Lambda_E -> Lambda_Q
```

is an additional operation and changes the cohomological correspondence.

The vanishing-cycle complex is correctly described by

```text
R Phi_P(M)_q=Cone(delta_q:Lambda->V_q).               (6)
```

The off-diagonal identity

```text
bucket=dot+off_diagonal
```

is a valid identity of bilinear forms.  It does not by itself put the two
endpoint maps on the cone (6).  Such an induced action requires
compatibility with the line `delta_q(Lambda)`.  A diagonal selector
`diag(alpha)` does not preserve that line unless `alpha` is constant on the
bucket.

## 4. The genuine generic connector

On a split pointed stable `M23`-model, let `B` be the normalized finite
`2A` branch-orbit scheme and let `I` be the normalized graph-orbit
incidence scheme.  For the two endpoint groups, take

```text
X=B,              L=Lambda_(B_eta),
Y=I,              M=Lambda_(I_eta),
C=L_Dy,           D=L_Dz,
```

where `C` maps `b` to `(w,b')` and `D` is the reverse correspondence.
The maps `u` and `v` are the canonical identity/fundamental-class maps on
the finite etale generic components.  Then

```text
F_eta=C_eta fibre_product_(X_eta times Y_eta) D_eta
```

is literally the lifted equality incidence.  The Fano--affine odd
fixed-point lemma gives

```text
degree <u,v>_eta=epsilon(Theta).                     (7)
```

This supplies the requested objects `X,Y,C,D,L,M,u,v` on the generic
fibre.  What is not presently constructed is a relative model whose
specialized maps can be computed as the claimed returned wild term plus a
universal finite unit.

In particular, replacing `I` by a special pinching that forgets the tag
does not solve the problem: (2)--(5) show that nearby cycles retain the
normalization labels.

## 5. What remains valid in the mod-four calculation

There is a correct secondary local calculation.  At a split tame node

```text
W=Spec R[[u,v]]/((u-v)(u+v))
```

the normalization quotient sends branch coefficients `(a_+,a_-)` to
`a_+-a_-`.  If the distinguished integral coefficient agrees on the two
branches modulo two, its first mod-four gluing obstruction is

```text
kappa(a_+,a_-)=(a_+-a_-)/2 mod 2.                   (8)
```

For the pushed augmented Mathieu packet the two coefficients are

```text
a_+=1078, a_-=112,
kappa=(1078-112)/2=483=1 mod 2.                      (9)
```

Equations (8)--(9) are a genuine Bockstein in the normalization sequence.
They are not, without further construction, the ordinary nearby-cycle
pairing of the lifted equality incidence in (7).

A viable relative proof would have to construct a self-dual two-level
complex retaining both graph provenance and quotient coefficients.  One
possible algebraic shape is

```text
K=[A_graph -> B_quotient],
```

where the augmented packet is a mod-two cycle because its integral image
is even, and (9) is its Bockstein at the tame node.  To finish the theorem
one must still:

1. define the two cohomological packet maps on `K_eta`;
2. construct the duality pairing whose generic value is
   `1+epsilon(Theta)`;
3. show that nearby cycles carry this same pairing to the returned
   Klein-four/Fano conductor class; and
4. prove that the latter has value `1+q(n)` under the pointed wild return.

The existing finite certificates prove the numerical values in items 3--4
once that comparison exists.  They do not construct the comparison.

## 6. Diagnostic using honest common tags

Intersecting the endpoint tag sets before pushing to sheet-pair branches
does not reproduce the universal Hadamard unit.  For the seven stored
Nielsen representatives, the parity of the elementary common-tag branch
difference is

```text
[1,1,0,1,0,1,0],
```

whereas multiplying the already pushed branch multiplicities gives the
constant value

```text
77^2-8^2=1 mod 2.
```

This is only a diagnostic of the sheet-graph model, not a replacement for
the missing relative calculation, but it confirms that the constant unit
comes from tag-forgetting rather than from the tagged pairing.

Run

```text
sage -gap -A -q notes/explore_tagged_tame_boundary.g
```

to reproduce the seven common-tag profiles.

## Conclusion

The correct local nearby-cycle computation preserves the tagged dot
product.  The bucket product is obtained only after an extra counit and
therefore belongs to a different cohomological correspondence.  The
mod-four fixed-node Bockstein remains a serious candidate for a corrected
proof, but the relative self-dual complex and its identification with the
wild returned class have not yet been constructed.  Until those maps are
given, the equality `epsilon(Theta)=e_sing` is established only by the
independent exact arithmetic/branch-cycle computation, not by the current
pointed nearby-cycle argument.
