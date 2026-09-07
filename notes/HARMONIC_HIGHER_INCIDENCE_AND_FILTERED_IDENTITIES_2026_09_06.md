# Harmonic incidence through eight orders, and the remaining induction identity

Private research, 6 September 2026. No manuscript or public files have
been changed. These calculations read no reconstructed generic cover.

Follow-up: HARMONIC_FIXED_HENSEL_SYSTEM_AND_INCIDENCE_IDEAL_2026_09_06.md
constructs a fixed eleven-variable analytic system with unit Jacobian.
This proves all-orders uniqueness for the selected necessary equations
and the stability of their correction matrix. Exact incidence remains
open; the follow-up isolates two values at the unique solution and
exhibits a defect in the naive constant-coefficient identity off the
solution locus.

## Status

Starting with the M23-compatible first jet selected in
HARMONIC_FIRST_JET_SELECTION_BY_WITT_DESIGNS_2026_09_06.md, necessary
degree-23-function and effective branch-fibre equations force the
osculating remainder to vanish through EIGHT orders:

    r0=r1=0 mod pi^9,       pi^2=-23.

Incidence is never imposed as an equation. It is evaluated independently
after solving the function and branch-fibre equations. At each stage its
value is also checked on the full affine solution space, including the
remaining undetermined next curve coefficient.

This is a finite-order result under the earlier explicit reduction-model
hypotheses. It is not proof of an all-orders incidence relation, global
rationality, or existence of a generic cover. The useful further finding
is that the correction matrix and short linear identities are identical
in orders 3 through 8. Section 6 identifies precisely what is still
needed to turn this observation into an induction.

## 1. Retain the effective fibre, not just its leading differential

On the smooth central chart let beta have its two total ramification
points at A,B. Normalize beta by beta(0), for computation only. Its actual
third branch value is then some scalar u close to one. Let

    eta=(d beta/beta)/23,
    R=z^8+pi*r1+pi^2*r2+...

be the monic Weierstrass polynomial of eta's effective zero divisor.
The three-point condition implies

    beta-u=R^2 F.

The square is essential. Previously, seven coefficients of eta1 were
used as necessary consequences of this identity. Here additional
coefficients of the identity itself are imposed.

The function is represented by a ratio N/D of canonical quintic sections,
with D a unit at z=0. Thus it is enough to require

    N-uD divisible by R^2.

One must not assume u=1 at every order: beta(0) first differs from the
critical value at order pi^5. The new calculation eliminates u by the
constant remainder. If P and Q are the remainders of N-D and D modulo
R^2, respectively, put a=P(0)/Q(0) and use P-aQ. The resulting scalar
u=1+a is the unique possible value for this truncated divisibility test.
All nonconstant remainder coefficients are retained as equations.

## 2. A precision lemma for the unknown ramification coefficients

This general algebraic observation justifies using a shorter expansion
of R than of beta. Suppose R0=z^8, r1 is divisible by z, and the special
quotient F0 in N-uD=R^2 F has order seven at z=0. Work modulo pi^(m+3),
where m>=2, but suppose R is known only modulo pi^(m+1).

Replacing R by

    R+pi^(m+1)h+pi^(m+2)j

cannot change the following selected remainder coefficients:

* order pi^m: all powers z^0,...,z^15;
* order pi^(m+1): powers z^0,...,z^14;
* order pi^(m+2): powers z^0,...,z^7.

To prove this, linearize division by the monic polynomial R^2. The first
change is minus 2R*h*F times pi^(m+1). Its reduction starts with
z^8*z^7=z^15. At the next order, the possible new terms involve r1*F0
or z^8*F1; both have order at least eight, since r1 is divisible by z.
Reducing terms of degree at least 16 modulo the deformed monic divisor
does not lower these bounds: its first variation is 2z^8*r1, of order
at least nine. The j term again starts with z^15. Quadratic changes in h
have order at least 2m+2, outside the retained precision.

The same bounds hold when the lower remainder vanishes only to the
preceding induction precision; its error terms enter beyond the selected
orders. Eliminating the scalar u does not spoil the bounds, because no
unavailable R coefficient changes the constant remainder at these orders.

The function modulo pi^(m+3) determines eta modulo pi^(m+1), hence the
needed R modulo pi^(m+1). No unknown higher R coefficient is set to zero
as a mathematical assumption. The checker additionally perturbs those
unavailable coefficients explicitly and verifies that every selected
coefficient is unchanged.

Weierstrass preparation is performed by successive monic polynomial
division, with a unit inverse modulo z^8. A series of length at least
8(m+1) suffices to determine this distinguished factor: omitted terms
cannot reduce to degree below eight within m positive pi-orders.

## 3. The second-order incidence is forced

Use the cubic coefficient complement

    H^3, H^2U, HU^2, U^3, H^2V, HUV, U^2V, HV^2, UV^2

and fix q=q0=HW-U^2-UV-V^2. The previous function/differential calculation
gave

    C1=H^3-5HU^2-7HUV+11HV^2,
    C2=6H^2U+5U^3+3H^2V+4U^2V+t*UV^2.

Allow all nine coefficients of C3, and mark motions

    z_A=1+pi+pi^2*a,
    z_B=-1+pi+pi^2*b.

The function equations modulo pi^5 and the selected effective-fibre
equations force

    t=3,          a=9, b=-9,
    C3=-H^3-10HU^2+6HUV+2HV^2+v*UV^2.

Here v remains free until the next stage. The exact Groebner basis,
with C3 coefficients v0,...,v8 and marks v9,v10, is

    t-3, v0+1, v1, v2+10, v3, v4, v5-6, v6, v7-2,
    v9-9, v10+9.

One last fibre coefficient is quadratic off this locus; it vanishes
on the resulting solution space. It is not incorrectly included in a
linear system. For the other equations the dependence on the C3 and
mark variables is affine, and the t-dependence has degree at most four.
Their coefficient polynomials are reconstructed with five t-values and
each affine basis direction, with three further mixed-point checks.

Independently computing the actual canonical osculating planes over
O/(pi^3) gives the order-two remainder coefficients

    i0=3(a+b),
    i1=16+14t+7(a-b).

Both vanish on the displayed ideal. More explicitly,

    i0=3[(a-9)+(b+9)],
    i1=14(t-3)+7(a-9)-7(b+9).

This is a genuine implication of the necessary local map equations,
after the first-jet Mathieu selection, not an osculation constraint added
to discover a preferred solution.

Certificate: investigate_harmonic_second_incidence_order.py.

## 4. Orders 3 through 8

At stage m there are twelve unknown residues: the remaining UV^2
coefficient of C_m, nine coefficients of C_(m+1), and the two mark
coefficients at order pi^m. Impose the degree-23-function equations
modulo pi^(m+3), then use the selected fibre coefficients from Section 2.

For every m=3,...,8 the resulting affine system has rank eleven, is
consistent, and has a one-dimensional kernel. That kernel changes only
the UV^2 coefficient of C_(m+1). Thus C_m and the order-m mark motions
are determined; the remaining next curve coefficient does not affect
the order-m incidence.

The tests of affineness have a finite degree bound. A change of C_m
first changes beta at order pi^(m+1), because the special linear
function variation vanishes. It changes eta and R at order pi^(m-1).
For m>=5, products of two such R variations have order at least
2m-2>=m+3 and disappear. For m=3,4 the selected residuals have degree
at most two in the twelve variables. The optional full quadratic audit
checks every negative basis vector and every sum of two distinct basis
vectors, in addition to the base and positive basis vectors. This tests
all possible quadratic coefficients, not just one random direction.

The resulting marked sections are

    z_A = 1+pi+9pi^2-7pi^3-4pi^4-6pi^5-2pi^6+7pi^7-6pi^8
          mod pi^9,
    z_B = -1+pi-9pi^2-7pi^3+4pi^4-6pi^5+2pi^6+7pi^7+6pi^8
          mod pi^9.

In the nine-monomial order of Section 3, the determined curve rows are:

| Order | Coefficients of C_m |
| --- | --- |
| 1 | [1,0,-5,0,0,-7,0,11,0] |
| 2 | [0,6,0,5,3,0,4,0,3] |
| 3 | [-1,0,-10,0,0,6,0,2,0] |
| 4 | [0,-7,0,-1,-2,0,10,0,6] |
| 5 | [-7,0,-9,0,0,0,0,7,0] |
| 6 | [0,-8,0,1,7,0,0,0,8] |
| 7 | [-3,0,2,0,0,-4,0,-4,0] |
| 8 | [0,3,0,-9,-8,0,-2,0,-10] |

These are centered integral representatives in a fixed local projective
gauge, not characteristic-zero equations for a completed cover.

At each order the restricted cubic on the line cut out by the two
osculating planes is divided by the restricted quadric. Its remainder
is zero modulo pi^(m+1). The actual plane equations are obtained from
the first three canonical jets at the moving marks; a first-order
formula is NOT reused in place of higher-order osculation.

Certificate: investigate_harmonic_higher_incidence_orders.py.

## 5. The same short row identities appear at every tested order

The correction matrix, including its labelled coefficient rows, is
identical in orders 3 through 8. Let E_(a,j) denote the coefficient of
pi^a*z^j in the scalar-adjusted fibre remainder. Let C_j denote the
coefficient of z^j in the order-pi^(m+2) function obstruction for the
first pair of quintic sections. Its sign and normalization are those
in the multiplier solver: it is minus the remaining cross-product
coefficient after solving the fifteen selected multiplier equations.

The two order-m osculating coefficients satisfy, on the full affine
space of extensions of the computed lower jet,

    i0 = 5E_(m+1,9) - E_(m+1,11) + 4E_(m+1,13) + 6C_7,

    i1 = 20E_(m,15) + 18E_(m+1,10) + 16E_(m+1,12)
         +13E_(m+1,14) + 14C_8 + 14C_10.

All these equalities are over F23. They are checked as affine identities,
not merely evaluated at the chosen solution. If E=M*x+b and
i=B*x+i_base, the script verifies BOTH

    W*M=B,           W*b=i_base.

The second equality matters. The row-space inclusion alone would prove
that incidence has a constant value on the solution space, not that
this constant is zero.

## 6. The exact remaining induction problem

The displayed row identities are proved for the computed lower jets in
orders 3 through 8. They have not been proved for every higher order.
In particular, an unchanged Jacobian is not proof that its nonlinear
constant terms continue to satisfy the required relation.

A sharply formulated next target is:

> For every m>=3, assuming the preceding function and effective-fibre
> equations and incidence through order m-1, prove the two affine
> identities of Section 5 for arbitrary extensions to order m.

An identity of this form, allowing coefficient functions with the same
reductions and terms in the earlier vanishing ideals, would complete
the induction. The first-order Mathieu selection and the separately
verified second-order calculation would supply the initial conditions.
Separatedness of the complete DVR would then give exact incidence.

This is now a problem about explicit function, divisor, and canonical
jet equations. It does not require introducing a hypothetical
cohomological correspondence. But it remains a mathematical problem:
eight successful orders cannot stand in for the universal identity.

## Checks and reproduction

All arithmetic keeps pi^2=-23. Quotients and inverses use only units at
23. Helpers are in harmonic_truncated_geometry.py; the existing mixed
function engine now also exposes its affine ratio and divided
differential at the retained precision. Off the global function-lifting
locus those outputs are labelled candidates and are used only to form
equations. Final solutions are re-evaluated with every cross product.

The full quadratic audits passed in orders 3 and 4 (78 additional
parameter points each). An independent helper check compares all 32
first-order directions with the older osculating calculation, checks
Weierstrass preparation in precisions 2 through 6, and detects the
nonzero mixed-characteristic base obstruction 19*pi^2*z. This last test
would fail if the arithmetic accidentally set 23 equal to zero beyond
the first infinitesimal order. The earlier pure-H^3 obstruction and
third-order function tests also still pass after the engine changes.

    env DOT_SAGE=/private/tmp/m23-cover-investigation-sage sage -python notes/investigate_harmonic_second_incidence_order.py --solve-parameter --workers 6
    env DOT_SAGE=/private/tmp/m23-cover-investigation-sage sage -python notes/investigate_harmonic_higher_incidence_orders.py --max-order 8 --workers 6
    env DOT_SAGE=/private/tmp/m23-cover-investigation-sage sage -python notes/investigate_harmonic_higher_incidence_orders.py --max-order 4 --workers 6 --audit-quadratic
    env DOT_SAGE=/private/tmp/m23-cover-investigation-sage sage -python notes/certify_harmonic_truncated_geometry.py
