# A fixed Hensel system for the harmonic lift

Private investigation, 6 September 2026. No manuscript, public repository,
or website changes. No reconstructed characteristic-zero cover was used
in this calculation.

Follow-up (same date):
HARMONIC_COMPLETE_LOCAL_THREE_POINT_CONSTRUCTION_2026_09_06.md proves
that the selected function and fibre equations suffice for the complete
three-point map, and proves generic smoothness. It supersedes the caution
below about omitted MAP equations. It does not remove the osculating
incidence gap or supply global Galois fixedness from local uniqueness.

Further follow-up (same date):
HARMONIC_EXACT_INCIDENCE_FROM_LOCAL_RECONSTRUCTION_2026_09_06.md now
proves the two local incidence constants zero by independently
constructing and exactly verifying a global map in this neighbourhood.
The proof uses the necessary-equation uniqueness established here; it
does not infer exact vanishing from the earlier finite precision tests.
The reduction/model hypotheses for a branch-cycle-only application
remain explicit. The historical account below is retained.

## What is proved, and what is not

The varying-order calculation in
HARMONIC_HIGHER_INCIDENCE_AND_FILTERED_IDENTITIES_2026_09_06.md can be
replaced by ELEVEN FIXED analytic equations in ELEVEN variables. After
the earlier first- and second-jet restrictions, their reduction is an
affine linear system with determinant 16 in F23.

Consequently these necessary equations have exactly one integral solution
to ALL orders in this normalized neighbourhood. Any actual three-point
cover satisfying the specified model and jet hypotheses must be that
solution. This also proves, at every higher order, the stability of the
linear correction operator for these equations; it is no longer just
an observation in orders 3 through 8.

This does NOT prove that the osculating incidence persists. Two values
at the unique solution still have to be shown to be exactly zero. The
previous eight-order calculation only places them in a high power of
the maximal ideal. Existence of a solution of our eleven selected
equations alone is also NOT asserted to construct a complete Belyi map:
there are additional function and branch-fibre equations.

The conditional identification with a genuine cover uses existence of a
cover in this neighbourhood. Neither such existence nor the compatible
integral-model hypotheses are supplied by an invertible Jacobian.

## 1. The normalized neighbourhood

Let O be the completion of Z23[pi], with pi^2=-23. The construction works
over unramified complete extensions as well. Use

    q0=HW-U^2-UV-V^2,       c0=UW^2-V^3,

and the cubic monomials, in this order,

    m0,...,m8 = H^3,H^2U,HU^2,U^3,H^2V,HUV,U^2V,HV^2,UV^2.

Set

    C1 = H^3-5HU^2-7HUV+11HV^2,
    C2 = 6H^2U+5U^3+3H^2V+4U^2V+3UV^2,
    C3o = -H^3-10HU^2+6HUV+2HV^2.

For eleven integral variables x0,...,x10, take the curve

    q=q0,
    c=c0+pi*C1+pi^2*C2+pi^3*C3o
         +pi^3*x0*UV^2+pi^4*sum(x_(i+1)*m_i, i=0,...,7),

with marks of z=V/W coordinates

    a=1+pi+9pi^2+pi^3*x9,
    b=-1+pi-9pi^2+pi^3*x10.

All further curve and mark coefficients are included here. The apparently
free next UV^2 coefficient in the earlier twelve-variable calculations
is simply the next digit of x0; it is not an additional formal parameter.

The good C1 was selected by the common-sheet Mathieu/Witt-design
calculation under its stated compatible-tail hypotheses. The preceding
function and effective-fibre calculations give C2, C3o and the displayed
mark congruences. They remain inputs to the present theorem, not
consequences of the Hensel argument itself.

### Check that these are actual curve coordinates

The special quadric is nonsingular and 2 is invertible. Formal changes
of projective coordinates trivialize its deformations. Once q=q0 is
fixed, the infinitesimal orthogonal group has dimension six. Its action
on c0, the cubic rescaling direction c0, and the four directions q0*H,
q0*U,q0*V,q0*W span an eleven-dimensional space of cubics. Adding
m0,...,m8 gives all twenty cubic monomials, with nonzero determinant
(10 with the basis ordering in the certificate).

The formal inverse-function theorem therefore supplies a slice with
these nine cubic coordinates, after choosing a projective trivialization
reducing to the identity. This is a local statement, not a claim about
global coordinates on the moduli space. A formal parametrization of the
orthogonal group may be obtained by the Cayley transform; its linear
term is invertible because 2 is a unit.

## 2. Define the candidate function without assuming that it glues

On W=1, put V=z. There are unique restricted power series H(z),U(z)
reducing to z^6+z^4+z^2,z^3 and satisfying q=c=0. The implicit Jacobian
in H,U is a unit on this central chart. These series can be evaluated
at a,b, which are integral and reduce to 1,-1.

Use the fixed 27 canonical quintic monomials in the earlier certificate.
The order-23 jet evaluation map at either mark has rank 23 modulo pi.
Choose the same four-dimensional kernel bases S_A and S_B by the fixed
invertible 23-by-23 minor. Their coefficients are integral analytic
functions of x.

Let D=M*S_B. After fixing one scalar coordinate of M, impose the same
fifteen independent coefficient equations from

    S_A,i D_j - S_A,j D_i.

Their matrix in the fifteen remaining entries of M is invertible modulo
pi. Thus they define M analytically for EVERY x in this neighbourhood,
including points for which other cross-products do not vanish. A ratio
of a chosen pair, normalized to value one at z=0, defines an analytic
LOCAL candidate beta=N/D. It must not be called a globally defined
degree-23 function unless the remaining relations vanish.

Write

    C_j = [z^j] (-(S_A,0 D_1-S_A,1 D_0)).

These are actual integral analytic functions, with a fixed sign and
normalization, not successive obstruction digits chosen afresh.

The diagnostic mode in the new calculation solves only the fifteen
selected rows at each digit. It never divides a previously nonvanishing
unselected row by a higher power of pi. The original global-lifting
test retains its stricter behaviour.

## 3. Retain the effective ramification factor

The differential

    eta=(d beta/beta)/23

is integral on this neighbourhood and reduces to

    2*z^8/(z^2-1) dz.

Here this is initially the differential of the local candidate only.
Let R(z), monic of degree eight and congruent to z^8 modulo pi, be its
Weierstrass factor. This factor depends analytically on x. For the
actual cover it represents the effective eight-point ramification
divisor on the central chart.

Divide N-D and D by R^2, obtaining remainders P,Q. Since Q(0) is a unit,
set

    E(z)=P(z)-P(0)/Q(0)*Q(z),       E_j=[z^j]E(z).

This eliminates the unknown third branch value. A genuine three-point
cover has E=0. As before, neither the effective square R^2 nor the
scalar branch-value adjustment can be omitted.

## 4. Eleven fixed equations

Consider, in the following order,

    F = (E15/pi^3,
         E9/pi^4,E10/pi^4,E11/pi^4,E12/pi^4,E13/pi^4,E14/pi^4,
         C7/pi^5,C8/pi^5,C9/pi^5,C10/pi^5).

All these quotients are integral analytic functions on the specified
neighbourhood. They are divided equations, not assertions that pi is
a unit in O.

Their reduction is

    F(x) = J*x+b mod pi,

with

    b = [0,3,0,5,0,11,0,21,0,20,0]

and

    J =
    [20,13, 0, 9, 0, 0, 5, 0, 4, 5,18]
    [ 0, 0,12, 0, 3, 7, 0,17, 0,15,15]
    [18,17, 0,14, 0, 0, 2, 0,14,14, 9]
    [ 0, 0,11, 0,21,15, 0,22, 0, 2, 2]
    [ 7,10, 0,20, 0, 0,12, 0, 4,22, 1]
    [ 0, 0,18, 0,19,13, 0, 6, 0, 6, 6]
    [14,16, 0, 9, 0, 0,19, 0,14, 3,20]
    [ 0, 0,22, 0,19,11, 0,20, 0,15,15]
    [19, 6, 0, 8, 0, 0,13, 0, 8,17, 6]
    [ 0, 0,20, 0, 7,13, 0,17, 0,16,16]
    [21,16, 0,17, 0, 0,10, 0, 0, 6,17].

In F23,

    det(J)=16,
    -J^(-1)b = [0,0,-7,0,-1,-2,0,10,0,-7,-7].

### Why the finite calculation certifies this reduction

It is important to distinguish a finite proof of this Jacobian statement
from finite evidence for the desired incidence theorem.

The curve and mark variables first enter at orders pi^3 or pi^4. The
first variation of the normalized function at the special curve is
zero; therefore a variable first affects beta at order at least pi^4,
and eta and R at order at least pi^2. Quadratic contributions to beta
start no earlier than pi^6; after division by pi^2 they start at pi^4
in eta. Cubic contributions to R start at pi^6. Consequently all
expressions needed modulo pi^6 have degree at most two in x, including
the lower-order coefficients whose vanishing proves the integrality of
the displayed quotients.

For eleven variables, values at zero, the eleven positive and eleven
negative coordinate vectors, and the 55 pairwise sums of distinct
coordinate vectors determine every polynomial of degree at most two.
The certificate checks all 78 values. All quadratic coefficients of
the normalized reductions vanish. Integrality of every division is
checked in these tests as well. Since this is a polynomial identity
test with a degree bound, it remains valid after extending F23; it is
not a search through a few possible residue values.

Weierstrass preparation is carried out modulo pi^4 from a sufficiently
long z series. The earlier z^15/z^8 precision lemma ensures that this
determines exactly E15 modulo pi^4 and E9,...,E14 modulo pi^5, as needed.
The checker perturbs both unavailable digits of R and verifies that
the retained normalized values do not change.

## 5. An all-orders theorem

**Proposition.** In the normalized neighbourhood of Section 1, the
equations F=0 have a unique integral solution x*. Formation of x*
commutes with complete extensions of O for which the displayed
congruences hold. In particular, any genuine cover satisfying those
congruences has the same canonical coefficients and marked points.

**Proof.** The constructions above are integral analytic: implicit
equations, the jet kernels, and M use only inverses of units. Division
by the indicated powers of pi is integral by Section 4. Weierstrass
preparation gives convergent coefficients in the pi-adically complete
parameter algebra. We may work with restricted power series over O.

Choose integral lifts Jt,bt of J,b. Then

    F(x)=Jt*x+bt+pi*G(x),

where G has integral restricted power-series coefficients. The map

    x |-> -Jt^(-1)bt-pi*Jt^(-1)G(x)

is a strict pi-adic contraction on integral points. It has a unique
fixed point. Alternatively, apply the multivariable form of Hensel's
lemma to the unique residual zero and the unit Jacobian determinant.
The same contraction proves uniqueness after complete base extension.
Every genuine cover obeys the eleven necessary equations, so its
coordinate vector equals x*. QED.

This is the ordinary Hensel/implicit-function argument, not a new
general lifting principle. See also the Stacks Project discussion of
complete local and henselian rings:
https://stacks.math.columbia.edu/tag/03QG and
https://stacks.math.columbia.edu/tag/04GE.

The map x |-> F(x) is itself an analytic change of coordinates on this
unit polydisc, by the same contraction with an arbitrary target y.
Thus the completed algebra defined by these SELECTED equations is O.
This statement does not identify it with the full Hurwitz algebra
without an existence input. If an actual cover in this neighbourhood
is known to exist, all its further necessary equations vanish at x*
by uniqueness, not because they were silently discarded.

Nor does uniqueness in this 23-adic neighbourhood imply uniqueness under
the GLOBAL arithmetic Galois action: a conjugate can lie in a different
neighbourhood. The counterexample in
LOCAL_PATTERN_DOES_NOT_FORCE_GLOBAL_FIXEDNESS_2026_09_06.md remains
relevant. Exact osculating incidence would supply the missing global,
place-independent condition; local rigidity by itself does not.

### Consequence for every higher correction step

For n>=1 and integral h,

    F(x+pi^n h)-F(x) = pi^n J h mod pi^(n+1).

This proves the uniform correction matrix for the fixed equations in
ALL orders. It does not identify the value of the incidence at their
zero. A constant term can remain even when this formula holds.

## 6. Put the remaining incidence problem in exact form

Compute the two genuine osculating planes and divide the cubic on their
intersection line by the quadric there. Let r0+r1*z be the remainder.
The preceding first- and second-order results imply that

    g=(r0/pi^3,r1/pi^3)

is integral analytic throughout the chosen neighbourhood. Let

    delta=(delta0,delta1)=g(x*).

The full fixed-mark incidence statement is now exactly

    delta0=delta1=0.

Since F gives analytic coordinates, there exist integral analytic
coefficient functions H_ij such that

    g_i(x)=delta_i+sum_j H_ij(x) F_j(x).

For example, express g in the coordinates y=F(x) and group every
nonconstant power-series monomial by one of its y factors. No division
by a possibly nonunit integer is involved. This gives an exact
decomposition, but does NOT prove that its constant delta_i is zero.
Equivalently, the missing assertion is that g0,g1 belong to the closed
ideal (F1,...,F11), not merely that their derivatives lie in its span.

The constant matrix of the earlier row witnesses is the reduction of
H:

    W =
    [ 0,5, 0,-1, 0,4, 0,6, 0,0, 0]
    [20,0,18, 0,16,0,13,0,14,0,14].

It satisfies both W*J=(the reduction of Dg) and W*b=g(0) mod pi.
The eight-order incidence computation gives delta in pi^6*O^2. It does
not show delta=0.

### A concrete warning against the constant-coefficient shortcut

Take the integral lift

    x=[0,0,16,0,22,21,0,10,0,16,16]+pi*[1,0,0,0,0,0,0,0,0,0,0].

This point has the correct residual solution but is not an actual
global function lift. The fixed candidate is nevertheless well defined.
Direct computation gives

    g(x)-W*F(x) = (5*pi^2,0) mod pi^3.

The stored centered representative is (-115,0), since pi^2=-23.
Thus the constant W is NOT an exact analytic syzygy on this
neighbourhood. This does not refute the proposed filtered induction:
correction terms involving earlier equations may remove this defect.
It identifies why higher coefficient functions, or such earlier-ideal
terms, cannot simply be omitted.

## 7. A finite algebraic presentation retaining the double fibre

There is a useful exact replacement for the long central power series
in the branch-fibre equations.

**Lemma.** Let A be pi-adically complete, let R be monic of degree eight
with R=z^8 mod pi, and suppose

    f(U,z)=U-z^3+pi*t(U,z)

is the affine curve equation obtained by eliminating H. In the finite
free rank-sixteen algebra B=A[z]/(R^2), there is a unique element Ubar
congruent to z^3 modulo pi with f(Ubar,z)=0. Put

    Hbar=Ubar^2+Ubar*z+z^2.

For the chosen homogeneous numerator and denominator sections N,D,
the local double-fibre condition is exactly

    N(Hbar,Ubar,z,1)-u*D(Hbar,Ubar,z,1)=0 in B.

**Proof.** The map U |-> z^3-pi*t(U,z) is a contraction on B. Since
B is finite free over A, it is pi-adically complete, and the unique
fixed point gives Ubar. Equivalently, the sixteen coefficient equations
f(Ubar,z)=0 modulo R^2 have derivative the identity modulo pi in the
sixteen coefficients of Ubar. Reduction of the implicit central series
U(z) into B satisfies the same equation and the same congruence, so it
equals Ubar. The denominator is a unit in B: modulo pi it has nonzero
constant coefficient in A/(pi)[z]/(z^16). Thus vanishing of the displayed
numerator is equivalent to beta=u on the length-sixteen subscheme
defined by R^2. This proves the effective, scheme-theoretic condition,
not just its values at the reduced support. QED.

In particular, one may introduce eight coefficients for R, sixteen for
Ubar, and one scalar u, with the TWO sets of sixteen polynomial remainder
equations

    f(Ubar,z)=0 mod R^2,
    N(Hbar,Ubar,z,1)-u*D(Hbar,Ubar,z,1)=0 mod R^2.

This is an actual finite presentation of the double-fibre condition.
The new certificate independently computes Ubar by contraction inside
B, evaluates the homogeneous quintics there, and compares the result
with the earlier long-series/remainder calculation at every sampled
parameter value, including an off-global-function value. They agree.

## 8. What this changes about the next calculation

There is no longer a need to prove the stability of an ever-changing
Jacobian. That part is supplied by a fixed analytic system and its
verified reduction. The useful next target is an exact ideal-membership
certificate for the two incidence expressions, or an independent exact
construction of a zero of F with incidence. Repeating the same linear
calculation at many more orders cannot replace this check.

Section 7 constructs the finite branch-fibre presentation. The canonical
jet kernels and global cross-product conditions are also finite
algebraic equations. Clearing only specified unit denominators, and
selecting the harmonic local component after the required pi
saturations, gives an exact algebraic ideal-membership problem.

This is a formulation of the next computation, NOT a claim that its
elimination or saturation has already been carried out.

## Reproduction and regression checks

    env DOT_SAGE=/private/tmp/m23-cover-investigation-sage sage -python \
      notes/certify_harmonic_fixed_hensel_system.py \
      --workers 6 --audit-quadratic --check-off-locus

This checks the slice, all normalized residual coefficients, the full
quadratic audit, the determinant, the unique residual solution, both
residual incidence identities, the finite double-fibre presentation,
and the explicit off-locus defect.

A separate run with --lift-digits 4 used only the SAME ELEVEN equations
at each step. The unselected function equations and incidence were then
checked independently at the resulting finite precision; they vanished.
This is a regression against the earlier varying-order calculation,
not the proof of its infinite-order incidence extension.

The pure-H^3 unit-obstruction regression still passes after adding the
analytic-candidate mode. Default failed lifts continue to be rejected.
