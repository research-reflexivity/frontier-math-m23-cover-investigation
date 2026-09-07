# Exact harmonic incidence from an independent local reconstruction

Private research, 6 September 2026. No manuscript, public repository,
or website changes. This calculation did not read the previously stored
characteristic-zero models, maps, or their osculating certificates.

Follow-up audit, same date: the private manuscript now incorporates the
conditional local theorem. The exact map, descent, local comparison and
full 78-point residual audit were rerun successfully. The separate
`verification/verify_harmonic_model_identification.py` gives an exact
isomorphism to the earlier certified map and thereby certifies M23
monodromy. That comparison deliberately reads the old data; the
independent construction below still does not. The historical monodromy
qualification below applies to the construction certificate in isolation.
See `HARMONIC_RECONSTRUCTION.md` for current scope and Magma run status.

## Status: the two local incidence constants are now zero

The fixed-mark incidence gap in the specified harmonic neighbourhood is
resolved. An intrinsic coordinate frame makes it possible to recognize a
small exact model from the independently computed local lift. Exact
calculations then prove that this model is smooth and carries the required
degree-23 three-point map. It belongs to the same normalized local
neighbourhood. The previously proved Hensel uniqueness therefore identifies
it with the local lift. A short exact factorization proves incidence.

This is an arithmetic proof using explicit equations and finite exact
certificates. It is NOT an equation-free geometric lifting principle.
Rational reconstruction supplied a candidate, not a proof; its subsequent
exact verification is essential.

In particular, the constants delta in
HARMONIC_FIXED_HENSEL_SYSTEM_AND_INCIDENCE_IDEAL_2026_09_06.md are exactly
zero, not merely divisible by the tested power of pi.

An additional exact semilinear symmetry gives descent of this map to Q.
In a rational target coordinate its branch locus is t=0 together with
23*t^2+1=0. Thus the new construction supplies a global map, not only a
23-adic object. The monodromy qualification below remains necessary.

There is still a distinction between this local theorem and a theorem
beginning with abstract Nielsen data alone. To apply it to that problem,
one must justify the harmonic integral model, its compatible component
models, and the choice of first jet by the common-sheet Witt-design
argument. Those hypotheses have NOT been silently discharged. The exact
ramification passport certified here does not by itself certify M23
monodromy.

## 1. The local theorem and its hypotheses

Put K0=Q(pi), pi^2=-23, and let O be the completion of its ring of integers
at the prime above 23.
The theorem also holds after a complete extension with the same stipulated
congruences. Use canonical coordinates [H:U:V:W], with

    q0=HW-U^2-UV-V^2,       c0=UW^2-V^3.

Let m0,...,m8 be

    H^3,H^2U,HU^2,U^3,H^2V,HUV,U^2V,HV^2,UV^2.

Retain the three rows C1,C2,C3o from the fixed Hensel construction:

    C1 = H^3-5HU^2-7HUV+11HV^2,
    C2 = 6H^2U+5U^3+3H^2V+4U^2V+3UV^2,
    C3o = -H^3-10HU^2+6HUV+2HV^2.

**Local incidence theorem.** Suppose a smooth degree-23 three-point map
with ramification types (23),(23),(2^8,1^7) has a canonical model and
ordered totally ramified sections A,B in the neighbourhood

    q=q0,
    c=c0+pi*C1+pi^2*C2+pi^3*C3o
         +pi^3*x0*UV^2+pi^4*sum(x_(i+1)*m_i, i=0,...,7),
    z_A=1+pi+9pi^2+pi^3*x9,
    z_B=-1+pi-9pi^2+pi^3*x10,

where all x_i are integral. Assume the normalized local map and its
effective critical fibre satisfy the hypotheses used to define the
eleven necessary Hensel equations: in particular, the primitive divided
logarithmic differential reduces to 2z^8/(z^2-1) dz, and its degree-eight
zero divisor is the critical cluster under consideration.

Then the order-three osculating planes at A and B have a common divisor
of degree two, consisting of two distinct geometric points. Equivalently,
the two osculating remainders r0,r1 vanish exactly.

The proof below uses NECESSITY and UNIQUENESS of the fixed eleven-equation
system, together with an independent exact existence certificate. It does
not require the subsequent Picard-functor argument establishing sufficiency
of the eleven equations. This separates the new proof from a possible
concern about redundant global equations.

## 2. Remove the coordinate ambiguity using two canonical flags

For i=0,1,2,3, take a generator s_i of

    H0(X,omega_X(-i A-(3-i)B)).

In this neighbourhood these four spaces have dimension one, and their
generators form a basis of H0(X,omega_X). Both assertions hold by unit
minors on the special fibre. Thus their projective coordinate lines are
intrinsic to the ordered marked curve, even though an initial generator
of each line is arbitrary.

In coordinates [x0:x1:x2:x3]=[s0:s1:s2:s3], the marks are

    A=[1:0:0:0],       B=[0:0:0:1].

Their vanishing sequences are (0,1,2,3) and (3,2,1,0), respectively.
Consequently the canonical quadric has the form

    a*x0*x3+b*x0*x2+d*x1*x2+e*x1^2+f*x1*x3+g*x2^2.

The coefficients a,b,d,e are units in the local neighbourhood. Rescale
the coordinates by x_i=d_i*y_i, fixing d0=1, and multiply the equation
by kappa. The unique normalization setting the first four displayed
coefficients to one is

    d1=b/d,        kappa=1/(e*d1^2),
    d2=1/(kappa*b),        d3=1/(kappa*a).

No square root is taken. These normalizations have determinant one on
the character lattice; they eliminate the remaining diagonal projective
ambiguity over the field of definition of the marked curve.

Finally reduce the cubic modulo the monic quadric in lexicographic order
x0>x1>x2>x3 and normalize its x0^2*x3 coefficient to one, another unit.
This gives an intrinsic normalized pair of equations. In particular,
unnecessary algebraic coordinate changes cannot disguise their field
of definition, as happened in the original fixed-quadric slice.

The finite local calculation now recognizes small coefficients. The two
remaining quadric coefficients are -127/229 and 12192/229^2. Recognition
of all cubic coefficients succeeds at precision pi^33. Nothing in the
proof relies on assuming that a sufficiently convincing approximation is
an exact value.

## 3. The exact candidate and its incidence

The recognized quadric is

    q=x0*x2+x0*x3+x1^2+x1*x2
       -(127/229)*x1*x3+(12192/52441)*x2^2.

The cubic coefficients, with all omitted monomials zero, are:

| Monomial | Coefficient |
| --- | --- |
| x0^2*x3 | 1 |
| x0*x1*x3 | (69+41*pi)/916 |
| x0*x3^2 | -(288+3552*pi)/52441 |
| x1^3 | -(5+pi)/458 |
| x1^2*x2 | -(3061+569*pi)/209764 |
| x1^2*x3 | -(8649+1515*pi)/104882 |
| x1*x2^2 | -(298479+49803*pi)/48035956 |
| x1*x2*x3 | -(1041970+181274*pi)/12008989 |
| x1*x3^2 | (2318419+389759*pi)/48035956 |
| x2^3 | -(2350008+338328*pi)/2750058481 |
| x2^2*x3 | -(55642056+9354216*pi)/2750058481 |

These are elements of K0, not just its completion. The useful denominator
identities are 52441=229^2, 48035956=4*229^3 and 2750058481=229^4.

The exact Jacobian-minor ideal is the unit ideal on each of the four
affine projective charts. Thus this is a smooth canonical (2,3) complete
intersection, of genus four. The determinant of twice the quadric Gram
matrix is

    77968/52441 = 4873*(4/229)^2,

so the quadric is nonsingular as well. The marked jet computations give
the exact vanishing sequences stated above. In particular, the two
osculating planes are x3=0 at A and x0=0 at B.

On their common line x0=x3=0 one has the polynomial identity

    c = (x1^2+x1*x2+(12192/52441)*x2^2)
         * (-(5+pi)/458*x1-(771+111*pi)/209764*x2).

The first factor is exactly the restricted quadric. Its discriminant is

    3673/52441 != 0.

It therefore cuts out precisely two reduced geometric points, neither
of which is A or B. This is the exact incidence certificate.

## 4. Verify the degree-23 map, not just the source curve

All the following calculations take place over K0 with exact arithmetic.

1. Compute the 27-dimensional degree-five canonical-ring piece. At A and
   B, take the kernels of the 23-jet evaluation matrices. Both kernels
   have dimension four. Denote their quintic section bases by S_A,S_B.

2. Check residual basepoint-freeness. The common ideal of q,c and the
   four sections, saturated by the corresponding marked-point ideal,
   is the unit ideal. At the mark at least one section has vanishing
   order exactly 23. Thus the common divisor of the four sections is
   exactly 23A or 23B, not that divisor plus an unobserved base point.

3. Reduce products in the 57-dimensional degree-ten canonical ring and
   solve for a multiplier matrix M. The coefficient matrix has rank 15
   and one-dimensional kernel; the resulting M is invertible. Verify
   ALL SIX cross-product identities

       S_A,i*(M*S_B)_j-S_A,j*(M*S_B)_i = 0 mod (q,c).

   A numerator/denominator pair N,D therefore gives a global function
   beta with div(beta)=23A-23B, hence degree exactly 23.

4. In the chart x1=1, compute the differential critical ideal using the
   determinant with rows grad(q), grad(c), and D*grad(N)-N*grad(D).
   Saturate away from N*D=0. The resulting quotient algebra has dimension
   eight. Its multiplication operator by x2 has squarefree characteristic
   polynomial of degree eight, proving that these are eight distinct
   geometric critical points.

5. In that algebra, check N=u*D for a nonzero scalar u in K0. All eight
   critical points therefore have the same value. The exact scalar and
   degree-eight polynomial are retained in the generated certificate.

Riemann--Hurwitz gives total ramification 2*4-2+2*23=52. The points A and B
contribute 44. The eight distinct critical points already contribute at
least eight, so each is simple and there can be no others, including in
an omitted affine chart. After dividing beta by u, its branch values are
0,1,infinity and its passport is (23),(23),(2^8,1^7).

This establishes existence of the required exact three-point map without
assuming the existence or equations of the formerly stored generic model.

## 5. Why this proves incidence for the local lift

The intrinsic framed equations of the local Hensel approximation agree
with the exact candidate modulo pi^33, coefficient by coefficient.
The comparison also checks the actual marked coordinate axes.

This high precision is much more than is needed for the logical argument.
An integral lift of the inverse frame gives a projective coordinate change
with unit determinant. Pulling back the exact candidate gives a canonical
deformation congruent to the specified harmonic model and marked jets.
The formal slice theorem of the fixed-Hensel note puts it into q=q0 and
the nine-monomial cubic slice; since all discrepancies are of higher
order, that correction preserves the stipulated low-order jets.

The global degree-23 function is determined up to scalar by its divisor.
Its normalized section construction and divided logarithmic differential
are therefore the same analytic constructions used in the Hensel system.
The exact map satisfies all eleven NECESSARY equations. Their unit
Jacobian and unique residual root imply that its slice parameters equal
the unique Hensel solution, to all orders.

Incidence is invariant under projective coordinate changes. The exact
factorization in Section 3 consequently proves r0=r1=0 on that solution.
Equivalently, the two constants delta in the local incidence ideal are
zero. This implication uses exact existence plus uniqueness, not a
precision or height extrapolation.

There is also an exact ideal-theoretic consequence. In the completed
local parameter ring, the eleven functions F are analytic coordinates.
The divided incidence functions g=(r0/pi^3,r1/pi^3) vanish at their
unique zero. Consequently

    (g0,g1) is contained in (F1,...,F11).

Indeed, write g in the F-coordinates. Its constant term is zero, so each
component lies in the ideal of those coordinates. This proves the
all-orders incidence identity whose constant terms had been missing;
it does not assert that the earlier constant multiplier matrix alone
was an exact identity, or provide short closed-form higher multipliers.

## 6. Explicit descent of the constructed map to Q

Let sigma send pi to -pi. Write a=-127/229 and b=12192/52441, and define
the rational projective transformation

    T(x0,x1,x2,x3)=(a*b*x3, b*x2, x1, x0/a).

There are exact identities

    T^2=b*Id,              q(Tx)=b*q(x),
    sigma(c)(Tx)=mu*c(x) mod q,
    mu=-(52130760+33951768*pi)/1525141603,
    mu*sigma(mu)=b^3.

The transformation exchanges the two marked coordinate axes. If
f=N/(u*D) is the map with third branch value one, the exact canonical-ring
calculation also gives

    f(x)*sigma(f)(Tx)=1.

This is a polynomial reduction identity after clearing denominators,
not an inference from numerical symmetry.

Set R=(b/mu)*T. Then R*sigma(R)=Id. The checker solves v=R*sigma(v) as
eight rational linear equations on the rational and pi coordinates of v.
Its solution space has dimension four over Q and supplies an invertible
matrix B of invariant columns. This constructs the descent explicitly,
without an assumed rational point on an unspecified twisted variety.

Substitute x=B*y. Normalize q(B*y); its coefficients are rational.
Reduce c(B*y) modulo that rational quadric and normalize its scalar;
the resulting cubic also has rational coefficients. These assertions
are checked coefficient by coefficient.

The rational target coordinate is

    t=(f-1)/(pi*(f+1)).

It is invariant under the semilinear descent because sigma exchanges f
with 1/f and sends pi to -pi. The certificate constructs rational
numerator and denominator forms for t: multiply its denominator by its
conjugate, and reduce the corresponding numerator in the descended
canonical ring. The coefficients of both reductions are checked to
belong to Q.

The branch values f=0,1,infinity become -1/pi,0,1/pi. Hence the descended
degree-23 map is defined over Q with branch locus

    t=0,                 23*t^2+1=0.

Its two order-23 branch points are conjugate, not individually rational.
This proves global arithmetic descent of this map. It does not identify
its monodromy group from its passport.

## 7. Scope for the global fixed-point question

There are now two separate verified routes to the same geometric incidence:
the earlier test on the seven reconstructed generic covers, and this
independent local reconstruction with exact existence and uniqueness.
The latter connects the stipulated harmonic deformation to the generic
incidence without using the former model as an input.

For a branch-cycle-only theorem, the remaining task is to establish all
the reduction and compatible-model hypotheses needed to place the actual
M23 harmonic lift in this neighbourhood. The common-sheet Witt-design
argument explains why M23 selects the positive first jet, subject to its
explicit geometric transport hypotheses. Those must still be audited.

Together with the previously proved empty special incidence for the other
six reduction types, the present theorem would give a place-independent
geometric separator under those hypotheses. Its uniqueness is then a
Galois-invariant condition. Local uniqueness alone is still not the
reason for global fixedness.

Nor is this a newly discovered cover: Huang, Jackson, Lee, Poonen, Pries,
and Zhang constructed the distinguished cover before this investigation.
The present contribution is an independent reconstruction from the
harmonic local data and the resulting exact implication for incidence.

## 8. Reproduction and records

Run from the private repository root:

    env DOT_SAGE=/private/tmp/m23-cover-investigation-sage sage -python \
      notes/certify_harmonic_recognized_model.py --ramification --descent \
      --frame notes/harmonic_intrinsic_frame_residues.json \
      --checkpoint notes/harmonic_local_recognition_checkpoint.json \
      --output notes/harmonic_independent_exact_model_certificate.json

The recorded intrinsic residues can be regenerated from
harmonic_local_recognition_checkpoint.json using
investigate_harmonic_intrinsic_frame.py. The checkpoint itself is a finite
local approximation and is explicitly labelled as such. Its digits can
be recomputed from the harmonic jets by
investigate_harmonic_exact_recognition.py; the exact-model certificate
does not assume their rational reconstructions are correct.

The exact certificate is a SageMath computation. No independent Magma
repetition of these NEW calculations or fresh branch-cycle continuation
for this NEW coordinate model is claimed.
