#!/usr/bin/env sage-python
"""Certify the weighted ADE target-gluing marker at p=23.

The coarse E8 tails become isomorphic after arbitrary target scaling.  A
stable weighted chart retains more: if pi=rho^m and the residual target has
rho-weight n, changing rho changes the target only by mu_m (because
gcd(m,n)=1 in the three cases below).  We compare Frobenius-conjugate tail
models and record whether their target transport is trivial modulo this
allowed root ambiguity.

The resulting values on the rational E8, unramified E8, and ramified A2+A6
closed points are 0,1,1.  This constructs the same Boolean function as the
singular-position idempotent on the reduced local Hurwitz scheme.  The
separate pinched nearby-cycle connector identifies that function with the
M23 normalizer-trace augmentation.
"""

from math import gcd

from sage.all import GF, PolynomialRing


F23 = GF(23)


def nontrivial_mod_roots_of_unity(value, root_degree):
    """Test membership in mu_m by the exact equation value^m=1."""

    return value**root_degree != 1


def norm_one_parity(value):
    """The unique quadratic character of the order-24 norm-one group."""

    assert value * value**23 == 1
    sign = value**12
    assert sign in (value.parent().one(), -value.parent().one())
    return int(sign == -value.parent().one())


def certify_hilbert90_orientation(field):
    """Identify U/U^2 with the canonical semilinear return parity.

    For tau in U choose a with a^(p-1)=tau.  The 22 choices differ by
    F_23^*, and every element of F_23^* is a square in F_(23^2).  Hence the
    square class of a is well defined.  More intrinsically, the semilinear
    map S_a(v)=a*v^23 has return S_a^2(v)=N(a)*v.  Replacing a by a base
    scalar b changes N(a) by b^2, so the Legendre symbol of the affine
    F_23 multiplier N(a) is independent of every Hilbert--90 choice and is
    exactly the norm-one character of tau.
    """

    nonzero = list(field)[1:]
    base_units = [field(a) for a in F23 if a]
    assert all(a.is_square() for a in base_units)
    # As permutations of the affine F_23 point set, multiplication by a
    # base unit has sign equal to its Legendre symbol.  Thus the same gauge
    # group contains both point-set orientations, even though every one of
    # its elements is a square after extension to F_(23^2).
    frame_orientation_changes = {
        int(a**((23 - 1) // 2) == -field.one()) for a in base_units
    }
    assert frame_orientation_changes == {0, 1}
    norm_one = [tau for tau in nonzero if tau * tau**23 == 1]
    assert len(norm_one) == 24
    for tau in norm_one:
        lifts = [a for a in nonzero if a**22 == tau]
        assert len(lifts) == 22
        lift_parities = {int(not a.is_square()) for a in lifts}
        assert lift_parities == {norm_one_parity(tau)}
        return_multipliers = {a**24 for a in lifts}
        assert all(
            multiplier**23 == multiplier for multiplier in return_multipliers
        )
        return_parities = {
            int(multiplier**11 == -field.one())
            for multiplier in return_multipliers
        }
        assert return_parities == {norm_one_parity(tau)}


# Rational E8: the residue field and every weighted coefficient are already
# F_23-rational, so Frobenius transport is the identity.
rational_e8_transport = F23.one()
assert not nontrivial_mod_roots_of_unity(rational_e8_transport, 15)


# Unramified E8 over F_23[r]/(r^2+18r+1).  This is the target component of
# the exact Frobenius transport already certified in
# certify_e8_tail_isomorphism.py.
Ru = PolynomialRing(F23, "u")
u = Ru.gen()
ku = F23.extension(u**2 + 18 * u + 1, "ru")
ru = ku.gen()
certify_hilbert90_orientation(ku)
unramified_e8_transport = 22 * ru + 5
assert unramified_e8_transport.multiplicative_order() == 8
assert unramified_e8_transport * unramified_e8_transport**23 == 1
assert gcd(23, 15) == gcd(22, 15) == 1
assert nontrivial_mod_roots_of_unity(unramified_e8_transport, 15)
assert norm_one_parity(unramified_e8_transport) == 1


# Ramified A2+A6 point over F_23[r]/(r^2+r+1).
Rr = PolynomialRing(F23, "v")
v = Rr.gen()
k = F23.extension(v**2 + v + 1, "r")
r = k.gen()
certify_hilbert90_orientation(k)
rp = r**23
assert rp == 22 * r + 22

RX = PolynomialRing(k, "X")
X = RX.gen()

# A2 tail:
#   Y^2=A*X^3+B*X,     phi=Y*h(X).
# The displayed alpha,beta give an isomorphism from the r-model to its
# Frobenius conjugate with target transport exactly 1.
A = 14 * r + 15
B = 14 * r + 11
Ap = 14 * rp + 15
Bp = 14 * rp + 11
h = (
    (18 * r + 7) * X**10
    + (8 * r + 11) * X**8
    + (3 * r + 15) * X**6
    + (21 * r + 12) * X**4
    + (19 * r + 6) * X**2
    + r
    + 5
)
hp = (
    (18 * rp + 7) * X**10
    + (8 * rp + 11) * X**8
    + (3 * rp + 15) * X**6
    + (21 * rp + 12) * X**4
    + (19 * rp + 6) * X**2
    + rp
    + 5
)
a2_alpha = 22 * r + 22
a2_beta = 7 * r + 11
a2_transport = k.one()
assert a2_beta**2 * A == Ap * a2_alpha**3
assert a2_beta**2 * B == Bp * a2_alpha
assert a2_beta * hp(a2_alpha * X) == a2_transport * h
assert a2_alpha * a2_alpha**23 == 1
assert a2_beta * a2_beta**23 == 1
assert gcd(23, 4) == 1
# A root change alters the Frobenius transport by
# zeta^(-23*(23-1))=zeta^2 in mu_4; surjectivity of this exponent on mu_4
# is not needed for the coset.
assert not nontrivial_mod_roots_of_unity(a2_transport, 4)
assert norm_one_parity(a2_transport) == 0


# A6 inner tail in the exact normalization returned by
# explore_p23_ade_deformation.py:
#   A6*X^7+B6*Y^2+C6=0,     phi=d6*X^4.
# Scaling the equation itself is harmless.  The displayed alpha,beta send
# the r-model to its Frobenius conjugate; gamma is the induced target
# transport.
A6 = 6 * r - 2
B6 = 2 * r + 9
C6 = -4 * r + 5
d6 = -4 * r - 9
A6p = 6 * rp - 2
B6p = 2 * rp + 9
C6p = -4 * rp + 5
d6p = -4 * rp - 9
a6_equation_scale = C6p / C6
a6_alpha = 12 * r + 19
a6_beta = k.one()
a6_transport = 11 * r + 7
assert A6p * a6_alpha**7 == a6_equation_scale * A6
assert B6p * a6_beta**2 == a6_equation_scale * B6
assert d6p * a6_alpha**4 == a6_transport * d6
assert a6_alpha * a6_alpha**23 == 1
assert a6_transport * a6_transport**23 == 1
assert a6_transport.multiplicative_order() == 24
assert gcd(23, 7) == gcd(22, 7) == 1
assert nontrivial_mod_roots_of_unity(a6_transport, 7)
assert norm_one_parity(a6_transport) == 1

# Root-choice independence of the binary character.  If the target has
# rho-weight n and rho is replaced by zeta*rho, Frobenius transport changes
# by zeta^(-n*(23-1)).  For (m,n)=(15,23) and (7,15) these factors have odd
# order, hence are squares in the norm-one torus.  For (m,n)=(4,23), the
# factor is in mu_2 and its twelfth power is also 1.
assert (-23 * 22) % 15 == 4
assert (-15 * 22) % 7 == 6
assert (-23 * 22) % 4 == 2


# A reducible ADE configuration is marked nontrivial when at least one of
# its attached weighted tails has nontrivial target-gluing transport.
gluing_marker = {
    "rational_E8": norm_one_parity(rational_e8_transport),
    "unramified_E8": norm_one_parity(unramified_e8_transport),
    "ramified_A2_A6": (
        norm_one_parity(a2_transport) + norm_one_parity(a6_transport)
    )
    % 2,
}
assert gluing_marker == {
    "rational_E8": 0,
    "unramified_E8": 1,
    "ramified_A2_A6": 1,
}

# Compare directly with the singular-position idempotent on
# F_23[u]/((u-16)(u^2+1)(u^2+u+1)).
RB = PolynomialRing(F23, "s")
s = RB.gen()
e_sing = 2 * s**4 + 2 * s**3 + 4 * s**2 + 2 * s + 3
assert e_sing.mod(s - 16) == 0
assert e_sing.mod(s**2 + 1) == 1
assert e_sing.mod(s**2 + s + 1) == 1
assert list(gluing_marker.values()) == [0, 1, 1]

print("rational_E8_target_transport=identity_mod_mu15")
print("unramified_E8_target_transport=order8_nontrivial_mod_mu15")
print("ramified_A2_target_transport=identity_mod_mu4")
print("ramified_A6_target_transport=order24_nontrivial_mod_mu7")
print("norm_one_quadratic_character=tau^12_in_{1,-1}")
print("hilbert90_lift_square_class=q(tau)_via_a^22=tau")
print("F23_units_are_squares_in_F23_squared=true")
print("Hilbert90_square_class_alone_does_not_orient_affine_F23_frame=true")
print("semilinear_return_multiplier=Norm_F23_squared_over_F23(a)=a^24")
print("semilinear_return_legendre_symbol_equals_q(tau)=true")
print("semilinear_return_parity_is_Hilbert90_gauge_independent=true")
print("nontrivial_transports_have_tau^12=-1")
print("ADE_gluing_marker=0,1,1")
print("ADE_gluing_marker_equals_singular_position_idempotent=true")
print("SCOPE=local_marker_only;normalizer_trace_comparison_is_separate_connector")
print("PASS weighted ADE target-gluing marker")
