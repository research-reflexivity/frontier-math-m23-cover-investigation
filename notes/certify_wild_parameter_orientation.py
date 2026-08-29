#!/usr/bin/env sage-python
"""Certify target transport = normalized wild-parameter orientation.

For a degree-p map with a totally ramified pole, normalize a local source
parameter u so that the target function starts u^(-p).  If a conjugate-tail
isomorphism transports u by eta and the target by tau, comparison of the
leading terms gives eta^(-p)=tau.  For p=23 on the norm-one group of order
24, eta^(-23)=eta, so eta=tau.

The E8 and A2 computations below derive eta from the exact curve parameters
and leading coefficients.  For the A6 outer chart, elimination of its
15th-root coordinate scaling proves that its target transport differs from
the inner A6 transport by an odd-order factor.  Hence their binary
orientations agree.  The raw outer chart has S23 monodromy, however, so
identification with a pointed M23 node is checked separately by the focused
A6 annular-normalization target.
"""

from sage.all import GF, PolynomialRing


p = 23
F23 = GF(p)


def norm_one_parity(value):
    assert value * value**p == 1
    sign = value**12
    assert sign in (value.parent().one(), -value.parent().one())
    return int(sign == -value.parent().one())


def normalized_parameter_transport(raw_transport, leading_coefficient):
    """Transport of u=s/delta when delta^23=c.

    In F_(23^2), the unique 23rd root of c is c^23.  Its Frobenius
    conjugate is c, giving eta=lambda*c^22.
    """

    field = leading_coefficient.parent()
    delta = leading_coefficient**p
    assert delta**p == leading_coefficient
    delta_conjugate = delta**p
    eta = raw_transport * delta / delta_conjugate
    assert eta in field
    return eta


# Unramified E8 tail:
#   Y^3=a*X^5+b,  phi=X*Y*(X^5+c0).
# At infinity s=Y/X^2 is a uniformizer.  Since
# X~a*s^-3 and Y~a^2*s^-5, phi~a^8*s^-23.
RE = PolynomialRing(F23, "e")
e = RE.gen()
kE = F23.extension(e**2 + 18 * e + 1, "rE")
rE = kE.gen()
e8_curve_a = 18 * rE + 5
e8_leading = e8_curve_a**8
e8_x_transport = 18 * rE + 1
e8_y_transport = rE + 18
e8_raw_parameter_transport = e8_y_transport / e8_x_transport**2
e8_target_transport = 22 * rE + 5
e8_parameter_transport = normalized_parameter_transport(
    e8_raw_parameter_transport, e8_leading
)
assert e8_leading == 15
# This checks the direction of the semilinear comparison before
# normalization: c^sigma*lambda^(-23)=tau*c.
assert (
    e8_leading**p * e8_raw_parameter_transport**(-p)
    == e8_target_transport * e8_leading
)
assert e8_parameter_transport == e8_target_transport
assert e8_parameter_transport**(-p) == e8_target_transport
assert norm_one_parity(e8_parameter_transport) == 1


# A2 tail:
#   Y^2=A*X^3+B*X,  phi=Y*(h10*X^10+...).
# At infinity s=X/Y is a uniformizer, with
# X~A^-1*s^-2 and Y~A^-1*s^-3.  Hence
# phi~h10*A^-11*s^-23.  The raw parameter transport is
# alpha/beta; normalization by the leading coefficient cancels its
# nonsquare class exactly.
RA = PolynomialRing(F23, "a")
a = RA.gen()
kA = F23.extension(a**2 + a + 1, "rA")
rA = kA.gen()
a2_curve_A = 14 * rA + 15
a2_h10 = 18 * rA + 7
a2_leading = a2_h10 * a2_curve_A**(-11)
a2_x_transport = 22 * rA + 22
a2_y_transport = 7 * rA + 11
a2_raw_parameter_transport = a2_x_transport / a2_y_transport
a2_target_transport = kA.one()
a2_parameter_transport = normalized_parameter_transport(
    a2_raw_parameter_transport, a2_leading
)
assert (
    a2_leading**p * a2_raw_parameter_transport**(-p)
    == a2_target_transport * a2_leading
)
assert a2_parameter_transport == a2_target_transport
assert a2_parameter_transport**(-p) == a2_target_transport
assert norm_one_parity(a2_parameter_transport) == 0
assert norm_one_parity(a2_raw_parameter_transport) == 1
assert norm_one_parity(a2_leading**22) == 1


# The A6 outer chart is A*X^23+B*X^8.  If X -> lambda*X identifies it with
# its Frobenius conjugate and gamma is the target transport, then
#
#   lambda^15=A*B^sigma/(B*A^sigma),
#   gamma=(B^sigma/B)*lambda^8.
#
# Eliminate lambda.  Relative to the inner A6 target transport tau, the
# exact coefficients give (gamma/tau)^15=r, hence (gamma/tau)^45=1.  The
# discrepancy has odd order, so it is invisible to the quadratic
# orientation.  This proves equality of the inner and outer binary target
# classes without selecting a 15th root.
a6_target_transport = 11 * rA + 7
assert a6_target_transport.multiplicative_order() == 24
assert (-p) % 24 == 1
a6_parameter_transport = a6_target_transport
assert a6_parameter_transport**(-p) == a6_target_transport
assert norm_one_parity(a6_parameter_transport) == 1

Rw = PolynomialRing(kA, "v")
v = Rw.gen()
kA4 = kA.extension(v**2 - (6 * rA + 22), "w")
w = kA4.gen()
r4 = kA4(rA)
outer_A = 15 * r4 * w
outer_B = -4 * r4 - 9
outer_A_sigma = outer_A**p
outer_B_sigma = outer_B**p
lambda_power_15 = (
    outer_A * outer_B_sigma / (outer_B * outer_A_sigma)
)
outer_over_inner_power_15 = (
    (outer_B_sigma / outer_B) ** 15
    * lambda_power_15**8
    / kA4(a6_target_transport) ** 15
)
assert outer_over_inner_power_15 == r4
assert r4**3 == 1
# For every possible root lambda, (gamma/tau)^45=1.  Thus the discrepancy
# belongs to an odd-order subgroup; its image in the quadratic quotient is
# trivial.


# Rational E8 is already defined over F23, so both transports are trivial.
rational_e8_target_transport = F23.one()
rational_e8_parameter_transport = F23.one()
assert rational_e8_parameter_transport**(-p) == rational_e8_target_transport


print("local_degree23_orientation_identity=eta^(-23)=tau")
print("norm_one_order=24_and_-23_congruent_1")
print("unramified_E8_normalized_parameter_transport=target_transport_order8")
print("A2_raw_parameter_nonsquare_cancelled_by_leading_coefficient=true")
print("A2_normalized_parameter_transport=target_transport_identity")
print("A6_outer_over_inner_transport_power15=r_of_order3")
print("A6_outer_and_inner_quadratic_orientations_agree=true")
print("orientation_values_rationalE8_unramifiedE8_A2_A6=0,1,0,1")
print("SCOPE=A6_pointed_node_realization_checked_by_verify_hurwitz_connector_a6")
print("PASS_WILD_PARAMETER_ORIENTATION")
