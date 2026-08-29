#!/usr/bin/env sage-python
"""Certify the explicit isomorphism between the two characteristic-23 E8 tails."""

from sage.all import GF, PolynomialRing
from math import gcd


F23 = GF(23)
R0 = PolynomialRing(F23, "u")
u = R0.gen()
k = GF(23**2, name="r", modulus=u**2 + 18 * u + 1)
r = k.gen()

alpha = 16 * r
beta = 18 * r + 5

# Rational E8 tail: y^3 = 22*x^5+6 and phi=3*x*y*(x^5+2).
# Unramified E8 tail: Y^3=(18*r+5)*X^5+(16*r+19) and
# phi_u=X*Y*(X^5+11*r).
assert alpha**5 == 17 * r
assert beta**3 == 18 * r + 7
assert (16 * r + 19) / 6 == beta**3
assert (18 * r + 5) * alpha**5 / beta**3 == 22
assert 11 * r == 2 * alpha**5

target_scale = alpha * beta * alpha**5 / 3
assert target_scale == 4 * r + 4

# The unramified quadratic closed point has two Frobenius-conjugate
# geometric tails.  Trivializing each by the rational tail gives the exact
# coarse Frobenius transport between them.
r_frobenius = r**23
source_x_transport = (16 * r_frobenius) / alpha
source_y_transport = (18 * r_frobenius + 5) / beta
target_transport = (4 * r_frobenius + 4) / target_scale
assert r_frobenius == 22 * r + 5
assert source_x_transport == 18 * r + 1
assert source_y_transport == r + 18
assert target_transport == 22 * r + 5
assert source_x_transport == alpha**23 / alpha
assert source_y_transport == beta**23 / beta
assert target_transport == target_scale**23 / target_scale
assert [
    value.multiplicative_order()
    for value in (source_x_transport, source_y_transport, target_transport)
] == [4, 8, 8]
assert all(
    value * value**23 == 1
    for value in (source_x_transport, source_y_transport, target_transport)
)

# The coarse isomorphism is allowed to rescale the target arbitrarily.  In
# the stable E8 chart, however, the smoothing base change is pi=rho^15 and
# the target has weight 23.  Replacing rho by zeta*rho changes the target
# coordinate by zeta^23; since gcd(23,15)=1, the residual ambiguity is
# exactly mu_15, not all scalar multiplications.  The Frobenius transition
# below has order 8, and hence survives this gluing ambiguity.
assert gcd(23, 15) == 1
assert target_scale.multiplicative_order() == 176
assert target_scale**15 != 1
assert target_transport.multiplicative_order() == 8
assert gcd(target_transport.multiplicative_order(), 15) == 1
assert target_transport**15 != 1
# Replacing rho by zeta*rho changes the target coordinate by zeta^(-23),
# hence changes its Frobenius transition by zeta^(-23*(23-1)).  Modulo 15
# this exponent is 4, which permutes mu_15.  In particular the transition
# coset modulo mu_15 is independent of the chosen fifteenth root.
root_change_exponent = (-23 * (23 - 1)) % 15
assert root_change_exponent == 4
assert gcd(root_change_exponent, 15) == 1

coefficient_x = 18 * r + 5
coefficient_constant = 16 * r + 19
coefficient_x_frobenius = 18 * r_frobenius + 5
coefficient_constant_frobenius = 16 * r_frobenius + 19
assert (
    source_y_transport**3 * coefficient_x
    == coefficient_x_frobenius * source_x_transport**5
)
assert (
    source_y_transport**3 * coefficient_constant
    == coefficient_constant_frobenius
)
assert (
    source_x_transport * source_y_transport * source_x_transport**5
    == target_transport
)
assert (
    source_x_transport * source_y_transport * 11 * r_frobenius
    == target_transport * 11 * r
)

print("source_change=X_u=(16*r)*X, Y_u=(18*r+5)*Y")
print("target_change=phi_u=(4*r+4)*phi_rational")
print("frobenius_transport_X=18*r+1_order4_norm1")
print("frobenius_transport_Y=r+18_order8_norm1")
print("frobenius_transport_target=22*r+5_order8_norm1")
print("coarse_frobenius_transport=isomorphism_coboundary")
print("stable_E8_target_root_ambiguity=mu15")
print("target_gluing_transport_not_in_mu15=true")
print("PASS the rational and unramified E8 tail covers are geometrically isomorphic")
