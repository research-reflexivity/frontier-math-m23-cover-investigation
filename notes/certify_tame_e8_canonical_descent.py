#!/usr/bin/env sage-python
"""Exact algebra for tame canonical descent from a (3,5) genus-four tail.

No reconstructed generic curve, map, or Hurwitz polynomial is read.
This checks the character arithmetic, the universal CLOSED normal form,
and elementary curve identities. It does not certify that a particular
patched G-cover has the specified inertia action or marked annuli.
The relative construction is proved in the accompanying research note.
"""
from sage.all import GF, QQ, PolynomialRing, PowerSeriesRing, matrix


k = GF(23)
weights = (0, 3, 5, 6)
names = ("H", "U", "V", "W")


def exponents(degree, slots=4):
    if slots == 1:
        yield (degree,)
        return
    for first in range(degree, -1, -1):
        for tail in exponents(degree - first, slots - 1):
            yield (first, *tail)


def label(exponent):
    return "*".join(
        name if power == 1 else f"{name}^{power}"
        for name, power in zip(names, exponent) if power
    )


closed = {}
first_allowed = {}
for degree, equation_weight in ((2, 6), (3, 15)):
    closed[degree] = set()
    first_allowed[degree] = {}
    print(f"DEGREE_{degree}_COEFFICIENT_ORDERS")
    for exponent in exponents(degree):
        weight = sum(e * w for e, w in zip(exponent, weights))
        good_order = (weight - equation_weight) % 15
        descended_order = good_order + equation_weight - weight
        assert descended_order >= 0 and descended_order % 15 == 0
        pi_order = descended_order // 15
        first_allowed[degree][label(exponent)] = pi_order
        if pi_order == 0:
            closed[degree].add(label(exponent))
        print(label(exponent), "weight", weight,
              "rho_order_good", good_order, "pi_order_descended", pi_order)

assert closed[2] == {"H*W", "U^2", "U*V", "U*W", "V^2", "V*W", "W^2"}
assert closed[3] == {"U*W^2", "V^3", "V^2*W", "V*W^2", "W^3"}
assert first_allowed[3]["H^3"] == 1
assert first_allowed[2]["H^2"] == 1
print("PASS_ALL_30_CHARACTER_AND_INTEGRALITY_CHECKS")

# Universal triangular normalization. The eight symbols are arbitrary
# closed-fibre coefficients, NOT values fitted to a known generic cover.
B = PolynomialRing(k, "a,b,c,d,e,f,g,h")
a, b, c, d, e, f, g, h = B.gens()
S = PolynomialRing(B, "z")
z = S.gen()
u = z**3 - f*z**2 - g*z - h
H = u**2 - a*u*z - b*u - c*z**2 - d*z - e
assert H - u**2 + a*u*z + b*u + c*z**2 + d*z + e == 0
assert u - z**3 + f*z**2 + g*z + h == 0

# This affine change fixes infinity. Cubic and lower canonical terms
# are subsequently removed by projective linear coordinate changes.
new_u = u(z + f/k(3))
new_H = H(z + f/k(3))
assert new_u[2] == 0 and new_u[3] == 1
assert new_H[5] == 0 and new_H[6] == 1
U_normal = new_u - new_u[1]*z - new_u[0]
H_normal = new_H - new_H[3]*z**3 - new_H[1]*z - new_H[0]
alpha, beta = H_normal[4], H_normal[2]
assert U_normal == z**3
assert H_normal == z**6 + alpha*z**4 + beta*z**2
assert H_normal - U_normal**2 - alpha*U_normal*z - beta*z**2 == 0
assert alpha(-1, 0, -1, 0, 0, 0, 0, 0) == 1
assert beta(-1, 0, -1, 0, 0, 0, 0, 0) == 1
print("PASS_UNIVERSAL_CLOSED_NORMALIZATION")
print("ALPHA", alpha)
print("BETA", beta)

semigroup = {3*i + 5*j for i in range(20) for j in range(20)}
gaps = [n for n in range(30) if n not in semigroup]
assert gaps == [1, 2, 4, 7]
assert all(n in semigroup for n in range(8, 30))
print("CUSP_GAPS", gaps, "CONDUCTOR", 8, "DELTA", len(gaps))

# Hilbert function of a (2,3) complete intersection in P^3.
T = PowerSeriesRing(QQ, "t", default_prec=12)
t = T.gen()
hilbert = (1-t**2)*(1-t**3)/(1-t)**4
assert [hilbert[i] for i in range(6)] == [1, 4, 9, 15, 21, 27]
assert all(hilbert[i] == 6*i-3 for i in range(2, 12))
print("CANONICAL_HILBERT_FUNCTION", [hilbert[i] for i in range(6)])

# Independently verify the good canonical fibre at lambda=1 is smooth.
# Any nonzero lambda is geometrically equivalent by scaling X and Y.
P = PolynomialRing(k, names)
H, U, V, W = P.gens()
q = H*W-U**2
c = U*W**2-V**3+H**3
jacobian = matrix(P, [[f.derivative(x) for x in P.gens()] for f in (q, c)])
singular = P.ideal([q, c, *jacobian.minors(2)])
for coordinate in P.gens():
    assert (singular + P.ideal(coordinate-1)).is_one()
print("PASS_GOOD_FIBRE_SMOOTH_ALL_FOUR_CHARTS")

# After descent, the only possible linear terms at [1:0:0:0] are
# W+a0*pi in the quadric and lambda*pi in the cubic. Their independent
# linear parts give a regular total space when lambda is a unit.
J = PolynomialRing(k, "a0,lambda")
a0, lam = J.gens()
linear_matrix = matrix(J, [[1, a0], [0, lam]])  # columns: W, pi
assert linear_matrix.det() == lam
print("REGULAR_TOTAL_SPACE_LINEAR_DETERMINANT", linear_matrix.det())

# Logarithmic annulus lengths, provided the source annuli map with
# degree 23 to the target annuli in the stated special-map construction.
for conductor in (7, 15):
    target_length = QQ(23*11)/((23-1)*conductor)  # v_23(23)=1
    source_length = target_length/23
    assert source_length == QQ(1)/(2*conductor)
    print("ANNULUS_LENGTH", conductor, "target", target_length,
          "source", source_length, "in_pi_units", 2*source_length)

# In an actual G-node xy=t^e the tame quotient has degree 11.
# For the mth r_c coefficient, v_t(r_c)=h*e, and dividing the
# displacement by the quotient smoothing scale subtracts 11*e.
# A surviving Laurent exponent must therefore be n=11-h*m.
M = PolynomialRing(QQ, "m")
m = M.gen()
for conductor in (7, 15):
    n = 11-conductor*m
    assert conductor*m+n-11 == 0
    assert (11-conductor)-n == conductor*(m-1)
    for order in range(2, 101):
        assert n(order) < n(1)
    print("SURVIVING_LAURENT_EXPONENT", conductor, n,
          "FIRST_EXPONENT", n(1))
assert 11-7 == 4 and 11-15 == -4
assert k(3)**11 == 1  # The X~xi^3 factor does not swap the two designs.
print("PASS_ANNULAR_EXPONENT_IDENTITIES")

print("SCOPE: algebra verified; inertia, patching, and marked-chart comparison require proofs")
