#!/usr/bin/env sage-python
"""Exact characteristic-23 special-deformation-datum calculation.

This certificate starts from the residual Kummer and refined-Swan formulas
already proved for the z-Gauss valuation.  After adjoining pi with pi^11=s,
the tame quotient is normalized by a^11=z/(z-2).  It verifies the resulting
differential, its C_11 character, its critical signature, and the lift of
branch exchange.
"""

from fractions import Fraction
from math import lcm

from sage.all import GF, PolynomialRing


k = GF(23)

# Verify the rational identities in the z-coordinate.
Rz = PolynomialRing(k, "z")
z0 = Rz.gen()
Kz = Rz.fraction_field()
z = Kz(z0)
E = Kz(
    8 * z0**9
    + 7 * z0**8
    + 11 * z0**7
    + 16 * z0**6
    + 3 * z0**5
    + 19 * z0**4
    + 21 * z0**3
    + 21 * z0**2
    + 17 * z0
    + 15
)
q = Kz(z0**2 + 13 * z0 + 10)
ybar = -E / (z - 1) ** 5
eta_coefficient = ybar.derivative() / ybar

assert eta_coefficient == 9 * z * (z - 2) ** 5 * q / E
assert q(z / (z - 1)) * (z - 1) ** 2 == q

# After pi^11=s, remove the displayed 11th power from the tame equation.
assert (
    (z * q / (z - 2) ** 3) ** 11 * (z / (z - 2)) ** 4
    == z**15 * q**11 / (z - 2) ** 37
)

# Work on Z_0=P^1_a, where z=2*a^11/(a^11-1).
Ra = PolynomialRing(k, "a")
a0 = Ra.gen()
Ka = Ra.fraction_field()
a = Ka(a0)
z_of_a = 2 * a**11 / (a**11 - 1)
assert z_of_a / (z_of_a - 2) == a**11


def pullback(function):
    """Substitute z=z(a) in an element of F_23(z)."""

    numerator = function.numerator()(z_of_a)
    denominator = function.denominator()(z_of_a)
    return Ka(numerator / denominator)


# Divide the refined-Swan tensor by its uniformizer degree.  Its residual
# coefficient is B, and B^{-1}*eta is the eigen-differential on Z_0.
B = (
    8
    * a**4
    * z_of_a
    * pullback(q)
    * (z_of_a - 2) ** 5
    * (z_of_a - 1)
    / pullback(E)
)
omega_coefficient = pullback(eta_coefficient) * z_of_a.derivative() / B
expected_omega_coefficient = 4 * a**6 / (a**22 - 1)

assert omega_coefficient == expected_omega_coefficient
assert omega_coefficient == 4 * a ** (-4) * z_of_a.derivative() / (z_of_a - 1)

# H=C_11 acts by a |-> lambda*a and omega |-> lambda^7*omega.
primitive = k.multiplicative_generator()
lam = primitive**2
assert lam.multiplicative_order() == 11
assert (
    expected_omega_coefficient(lam * a) * lam
    == lam**7 * expected_omega_coefficient
)

# Divisor of omega=(4*a^6/(a^22-1)) da.
ord_zero = 6
ord_infinity = 22 - 6 - 2
assert ord_infinity == 14
assert (a0**22 - 1).is_squarefree()
assert set((a0**22 - 1).roots(multiplicities=False)) == set(k) - {k.zero()}

# An explicit dlog witness.  At a=c in F_23^x the residue is -4*c^7.
# Taking the least nonnegative integer representative of every residue gives
# a rational function with precisely this logarithmic derivative.
logarithmic_function = Ka.one()
for element in k:
    if element != 0:
        residue = -4 * element**7
        logarithmic_function *= (a - element) ** int(residue)
assert logarithmic_function.derivative() / logarithmic_function == expected_omega_coefficient

signature = {
    "z=0": (11, ord_zero + 1, Fraction(ord_zero + 1, 11)),
    "z=1": (1, 0, Fraction(0, 1)),
    "z=infinity": (1, 0, Fraction(0, 1)),
    "z=2": (11, ord_infinity + 1, Fraction(ord_infinity + 1, 11)),
}
assert signature == {
    "z=0": (11, 7, Fraction(7, 11)),
    "z=1": (1, 0, Fraction(0, 1)),
    "z=infinity": (1, 0, Fraction(0, 1)),
    "z=2": (11, 15, Fraction(15, 11)),
}
assert sum(value[2] - 1 for value in signature.values()) == -2
assert Fraction(7, 11) + Fraction(4, 11) == 1

# For any special G-map realizing this deformation datum, Wewers's
# patching-data theorem gives
#
#   |P(fbar)|=(p-1)*prod_j h_j,
#   orbit length=(p-1)*lcm_j(h_j),
#
# over the non-wild primitive and new tails.  Here h_j=7,15 are coprime,
# so the *unquotiented* patching torsor is one orbit of size 2310.  The
# numerical factor 15 matches the number of C_G(x)\G/(23:11) pointing
# colors, but this calculation does not construct a morphism between the
# two sets; that incidence map is the remaining geometric problem.
nonwild_conductors = [7, 15]
patching_data_count = (23 - 1)
for conductor in nonwild_conductors:
    patching_data_count *= conductor
patching_orbit_length = (23 - 1) * lcm(*nonwild_conductors)
pointing_color_count = 15
putative_color_fiber_size = patching_data_count // pointing_color_count
assert patching_data_count == 2310
assert patching_orbit_length == patching_data_count
assert putative_color_fiber_size == 154 == 22 * 7

# The 22 poles split into square and nonsquare H-orbits, above infinity and 1.
squares = {element**2 for element in k if element != 0}
nonsquares = set(k) - {k.zero()} - squares
assert len(squares) == len(nonsquares) == 11
assert all(element**11 == 1 for element in squares)
assert all(element**11 == -1 for element in nonsquares)

# tau(z)=z/(z-1) lifts to a |-> -a and reverses omega.
tau_z = z_of_a / (z_of_a - 1)
assert z_of_a(-a) == tau_z
assert expected_omega_coefficient(-a) * (-1) == -expected_omega_coefficient

print("tame_quotient=a^11=z/(z-2)")
print("omega=4*a^6*da/(a^22-1)")
print("omega_is_logarithmic=true")
print("critical_signature=z0:(11,7),z1:(1,0),zinf:(1,0),z2:(11,15)")
print("tails=primitive_at_z0,new_at_z2")
print("branch_exchange=a_to_-a,omega_to_-omega")
print("wewers_patching_data_count=2310")
print("wewers_patching_galois_orbit_length=2310")
print("numerical_fifteen_color_fiber_size=154")
print("SCOPE_no_canonical_map_from_patching_data_to_pointing_colors_constructed")
print("PASS_P23_SPECIAL_DEFORMATION_DATUM")
