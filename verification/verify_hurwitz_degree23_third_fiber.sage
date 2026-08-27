#!/usr/bin/env sage
"""Long exact certificate for the third fibre on the sextic Hurwitz component.

This deliberately avoids numerical root finding.  It reconstructs the exact
canonical complete intersection and exact quintic numerator/denominator over
the degree-12 Hurwitz field, projects the curve birationally to a plane sextic,
and compares a generic fibre with the exact third fibre.  The latter acquires
eight additional, distinct double roots.

Expected runtime on the reference machine is tens of minutes.  It is therefore
an opt-in certificate rather than part of the quick ``make verify-all`` suite.
"""

import json
from pathlib import Path

from sage.all import PolynomialRing, QQ, prod


ROOT = Path(__file__).resolve().parents[1]
CANONICAL_PATH = ROOT / "data" / "hurwitz_canonical_models_candidate.json"
MAPS_PATH = ROOT / "data" / "hurwitz_degree23_maps_candidate.json"
BRANCH_PATH = ROOT / "data" / "hurwitz_degree23_branch_candidate.json"


def rational(record):
    return QQ(record["numerator"]) / QQ(record["denominator"])


canonical = json.loads(CANONICAL_PATH.read_text())
maps = json.loads(MAPS_PATH.read_text())
branch = json.loads(BRANCH_PATH.read_text())

field_ring = PolynomialRing(QQ, "a")
a_polynomial = field_ring([
    rational(value)
    for value in canonical["absolute_field_polynomial_coefficients_ascending"]
])
K = QQ.extension(a_polynomial, "a")
integral_basis = K.integral_basis()


def integral_basis_element(record):
    return sum(
        rational(value) * basis
        for value, basis in zip(
            record["sextic_component"]["integral_basis_coordinates"],
            integral_basis,
        )
    )


def power_basis_element(records):
    return K([rational(value) for value in records])


quadric_monomials = [tuple(value) for value in canonical["quadric"]["monomials"]]
quadric_coefficients = [
    integral_basis_element(value) for value in canonical["quadric"]["coefficients"]
]
cubic_monomials = [tuple(value) for value in canonical["petri_cubic"]["monomials"]]
cubic_coefficients = [
    integral_basis_element(value)
    for value in canonical["petri_cubic"]["coefficients"]
]

sextic_map = next(
    value for value in maps["components"] if value["component"] == "sextic"
)


def sparse_section(records):
    return [
        (
            tuple(record["monomial"]),
            power_basis_element(record["coefficient_power_basis"]),
        )
        for record in records
    ]


numerator_section = sparse_section(sextic_map["numerator"])
denominator_section = sparse_section(sextic_map["denominator"])
branch_value = K([
    QQ(value) / QQ(branch["lambda"]["denominator"])
    for value in branch["lambda"]["power_basis_numerators"]
])

plane_ring = PolynomialRing(K, 2, names=("x", "y"))
x, y = plane_ring.gens()
plane_fraction = plane_ring.fraction_field()
z_ring = PolynomialRing(plane_fraction, "z")
z = z_ring.gen()
coordinates = [z_ring(1), z_ring(x), z_ring(y), z]


def indexed_form(monomials, coefficients):
    return sum(
        z_ring(coefficient) * prod(coordinates[index] for index in monomial)
        for monomial, coefficient in zip(monomials, coefficients)
    )


def exponent_form(section):
    return sum(
        z_ring(coefficient)
        * prod(
            coordinates[index] ** exponent
            for index, exponent in enumerate(monomial)
        )
        for monomial, coefficient in section
    )


quadric = indexed_form(quadric_monomials, quadric_coefficients)
cubic = indexed_form(cubic_monomials, cubic_coefficients)
plane_curve = quadric.resultant(cubic).numerator()
plane_curve /= plane_curve.content()
linear_subresultant = next(
    value for value in reversed(quadric.subresultants(cubic))
    if value.degree() == 1
)
z_expression = -linear_subresultant[0] / linear_subresultant[1]

map_fraction = plane_fraction(
    exponent_form(numerator_section)(z_expression)
    / exponent_form(denominator_section)(z_expression)
)
plane_numerator = map_fraction.numerator()
plane_denominator = map_fraction.denominator()

assert plane_curve.total_degree() == 6
assert plane_numerator.total_degree() == 7
assert plane_denominator.total_degree() == 7

x_ring = PolynomialRing(K, "x_fiber")
x_fiber = x_ring.gen()
y_ring = PolynomialRing(x_ring, "y_fiber")
y_fiber = y_ring.gen()


def lift(polynomial):
    return y_ring(sum(
        coefficient
        * x_fiber ** exponents[0]
        * y_fiber ** exponents[1]
        for exponents, coefficient in polynomial.dict().items()
    ))


lifted_curve = lift(plane_curve)
lifted_numerator = lift(plane_numerator)
lifted_denominator = lift(plane_denominator)


def fibre_resultant(value):
    resultant = x_ring(
        lifted_curve.resultant(lifted_numerator - value * lifted_denominator)
    )
    return resultant, resultant.gcd(resultant.derivative())


generic_resultant, generic_gcd = fibre_resultant(K(2))
branch_resultant, branch_gcd = fibre_resultant(branch_value)
common_gcd = generic_gcd.gcd(branch_gcd)
extra_gcd = branch_gcd // common_gcd

assert generic_resultant.degree() == 42
assert branch_resultant.degree() == 42
assert generic_gcd.degree() == 6
assert branch_gcd.degree() == 14
assert common_gcd.degree() == 6
assert extra_gcd.degree() == 8
assert extra_gcd.gcd(extra_gcd.derivative()).degree() == 0

print("PASS exact characteristic-zero third fibre")
print("PASS generic/third gcd degrees 6/14 with squarefree extra degree 8")
