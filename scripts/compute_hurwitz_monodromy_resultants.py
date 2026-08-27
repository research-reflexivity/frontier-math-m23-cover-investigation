#!/usr/bin/env sage-python
"""Compute exact fibre resultants for the degree-23 Hurwitz eliminant.

This uses the fact that the reconstructed quintics are linear in the fourth
canonical coordinate.  It avoids generic fraction-field simplification.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from sage.all import PolynomialRing, QQ, prod, save


ROOT = Path(__file__).resolve().parents[1]
CANONICAL_PATH = ROOT / "data" / "hurwitz_canonical_models_candidate.json"
MAPS_PATH = ROOT / "data" / "hurwitz_degree23_maps_candidate.json"


def rational(record):
    return QQ(record["numerator"]) / QQ(record["denominator"])


parser = argparse.ArgumentParser()
parser.add_argument("--target", type=int, action="append")
parser.add_argument("--output-dir", type=Path)
arguments = parser.parse_args()
targets = arguments.target or [2, 3]

started = time.time()
canonical = json.loads(CANONICAL_PATH.read_text())
maps = json.loads(MAPS_PATH.read_text())
field_ring = PolynomialRing(QQ, "a")
K = QQ.extension(
    field_ring([
        rational(value)
        for value in canonical["absolute_field_polynomial_coefficients_ascending"]
    ]),
    "a",
)
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


plane_ring = PolynomialRing(K, 2, names=("x", "y"))
x, y = plane_ring.gens()
z_ring = PolynomialRing(plane_ring, "z")
z = z_ring.gen()
coordinates = [z_ring(1), z_ring(x), z_ring(y), z]


def indexed_form(monomials, coefficients):
    return sum(
        z_ring(coefficient) * prod(coordinates[index] for index in monomial)
        for monomial, coefficient in zip(monomials, coefficients)
    )


quadric = indexed_form(
    [tuple(value) for value in canonical["quadric"]["monomials"]],
    [integral_basis_element(value) for value in canonical["quadric"]["coefficients"]],
)
cubic = indexed_form(
    [tuple(value) for value in canonical["petri_cubic"]["monomials"]],
    [
        integral_basis_element(value)
        for value in canonical["petri_cubic"]["coefficients"]
    ],
)

plane_curve = plane_ring(quadric.resultant(cubic))
plane_curve /= plane_curve.content()

q0, q1, q2 = (quadric[index] for index in range(3))
c0, c1, c2, c3 = (cubic[index] for index in range(4))
auxiliary = c2 * q2 - c3 * q1
linear_coefficient = q2 * (c1 * q2 - c3 * q0) - auxiliary * q1
constant_coefficient = c0 * q2**2 - auxiliary * q0
common = linear_coefficient.gcd(constant_coefficient)
linear_coefficient //= common
constant_coefficient //= common

sextic_map = next(
    value for value in maps["components"] if value["component"] == "sextic"
)


def plane_section(records):
    value = plane_ring.zero()
    for record in records:
        exponents = record["monomial"]
        coefficient = power_basis_element(record["coefficient_power_basis"])
        if exponents[3] == 0:
            multiplier = linear_coefficient
        elif exponents[3] == 1:
            multiplier = -constant_coefficient
        else:
            raise AssertionError("section is not linear in the eliminated coordinate")
        value += coefficient * x ** exponents[1] * y ** exponents[2] * multiplier
    return value


plane_numerator = plane_section(sextic_map["numerator"])
plane_denominator = plane_section(sextic_map["normalized_denominator"])
section_common = plane_numerator.gcd(plane_denominator)
plane_numerator //= section_common
plane_denominator //= section_common
print(
    "plane",
    plane_curve.total_degree(),
    plane_numerator.total_degree(),
    plane_denominator.total_degree(),
    "linear remainder degrees",
    linear_coefficient.total_degree(),
    constant_coefficient.total_degree(),
    "seconds",
    time.time() - started,
    flush=True,
)

x_ring = PolynomialRing(K, "x_fiber")
x_fiber = x_ring.gen()
y_ring = PolynomialRing(x_ring, "y_fiber")
y_fiber = y_ring.gen()


def lift(polynomial):
    return y_ring(
        sum(
            coefficient * x_fiber ** exponents[0] * y_fiber ** exponents[1]
            for exponents, coefficient in polynomial.dict().items()
        )
    )


lifted_curve = lift(plane_curve)
lifted_numerator = lift(plane_numerator)
lifted_denominator = lift(plane_denominator)


def fibre_resultant(value):
    local_started = time.time()
    result = x_ring(
        lifted_curve.resultant(lifted_numerator - K(value) * lifted_denominator)
    )
    print(
        "resultant",
        value,
        "degree",
        result.degree(),
        "seconds",
        time.time() - local_started,
        flush=True,
    )
    return result


resultants = {target: fibre_resultant(target) for target in targets}
if arguments.output_dir is not None:
    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    for target, resultant in resultants.items():
        output = arguments.output_dir / f"hurwitz_resultant_{target}.sobj"
        save(resultant, str(output))
        print("saved", output, flush=True)
if len(resultants) >= 2:
    first, second = list(resultants.values())[:2]
    base_factor = first.gcd(second)
    print(
        "base factor degree",
        base_factor.degree(),
        "quotient degrees",
        *(resultant.quo_rem(base_factor)[0].degree() for resultant in resultants.values()),
        "total seconds",
        time.time() - started,
        flush=True,
    )
