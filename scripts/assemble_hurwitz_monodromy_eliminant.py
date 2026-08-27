#!/usr/bin/env sage-python
"""Assemble the exact degree-23 target eliminant from parallel resultants."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ, load


ROOT = Path(__file__).resolve().parents[1]
CANONICAL_PATH = ROOT / "data" / "hurwitz_canonical_models_candidate.json"
MAPS_PATH = ROOT / "data" / "hurwitz_degree23_maps_candidate.json"
BRANCH_PATH = ROOT / "data" / "hurwitz_degree23_branch_candidate.json"
DEFAULT_OUTPUT = ROOT / "data" / "hurwitz_monodromy_eliminant_candidate.json"


def encode_rational(value):
    value = QQ(value)
    return {"numerator": int(value.numerator()), "denominator": int(value.denominator())}


def encode_element(value, degree):
    coefficients = list(value) + [QQ.zero()] * degree
    return [encode_rational(coefficient) for coefficient in coefficients[:degree]]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()

    loaded = {
        target: load(
            str(arguments.input_dir / f"hurwitz_resultant_{target}.sobj")
        )
        for target in range(2, 10)
    }
    parent = loaded[2].parent()
    resultants = {target: parent(value) for target, value in loaded.items()}
    field = parent.base_ring()
    base_factor = resultants[2].gcd(resultants[3]).monic()
    if base_factor.degree() != 19:
        raise AssertionError("target-independent base factor is not degree 19")

    quotients = {}
    for target, resultant in resultants.items():
        quotient, remainder = resultant.quo_rem(base_factor)
        if remainder:
            raise AssertionError(f"nonzero base-factor remainder at target {target}")
        if quotient.degree() != 23:
            raise AssertionError(f"quotient at target {target} is not degree 23")
        quotients[target] = quotient

    target_ring = PolynomialRing(field, "target")
    target = target_ring.gen()
    coefficient_polynomials = []
    for root_degree in range(24):
        values = [
            (field(sample), quotients[sample][root_degree])
            for sample in range(2, 9)
        ]
        coefficient = target_ring.lagrange_polynomial(values)
        if coefficient.degree() > 6:
            raise AssertionError("target degree exceeds resultant bound")
        coefficient_polynomials.append(coefficient)

    def evaluate_eliminant(value):
        return parent([
            coefficient(field(value)) for coefficient in coefficient_polynomials
        ])

    for sample in range(2, 9):
        if evaluate_eliminant(sample) != quotients[sample]:
            raise AssertionError(f"interpolation failed at target {sample}")
    if evaluate_eliminant(9) != quotients[9]:
        raise AssertionError("unused target-9 holdout failed")

    zero_fibre = evaluate_eliminant(0)
    one_fibre = evaluate_eliminant(1)
    infinity_fibre = parent([
        coefficient[6] for coefficient in coefficient_polynomials
    ])
    zero_gcd_degree = int(zero_fibre.gcd(zero_fibre.derivative()).degree())
    one_gcd = one_fibre.gcd(one_fibre.derivative())
    infinity_gcd_degree = int(
        infinity_fibre.gcd(infinity_fibre.derivative()).degree()
    )
    if zero_gcd_degree != 22:
        raise AssertionError("zero fibre is not totally ramified")
    if one_gcd.degree() != 8 or one_gcd.gcd(one_gcd.derivative()).degree() != 0:
        raise AssertionError("one fibre does not have eight distinct double roots")
    if infinity_gcd_degree != 22:
        raise AssertionError("infinity fibre is not totally ramified")

    field_degree = int(field.degree())
    coefficients = []
    for target_degree in range(7):
        row = []
        for root_degree in range(24):
            row.append(
                encode_element(
                    coefficient_polynomials[root_degree][target_degree], field_degree
                )
            )
        coefficients.append(row)

    sources = [CANONICAL_PATH, MAPS_PATH, BRANCH_PATH]
    payload = {
        "schema": "m23.cover-investigation.monodromy-eliminant-candidate.v1",
        "status": "exact_degree23_eliminant_for_reconstructed_sextic_component",
        "field_degree": field_degree,
        "field_defining_polynomial_coefficients_ascending": [
            encode_rational(value) for value in field.polynomial()
        ],
        "coordinate": "x=X1/X0 in the canonical plane projection X0=1",
        "normalized_target": "beta=N/(lambda*D)",
        "coefficient_convention": (
            "coefficients_target_then_x[i][j] is the power-basis encoding "
            "of [target^i*x^j] P(target,x)"
        ),
        "target_degree": 6,
        "root_degree": 23,
        "coefficients_target_then_x": coefficients,
        "construction": {
            "plane_curve_degree": 6,
            "plane_section_degrees": [7, 7],
            "raw_resultant_degree": 42,
            "target_independent_base_factor_degree": 19,
            "interpolation_targets": list(range(2, 9)),
            "holdout_target": 9,
            "holdout_exact": True,
        },
        "exact_branch_fibres": {
            "zero_gcd_degree": zero_gcd_degree,
            "one_gcd_degree": int(one_gcd.degree()),
            "one_gcd_squarefree": True,
            "infinity_gcd_degree": infinity_gcd_degree,
            "passport": ["23", "2^8 1^7", "23"],
        },
        "input_sha256": {
            str(path.relative_to(ROOT)): sha256(path) for path in sources
        },
    }
    arguments.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print("PASS exact degree-19 base factor and degree-23 quotient")
    print("PASS seven-point interpolation and unused target-9 holdout")
    print("PASS exact eliminant passport (23), (2^8 1^7), (23)")
    print(f"wrote {arguments.output}")


if __name__ == "__main__":
    main()
