#!/usr/bin/env sage-python
"""Recognize a numerical seven-root fingerprint over Q(sqrt(-23)).

This is deliberately a reconstruction tool, not an exact certificate.  It
reads a scale-free canonical-quadric invariant from seven Acb model JSON
files, forms its monic root polynomial, and finds bounded-denominator rational
approximations to the real and ``sqrt(-23)`` parts of every coefficient.
Exact cover substitution is still required before accepting the output as
the Hurwitz algebra.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from sage.all import ComplexField, PolynomialRing, QQ, RealField, I


def rational_record(value, maximum_denominator: int) -> dict[str, object]:
    fraction = Fraction(value.str(truncate=False)).limit_denominator(
        maximum_denominator
    )
    approximation = value.parent()(fraction.numerator) / fraction.denominator
    return {
        "numerator": fraction.numerator,
        "denominator": fraction.denominator,
        "absolute_error": str(abs(value - approximation)),
    }


def invariant_value(payload: dict[str, object], name: str, field):
    quadric = payload["canonical_quadric"]
    coefficients = {
        tuple(monomial): field(coefficient)
        for monomial, coefficient in zip(
            quadric["monomials"], quadric["coefficients"]
        )
    }
    A = coefficients[(0, 3)]
    invariants = {
        "A": A,
        "q13_over_A2": coefficients[(1, 3)] / A**2,
        "q22_over_A2": coefficients[(2, 2)] / A**2,
        "q23_over_A3": coefficients[(2, 3)] / A**3,
        "q33_over_A4": coefficients[(3, 3)] / A**4,
        "q13_over_q22": coefficients[(1, 3)] / coefficients[(2, 2)],
    }
    return invariants[name]


def rational_pair(record: dict[str, object]):
    return QQ(record["numerator"]) / QQ(record["denominator"])


def encode_rational(value) -> dict[str, int]:
    return {"numerator": int(value.numerator()), "denominator": int(value.denominator())}


def encode_quadratic(value) -> dict[str, dict[str, int]]:
    coordinates = value.list()
    coordinates += [QQ(0)] * (2 - len(coordinates))
    return {
        "rational_part": encode_rational(coordinates[0]),
        "sqrt_minus_23_part": encode_rational(coordinates[1]),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("models", nargs=7, type=Path)
    parser.add_argument("--precision", type=int, default=1024)
    parser.add_argument("--max-denominator", type=int, default=10**40)
    parser.add_argument(
        "--invariant",
        choices=[
            "A",
            "q13_over_A2",
            "q22_over_A2",
            "q23_over_A3",
            "q33_over_A4",
            "q13_over_q22",
        ],
        default="q13_over_A2",
    )
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()

    complex_field = ComplexField(arguments.precision)
    real_field = RealField(arguments.precision)
    polynomial_ring = PolynomialRing(complex_field, "X")
    X = polynomial_ring.gen()
    roots = []
    class_ids = []
    payloads = []
    for model_path in arguments.models:
        payload = json.loads(model_path.read_text())
        payloads.append(payload)
        roots.append(invariant_value(payload, arguments.invariant, complex_field))
        class_ids.append(int(payload["class_id"]))
    if sorted(class_ids) != list(range(1, 8)):
        raise SystemExit("the seven inputs must contain class IDs 1 through 7")

    polynomial = polynomial_ring.one()
    for root in roots:
        polynomial *= X - root

    sqrt_23 = real_field(23).sqrt()
    records = []
    reconstructed = []
    for coefficient in polynomial.list():
        real_record = rational_record(
            real_field(coefficient.real()), arguments.max_denominator
        )
        sqrt_minus_23_record = rational_record(
            real_field(coefficient.imag()) / sqrt_23,
            arguments.max_denominator,
        )
        real_part = (
            complex_field(real_record["numerator"])
            / real_record["denominator"]
        )
        imaginary_part = (
            complex_field(sqrt_minus_23_record["numerator"])
            / sqrt_minus_23_record["denominator"]
            * complex_field(I)
            * complex_field(23).sqrt()
        )
        candidate = real_part + imaginary_part
        reconstructed.append(candidate)
        records.append(
            {
                "numerical": str(coefficient),
                "rational_part": real_record,
                "sqrt_minus_23_part": sqrt_minus_23_record,
                "reconstruction_error": str(abs(coefficient - candidate)),
            }
        )

    candidate_polynomial = polynomial_ring(reconstructed)
    root_residuals = [abs(candidate_polynomial(root)) for root in roots]

    rational_polynomial_ring = PolynomialRing(QQ, "s")
    s = rational_polynomial_ring.gen()
    quadratic_field = QQ.extension(s**2 + 23, "sqrt_minus_23")
    sqrt_minus_23 = quadratic_field.gen()
    exact_polynomial_ring = PolynomialRing(quadratic_field, "X")
    exact_coefficients = [
        quadratic_field(rational_pair(record["rational_part"]))
        + quadratic_field(rational_pair(record["sqrt_minus_23_part"]))
        * sqrt_minus_23
        for record in records
    ]
    exact_polynomial = exact_polynomial_ring(exact_coefficients)
    class_6_index = class_ids.index(6)
    class_6_root = roots[class_6_index]
    class_6_real = rational_record(real_field(class_6_root.real()), 10**12)
    class_6_sqrt = rational_record(
        real_field(class_6_root.imag()) / sqrt_23, 10**12
    )
    degree_one_exact = (
        quadratic_field(rational_pair(class_6_real))
        + quadratic_field(rational_pair(class_6_sqrt)) * sqrt_minus_23
    )
    exact_factors = list(exact_polynomial.factor())
    factor_degrees = [int(factor.degree()) for factor, _ in exact_factors]
    exact_linear_factor_check = exact_polynomial(degree_one_exact) == 0
    sextic_irreducible = sorted(factor_degrees) == [1, 6]
    result = {
        "status": "NUMERICAL_RECONSTRUCTION_CANDIDATE_ONLY",
        "scope": (
            "bounded-denominator recognition of the A-root polynomial; "
            "not an exact Hurwitz-algebra or cover certificate"
        ),
        "class_ids": class_ids,
        "invariant": arguments.invariant,
        "precision_bits": arguments.precision,
        "maximum_denominator": arguments.max_denominator,
        "coefficient_order": "ascending",
        "coefficients": records,
        "numerical_roots_by_class": {
            str(class_id): str(root) for class_id, root in zip(class_ids, roots)
        },
        "model_summaries": [
            {
                "class_id": int(payload["class_id"]),
                "terms": int(payload["terms"]),
                "samples": int(payload["samples"]),
                "precision_bits": int(payload["precision_bits"]),
                "maximum_basis_residual": max(
                    item["final_residual_norm"]
                    for item in payload["basis_refinement"]
                ),
                "quadric_balanced_jet_residual": payload["canonical_quadric"].get(
                    "balanced_jet_residual_norm"
                ),
                "cubic_balanced_jet_residual": payload.get("petri_cubic", {}).get(
                    "balanced_jet_residual_norm"
                ),
            }
            for payload in payloads
        ],
        "class_6_root": {
            "rational_part": class_6_real,
            "sqrt_minus_23_part": class_6_sqrt,
            "exact_linear_factor_check": exact_linear_factor_check,
        },
        "exact_factor_degrees": factor_degrees,
        "exact_sextic_irreducible_over_Q_sqrt_minus_23": sextic_irreducible,
        "exact_factors": [
            {
                "multiplicity": int(multiplicity),
                "coefficients_ascending": [
                    encode_quadratic(coefficient) for coefficient in factor.list()
                ],
            }
            for factor, multiplicity in exact_factors
        ],
        "maximum_root_residual": str(max(root_residuals)),
        "root_residuals": [str(value) for value in root_residuals],
    }
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    print(rendered, end="")
    if arguments.output:
        arguments.output.write_text(rendered)


if __name__ == "__main__":
    main()
