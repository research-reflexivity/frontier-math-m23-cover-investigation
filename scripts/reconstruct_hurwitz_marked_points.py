#!/usr/bin/env sage-python
"""Reconstruct the two totally ramified points on the canonical models.

The b-chart normalization gives b=[1:0:0:0].  Models computed with the
optional c chart give the opposite point directly from its constant terms.
This script reconstructs c over both components of the Hurwitz algebra and
requires exact substitution into the reconstructed quadric and cubic.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

from sage.all import ComplexField, PolynomialRing, QQ, RealField, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ALGEBRA = ROOT / "data" / "hurwitz_algebra_candidate.json"
DEFAULT_MODELS = ROOT / "data" / "hurwitz_canonical_models_candidate.json"


def rational(record):
    return QQ(record["numerator"]) / QQ(record["denominator"])


def encode_rational(value):
    return {"numerator": int(value.numerator()), "denominator": int(value.denominator())}


def decode_quadratic(record, field, generator):
    return field(rational(record["rational_part"])) + field(
        rational(record["sqrt_minus_23_part"])
    ) * generator


def bounded_rational(value, bound):
    return Fraction(value.str(truncate=False)).limit_denominator(bound)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("models", nargs=7, type=Path)
    parser.add_argument("--algebra", type=Path, default=DEFAULT_ALGEBRA)
    parser.add_argument("--canonical-models", type=Path, default=DEFAULT_MODELS)
    parser.add_argument("--precision", type=int, default=768)
    parser.add_argument("--sextic-low", type=int, default=10**16)
    parser.add_argument("--sextic-high", type=int, default=10**20)
    parser.add_argument("--degree-one-low", type=int, default=10**9)
    parser.add_argument("--degree-one-high", type=int, default=10**13)
    parser.add_argument("--simultaneous-low-digits", type=int, default=34)
    parser.add_argument("--simultaneous-high-digits", type=int, default=38)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()

    complex_field = ComplexField(arguments.precision)
    real_field = RealField(arguments.precision)
    algebra = json.loads(arguments.algebra.read_text())
    canonical = json.loads(arguments.canonical_models.read_text())

    base_polynomial_ring = PolynomialRing(QQ, "s")
    s = base_polynomial_ring.gen()
    base_field = QQ.extension(s**2 + 23, "sqrt_minus_23")
    sqrt_minus_23 = base_field.gen()
    hurwitz_polynomial_ring = PolynomialRing(base_field, "J")
    factors = [
        hurwitz_polynomial_ring(
            [
                decode_quadratic(coefficient, base_field, sqrt_minus_23)
                for coefficient in factor["coefficients_ascending"]
            ]
        )
        for factor in algebra["exact_factors"]
    ]
    sextic = next(factor for factor in factors if factor.degree() == 6)
    relative_field = base_field.extension(sextic, "j")
    j_relative = relative_field.gen()
    absolute_field = relative_field.absolute_field("a")
    _, to_absolute = absolute_field.structure()
    j_absolute = to_absolute(j_relative)
    sqrt_minus_23_absolute = to_absolute(relative_field(sqrt_minus_23))
    integral_basis = absolute_field.integral_basis()

    payloads = {}
    numerical_points = {}
    numerical_J = {}
    sources = []
    for path in arguments.models:
        raw = path.read_bytes()
        payload = json.loads(raw)
        class_id = int(payload["class_id"])
        point = payload.get("opposite_branch_point")
        if point is None:
            raise SystemExit(f"class {class_id} has no c-chart point")
        numerical_points[class_id] = [
            complex_field(value)
            for value in point["scale_free_coordinates_y0_equal_1"]
        ]
        monomials = [tuple(value) for value in payload["canonical_quadric"]["monomials"]]
        coefficients = {
            monomial: complex_field(value)
            for monomial, value in zip(
                monomials, payload["canonical_quadric"]["coefficients"]
            )
        }
        A = coefficients[(0, 3)]
        numerical_J[class_id] = coefficients[(1, 3)] / A**2
        payloads[class_id] = payload
        sources.append(
            {
                "class_id": class_id,
                "sha256": hashlib.sha256(raw).hexdigest(),
                "terms": payload["terms"],
                "samples": payload["samples"],
                "precision_bits": payload["precision_bits"],
                "patch_count": payload["patch_count"],
                "c_patch": point["patch"],
                "refinement_rounds": len(payload["basis_refinement"][0]["rounds"]),
                "maximum_basis_residual": max(
                    item["final_residual_norm"] for item in payload["basis_refinement"]
                ),
            }
        )
    if sorted(payloads) != list(range(1, 8)):
        raise SystemExit("inputs must contain class IDs 1 through 7")

    embeddings = [
        embedding
        for embedding in absolute_field.embeddings(complex_field)
        if embedding(sqrt_minus_23_absolute).imag() > 0
    ]
    unused = {1, 2, 3, 4, 5, 7}
    matched = []
    for embedding in embeddings:
        class_id = min(unused, key=lambda value: abs(embedding(j_absolute) - numerical_J[value]))
        unused.remove(class_id)
        matched.append((embedding, class_id))
    embedding_rows = []
    for embedding, _ in matched:
        embedding_rows.append([real_field(embedding(value).real()) for value in integral_basis])
        embedding_rows.append([real_field(embedding(value).imag()) for value in integral_basis])
    embedding_matrix = matrix(real_field, embedding_rows)

    def simultaneous_rationals(values, scale_digits):
        scale = ZZ(10) ** scale_digits
        rounded = [ZZ((real_field(scale) * value).round()) for value in values]
        rows = [[ZZ(1), *rounded]]
        for index in range(len(values)):
            row = [ZZ(0)] * (len(values) + 1)
            row[index + 1] = scale
            rows.append(row)
        reduced = matrix(ZZ, rows).LLL(delta=0.99)
        candidates = []
        for row in reduced.rows():
            denominator = abs(ZZ(row[0]))
            if denominator == 0:
                continue
            numerators = [ZZ((real_field(denominator) * value).round()) for value in values]
            common = denominator
            for numerator in numerators:
                common = common.gcd(abs(numerator))
            denominator //= common
            numerators = [numerator // common for numerator in numerators]
            residual = max(
                abs(real_field(denominator) * value - real_field(numerator))
                for value, numerator in zip(values, numerators)
            )
            candidates.append((residual, denominator, numerators))
        residual, denominator, numerators = min(candidates, key=lambda item: item[0])
        return (
            [Fraction(int(numerator), int(denominator)) for numerator in numerators],
            residual,
            denominator,
        )

    def reconstruct_degree_one(value):
        low = (
            bounded_rational(real_field(value.real()), arguments.degree_one_low),
            bounded_rational(real_field(value.imag()) / real_field(23).sqrt(), arguments.degree_one_low),
        )
        high = (
            bounded_rational(real_field(value.real()), arguments.degree_one_high),
            bounded_rational(real_field(value.imag()) / real_field(23).sqrt(), arguments.degree_one_high),
        )
        exact = base_field(QQ(high[0].numerator) / high[0].denominator) + base_field(
            QQ(high[1].numerator) / high[1].denominator
        ) * sqrt_minus_23
        numerical = complex_field(exact[0]) + complex_field(exact[1]) * complex_field(23).sqrt() * complex_field.gen()
        return exact, {
            "plateau": low == high,
            "rational_part": encode_rational(QQ(high[0].numerator) / high[0].denominator),
            "sqrt_minus_23_part": encode_rational(QQ(high[1].numerator) / high[1].denominator),
            "numerical_error": str(abs(value - numerical)),
        }

    def reconstruct_sextic(values):
        right_hand_side = []
        for _, class_id in matched:
            value = values[class_id]
            right_hand_side.extend([real_field(value.real()), real_field(value.imag())])
        coordinates = embedding_matrix.solve_right(vector(real_field, right_hand_side))
        low = [bounded_rational(value, arguments.sextic_low) for value in coordinates]
        high = [bounded_rational(value, arguments.sextic_high) for value in coordinates]
        simultaneous_low = simultaneous_rationals(coordinates, arguments.simultaneous_low_digits)
        simultaneous_high = simultaneous_rationals(coordinates, arguments.simultaneous_high_digits)
        coordinate_plateau = low == high
        simultaneous_plateau = simultaneous_low[0] == simultaneous_high[0]
        if simultaneous_plateau and not coordinate_plateau:
            selected = simultaneous_high[0]
            method = "common_denominator_lll"
        else:
            selected = high
            method = "coordinatewise_continued_fractions"
        exact = sum(
            QQ(value.numerator) / value.denominator * basis_value
            for value, basis_value in zip(selected, integral_basis)
        )
        errors = [abs(embedding(exact) - values[class_id]) for embedding, class_id in matched]
        return exact, {
            "plateau": coordinate_plateau or simultaneous_plateau,
            "recognition_method": method,
            "coordinatewise_plateau": coordinate_plateau,
            "simultaneous_plateau": simultaneous_plateau,
            "simultaneous_common_denominator": int(simultaneous_high[2]),
            "simultaneous_maximum_scaled_residual": str(simultaneous_high[1]),
            "integral_basis_coordinates": [
                encode_rational(QQ(value.numerator) / value.denominator) for value in selected
            ],
            "maximum_numerical_error": str(max(errors)),
        }

    degree_one_exact = []
    degree_one_records = []
    sextic_exact = []
    sextic_records = []
    for coordinate in range(4):
        exact, record = reconstruct_degree_one(numerical_points[6][coordinate])
        degree_one_exact.append(exact)
        degree_one_records.append(record)
        exact, record = reconstruct_sextic(
            {class_id: numerical_points[class_id][coordinate] for class_id in numerical_points}
        )
        sextic_exact.append(exact)
        sextic_records.append(record)

    def degree_one_coefficient(record):
        return decode_quadratic(record["degree_one_component"], base_field, sqrt_minus_23)

    def sextic_coefficient(record):
        return sum(
            rational(coordinate) * basis_value
            for coordinate, basis_value in zip(
                record["sextic_component"]["integral_basis_coordinates"], integral_basis
            )
        )

    quadric_monomials = [tuple(value) for value in canonical["quadric"]["monomials"]]
    cubic_monomials = [tuple(value) for value in canonical["petri_cubic"]["monomials"]]
    degree_one_quadric = [degree_one_coefficient(value) for value in canonical["quadric"]["coefficients"]]
    sextic_quadric = [sextic_coefficient(value) for value in canonical["quadric"]["coefficients"]]
    degree_one_cubic = [degree_one_coefficient(value) for value in canonical["petri_cubic"]["coefficients"]]
    sextic_cubic = [sextic_coefficient(value) for value in canonical["petri_cubic"]["coefficients"]]

    def evaluate_quadric(coefficients, point):
        return sum(
            coefficient * point[first] * point[second]
            for (first, second), coefficient in zip(quadric_monomials, coefficients)
        )

    def evaluate_cubic(coefficients, point):
        return sum(
            coefficient * point[first] * point[second] * point[third]
            for (first, second, third), coefficient in zip(cubic_monomials, coefficients)
        )

    # The third nontrivial coordinate has substantially larger height than
    # the first two.  Once y_1 and y_2 are recognized, recover y_3 exactly as
    # the common root of Q(1,y_1,y_2,Z) and C(1,y_1,y_2,Z).
    point_polynomial_ring = PolynomialRing(absolute_field, "Z")
    Z = point_polynomial_ring.gen()
    partial_point = [absolute_field(1), sextic_exact[1], sextic_exact[2], Z]
    quadric_in_Z = point_polynomial_ring(
        evaluate_quadric(sextic_quadric, partial_point)
    )
    cubic_in_Z = point_polynomial_ring(evaluate_cubic(sextic_cubic, partial_point))
    common_factor = quadric_in_Z.gcd(cubic_in_Z).monic()
    if common_factor.degree() != 1:
        raise SystemExit(
            "recognized c_1,c_2 do not give a unique exact common curve point"
        )
    derived_c3 = -common_factor[0]
    preliminary_c3 = sextic_records[3]
    sextic_exact[3] = derived_c3
    degree = absolute_field.degree()
    basis_matrix = matrix(
        QQ,
        degree,
        degree,
        lambda row, column: (
            integral_basis[column].list()[row]
            if row < len(integral_basis[column].list())
            else 0
        ),
    )
    derived_power_coordinates = derived_c3.list() + [QQ(0)] * (
        degree - len(derived_c3.list())
    )
    derived_integral_coordinates = basis_matrix.solve_right(
        vector(QQ, derived_power_coordinates)
    )
    derived_errors = [
        abs(embedding(derived_c3) - numerical_points[class_id][3])
        for embedding, class_id in matched
    ]
    sextic_records[3] = {
        "plateau": True,
        "recognition_method": "exact_curve_gcd_from_c1_c2",
        "coordinatewise_plateau": False,
        "simultaneous_plateau": False,
        "integral_basis_coordinates": [
            encode_rational(value) for value in derived_integral_coordinates
        ],
        "maximum_numerical_error": str(max(derived_errors)),
        "exact_common_factor_coefficients_ascending": [
            str(value) for value in common_factor.list()
        ],
        "preliminary_numerical_reconstruction": preliminary_c3,
    }

    substitutions = {
        "degree_one_quadric": str(evaluate_quadric(degree_one_quadric, degree_one_exact)),
        "degree_one_cubic": str(evaluate_cubic(degree_one_cubic, degree_one_exact)),
        "sextic_quadric": str(evaluate_quadric(sextic_quadric, sextic_exact)),
        "sextic_cubic": str(evaluate_cubic(sextic_cubic, sextic_exact)),
    }
    substitution_failures = sum(value != "0" for value in substitutions.values())
    failures = sum(not value["plateau"] for value in degree_one_records + sextic_records)
    result = {
        "status": "EXACT_MARKED_POINT_RECONSTRUCTION_CANDIDATE_ONLY",
        "scope": "b and c on the reconstructed canonical curves; the divisor relation 23(b-c) and Belyi functions remain outstanding",
        "b_coordinates": [1, 0, 0, 0],
        "c_coordinates": {
            "degree_one_component": degree_one_records,
            "sextic_component": sextic_records,
        },
        "plateau_failure_count": failures,
        "exact_substitution": substitutions,
        "exact_substitution_failure_count": substitution_failures,
        "embedding_matches": [
            {"class_id": class_id, "J_error": str(abs(embedding(j_absolute) - numerical_J[class_id]))}
            for embedding, class_id in matched
        ],
        "sources": sorted(sources, key=lambda item: item["class_id"]),
    }
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    print(rendered, end="")
    if arguments.output:
        arguments.output.write_text(rendered)
    if failures:
        raise SystemExit(f"{failures} marked-point coordinates lack a plateau")
    if substitution_failures:
        raise SystemExit(
            f"{substitution_failures} exact marked-point substitutions failed"
        )


if __name__ == "__main__":
    main()
