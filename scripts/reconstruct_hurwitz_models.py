#!/usr/bin/env sage-python
"""Reconstruct scale-free canonical models in the Hurwitz number field.

The six remaining embeddings are reconstructed in an integral basis of
the degree-12 absolute field underlying the irreducible sextic component.
The embedding on the degree-one component is reconstructed separately over
Q(sqrt(-23)).
Output is still a numerical reconstruction candidate until an exact cover
map is recovered and substituted.
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


def rational(record):
    return QQ(record["numerator"]) / QQ(record["denominator"])


def encode_rational(value) -> dict[str, int]:
    return {
        "numerator": int(value.numerator()),
        "denominator": int(value.denominator()),
    }


def encode_power_basis(value, degree: int) -> list[dict[str, int]]:
    coefficients = value.list()
    coefficients += [QQ(0)] * (degree - len(coefficients))
    return [encode_rational(coefficient) for coefficient in coefficients]


def bounded_rational(value, bound: int) -> Fraction:
    return Fraction(value.str(truncate=False)).limit_denominator(bound)


def decode_quadratic(record, field, generator):
    return field(rational(record["rational_part"])) + field(
        rational(record["sqrt_minus_23_part"])
    ) * generator


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("models", nargs=7, type=Path)
    parser.add_argument(
        "--cross-check-models",
        nargs=6,
        type=Path,
        help="independent models for classes 1,2,3,4,5,7",
    )
    parser.add_argument("--algebra", type=Path, default=DEFAULT_ALGEBRA)
    parser.add_argument("--precision", type=int, default=1536)
    parser.add_argument("--plateau-low", type=int, default=10**40)
    parser.add_argument("--plateau-high", type=int, default=10**44)
    parser.add_argument("--degree-one-low", type=int, default=10**20)
    parser.add_argument("--degree-one-high", type=int, default=10**30)
    parser.add_argument("--simultaneous-scale-low-digits", type=int, default=96)
    parser.add_argument("--simultaneous-scale-high-digits", type=int, default=100)
    parser.add_argument(
        "--cubic-normalized-monomial",
        default="2,2,2",
        help="comma-separated cubic monomial whose scale-free coefficient is set to one",
    )
    parser.add_argument(
        "--cubic-gauge-monomials",
        default="0,0,2;0,1,2;0,2,2;0,2,3",
        help="semicolon-separated cubic monomials set to zero modulo Q times a linear form",
    )
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    cubic_normalized_monomial = tuple(
        int(value) for value in arguments.cubic_normalized_monomial.split(",")
    )
    if len(cubic_normalized_monomial) != 3:
        raise SystemExit("the cubic normalization must have three indices")
    cubic_gauge_monomials = [
        tuple(int(value) for value in monomial.split(","))
        for monomial in arguments.cubic_gauge_monomials.split(";")
    ]
    if len(cubic_gauge_monomials) != 4 or any(
        len(monomial) != 3 for monomial in cubic_gauge_monomials
    ):
        raise SystemExit("the cubic gauge must contain four cubic monomials")

    complex_field = ComplexField(arguments.precision)
    real_field = RealField(arguments.precision)
    algebra = json.loads(arguments.algebra.read_text())
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
    models = {}
    for path in arguments.models:
        payload = json.loads(path.read_text())
        class_id = int(payload["class_id"])
        payloads[class_id] = payload
        quadric = {
            tuple(monomial): complex_field(coefficient)
            for monomial, coefficient in zip(
                payload["canonical_quadric"]["monomials"],
                payload["canonical_quadric"]["coefficients"],
            )
        }
        cubic = {
            tuple(monomial): complex_field(coefficient)
            for monomial, coefficient in zip(
                payload["petri_cubic"]["monomials"],
                payload["petri_cubic"]["coefficients"],
            )
        }
        A = quadric[(0, 3)]
        models[class_id] = {
            "A": A,
            "J": quadric[(1, 3)] / A**2,
            "quadric": quadric,
            "cubic": cubic,
        }
    if sorted(models) != list(range(1, 8)):
        raise SystemExit("the model inputs must contain class IDs 1 through 7")

    quadric_monomials = list(models[1]["quadric"])
    cubic_monomials = list(models[1]["cubic"])
    cubic_index = {monomial: index for index, monomial in enumerate(cubic_monomials)}
    for model in models.values():
        A = model["A"]
        scaled_quadric = {
            monomial: coefficient * A ** (2 - sum(monomial))
            for monomial, coefficient in model["quadric"].items()
        }
        scaled_cubic = {
            monomial: coefficient * A ** (6 - sum(monomial))
            for monomial, coefficient in model["cubic"].items()
        }
        multipliers = [[complex_field(0) for _ in range(4)] for _ in cubic_monomials]
        for variable in range(4):
            for monomial, coefficient in scaled_quadric.items():
                triple = tuple(sorted((*monomial, variable)))
                multipliers[cubic_index[triple]][variable] += coefficient
        gauge_matrix = matrix(
            complex_field,
            [
                multipliers[cubic_index[monomial]]
                for monomial in cubic_gauge_monomials
            ],
        )
        if abs(gauge_matrix.det()) < real_field("1e-30"):
            raise SystemExit("the requested cubic gauge is numerically singular")
        correction = gauge_matrix.solve_right(
            vector(
                complex_field,
                [-scaled_cubic[monomial] for monomial in cubic_gauge_monomials],
            )
        )
        gauged = {
            monomial: scaled_cubic[monomial]
            + sum(
                multipliers[cubic_index[monomial]][variable] * correction[variable]
                for variable in range(4)
            )
            for monomial in cubic_monomials
        }
        normalizer = gauged[cubic_normalized_monomial]
        if abs(normalizer) < real_field("1e-20"):
            raise SystemExit("the requested cubic normalization is numerically zero")
        model["scaled_cubic"] = {
            monomial: coefficient / normalizer
            for monomial, coefficient in gauged.items()
        }
        symmetric_quadric = matrix(complex_field, 4, 4)
        for (first, second), coefficient in scaled_quadric.items():
            if first == second:
                symmetric_quadric[first, second] = coefficient
            else:
                symmetric_quadric[first, second] = coefficient / 2
                symmetric_quadric[second, first] = coefficient / 2
        model["quadric_determinant"] = symmetric_quadric.det()

    embeddings = [
        embedding
        for embedding in absolute_field.embeddings(complex_field)
        if embedding(sqrt_minus_23_absolute).imag() > 0
    ]
    unused = {1, 2, 3, 4, 5, 7}
    matched = []
    for embedding in embeddings:
        class_id = min(
            unused,
            key=lambda value: abs(embedding(j_absolute) - models[value]["J"]),
        )
        unused.remove(class_id)
        matched.append((embedding, class_id))

    embedding_rows = []
    for embedding, _ in matched:
        embedding_rows.append(
            [real_field(embedding(value).real()) for value in integral_basis]
        )
        embedding_rows.append(
            [real_field(embedding(value).imag()) for value in integral_basis]
        )
    embedding_matrix = matrix(real_field, embedding_rows)

    def simultaneous_rationals(values, scale_digits):
        """Recognize rational coordinates with one denominator using LLL."""

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
            numerators = [
                ZZ((real_field(denominator) * value).round()) for value in values
            ]
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
        if not candidates:
            raise ArithmeticError("LLL produced no nonzero candidate denominator")
        residual, denominator, numerators = min(candidates, key=lambda item: item[0])
        return (
            [Fraction(int(numerator), int(denominator)) for numerator in numerators],
            residual,
            denominator,
        )

    def reconstruct_degree_one(value):
        low = (
            bounded_rational(real_field(value.real()), arguments.degree_one_low),
            bounded_rational(
                real_field(value.imag()) / real_field(23).sqrt(),
                arguments.degree_one_low,
            ),
        )
        high = (
            bounded_rational(real_field(value.real()), arguments.degree_one_high),
            bounded_rational(
                real_field(value.imag()) / real_field(23).sqrt(),
                arguments.degree_one_high,
            ),
        )
        exact = base_field(QQ(high[0].numerator) / high[0].denominator) + base_field(
            QQ(high[1].numerator) / high[1].denominator
        ) * sqrt_minus_23
        numerical = complex_field(exact[0]) + complex_field(exact[1]) * complex_field(
            23
        ).sqrt() * complex_field.gen()
        return {
            "plateau": low == high,
            "rational_part": encode_rational(
                QQ(high[0].numerator) / high[0].denominator
            ),
            "sqrt_minus_23_part": encode_rational(
                QQ(high[1].numerator) / high[1].denominator
            ),
            "numerical_error": str(abs(value - numerical)),
        }

    def reconstruct_sextic(values):
        right_hand_side = []
        for _, class_id in matched:
            value = values[class_id]
            right_hand_side.extend(
                [real_field(value.real()), real_field(value.imag())]
            )
        numerical_coordinates = embedding_matrix.solve_right(
            vector(real_field, right_hand_side)
        )
        low = [
            bounded_rational(value, arguments.plateau_low)
            for value in numerical_coordinates
        ]
        high = [
            bounded_rational(value, arguments.plateau_high)
            for value in numerical_coordinates
        ]
        simultaneous_low = simultaneous_rationals(
            numerical_coordinates, arguments.simultaneous_scale_low_digits
        )
        simultaneous_high = simultaneous_rationals(
            numerical_coordinates, arguments.simultaneous_scale_high_digits
        )
        coordinate_plateau = low == high
        simultaneous_plateau = simultaneous_low[0] == simultaneous_high[0]
        if simultaneous_plateau and not coordinate_plateau:
            selected = simultaneous_high[0]
            recognition_method = "common_denominator_lll"
        else:
            selected = high
            recognition_method = "coordinatewise_continued_fractions"
        exact = sum(
            QQ(coordinate.numerator) / coordinate.denominator * basis_value
            for coordinate, basis_value in zip(selected, integral_basis)
        )
        errors = [
            abs(embedding(exact) - values[class_id])
            for embedding, class_id in matched
        ]
        return {
            "plateau": coordinate_plateau or simultaneous_plateau,
            "recognition_method": recognition_method,
            "coordinatewise_plateau": coordinate_plateau,
            "simultaneous_plateau": simultaneous_plateau,
            "simultaneous_common_denominator": int(simultaneous_high[2]),
            "simultaneous_maximum_scaled_residual": str(simultaneous_high[1]),
            "integral_basis_coordinates": [
                encode_rational(QQ(value.numerator) / value.denominator)
                for value in selected
            ],
            "maximum_numerical_error": str(max(errors)),
            "numerical_errors_by_class": {
                str(class_id): str(error)
                for (_, class_id), error in zip(matched, errors)
            },
        }

    def reconstruct_coefficient(values):
        return {
            "degree_one_component": reconstruct_degree_one(values[6]),
            "sextic_component": reconstruct_sextic(values),
        }

    quadric_records = []
    for monomial in quadric_monomials:
        values = {
            class_id: model["quadric"][monomial]
            * model["A"] ** (2 - sum(monomial))
            for class_id, model in models.items()
        }
        quadric_records.append(reconstruct_coefficient(values))
    cubic_records = []
    for monomial in cubic_monomials:
        values = {
            class_id: model["scaled_cubic"][monomial]
            for class_id, model in models.items()
        }
        cubic_records.append(reconstruct_coefficient(values))

    quadric_determinant_record = reconstruct_coefficient(
        {
            class_id: model["quadric_determinant"]
            for class_id, model in models.items()
        }
    )

    # The determinant is useful diagnostically, but it is not an additional
    # coefficient needed to define the complete intersection, so do not let
    # its separate height control model acceptance.
    all_records = quadric_records + cubic_records
    plateau_failures = sum(
        not component["plateau"]
        for record in all_records
        for component in (record["degree_one_component"], record["sextic_component"])
    )
    maximum_error = max(
        real_field(record["sextic_component"]["maximum_numerical_error"])
        for record in all_records
    )
    independent_cross_check = None
    if arguments.cross_check_models:
        cross_models = {}
        cross_sources = []
        for path in arguments.cross_check_models:
            raw = path.read_bytes()
            payload = json.loads(raw)
            class_id = int(payload["class_id"])
            quadric = {
                tuple(monomial): complex_field(coefficient)
                for monomial, coefficient in zip(
                    payload["canonical_quadric"]["monomials"],
                    payload["canonical_quadric"]["coefficients"],
                )
            }
            cubic = {
                tuple(monomial): complex_field(coefficient)
                for monomial, coefficient in zip(
                    payload["petri_cubic"]["monomials"],
                    payload["petri_cubic"]["coefficients"],
                )
            }
            A = quadric[(0, 3)]
            scaled_quadric = {
                monomial: coefficient * A ** (2 - sum(monomial))
                for monomial, coefficient in quadric.items()
            }
            scaled_cubic = {
                monomial: coefficient * A ** (6 - sum(monomial))
                for monomial, coefficient in cubic.items()
            }
            multipliers = [
                [complex_field(0) for _ in range(4)] for _ in cubic_monomials
            ]
            for variable in range(4):
                for monomial, coefficient in scaled_quadric.items():
                    triple = tuple(sorted((*monomial, variable)))
                    multipliers[cubic_index[triple]][variable] += coefficient
            gauge_matrix = matrix(
                complex_field,
                [
                    multipliers[cubic_index[monomial]]
                    for monomial in cubic_gauge_monomials
                ],
            )
            correction = gauge_matrix.solve_right(
                vector(
                    complex_field,
                    [
                        -scaled_cubic[monomial]
                        for monomial in cubic_gauge_monomials
                    ],
                )
            )
            gauged = {
                monomial: scaled_cubic[monomial]
                + sum(
                    multipliers[cubic_index[monomial]][variable]
                    * correction[variable]
                    for variable in range(4)
                )
                for monomial in cubic_monomials
            }
            normalizer = gauged[cubic_normalized_monomial]
            cross_models[class_id] = {
                "quadric": scaled_quadric,
                "cubic": {
                    monomial: coefficient / normalizer
                    for monomial, coefficient in gauged.items()
                },
            }
            cross_sources.append(
                {
                    "class_id": class_id,
                    "sha256": hashlib.sha256(raw).hexdigest(),
                    "terms": payload["terms"],
                    "samples": payload["samples"],
                    "precision_bits": payload["precision_bits"],
                    "refinement_rounds": len(
                        payload["basis_refinement"][0]["rounds"]
                    ),
                }
            )
        if sorted(cross_models) != [1, 2, 3, 4, 5, 7]:
            raise SystemExit("cross-check models must have class IDs 1,2,3,4,5,7")

        def exact_sextic_value(record):
            return sum(
                rational(coordinate) * basis_value
                for coordinate, basis_value in zip(
                    record["sextic_component"]["integral_basis_coordinates"],
                    integral_basis,
                )
            )

        cross_check_coefficients = []
        for kind, monomials, records in (
            ("quadric", quadric_monomials, quadric_records),
            ("petri_cubic", cubic_monomials, cubic_records),
        ):
            for monomial, record in zip(monomials, records):
                exact = exact_sextic_value(record)
                errors = [
                    abs(
                        embedding(exact)
                        - cross_models[class_id][
                            "quadric" if kind == "quadric" else "cubic"
                        ][monomial]
                    )
                    for embedding, class_id in matched
                ]
                cross_check_coefficients.append(
                    {
                        "kind": kind,
                        "monomial": list(monomial),
                        "maximum_error": str(max(errors)),
                    }
                )
        independent_cross_check = {
            "description": (
                "exact reconstructed sextic-component coefficients evaluated "
                "against an independently truncated N=480, five-round batch"
            ),
            "sources": sorted(cross_sources, key=lambda item: item["class_id"]),
            "coefficients": cross_check_coefficients,
            "maximum_error": str(
                max(
                    real_field(item["maximum_error"])
                    for item in cross_check_coefficients
                )
            ),
        }
    result = {
        "status": "NUMERICAL_RECONSTRUCTION_CANDIDATE_ONLY",
        "scope": (
            "scale-free canonical quadrics and fixed-gauge Petri cubics over "
            "the degree-one component and the absolute sextic field; exact "
            "degree-23 maps and substitution certificates remain outstanding"
        ),
        "coordinate_change": "y_i = A^i x_i",
        "quadric_scaling": "c_ij -> c_ij A^(2-i-j)",
        "cubic_scaling": "c_ijk -> c_ijk A^(6-i-j-k)",
        "plateau_bounds": {
            "sextic_low": arguments.plateau_low,
            "sextic_high": arguments.plateau_high,
            "degree_one_low": arguments.degree_one_low,
            "degree_one_high": arguments.degree_one_high,
            "simultaneous_scale_low_digits": (
                arguments.simultaneous_scale_low_digits
            ),
            "simultaneous_scale_high_digits": (
                arguments.simultaneous_scale_high_digits
            ),
        },
        "plateau_failure_count": plateau_failures,
        "maximum_sextic_embedding_error": str(maximum_error),
        "independent_cross_check": independent_cross_check,
        "absolute_field_polynomial_coefficients_ascending": [
            encode_rational(coefficient)
            for coefficient in absolute_field.polynomial().list()
        ],
        "integral_basis_power_coefficients": [
            encode_power_basis(value, absolute_field.degree())
            for value in integral_basis
        ],
        "embedding_matches": [
            {
                "class_id": class_id,
                "J_error": str(abs(embedding(j_absolute) - models[class_id]["J"])),
            }
            for embedding, class_id in matched
        ],
        "quadric": {
            "monomials": [list(monomial) for monomial in quadric_monomials],
            "coefficients": quadric_records,
            "symmetric_matrix_determinant": quadric_determinant_record,
        },
        "petri_cubic": {
            "monomials": [list(monomial) for monomial in cubic_monomials],
            "gauge_zero_monomials": [
                list(monomial) for monomial in cubic_gauge_monomials
            ],
            "normalized_monomial": list(cubic_normalized_monomial),
            "coefficients": cubic_records,
        },
        "model_summaries": [
            {
                "class_id": class_id,
                "terms": payloads[class_id]["terms"],
                "samples": payloads[class_id]["samples"],
                "precision_bits": payloads[class_id]["precision_bits"],
                "maximum_basis_residual": max(
                    item["final_residual_norm"]
                    for item in payloads[class_id]["basis_refinement"]
                ),
            }
            for class_id in range(1, 8)
        ],
    }
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    print(rendered, end="")
    if arguments.output:
        arguments.output.write_text(rendered)
    if plateau_failures:
        raise SystemExit(f"{plateau_failures} rational reconstructions lack a plateau")


if __name__ == "__main__":
    main()
