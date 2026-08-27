#!/usr/bin/env sage-python
"""Verify the exact internal structure of the reconstructed canonical models.

This checks the degree-12 coefficient field, all rational reconstructions,
the fixed coordinate gauge, and nonsingularity of the unique canonical
quadric on both Hurwitz-algebra components.  It does not prove that the
recognized algebraic coefficients equal the exact automorphic invariants;
the data therefore remains explicitly labelled a candidate.
"""

from __future__ import annotations

import json
from pathlib import Path

from sage.all import ComplexField, PolynomialRing, QQ, matrix


ROOT = Path(__file__).resolve().parents[1]
ALGEBRA_PATH = ROOT / "data" / "hurwitz_algebra_candidate.json"
MODELS_PATH = ROOT / "data" / "hurwitz_canonical_models_candidate.json"


def rational(record):
    return QQ(record["numerator"]) / QQ(record["denominator"])


def quadratic_element(record, field, generator):
    return field(rational(record["rational_part"])) + field(
        rational(record["sqrt_minus_23_part"])
    ) * generator


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def main() -> None:
    algebra = json.loads(ALGEBRA_PATH.read_text())
    candidate = json.loads(MODELS_PATH.read_text())

    base_polynomial_ring = PolynomialRing(QQ, "s")
    s = base_polynomial_ring.gen()
    base_field = QQ.extension(s**2 + 23, "sqrt_minus_23")
    sqrt_minus_23 = base_field.gen()
    hurwitz_polynomial_ring = PolynomialRing(base_field, "J")
    factors = [
        hurwitz_polynomial_ring(
            [
                quadratic_element(coefficient, base_field, sqrt_minus_23)
                for coefficient in factor["coefficients_ascending"]
            ]
        )
        for factor in algebra["exact_factors"]
    ]
    sextic = next(factor for factor in factors if factor.degree() == 6)
    relative_field = base_field.extension(sextic, "j")
    absolute_field = relative_field.absolute_field("a")
    integral_basis = absolute_field.integral_basis()

    recorded_absolute_polynomial = [
        rational(coefficient)
        for coefficient in candidate[
            "absolute_field_polynomial_coefficients_ascending"
        ]
    ]
    require(
        absolute_field.polynomial().list() == recorded_absolute_polynomial,
        "absolute coefficient field polynomial changed",
    )
    require(len(integral_basis) == 12, "absolute field must have degree 12")
    recorded_integral_basis = [
        absolute_field([rational(coefficient) for coefficient in basis_value])
        for basis_value in candidate["integral_basis_power_coefficients"]
    ]
    require(
        recorded_integral_basis == integral_basis,
        "recorded integral basis changed",
    )
    require(candidate["plateau_failure_count"] == 0, "a coefficient lacks a plateau")
    require(
        candidate["status"] == "NUMERICAL_RECONSTRUCTION_CANDIDATE_ONLY",
        "candidate status must not overclaim an exact cover certificate",
    )

    def degree_one_element(record):
        return quadratic_element(record["degree_one_component"], base_field, sqrt_minus_23)

    def sextic_element(record):
        coordinates = record["sextic_component"]["integral_basis_coordinates"]
        return sum(
            rational(coordinate) * basis_value
            for coordinate, basis_value in zip(coordinates, integral_basis)
        )

    quadric_monomials = [tuple(value) for value in candidate["quadric"]["monomials"]]
    cubic_monomials = [
        tuple(value) for value in candidate["petri_cubic"]["monomials"]
    ]
    require(len(quadric_monomials) == 10, "quadric must have ten monomials")
    require(len(cubic_monomials) == 20, "cubic must have twenty monomials")

    all_records = (
        candidate["quadric"]["coefficients"]
        + candidate["petri_cubic"]["coefficients"]
    )
    for record in all_records:
        require(record["degree_one_component"]["plateau"], "degree-one coefficient failed")
        require(record["sextic_component"]["plateau"], "sextic coefficient failed")

    def coefficient_dictionary(section, decoder):
        return {
            tuple(monomial): decoder(record)
            for monomial, record in zip(
                section["monomials"], section["coefficients"]
            )
        }

    degree_one_quadric = coefficient_dictionary(candidate["quadric"], degree_one_element)
    sextic_quadric = coefficient_dictionary(candidate["quadric"], sextic_element)
    degree_one_cubic = coefficient_dictionary(candidate["petri_cubic"], degree_one_element)
    sextic_cubic = coefficient_dictionary(candidate["petri_cubic"], sextic_element)

    fixed_quadric = {
        (0, 0): 0,
        (0, 1): 0,
        (0, 2): 1,
        (0, 3): 1,
        (1, 1): -1,
        (1, 2): -1,
    }
    fixed_cubic = {
        (0, 0, 0): 0,
        (0, 0, 1): 0,
        (0, 0, 2): 0,
        (0, 1, 1): 0,
        (0, 1, 2): 0,
        (0, 2, 2): 0,
        (0, 2, 3): 0,
        (2, 2, 2): 1,
    }
    for coefficients in (degree_one_quadric, sextic_quadric):
        for monomial, expected in fixed_quadric.items():
            require(coefficients[monomial] == expected, f"bad quadric gauge at {monomial}")
    for coefficients in (degree_one_cubic, sextic_cubic):
        for monomial, expected in fixed_cubic.items():
            require(coefficients[monomial] == expected, f"bad cubic gauge at {monomial}")

    def symmetric_quadric(coefficients, field):
        result = matrix(field, 4, 4)
        for (first, second), coefficient in coefficients.items():
            if first == second:
                result[first, second] = coefficient
            else:
                result[first, second] = coefficient / 2
                result[second, first] = coefficient / 2
        return result

    require(
        symmetric_quadric(degree_one_quadric, base_field).det() != 0,
        "degree-one canonical quadric is singular",
    )
    require(
        symmetric_quadric(sextic_quadric, absolute_field).det() != 0,
        "sextic-component canonical quadric is singular",
    )

    def verify_smooth_complete_intersection(quadric, cubic, field, label):
        coordinate_ring = PolynomialRing(field, 4, names=("x0", "x1", "x2", "x3"))
        coordinates = coordinate_ring.gens()
        quadric_polynomial = sum(
            field(coefficient) * coordinates[first] * coordinates[second]
            for (first, second), coefficient in quadric.items()
        )
        cubic_polynomial = sum(
            field(coefficient)
            * coordinates[first]
            * coordinates[second]
            * coordinates[third]
            for (first, second, third), coefficient in cubic.items()
        )
        quadric_gradient = [quadric_polynomial.derivative(value) for value in coordinates]
        cubic_gradient = [cubic_polynomial.derivative(value) for value in coordinates]
        jacobian_minors = [
            quadric_gradient[first] * cubic_gradient[second]
            - quadric_gradient[second] * cubic_gradient[first]
            for first in range(4)
            for second in range(first + 1, 4)
        ]
        for patch_index, coordinate in enumerate(coordinates):
            singular_patch = coordinate_ring.ideal(
                [
                    quadric_polynomial,
                    cubic_polynomial,
                    coordinate - 1,
                    *jacobian_minors,
                ]
            )
            require(
                coordinate_ring(1) in singular_patch,
                f"{label} curve has a singular point in patch {patch_index}",
            )

    verify_smooth_complete_intersection(
        degree_one_quadric, degree_one_cubic, base_field, "degree-one-component"
    )
    verify_smooth_complete_intersection(
        sextic_quadric, sextic_cubic, absolute_field, "sextic-component"
    )

    q33_index = quadric_monomials.index((3, 3))
    c333_index = cubic_monomials.index((3, 3, 3))
    for record in (
        candidate["quadric"]["coefficients"][q33_index],
        candidate["petri_cubic"]["coefficients"][c333_index],
    ):
        sextic_record = record["sextic_component"]
        require(sextic_record["simultaneous_plateau"], "LLL plateau missing")
        require(
            sextic_record["recognition_method"] == "common_denominator_lll",
            "hard coefficient was not selected from simultaneous reconstruction",
        )

    require(
        float(candidate["maximum_sextic_embedding_error"]) < 2e-105,
        "reconstruction embedding error is too large",
    )
    cross_check = candidate["independent_cross_check"]
    require(cross_check is not None, "independent N=480 cross-check is missing")
    require(len(cross_check["coefficients"]) == 30, "cross-check must cover 30 coefficients")
    require(
        float(cross_check["maximum_error"]) < 1e-89,
        "independent N=480 coefficient discrepancy is too large",
    )
    cross_sources = {item["class_id"]: item for item in cross_check["sources"]}
    require(sorted(cross_sources) == [1, 2, 3, 4, 5, 7], "bad cross-check classes")
    for class_id, source in cross_sources.items():
        require(source["terms"] == 480, f"class {class_id} cross-check has wrong N")
        require(source["samples"] == 1280, f"class {class_id} cross-check has wrong Q")
        require(source["precision_bits"] == 1024, "cross-check precision changed")
        require(source["refinement_rounds"] == 5, "cross-check round count changed")
        require(len(source["sha256"]) == 64, "cross-check source digest is malformed")
    summaries = {item["class_id"]: item for item in candidate["model_summaries"]}
    require(sorted(summaries) == list(range(1, 8)), "missing numerical model summary")
    for class_id in (1, 2, 3, 4, 5, 7):
        require(
            summaries[class_id]["maximum_basis_residual"] < 2e-108,
            f"class {class_id} finite-equation residual is too large",
        )
    require(
        summaries[6]["maximum_basis_residual"] < 2e-77,
        "class-6 finite-equation residual is too large",
    )

    print("PASS exact degree-12 coefficient field and integral-basis data")
    print("PASS all 30 canonical coefficients have stable rational reconstructions")
    print("PASS fixed gauge and exact smoothness of both canonical complete intersections")
    print("PASS finite-equation metadata and independent N=480 cross-check")
    print(
        "SCOPE reconstruction checks alone do not identify the Hurwitz algebra; "
        "the independent branch-cycle certificate does"
    )


if __name__ == "__main__":
    main()
