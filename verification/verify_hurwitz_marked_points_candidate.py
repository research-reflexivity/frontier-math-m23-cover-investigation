#!/usr/bin/env sage-python
"""Verify the reconstructed b and c ramification points on both components."""

from __future__ import annotations

import json
from pathlib import Path

from sage.all import PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[1]


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
    algebra = json.loads((ROOT / "data/hurwitz_algebra_candidate.json").read_text())
    canonical = json.loads(
        (ROOT / "data/hurwitz_canonical_models_candidate.json").read_text()
    )
    marked = json.loads(
        (ROOT / "data/hurwitz_marked_points_candidate.json").read_text()
    )

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

    def degree_one_curve_coefficient(record):
        return quadratic_element(record["degree_one_component"], base_field, sqrt_minus_23)

    def sextic_curve_coefficient(record):
        return sum(
            rational(coordinate) * basis_value
            for coordinate, basis_value in zip(
                record["sextic_component"]["integral_basis_coordinates"],
                integral_basis,
            )
        )

    def sextic_point_coordinate(record):
        return sum(
            rational(coordinate) * basis_value
            for coordinate, basis_value in zip(
                record["integral_basis_coordinates"], integral_basis
            )
        )

    b_degree_one = [base_field(value) for value in marked["b_coordinates"]]
    b_sextic = [absolute_field(value) for value in marked["b_coordinates"]]
    c_degree_one = [
        quadratic_element(record, base_field, sqrt_minus_23)
        for record in marked["c_coordinates"]["degree_one_component"]
    ]
    c_sextic = [
        sextic_point_coordinate(record)
        for record in marked["c_coordinates"]["sextic_component"]
    ]

    quadric_monomials = [tuple(value) for value in canonical["quadric"]["monomials"]]
    cubic_monomials = [tuple(value) for value in canonical["petri_cubic"]["monomials"]]
    degree_one_quadric = [
        degree_one_curve_coefficient(record)
        for record in canonical["quadric"]["coefficients"]
    ]
    degree_one_cubic = [
        degree_one_curve_coefficient(record)
        for record in canonical["petri_cubic"]["coefficients"]
    ]
    sextic_quadric = [
        sextic_curve_coefficient(record)
        for record in canonical["quadric"]["coefficients"]
    ]
    sextic_cubic = [
        sextic_curve_coefficient(record)
        for record in canonical["petri_cubic"]["coefficients"]
    ]

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

    require(marked["plateau_failure_count"] == 0, "a point coordinate lacks recognition")
    require(marked["exact_substitution_failure_count"] == 0, "stored substitution failed")
    for record in marked["c_coordinates"]["degree_one_component"]:
        require(record["plateau"], "degree-one c coordinate has no plateau")
    for record in marked["c_coordinates"]["sextic_component"]:
        require(record["plateau"], "sextic c coordinate was not recognized")

    for label, quadric, cubic, b_point, c_point in (
        ("degree-one", degree_one_quadric, degree_one_cubic, b_degree_one, c_degree_one),
        ("sextic", sextic_quadric, sextic_cubic, b_sextic, c_sextic),
    ):
        require(evaluate_quadric(quadric, b_point) == 0, f"{label} b not on Q")
        require(evaluate_cubic(cubic, b_point) == 0, f"{label} b not on C")
        require(evaluate_quadric(quadric, c_point) == 0, f"{label} c not on Q")
        require(evaluate_cubic(cubic, c_point) == 0, f"{label} c not on C")
        require(b_point != c_point, f"{label} marked points coincide")

    point_ring = PolynomialRing(absolute_field, "Z")
    Z = point_ring.gen()
    partial_point = [absolute_field(1), c_sextic[1], c_sextic[2], Z]
    quadric_in_Z = point_ring(evaluate_quadric(sextic_quadric, partial_point))
    cubic_in_Z = point_ring(evaluate_cubic(sextic_cubic, partial_point))
    common_factor = quadric_in_Z.gcd(cubic_in_Z).monic()
    require(common_factor.degree() == 1, "c_1,c_2 do not determine a unique c_3")
    require(-common_factor[0] == c_sextic[3], "stored c_3 is not the exact common root")

    sources = {item["class_id"]: item for item in marked["sources"]}
    require(sorted(sources) == list(range(1, 8)), "marked-point sources are incomplete")
    for class_id, source in sources.items():
        require(source["patch_count"] == 48, f"class {class_id} lacks the c chart")
        require(source["c_patch"] == 46, f"class {class_id} has the wrong c patch")
        require(source["maximum_basis_residual"] < 6e-43, "source residual too large")
        require(len(source["sha256"]) == 64, "source digest is malformed")

    require(
        marked["c_coordinates"]["sextic_component"][3]["recognition_method"]
        == "exact_curve_gcd_from_c1_c2",
        "c_3 derivation method changed",
    )
    print("PASS exact b and c coordinates on both Hurwitz-algebra components")
    print("PASS exact quadric/cubic substitution and linear recovery of c_3")
    print("PASS 48-chart finite-equation source metadata")
    print(
        "SCOPE reconstruction checks alone do not identify the Hurwitz algebra; "
        "the independent branch-cycle certificate does"
    )


if __name__ == "__main__":
    main()
