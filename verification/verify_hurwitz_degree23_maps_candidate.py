#!/usr/bin/env sage-python
"""Check the serialized exact degree-23 map sections and normalization."""

from __future__ import annotations

import json
from pathlib import Path

from sage.all import PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[1]
MAPS_PATH = ROOT / "data" / "hurwitz_degree23_maps_candidate.json"
BRANCH_PATH = ROOT / "data" / "hurwitz_degree23_branch_candidate.json"


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def rational(record):
    return QQ(record["numerator"]) / QQ(record["denominator"])


def decode_element(record, field):
    return field([rational(value) for value in record])


def decode_section(record, field):
    result = {}
    for term in record:
        monomial = tuple(term["monomial"])
        require(len(monomial) == 4 and sum(monomial) == 5, "non-quintic term")
        require(monomial not in result, "duplicate section monomial")
        coefficient = decode_element(term["coefficient_power_basis"], field)
        require(coefficient != 0, "serialized zero coefficient")
        result[monomial] = coefficient
    return result


def main():
    maps = json.loads(MAPS_PATH.read_text())
    branch = json.loads(BRANCH_PATH.read_text())
    require(
        maps["schema"] == "m23.cover-investigation.degree23-maps-candidate.v1",
        "unexpected map schema",
    )
    require(
        maps["status"] == "exact_map_sections_and_exact_third_branch_normalization",
        "map status is not exact",
    )
    require(
        maps["third_fiber_certificate"]
        == {
            "certificate": "verification/verify_hurwitz_degree23_third_fiber.sage",
            "extra_gcd_degree": 8,
            "extra_gcd_squarefree": True,
            "magma_collision_free_residue_embeddings": [2, 3, 4, 5, 9, 12],
            "magma_geometry_certificate": "verification/verify_hurwitz_degree23_geometry.m",
            "projected_point_collisions_excluded": True,
            "status": "verified_exact",
        },
        "third-fibre certificate metadata changed",
    )
    require(
        [item["component"] for item in maps["components"]]
        == ["degree_one", "sextic"],
        "map components changed",
    )

    base_ring = PolynomialRing(QQ, "s")
    s = base_ring.gen()
    base_field = QQ.extension(s**2 + 23, "sqrt_minus_23")
    absolute_ring = PolynomialRing(QQ, "a")
    absolute_polynomial = absolute_ring([
        rational(value)
        for value in branch["field"]["defining_polynomial_coefficients_ascending"]
    ])
    absolute_field = QQ.extension(absolute_polynomial, "a")
    fields = {"degree_one": base_field, "sextic": absolute_field}

    for component in maps["components"]:
        label = component["component"]
        field = fields[label]
        require(component["field_degree"] == field.degree(), f"{label} field degree")
        require(component["degree_5_quotient_dimension"] == 27, f"{label} degree-5 Hilbert value")
        require(component["degree_10_quotient_dimension"] == 57, f"{label} degree-10 Hilbert value")
        require(component["jet_ranks"] == {"b": 23, "c": 23}, f"{label} jet ranks")
        require(
            component["section_space_dimensions"] == {"b": 4, "c": 4},
            f"{label} section dimensions",
        )
        require(component["multiplier_system_rank"] == 15, f"{label} multiplier rank")
        require(component["multiplier_kernel_dimension"] == 1, f"{label} multiplier kernel")
        require(component["multiplier_determinant_nonzero"], f"{label} singular multiplier")
        require(
            component["orders"]
            == {"numerator_at_b": 23, "denominator_at_c": 23},
            f"{label} local orders",
        )
        numerator = decode_section(component["numerator"], field)
        denominator = decode_section(component["denominator"], field)
        require(numerator and denominator, f"{label} empty map section")
        first_monomial = next(iter(numerator))
        require(numerator[first_monomial] == 1, f"{label} numerator normalization")

        if label == "sextic":
            lambda_value = absolute_field([
                QQ(value) / QQ(branch["lambda"]["denominator"])
                for value in branch["lambda"]["power_basis_numerators"]
            ])
            serialized_lambda = decode_element(
                component["third_branch_value_power_basis"], absolute_field
            )
            require(serialized_lambda == lambda_value, "third branch values disagree")
            normalized = decode_section(component["normalized_denominator"], field)
            require(set(normalized) == set(denominator), "normalized denominator support")
            require(
                all(normalized[monomial] == lambda_value * value for monomial, value in denominator.items()),
                "normalized denominator is not lambda times denominator",
            )
            require(
                component["normalized_map"]
                == "beta=numerator/normalized_denominator",
                "normalized map convention changed",
            )

    print("PASS exact serialized map sections over both Hurwitz components")
    print("PASS degree-12 normalization beta=N/(lambda D)")
    print("PASS recorded ranks and local orders are internally consistent")
    print("PASS exact characteristic-zero third-fibre certificate metadata")


if __name__ == "__main__":
    main()
