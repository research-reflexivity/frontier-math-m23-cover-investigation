#!/usr/bin/env sage-python
"""Verify the serialized exact degree-23 monodromy eliminant."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = ROOT / "data" / "hurwitz_monodromy_eliminant_candidate.json"
CANONICAL_PATH = ROOT / "data" / "hurwitz_canonical_models_candidate.json"


def rational(record):
    return QQ(record["numerator"]) / QQ(record["denominator"])


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    payload = json.loads(DATA_PATH.read_text())
    canonical = json.loads(CANONICAL_PATH.read_text())
    assert payload["schema"] == "m23.cover-investigation.monodromy-eliminant-candidate.v1"
    assert payload["status"] == "exact_degree23_eliminant_for_reconstructed_sextic_component"
    assert payload["field_degree"] == 12
    field_coefficients = [
        rational(value)
        for value in payload["field_defining_polynomial_coefficients_ascending"]
    ]
    assert field_coefficients == [
        rational(value)
        for value in canonical["absolute_field_polynomial_coefficients_ascending"]
    ]
    field_ring = PolynomialRing(QQ, "a")
    field = QQ.extension(field_ring(field_coefficients), "a")
    assert field.polynomial().is_irreducible()

    def decode(records):
        assert len(records) == 12
        return field([rational(value) for value in records])

    encoded = payload["coefficients_target_then_x"]
    assert len(encoded) == 7 and all(len(row) == 24 for row in encoded)
    coefficients = [
        [decode(encoded[i][j]) for j in range(24)] for i in range(7)
    ]
    target_ring = PolynomialRing(field, "target")
    target = target_ring.gen()
    coefficient_polynomials = [
        target_ring([coefficients[i][j] for i in range(7)]) for j in range(24)
    ]
    assert max(value.degree() for value in coefficient_polynomials) == 6

    root_ring = PolynomialRing(field, "x")

    def fibre(value):
        return root_ring([coefficient(value) for coefficient in coefficient_polynomials])

    zero = fibre(0)
    one = fibre(1)
    infinity = root_ring([coefficient[6] for coefficient in coefficient_polynomials])
    assert zero.degree() == one.degree() == infinity.degree() == 23
    assert zero.gcd(zero.derivative()).degree() == 22
    one_gcd = one.gcd(one.derivative())
    assert one_gcd.degree() == 8
    assert one_gcd.gcd(one_gcd.derivative()).degree() == 0
    assert infinity.gcd(infinity.derivative()).degree() == 22
    assert payload["exact_branch_fibres"] == {
        "infinity_gcd_degree": 22,
        "one_gcd_degree": 8,
        "one_gcd_squarefree": True,
        "passport": ["23", "2^8 1^7", "23"],
        "zero_gcd_degree": 22,
    }
    assert payload["construction"] == {
        "holdout_exact": True,
        "holdout_target": 9,
        "interpolation_targets": [2, 3, 4, 5, 6, 7, 8],
        "plane_curve_degree": 6,
        "plane_section_degrees": [7, 7],
        "raw_resultant_degree": 42,
        "target_independent_base_factor_degree": 19,
    }
    for relative, expected in payload["input_sha256"].items():
        assert sha256(ROOT / relative) == expected

    print("PASS exact degree-23 eliminant over the degree-12 Hurwitz field")
    print("PASS exact eliminant passport (23), (2^8 1^7), (23)")
    print("PASS recorded seven-fibre interpolation and exact holdout metadata")


if __name__ == "__main__":
    main()
