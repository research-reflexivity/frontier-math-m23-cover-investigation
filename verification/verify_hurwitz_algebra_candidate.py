#!/usr/bin/env sage-python
"""Verify exact internal consistency of the reconstructed Hurwitz algebra.

This script proves statements about the stored rational reconstruction; it
does not by itself identify that reconstruction with the true Hurwitz
algebra.  The independent exact-eliminant and Arb branch-cycle certificates
in this repository supply that identification.
"""

from __future__ import annotations

import json
from pathlib import Path

from sage.all import PolynomialRing, QQ, RealField


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "hurwitz_algebra_candidate.json"


def rational(record):
    return QQ(record["numerator"]) / QQ(record["denominator"])


def quadratic(record, field, generator):
    return field(rational(record["rational_part"])) + field(
        rational(record["sqrt_minus_23_part"])
    ) * generator


def main() -> None:
    payload = json.loads(DATA.read_text())
    assert payload["status"] == "NUMERICAL_RECONSTRUCTION_CANDIDATE_ONLY"
    assert payload["invariant"] == "q13_over_A2"
    assert payload["class_ids"] == list(range(1, 8))

    base_ring = PolynomialRing(QQ, "s")
    s = base_ring.gen()
    field = QQ.extension(s**2 + 23, "sqrt_minus_23")
    sqrt_minus_23 = field.gen()
    polynomial_ring = PolynomialRing(field, "X")
    coefficients = [
        quadratic(record, field, sqrt_minus_23)
        for record in payload["coefficients"]
    ]
    polynomial = polynomial_ring(coefficients)
    assert polynomial.is_monic() and polynomial.degree() == 7

    class_6 = payload["class_6_root"]
    class_6_root = (
        field(rational(class_6["rational_part"]))
        + field(rational(class_6["sqrt_minus_23_part"])) * sqrt_minus_23
    )
    assert class_6_root == (
        field(QQ(148227) / 142129)
        - field(QQ(34830) / 142129) * sqrt_minus_23
    )
    assert class_6["exact_linear_factor_check"] is True
    assert polynomial(class_6_root) == 0

    factorization = list(polynomial.factor())
    assert [factor.degree() for factor, _ in factorization] == [1, 6]
    assert all(multiplicity == 1 for _, multiplicity in factorization)
    assert payload["exact_factor_degrees"] == [1, 6]
    assert payload["exact_sextic_irreducible_over_Q_sqrt_minus_23"] is True
    stored_factors = [
        polynomial_ring(
            [
                quadratic(record, field, sqrt_minus_23)
                for record in factor["coefficients_ascending"]
            ]
        )
        for factor in payload["exact_factors"]
    ]
    assert stored_factors == [factor for factor, _ in factorization]
    assert stored_factors[0] * stored_factors[1] == polynomial

    real = RealField(256)
    assert real(payload["maximum_root_residual"]) < real("5e-70")
    assert real(class_6["rational_part"]["absolute_error"]) < real("1e-90")
    assert real(class_6["sqrt_minus_23_part"]["absolute_error"]) < real("1e-90")
    summaries = payload["model_summaries"]
    assert [entry["class_id"] for entry in summaries] == list(range(1, 8))
    for entry in summaries:
        assert real(entry["maximum_basis_residual"]) < real("4e-77")
        if entry["class_id"] == 6:
            assert entry["terms"] == 480
        else:
            assert entry["terms"] == 360
            assert real(entry["quadric_balanced_jet_residual"]) < real("1e-74")
            assert real(entry["cubic_balanced_jet_residual"]) < real("1e-74")

    print("PASS reconstructed degree-seven polynomial over Q(sqrt(-23))")
    print("PASS exact class-6 linear factor and irreducible sextic factor")
    print("PASS numerical root residual below 5e-70")
    print(
        "SCOPE reconstruction checks alone do not identify the Hurwitz algebra; "
        "the independent branch-cycle certificate does"
    )


if __name__ == "__main__":
    main()
