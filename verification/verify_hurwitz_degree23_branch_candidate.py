#!/usr/bin/env sage-python
"""Verify the exact CRT/LLL reconstruction of the third branch value.

This checks the arithmetic certificate stored in
``data/hurwitz_degree23_branch_candidate.json``.  It does not recompute the
canonical degree-five pencil; that slower geometric check is kept separate.
"""

from __future__ import annotations

import json
from pathlib import Path

from sage.all import CRT, IntegerModRing, PolynomialRing, QQ, ZZ, matrix, prod


ROOT = Path(__file__).resolve().parents[1]
BRANCH_PATH = ROOT / "data" / "hurwitz_degree23_branch_candidate.json"
CANONICAL_PATH = ROOT / "data" / "hurwitz_canonical_models_candidate.json"


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def rational(record):
    return QQ(record["numerator"]) / QQ(record["denominator"])


def normalized_vector(row):
    if row[0] == 0:
        return None
    denominator = abs(ZZ(row[0]))
    sign = 1 if row[0] > 0 else -1
    numerators = [sign * ZZ(value) for value in row[1:]]
    common = denominator
    for numerator in numerators:
        common = common.gcd(abs(numerator))
    return denominator // common, [value // common for value in numerators]


def main():
    branch = json.loads(BRANCH_PATH.read_text())
    canonical = json.loads(CANONICAL_PATH.read_text())
    require(
        branch["schema"] == "m23.cover-investigation.degree23-branch-candidate.v1",
        "unexpected branch-certificate schema",
    )
    require(
        branch["status"]
        == "exact_reconstruction_and_characteristic_zero_fiber_verified",
        "third-fibre status is not exact",
    )
    require(
        branch["characteristic_zero_third_fiber"]
        == {
            "branch_gcd_degree": 14,
            "branch_resultant_degree": 42,
            "certificate": "verification/verify_hurwitz_degree23_third_fiber.sage",
            "common_gcd_degree": 6,
            "extra_gcd_degree": 8,
            "extra_gcd_squarefree": True,
            "generic_gcd_degree": 6,
            "generic_resultant_degree": 42,
            "magma_collision_free_residue_embeddings": [2, 3, 4, 5, 9, 12],
            "magma_geometry_certificate": "verification/verify_hurwitz_degree23_geometry.m",
            "magma_residue_embeddings": 12,
            "projected_point_collisions_excluded": True,
            "status": "verified_exact",
        },
        "third-fibre certificate metadata changed",
    )
    field_coefficients = [
        rational(value)
        for value in branch["field"]["defining_polynomial_coefficients_ascending"]
    ]
    canonical_coefficients = [
        rational(value)
        for value in canonical["absolute_field_polynomial_coefficients_ascending"]
    ]
    require(field_coefficients == canonical_coefficients, "number fields disagree")

    reconstruction = {
        ZZ(prime): [ZZ(value) for value in values]
        for prime, values in branch["crt_reconstruction"]["records"].items()
    }
    holdouts = {
        ZZ(prime): [ZZ(value) for value in values]
        for prime, values in branch["independent_holdouts"]["records"].items()
    }
    require(len(reconstruction) == 71, "reconstruction-prime count changed")
    require(len(holdouts) >= 20, "too few independent holdout primes")
    require(set(reconstruction).isdisjoint(holdouts), "holdouts leaked into CRT")

    polynomial_ring = PolynomialRing(QQ, "a")
    defining_polynomial = polynomial_ring(field_coefficients)
    require(defining_polynomial.is_monic(), "absolute polynomial is not monic")
    require(defining_polynomial.is_irreducible(), "absolute polynomial is reducible")
    for prime in [*reconstruction, *holdouts]:
        require(prime.is_prime(), f"{prime} is not prime")
        finite_ring = PolynomialRing(IntegerModRing(prime), "a")
        reduced = finite_ring([
            IntegerModRing(prime)(value.numerator())
            / IntegerModRing(prime)(value.denominator())
            for value in field_coefficients
        ])
        require(
            len(reduced.roots()) == 12,
            f"absolute field does not split completely at {prime}",
        )

    primes = list(reconstruction)
    modulus = prod(primes)
    require(
        len(str(modulus)) == branch["crt_reconstruction"]["modulus_decimal_digits"],
        "recorded CRT modulus size changed",
    )
    residues = [
        CRT([reconstruction[prime][coordinate] for prime in primes], primes)
        for coordinate in range(12)
    ]
    rows = [[ZZ(1), *map(ZZ, residues)]]
    for coordinate in range(12):
        row = [ZZ(0)] * 13
        row[coordinate + 1] = modulus
        rows.append(row)
    reduced_basis = matrix(ZZ, rows).LLL(delta=0.99)
    candidates = [
        candidate
        for candidate in (normalized_vector(row) for row in reduced_basis.rows())
        if candidate is not None
    ]
    candidates.sort(
        key=lambda item: max(
            [item[0].nbits(), *[abs(value).nbits() for value in item[1]]]
        )
    )
    denominator, numerators = candidates[0]
    stored_denominator = ZZ(branch["lambda"]["denominator"])
    stored_numerators = [ZZ(value) for value in branch["lambda"]["power_basis_numerators"]]
    require(denominator == stored_denominator, "LLL denominator changed")
    require(numerators == stored_numerators, "LLL numerator vector changed")
    require(
        max([denominator.nbits(), *[abs(value).nbits() for value in numerators]])
        == branch["lambda"]["lattice_vector_bits"],
        "recorded lattice-vector size changed",
    )
    if len(candidates) > 1:
        first_bits = max([denominator.nbits(), *[abs(value).nbits() for value in numerators]])
        second_denominator, second_numerators = candidates[1]
        second_bits = max([
            second_denominator.nbits(),
            *[abs(value).nbits() for value in second_numerators],
        ])
        require(second_bits - first_bits >= 100, "short-vector gap is too small")

    def check_records(records, label):
        for prime, expected in records.items():
            finite_field = IntegerModRing(prime)
            actual = [
                ZZ(finite_field(numerator) / finite_field(denominator))
                for numerator in numerators
            ]
            require(actual == expected, f"{label} mismatch at {prime}")

    check_records(reconstruction, "CRT")
    check_records(holdouts, "holdout")

    absolute_field = QQ.extension(defining_polynomial, "a")
    branch_value = absolute_field([
        QQ(value) / QQ(denominator) for value in numerators
    ])
    require(branch_value != 0, "third branch value is zero")
    require(branch_value.minpoly().degree() == 12, "third branch value is not primitive")

    print("PASS exact 71-prime CRT/LLL reconstruction of lambda")
    print(f"PASS {len(holdouts)} independent split-prime holdouts")
    print("PASS lambda has degree 12 in the reconstructed Hurwitz field")
    print("PASS exact characteristic-zero third-fibre metadata")


if __name__ == "__main__":
    main()
