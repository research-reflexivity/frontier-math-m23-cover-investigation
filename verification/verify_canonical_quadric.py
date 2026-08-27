#!/usr/bin/env python3
"""Verify the reconstructed canonical quadric and its ruling field exactly."""

from __future__ import annotations

import json
from pathlib import Path

from sage.all import GF, QQ, PolynomialRing, matrix

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"


def fraction(entry: dict[str, int]):
    return QQ(entry["numerator"]) / QQ(entry["denominator"])


def main() -> None:
    payload = json.loads((DATA / "canonical_quadric_Q.json").read_text())
    assert payload["status"] == "PASS_CANONICAL_QUADRIC_Q_HOLDOUT"
    assert payload["basis_modulus_bits"] == 604
    assert payload["quadric_modulus_bits"] == 141
    assert payload["holdout_prime"] == 137
    assert payload["basis_holdout_mismatch_count"] == 0
    assert payload["quadric_holdout_mismatch_count"] == 0

    vectors = [[fraction(x) for x in row] for row in payload["canonical_vectors"]]
    assert len(vectors) == 4 and all(len(row) == 88 for row in vectors)
    # The last four adjoint slots are the free columns in the canonical RREF.
    assert matrix(QQ, [row[84:88] for row in vectors]) == matrix.identity(QQ, 4)

    rt = PolynomialRing(QQ, "T")
    t = rt.gen()
    kt = rt.fraction_field()
    rv = PolynomialRing(kt, "V")
    v = rv.gen()
    zlines = [
        line for line in (DATA / "Fint_coefficients_Z.json").read_text().splitlines()
        if line
    ]
    source = json.loads(zlines[-1])
    f = rv.zero()
    for j, row in enumerate(source):
        f += rv(sum(QQ(c) * t**i for i, c in enumerate(row))) * v**j
    fhat = f.monic()
    assert fhat.degree() == 23

    slots = [tuple(slot) for slot in payload["adjoint_slots_j_m_i"]]
    assert len(slots) == 88
    d = t**2 + 23

    def numerator(vector):
        answer = rv.zero()
        for coefficient, (j, m, i) in zip(vector, slots):
            answer += coefficient * t**i * d**m * v**j
        return answer

    canonical = [numerator(vector) for vector in vectors]
    pairs = [tuple(pair) for pair in payload["quadric_pairs"]]
    assert pairs == [(i, j) for i in range(4) for j in range(i, 4)]
    products = [(canonical[i] * canonical[j]) % fhat for i, j in pairs]
    coefficients = [fraction(x) for x in payload["quadric_coefficients"]]
    assert coefficients[-1] == 1
    relation = sum((c * product for c, product in zip(coefficients, products)), rv.zero())
    assert relation == 0

    product_matrix = matrix(
        kt,
        23,
        10,
        lambda row, column: products[column][row],
    )
    assert product_matrix.rank() == 9
    assert product_matrix * matrix(kt, 10, 1, coefficients) == 0

    twice_symmetric = matrix(
        QQ,
        [[fraction(x) for x in row] for row in payload["twice_symmetric_matrix"]],
    )
    determinant = fraction(payload["twice_symmetric_determinant"])
    assert twice_symmetric.det() == determinant
    assert determinant == QQ(644454138716416151027888) / QQ(970299)
    expected_factorization = (
        QQ(2) ** 4
        * QQ(3) ** (-6)
        * QQ(11) ** (-3)
        * QQ(23) ** 2
        * QQ(443)
        * QQ(414578063) ** 2
    )
    assert determinant == expected_factorization
    assert payload["determinant_squarefree_part"] == 4873 == 11 * 443
    square_root = fraction(payload["determinant_over_squarefree_part_sqrt"])
    assert square_root == QQ(38141181796) / QQ(3267)
    assert determinant == 4873 * square_root**2
    assert not (QQ(4873) / QQ(-23)).is_square()
    assert payload["ruling_field"] == "Q(sqrt(4873)) = Q(sqrt(11*443))"
    assert payload["branch_orientation_field"] == "Q(sqrt(-23))"
    assert payload["fields_equal"] is False

    # Independent compatibility with the directly recomputed canonical RREF mod 31.
    mod31 = json.loads((DATA / "canonical_pencil_mod31.json").read_text())
    k31 = GF(31)
    reduced_vectors = [
        [k31(x.numerator()) / k31(x.denominator()) for x in row]
        for row in vectors
    ]
    assert reduced_vectors == [[k31(x) for x in row] for row in mod31["canonical_vectors"]]
    assert k31(determinant) == 12
    assert not k31(determinant).is_square()
    assert k31(-23).is_square()

    print("PASS exact canonical-quadric identity and uniqueness over Q(T)")
    print("PASS canonical ruling field Q(sqrt(4873)) differs from Q(sqrt(-23))")


if __name__ == "__main__":
    main()
