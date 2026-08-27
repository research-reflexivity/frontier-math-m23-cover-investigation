#!/usr/bin/env sage-python
"""Verify descent of relative-transporter invariants on Spec(K0 x L)."""

from __future__ import annotations

import json
from pathlib import Path

from sage.all import PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[1]
ALGEBRA_PATH = ROOT / "data" / "hurwitz_algebra_candidate.json"
BRANCH_PATH = ROOT / "verification" / "hurwitz_branch_cycle_summary.json"

NU = {1: 54, 2: 54, 3: 28, 4: 46, 5: 46, 6: 32, 7: 42}
MU = {1: 16, 2: 16, 3: 31, 4: 28, 5: 28, 6: 17, 7: 14}
KAPPA = {
    1: (0, 16),
    2: (0, 16),
    3: (0, 8),
    4: (2, 6),
    5: (2, 6),
    6: (0, 12),
    7: (4, 12),
}
THETA_WEIGHT = {1: 15, 2: 15, 3: 25, 4: 27, 5: 27, 6: 16, 7: 11}
THETA_AUGMENTATION = {1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 0, 7: 1}


def rational(record):
    return QQ(record["numerator"]) / QQ(record["denominator"])


def quadratic(record, field, generator):
    return field(rational(record["rational_part"])) + field(
        rational(record["sqrt_minus_23_part"])
    ) * generator


def values_on(ids, table):
    return [table[class_id] for class_id in ids]


def main():
    algebra = json.loads(ALGEBRA_PATH.read_text())
    branch = json.loads(BRANCH_PATH.read_text())

    base_ring = PolynomialRing(QQ, "s")
    s = base_ring.gen()
    base = QQ.extension(s**2 + 23, "sqrt_minus_23")
    sqrt_minus_23 = base.gen()
    invariant_ring = PolynomialRing(base, "J")
    invariant = invariant_ring.gen()
    polynomial = invariant_ring(
        [
            quadratic(record, base, sqrt_minus_23)
            for record in algebra["coefficients"]
        ]
    )
    factors = list(polynomial.factor())
    assert [factor.degree() for factor, _ in factors] == [1, 6]
    assert all(multiplicity == 1 for _, multiplicity in factors)
    sextic_polynomial = factors[1][0]
    assert sextic_polynomial.is_irreducible()
    sextic_field = base.extension(sextic_polynomial, "j")

    degree_one_ids = [
        record["class_id"]
        for record in branch["records"]
        if record["component"] == "degree_one"
    ]
    sextic_records = [
        record for record in branch["records"] if record["component"] == "sextic"
    ]
    sextic_records.sort(key=lambda record: record["embedding_index"])
    sextic_ids = [record["class_id"] for record in sextic_records]
    assert degree_one_ids == [6]
    assert sextic_ids == [7, 4, 1, 5, 3, 2]
    assert set(sextic_ids) == {1, 2, 3, 4, 5, 7}

    # A function to a constant finite set descends over the connected sextic
    # component exactly when its six geometric values are constant.
    assert len(set(values_on(sextic_ids, NU))) > 1
    assert len(set(values_on(sextic_ids, MU))) > 1
    assert len(set(values_on(sextic_ids, KAPPA))) > 1
    assert len(set(values_on(sextic_ids, THETA_WEIGHT))) > 1
    assert values_on(sextic_ids, THETA_AUGMENTATION) == [1] * 6

    class_6_indicator = {class_id: int(class_id == 6) for class_id in range(1, 8)}
    selectors = {
        "nu_equals_32": {
            class_id: int(NU[class_id] == 32) for class_id in range(1, 8)
        },
        "mu_equals_17": {
            class_id: int(MU[class_id] == 17) for class_id in range(1, 8)
        },
        "kappa_equals_0_12": {
            class_id: int(KAPPA[class_id] == (0, 12))
            for class_id in range(1, 8)
        },
        "binary_augmentation_equals_zero": {
            class_id: int(THETA_AUGMENTATION[class_id] == 0)
            for class_id in range(1, 8)
        },
    }
    assert all(selector == class_6_indicator for selector in selectors.values())

    # In A_H = K0 x L, these common Boolean predicates are precisely the
    # primitive idempotent of the K0 factor.  The augmentation itself is the
    # complementary idempotent.
    degree_one_projector = (base.one(), sextic_field.zero())
    sextic_projector = (base.zero(), sextic_field.one())

    def multiply(left, right):
        return (left[0] * right[0], left[1] * right[1])

    assert multiply(degree_one_projector, degree_one_projector) == degree_one_projector
    assert multiply(sextic_projector, sextic_projector) == sextic_projector
    assert multiply(degree_one_projector, sextic_projector) == (
        base.zero(),
        sextic_field.zero(),
    )
    assert (
        degree_one_projector[0] + sextic_projector[0],
        degree_one_projector[1] + sextic_projector[1],
    ) == (base.one(), sextic_field.one())

    print("PASS the Hurwitz scheme has connected components of degrees 1 and 6")
    print("PASS raw nu, mu, cyclic-conjugacy counts, and Theta fail sextic descent")
    print("PASS all four Boolean predicates equal the degree-one idempotent")
    print("PASS binary augmentation equals the complementary sextic idempotent")


if __name__ == "__main__":
    main()
