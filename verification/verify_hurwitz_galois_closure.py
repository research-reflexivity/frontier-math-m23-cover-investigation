#!/usr/bin/env sage-python
"""Certify the Galois closure of the sextic Hurwitz component."""

from __future__ import annotations

import json
from pathlib import Path

from sage.all import (
    GF,
    PolynomialRing,
    QQ,
    TransitiveGroups,
    ZZ,
)


ROOT = Path(__file__).resolve().parents[1]
ALGEBRA_PATH = ROOT / "data" / "hurwitz_algebra_candidate.json"
BRANCH_PATH = ROOT / "verification" / "hurwitz_branch_cycle_summary.json"

ABSOLUTE_COEFFICIENTS_ASCENDING = [
    1,
    -6,
    20,
    -32,
    44,
    -22,
    6,
    -22,
    44,
    -32,
    20,
    -6,
    1,
]
TRACE_COEFFICIENTS_ASCENDING = [-44, 44, -27, -2, 14, -6, 1]

# All three primes split in Q(sqrt(-23)).  Their factor degrees give a
# 6-cycle, a 5-cycle, and a transposition in the relative Galois group.
FROBENIUS_WITNESSES = {
    3: (6,),
    139: (1, 5),
    2671: (1, 1, 1, 1, 2),
}


def rational(record):
    return QQ(record["numerator"]) / QQ(record["denominator"])


def quadratic(record, field, generator):
    return field(rational(record["rational_part"])) + field(
        rational(record["sqrt_minus_23_part"])
    ) * generator


def factor_degrees_mod(polynomial, prime):
    reduction = polynomial.change_ring(GF(prime))
    assert reduction.is_squarefree()
    return tuple(
        sorted(
            factor.degree()
            for factor, multiplicity in reduction.factor()
            for _ in range(multiplicity)
        )
    )


def cycle_types(group):
    return {tuple(sorted(element.cycle_type())) for element in group}


def main():
    algebra = json.loads(ALGEBRA_PATH.read_text())
    branch = json.loads(BRANCH_PATH.read_text())

    rational_ring = PolynomialRing(QQ, "x")
    x = rational_ring.gen()
    absolute_polynomial = rational_ring(ABSOLUTE_COEFFICIENTS_ASCENDING)
    trace_polynomial = rational_ring(TRACE_COEFFICIENTS_ASCENDING)

    # Reconstruct the exact relative sextic directly from the Hurwitz data.
    base_ring = PolynomialRing(QQ, "s")
    s0 = base_ring.gen()
    base = QQ.extension(s0**2 + 23, "sqrt_minus_23")
    sqrt_minus_23 = base.gen()
    invariant_ring = PolynomialRing(base, "J")
    factors = [
        invariant_ring(
            [
                quadratic(coefficient, base, sqrt_minus_23)
                for coefficient in factor["coefficients_ascending"]
            ]
        )
        for factor in algebra["exact_factors"]
    ]
    relative_sextic = next(factor for factor in factors if factor.degree() == 6)
    assert relative_sextic.is_irreducible()
    relative_field = base.extension(relative_sextic, "j")

    # Tie its absolute field to the small reciprocal polynomial used below.
    absolute_field = relative_field.absolute_field("a")
    assert absolute_field.polynomial() == absolute_polynomial
    _, to_absolute = absolute_field.structure()
    a = absolute_field.gen()
    s_in_absolute = to_absolute(relative_field(sqrt_minus_23))
    y = a + 1 / a

    assert s_in_absolute**2 == -23
    assert s_in_absolute.minpoly() == x**2 + 23
    assert y.minpoly() == trace_polynomial
    assert (y + s_in_absolute).minpoly().degree() == 12
    assert x**6 * trace_polynomial(x + 1 / x) == absolute_polynomial

    # Maximal-order discriminants and the relative discriminant ideal.
    trace_field = QQ.extension(trace_polynomial, "y")
    assert trace_polynomial.discriminant() == 2**22 * 11 * 23**4
    assert trace_field.discriminant() == 2**4 * 11 * 23**4
    assert absolute_field.discriminant() == 2**8 * 11**2 * 23**8
    assert relative_field.relative_discriminant() == base.ideal(2**4 * 11 * 23)

    # Frobenius cycle types at primes unramified in the trace field and split
    # in the quadratic base field.
    for prime, expected_degrees in FROBENIUS_WITNESSES.items():
        assert ZZ(-23).kronecker(prime) == 1
        assert trace_polynomial.discriminant() % prime != 0
        assert factor_degrees_mod(trace_polynomial, prime) == expected_degrees

    # The mod-3 irreducibility proves transitivity.  An exhaustive census of
    # the degree-six transitive groups shows that the three certified cycle
    # types force 6T16 = S6.
    required_types = set(FROBENIUS_WITNESSES.values())
    compatible_groups = []
    for group in TransitiveGroups(6):
        if required_types <= cycle_types(group):
            compatible_groups.append(group)
    assert len(compatible_groups) == 1
    galois_group = compatible_groups[0]
    assert galois_group.transitive_number() == 16
    assert galois_group.order() == 720
    assert galois_group.structure_description() == "S6"

    # The unique quadratic subfield of an S6 splitting field is its
    # discriminant field.  Its square class is 11, not -23, so adjoining the
    # Hurwitz orientation field is linearly disjoint from the S6 closure.
    discriminant_square_class = ZZ(trace_polynomial.discriminant()).squarefree_part()
    assert discriminant_square_class == 11
    assert discriminant_square_class != -23

    sextic_records = sorted(
        (
            record
            for record in branch["records"]
            if record["component"] == "sextic"
        ),
        key=lambda record: record["embedding_index"],
    )
    nielsen_ids = tuple(record["class_id"] for record in sextic_records)
    assert nielsen_ids == (7, 4, 1, 5, 3, 2)

    print("PASS L is the compositum of Q(sqrt(-23)) and the trace sextic E")
    print("PASS disc(E)=2^4*11*23^4 and relative_disc(L/K0)=(2^4*11*23)")
    print("PASS split-prime Frobenius types 6, 1+5, and 1+1+1+1+2 force S6")
    print("PASS Gal(L^gal/K0)=S6 and Gal(L^gal/Q)=S6 x C2")
    print("PASS the natural six-point action has Nielsen order (7,4,1,5,3,2)")


if __name__ == "__main__":
    main()
