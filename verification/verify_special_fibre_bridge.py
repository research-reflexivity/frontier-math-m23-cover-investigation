#!/usr/bin/env sage-python
"""Exact bridge between the original and optimal T=0 fibres.

All polynomials are rebuilt from the canonical JSON tables.  The certificate
checks the factorization patterns and proves that the specialized rational
function W=J0/J1 generates the septic and octic factors of the optimal fibre
from the corresponding factors of the original fibre.
"""

from __future__ import annotations

import json
from pathlib import Path

from sage.all import QQ, ZZ, PolynomialRing, gcd, lcm


ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"


def load_json(name: str):
    return json.loads((DATA / name).read_text())


def primitive_part(poly):
    content = gcd([abs(c) for c in poly.list()])
    return poly / content, content


def primitive_integral(poly):
    denominator = lcm([coefficient.denominator() for coefficient in poly.list()])
    integral = (denominator * poly).change_ring(ZZ)
    content = gcd([abs(coefficient) for coefficient in integral.list()])
    integral /= content
    if integral.leading_coefficient() < 0:
        integral = -integral
    return integral


def main() -> None:
    source_table = load_json("Fint_coefficients_Z.json")
    optimal = load_json("optimal_23_4_Z.json")
    pencil = load_json("optimal_degree4_pencil.json")

    RV = PolynomialRing(QQ, "V")
    V = RV.gen()
    RW = PolynomialRing(QQ, "W")
    W = RW.gen()

    F0 = sum(QQ(row[0]) * V**v_degree
             for v_degree, row in enumerate(source_table))
    F0_primitive, F0_content = primitive_part(F0)
    source_factors = list(F0_primitive.factor())
    assert [(factor.degree(), multiplicity)
            for factor, multiplicity in source_factors] == [(7, 1), (8, 2)]
    H7, R8 = source_factors[0][0], source_factors[1][0]
    assert H7.is_irreducible() and R8.is_irreducible()
    assert F0_content == 23**4

    P0 = sum(QQ(coefficient) * W**w_degree
             for w_degree, coefficient
             in enumerate(optimal["coefficients_T_then_W"][0]))
    P0_primitive, P0_content = primitive_part(P0)
    optimal_factors = list(P0_primitive.factor())
    assert [(factor.degree(), multiplicity)
            for factor, multiplicity in optimal_factors] == [(7, 1), (8, 2)]
    Q7, S8 = optimal_factors[0][0], optimal_factors[1][0]
    assert Q7.is_irreducible() and S8.is_irreducible()
    Q7_integral = primitive_integral(Q7)
    S8_integral = primitive_integral(S8)
    assert P0 == -(23**3) * RW(Q7_integral) * RW(S8_integral)**2

    specialized_pencil = []
    for encoded in pencil["expanded_coefficients_V_then_T"]:
        value = RV.zero()
        for key, coefficient in encoded.items():
            v_degree, t_degree = map(int, key.split(","))
            if t_degree == 0:
                value += QQ(coefficient) * V**v_degree
        specialized_pencil.append(value)
    J0, J1 = specialized_pencil
    assert gcd(J1, H7 * R8) == 1

    transported = []
    for source_factor, target_factor in ((H7, Q7), (R8, S8)):
        K = QQ.extension(source_factor, "v")
        v = K.gen()
        w = J0(v) / J1(v)
        minimal = w.minpoly("W")
        minimal = RW(minimal)
        minimal = minimal / minimal.leading_coefficient()
        target_monic = target_factor / target_factor.leading_coefficient()
        assert minimal == target_monic
        assert target_factor(w) == 0
        assert minimal.degree() == source_factor.degree()
        transported.append({
            "source_degree": int(source_factor.degree()),
            "target_degree": int(target_factor.degree()),
            "generator_minimal_polynomial": str(minimal),
        })

    summary = {
        "status": "PASS_SPECIAL_FIBRE_BRIDGE",
        "source_special_fibre_content": int(F0_content),
        "source_factor_degrees_and_multiplicities": [[7, 1], [8, 2]],
        "optimal_special_fibre_content": int(P0_content),
        "optimal_factor_degrees_and_multiplicities": [[7, 1], [8, 2]],
        "optimal_integral_factorization_scalar": -(23**3),
        "optimal_septic": str(Q7_integral),
        "optimal_octic": str(S8_integral),
        "J1_nonzero_on_all_geometric_special_fibre_points": True,
        "specialized_W_generates_both_factor_fields": True,
        "transported_factors": transported,
        "conclusion": (
            "The degree-7 Fano point field and degree-8 affine-point field "
            "are carried isomorphically from the original V-model to the "
            "optimal W-model at T=0."
        ),
    }
    output = ROOT / "verification" / "bridge_summary.json"
    output.write_text(json.dumps(summary, indent=2) + "\n")

    print(summary["status"])
    print("source factors: 7^1, 8^2")
    print("optimal factors: 7^1, 8^2")
    print("J1 is invertible on both factors")
    print("W generates the septic and octic target fields")


if __name__ == "__main__":
    main()
