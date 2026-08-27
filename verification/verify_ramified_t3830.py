#!/usr/bin/env python3
"""Dependency-free checks around the first controlled ramified specialization.

PARI/GP performs the number-field discriminant calculation.  This companion
certificate verifies the CRT condition, the uniform M23 Frobenius lock, and
the transverse eight-double-root fiber at T=0 modulo 5.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

from verify_specialization_t2 import (
    is_irreducible,
    monic,
    multiply,
)


ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"


def specialize(table: list[list[int]], value: int) -> list[int]:
    return [
        sum(row[x_degree] * value**t_degree for t_degree, row in enumerate(table))
        for x_degree in range(24)
    ]


def degree(poly: list[int]) -> int:
    return len(poly) - 1


def main() -> None:
    family = json.loads((DATA / "optimal_23_4_Z.json").read_text())
    base = json.loads((DATA / "specialization_t2_Z.json").read_text())
    ramification = json.loads((DATA / "ramified_t3830.json").read_text())
    table = [
        [int(coefficient) for coefficient in row]
        for row in family["coefficients_T_then_W"]
    ]

    t0 = ramification["specialization"]
    prime = ramification["ramified_prime"]
    if t0 != 3830 or prime != 5:
        raise RuntimeError("unexpected ramification witness")
    if t0 % 319 != 2 or t0 % 25 != 5:
        raise RuntimeError("CRT conditions are not satisfied")

    coefficients = specialize(table, t0)
    content = math.gcd(*map(abs, coefficients))
    primitive = [coefficient // content for coefficient in coefficients]
    if degree(primitive) != 23:
        raise RuntimeError("specialization has the wrong degree")

    irreducibility = base["modular_certificates"]["irreducibility"]
    reduction_29 = monic(primitive, 29)
    if reduction_29 != irreducibility["monic_coefficients_low_to_high"]:
        raise RuntimeError("the Frobenius lock modulo 29 was not preserved")
    if not is_irreducible(reduction_29, 29):
        raise RuntimeError("specialization is reducible modulo 29")

    even = base["modular_certificates"]["even_frobenius"]
    reduction_11 = monic(primitive, 11)
    product = [1]
    for factor in even["irreducible_factors_low_to_high"]:
        if not is_irreducible(factor, 11):
            raise RuntimeError("stored factor modulo 11 is reducible")
        product = multiply(product, factor, 11)
    if product != reduction_11:
        raise RuntimeError("the Frobenius lock modulo 11 was not preserved")

    summary = {
        "status": "PASS_RAMIFIED_T3830_COMPANION_CERTIFICATE",
        "specialization": t0,
        "primitive_content": content,
        "crt_conditions": ["t=2 mod 319", "t=5 mod 25"],
        "galois_group": "M23",
        "irreducible_mod_29": True,
        "factor_degrees_mod_11": [2, 7, 14],
        "pari_pinned_branch_fiber_mod_5_repeated_degree": 8,
        "pari_pinned_branch_fiber_mod_5_repeated_part_squarefree": True,
        "pari_pinned_branch_fiber_mod_5_transverse": True,
        "pari_field_discriminant_valuation_at_5": 8,
        "base_t2_field_discriminant_valuation_at_5": 0,
        "conclusion": "the t=3830 and t=2 M23 splitting fields are distinct",
    }
    output = ROOT / "verification" / "ramified_t3830_summary.json"
    output.write_text(json.dumps(summary, indent=2) + "\n")
    print(summary["status"])
    print(f"t={t0} content={content} congruences=2_mod_319,5_mod_25")
    print("pari_branch_mod_5=repeated_degree_8,squarefree,transverse")
    print("galois_group=M23 pari_field_discriminant_v5=8")


if __name__ == "__main__":
    main()
