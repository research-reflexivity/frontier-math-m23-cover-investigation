#!/usr/bin/env python3
"""Exact certificate for the arithmetic progression T = 2 (mod 319).

The specialization at T=2 has two useful reductions: an irreducible
degree-23 reduction modulo 29 and squarefree factor degrees (2,7,14) modulo
11.  Since 319=11*29, both reductions are unchanged throughout the
progression.  This script verifies that universal congruence coefficient by
coefficient, reruns the finite-field tests, and checks the maximal-subgroup
exclusion.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

from verify_specialization_t2 import (
    derivative,
    gcd_poly,
    is_irreducible,
    monic,
    multiply,
)


ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RESIDUE = 2
MODULUS = 319
PRIMES = (11, 29)


def coefficient_after_affine_substitution(
    coefficients_in_t: list[int], power_of_n: int
) -> int:
    """Coefficient of n^power in a(RESIDUE + MODULUS*n)."""
    return sum(
        coefficient
        * math.comb(t_degree, power_of_n)
        * RESIDUE ** (t_degree - power_of_n)
        * MODULUS**power_of_n
        for t_degree, coefficient in enumerate(coefficients_in_t)
        if t_degree >= power_of_n
    )


def main() -> None:
    family = json.loads((DATA / "optimal_23_4_Z.json").read_text())
    base = json.loads((DATA / "specialization_t2_Z.json").read_text())
    progression = json.loads((DATA / "progression_319.json").read_text())
    table = [
        [int(coefficient) for coefficient in row]
        for row in family["coefficients_T_then_W"]
    ]
    base_coefficients = [int(value) for value in base["coefficients_low_to_high"]]

    if MODULUS != math.prod(PRIMES) or progression["modulus"] != MODULUS:
        raise RuntimeError("incorrect progression modulus")
    if progression["residue"] != RESIDUE or progression["witness_primes"] != list(PRIMES):
        raise RuntimeError("progression metadata are inconsistent")

    evaluated_at_base = [
        sum(table[t_degree][x_degree] * RESIDUE**t_degree for t_degree in range(5))
        for x_degree in range(24)
    ]
    if evaluated_at_base != base_coefficients:
        raise RuntimeError("stored base specialization is not P(2,x)")

    # Prove P(2+319*n,x)-P(2,x) is zero modulo 11 and modulo 29 as a
    # polynomial in the two indeterminates n and x.
    for prime in PRIMES:
        if MODULUS % prime:
            raise RuntimeError(f"{prime} does not divide the progression modulus")
        for x_degree in range(24):
            coefficients_in_t = [table[t_degree][x_degree] for t_degree in range(5)]
            constant = coefficient_after_affine_substitution(coefficients_in_t, 0)
            if constant != base_coefficients[x_degree]:
                raise RuntimeError("constant term after substitution is incorrect")
            for power_of_n in range(1, 5):
                transformed = coefficient_after_affine_substitution(
                    coefficients_in_t, power_of_n
                )
                if transformed % prime:
                    raise RuntimeError(
                        f"universal congruence failed modulo {prime} "
                        f"at x^{x_degree} n^{power_of_n}"
                    )

    irreducibility = base["modular_certificates"]["irreducibility"]
    prime_irreducible = irreducibility["prime"]
    reduction_irreducible = monic(base_coefficients, prime_irreducible)
    if reduction_irreducible != irreducibility["monic_coefficients_low_to_high"]:
        raise RuntimeError("stored reduction modulo 29 is incorrect")
    if not is_irreducible(reduction_irreducible, prime_irreducible):
        raise RuntimeError("Rabin irreducibility test failed modulo 29")
    if base_coefficients[-1] % prime_irreducible == 0:
        raise RuntimeError("leading coefficient vanishes modulo 29")

    even = base["modular_certificates"]["even_frobenius"]
    prime_even = even["prime"]
    reduction_even = monic(base_coefficients, prime_even)
    product = [1]
    for factor in even["irreducible_factors_low_to_high"]:
        if not is_irreducible(factor, prime_even):
            raise RuntimeError("stored factor modulo 11 is reducible")
        product = multiply(product, factor, prime_even)
    if product != reduction_even:
        raise RuntimeError("stored factors do not multiply to the reduction modulo 11")
    if len(gcd_poly(reduction_even, derivative(reduction_even, prime_even), prime_even)) != 1:
        raise RuntimeError("reduction modulo 11 is not squarefree")
    if base_coefficients[-1] % prime_even == 0:
        raise RuntimeError("leading coefficient vanishes modulo 11")
    factor_degrees = [len(factor) - 1 for factor in even["irreducible_factors_low_to_high"]]
    if factor_degrees != [2, 7, 14] or math.lcm(*factor_degrees) != 14:
        raise RuntimeError("unexpected Frobenius data modulo 11")

    maximal = base["maximal_subgroup_certificate"]
    maximal_with_23 = [
        entry for entry in maximal["maximal_subgroups"] if entry["order"] % 23 == 0
    ]
    if maximal_with_23 != [{"structure": "23:11", "order": 253}]:
        raise RuntimeError("maximal-subgroup exclusion data are inconsistent")

    summary = {
        "status": "PASS_PROGRESSION_319_M23_CERTIFICATE",
        "parameterization": "t=2+319*n, n in Z",
        "primitive_part": "prim_x P(t,x)",
        "universal_polynomial_congruence_mod_11": True,
        "universal_polynomial_congruence_mod_29": True,
        "full_degree_throughout_progression": True,
        "content_invertible_mod_11_and_29": True,
        "irreducible_mod_29_throughout_progression": True,
        "frobenius_cycle_degrees_mod_29": [23],
        "squarefree_mod_11_throughout_progression": True,
        "frobenius_cycle_degrees_mod_11": factor_degrees,
        "frobenius_order_mod_11": 14,
        "ambient_specialization_group": "subgroup of M23",
        "maximal_subgroup_exclusion": "only 23:11 contains order 23, but it has odd order",
        "conclusion": "Gal(prim_x P(t,x)/Q)=M23 for every t congruent to 2 modulo 319",
    }
    output = ROOT / "verification" / "progression_319_summary.json"
    output.write_text(json.dumps(summary, indent=2) + "\n")
    print(summary["status"])
    print("parameterization=t=2+319*n universal_congruences=1")
    print("irreducible_mod_29=1 factor_degrees_mod_11=[2,7,14]")
    print("galois_group_of_every_primitive_part=M23")


if __name__ == "__main__":
    main()
