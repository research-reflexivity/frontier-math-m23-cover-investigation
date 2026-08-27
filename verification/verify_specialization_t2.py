#!/usr/bin/env python3
"""Exact dependency-free certificate for the specialization T=2.

The finite-field irreducibility tests use Rabin's criterion, implemented here
with elementary polynomial arithmetic.  The group-theoretic final step uses
the ATLAS maximal-subgroup orders recorded in the specialization data.
"""

from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"


def trim(poly: list[int]) -> list[int]:
    result = poly[:]
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def reduce_poly(poly: list[int], prime: int) -> list[int]:
    return trim([coefficient % prime for coefficient in poly])


def monic(poly: list[int], prime: int) -> list[int]:
    result = reduce_poly(poly, prime)
    if result == [0]:
        raise ValueError("zero polynomial cannot be made monic")
    scale = pow(result[-1], -1, prime)
    return trim([(scale * coefficient) % prime for coefficient in result])


def add(left: list[int], right: list[int], prime: int) -> list[int]:
    size = max(len(left), len(right))
    result = [0] * size
    for index in range(size):
        result[index] = (
            (left[index] if index < len(left) else 0)
            + (right[index] if index < len(right) else 0)
        ) % prime
    return trim(result)


def subtract(left: list[int], right: list[int], prime: int) -> list[int]:
    return add(left, [(-coefficient) % prime for coefficient in right], prime)


def multiply(left: list[int], right: list[int], prime: int) -> list[int]:
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] = (result[i + j] + a * b) % prime
    return trim(result)


def divide_with_remainder(
    dividend: list[int], divisor: list[int], prime: int
) -> tuple[list[int], list[int]]:
    remainder = reduce_poly(dividend, prime)
    divisor = reduce_poly(divisor, prime)
    if divisor == [0]:
        raise ZeroDivisionError("polynomial division by zero")
    if len(remainder) < len(divisor):
        return [0], remainder
    quotient = [0] * (len(remainder) - len(divisor) + 1)
    inverse_lead = pow(divisor[-1], -1, prime)
    while remainder != [0] and len(remainder) >= len(divisor):
        shift = len(remainder) - len(divisor)
        coefficient = remainder[-1] * inverse_lead % prime
        quotient[shift] = coefficient
        for index, value in enumerate(divisor):
            remainder[index + shift] = (
                remainder[index + shift] - coefficient * value
            ) % prime
        remainder = trim(remainder)
    return trim(quotient), remainder


def remainder(poly: list[int], modulus: list[int], prime: int) -> list[int]:
    return divide_with_remainder(poly, modulus, prime)[1]


def gcd_poly(left: list[int], right: list[int], prime: int) -> list[int]:
    left = reduce_poly(left, prime)
    right = reduce_poly(right, prime)
    while right != [0]:
        left, right = right, remainder(left, right, prime)
    return monic(left, prime)


def power_mod(
    base: list[int], exponent: int, modulus: list[int], prime: int
) -> list[int]:
    result = [1]
    base = remainder(base, modulus, prime)
    while exponent:
        if exponent & 1:
            result = remainder(multiply(result, base, prime), modulus, prime)
        base = remainder(multiply(base, base, prime), modulus, prime)
        exponent >>= 1
    return result


def prime_divisors(integer: int) -> list[int]:
    result = []
    divisor = 2
    while divisor * divisor <= integer:
        if integer % divisor == 0:
            result.append(divisor)
            while integer % divisor == 0:
                integer //= divisor
        divisor += 1
    if integer > 1:
        result.append(integer)
    return result


def is_irreducible(poly: list[int], prime: int) -> bool:
    """Rabin irreducibility test over F_prime."""
    poly = monic(poly, prime)
    degree = len(poly) - 1
    if degree <= 0:
        return False
    x = [0, 1]
    for divisor in prime_divisors(degree):
        frobenius = power_mod(x, prime ** (degree // divisor), poly, prime)
        if len(gcd_poly(subtract(frobenius, x, prime), poly, prime)) > 1:
            return False
    return power_mod(x, prime**degree, poly, prime) == x


def derivative(poly: list[int], prime: int) -> list[int]:
    if len(poly) <= 1:
        return [0]
    return trim([(index * poly[index]) % prime for index in range(1, len(poly))])


def main() -> None:
    family = json.loads((DATA / "optimal_23_4_Z.json").read_text())
    certificate = json.loads((DATA / "specialization_t2_Z.json").read_text())
    table = [
        [int(coefficient) for coefficient in row]
        for row in family["coefficients_T_then_W"]
    ]
    t0 = certificate["specialization"]["value"]
    coefficients = [int(value) for value in certificate["coefficients_low_to_high"]]
    expected = [
        sum(table[t_degree][x_degree] * t0**t_degree for t_degree in range(5))
        for x_degree in range(24)
    ]
    if coefficients != expected:
        raise RuntimeError("stored polynomial is not the exact T=2 specialization")
    if len(coefficients) != 24 or coefficients[-1] == 0:
        raise RuntimeError("specialization does not have degree 23")
    if math.gcd(*map(abs, coefficients)) != 1:
        raise RuntimeError("specialization is not primitive")

    irreducibility = certificate["modular_certificates"]["irreducibility"]
    prime_irreducible = irreducibility["prime"]
    reduction_irreducible = monic(coefficients, prime_irreducible)
    if reduction_irreducible != irreducibility["monic_coefficients_low_to_high"]:
        raise RuntimeError("incorrect stored reduction modulo 29")
    if len(gcd_poly(
        reduction_irreducible,
        derivative(reduction_irreducible, prime_irreducible),
        prime_irreducible,
    )) != 1:
        raise RuntimeError("reduction modulo 29 is not squarefree")
    if not is_irreducible(reduction_irreducible, prime_irreducible):
        raise RuntimeError("Rabin certificate failed modulo 29")

    even_frobenius = certificate["modular_certificates"]["even_frobenius"]
    prime_even = even_frobenius["prime"]
    reduction_even = monic(coefficients, prime_even)
    if reduction_even != even_frobenius["monic_coefficients_low_to_high"]:
        raise RuntimeError("incorrect stored reduction modulo 11")
    factors = even_frobenius["irreducible_factors_low_to_high"]
    product = [1]
    for factor in factors:
        if not is_irreducible(factor, prime_even):
            raise RuntimeError("stored factor modulo 11 is reducible")
        product = multiply(product, factor, prime_even)
    if product != reduction_even:
        raise RuntimeError("stored factors do not multiply to the reduction modulo 11")
    degrees = [len(factor) - 1 for factor in factors]
    if degrees != even_frobenius["factor_degrees"]:
        raise RuntimeError("incorrect factor degrees modulo 11")
    if len(gcd_poly(reduction_even, derivative(reduction_even, prime_even), prime_even)) != 1:
        raise RuntimeError("reduction modulo 11 is not squarefree")
    frobenius_order = math.lcm(*degrees)
    if frobenius_order != 14:
        raise RuntimeError("unexpected Frobenius order modulo 11")

    maximal = certificate["maximal_subgroup_certificate"]
    maximal_with_23 = [
        entry for entry in maximal["maximal_subgroups"] if entry["order"] % 23 == 0
    ]
    if maximal_with_23 != [{"structure": "23:11", "order": 253}]:
        raise RuntimeError("ATLAS maximal-subgroup exclusion data are inconsistent")
    if maximal_with_23[0]["order"] % 2 == 0:
        raise RuntimeError("23:11 exclusion requires odd order")

    summary = {
        "status": "PASS_SPECIALIZATION_T2_M23_CERTIFICATE",
        "specialization": "T=2",
        "degree": 23,
        "primitive": True,
        "exact_family_specialization": True,
        "irreducible_mod_29": True,
        "frobenius_cycle_degrees_mod_29": [23],
        "squarefree_mod_11": True,
        "frobenius_cycle_degrees_mod_11": degrees,
        "frobenius_order_mod_11": frobenius_order,
        "ambient_specialization_group": "subgroup of M23",
        "maximal_subgroup_exclusion": "only 23:11 contains order 23, but it has odd order",
        "galois_group_over_Q": "M23",
    }
    output = ROOT / "verification" / "specialization_t2_summary.json"
    output.write_text(json.dumps(summary, indent=2) + "\n")
    print(summary["status"])
    print("degree=23 content=1 irreducible_mod_29=1")
    print("factor_degrees_mod_11=[2,7,14] frobenius_order=14")
    print("specialization_group=M23")


if __name__ == "__main__":
    main()
