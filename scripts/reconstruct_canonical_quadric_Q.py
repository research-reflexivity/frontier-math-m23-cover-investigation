#!/usr/bin/env python3
"""CRT-lift the canonical RREF basis and its unique quadric to Q."""

from __future__ import annotations

import json
import math
import subprocess
import sys
import tempfile
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RECOVER = ROOT / "scripts" / "recover_canonical_quadric_mod_prime.py"
OUTPUT = ROOT / "data" / "canonical_quadric_Q.json"

SMALL_PRIMES = [
    29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
    79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131,
]
LARGE_PRIMES = [
    500000003, 500001001, 500002003, 500003017,
    500004073, 500005039, 500006011, 500007007,
    500008007, 500009089, 500010029, 500011003,
    500012021, 500013001, 500014003, 500015017,
]
HOLDOUT_PRIME = 137


def rational_reconstruct(a: int, modulus: int) -> Fraction | None:
    a %= modulus
    bound = math.isqrt(modulus // 2)
    r0, r1 = modulus, a
    s0, s1 = 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound or math.gcd(r1, s1) != 1:
        return None
    numerator, denominator = r1, s1
    if denominator < 0:
        numerator, denominator = -numerator, -denominator
    if (a * denominator - numerator) % modulus:
        return None
    return Fraction(numerator, denominator)


def reduce_fraction(value: Fraction, p: int) -> int | None:
    if value.denominator % p == 0:
        return None
    return value.numerator * pow(value.denominator, -1, p) % p


def crt_tables(
    tables: list[tuple[int, list[list[int]]]], rows: int, columns: int
) -> tuple[int, list[list[Fraction]]]:
    residues = [[0] * columns for _ in range(rows)]
    modulus = 1
    for p, table in tables:
        old_modulus = modulus
        modulus *= p
        for i in range(rows):
            for j in range(columns):
                residues[i][j] = (
                    residues[i][j]
                    + (
                        (table[i][j] - residues[i][j])
                        * pow(old_modulus, -1, p)
                        % p
                    )
                    * old_modulus
                ) % modulus
    reconstructed = [
        [rational_reconstruct(entry, modulus) for entry in row]
        for row in residues
    ]
    if any(entry is None for row in reconstructed for entry in row):
        raise RuntimeError("CRT modulus is too small for rational reconstruction")
    return modulus, [[entry for entry in row if entry is not None] for row in reconstructed]


def determinant(matrix: list[list[Fraction]]) -> Fraction:
    a = [row[:] for row in matrix]
    answer = Fraction(1)
    for column in range(len(a)):
        pivot = next((row for row in range(column, len(a)) if a[row][column]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            a[column], a[pivot] = a[pivot], a[column]
            answer = -answer
        q = a[column][column]
        answer *= q
        for row in range(column + 1, len(a)):
            scale = a[row][column] / q
            for j in range(column, len(a)):
                a[row][j] -= scale * a[column][j]
    return answer


def encode(value: Fraction) -> dict[str, int]:
    return {"numerator": value.numerator, "denominator": value.denominator}


def main() -> None:
    train_primes = SMALL_PRIMES + LARGE_PRIMES
    with tempfile.TemporaryDirectory(prefix="m23_quadric_crt_") as dirname:
        work = Path(dirname)
        modular: dict[int, dict[str, object]] = {}
        for p in train_primes + [HOLDOUT_PRIME]:
            target = work / f"canonical_{p}.json"
            command = [sys.executable, str(RECOVER), str(p), "--output", str(target)]
            if p in LARGE_PRIMES:
                command.append("--basis-only")
            subprocess.run(command, cwd=ROOT, check=True)
            modular[p] = json.loads(target.read_text())

    basis_tables = [
        (p, modular[p]["canonical_vectors"])  # type: ignore[arg-type]
        for p in train_primes
    ]
    basis_modulus, basis = crt_tables(basis_tables, 4, 88)

    # The point-normalized quadric has its X_3^2 coefficient equal to one.
    quadric_tables = [
        (p, [modular[p]["quadric_coefficients"]])  # type: ignore[list-item]
        for p in SMALL_PRIMES
    ]
    quadric_modulus, quadric_rows = crt_tables(quadric_tables, 1, 10)
    quadric = quadric_rows[0]

    holdout_basis = modular[HOLDOUT_PRIME]["canonical_vectors"]
    holdout_quadric = modular[HOLDOUT_PRIME]["quadric_coefficients"]
    basis_mismatches = [
        [i, j]
        for i in range(4)
        for j in range(88)
        if reduce_fraction(basis[i][j], HOLDOUT_PRIME) != holdout_basis[i][j]
    ]
    quadric_mismatches = [
        i
        for i in range(10)
        if reduce_fraction(quadric[i], HOLDOUT_PRIME) != holdout_quadric[i]
    ]
    if basis_mismatches or quadric_mismatches:
        raise RuntimeError("independent holdout prime does not match reconstruction")

    pairs = [(i, j) for i in range(4) for j in range(i, 4)]
    twice_symmetric = [[Fraction(0) for _ in range(4)] for _ in range(4)]
    for coefficient, (i, j) in zip(quadric, pairs):
        if i == j:
            twice_symmetric[i][j] = 2 * coefficient
        else:
            twice_symmetric[i][j] = coefficient
            twice_symmetric[j][i] = coefficient
    det = determinant(twice_symmetric)
    squarefree_part = 4873
    quotient = det / squarefree_part
    square_root = Fraction(math.isqrt(quotient.numerator), math.isqrt(quotient.denominator))
    if square_root * square_root != quotient:
        raise RuntimeError("unexpected canonical-quadric discriminant square class")

    nonzero_basis = [entry for row in basis for entry in row if entry]
    payload = {
        "status": "PASS_CANONICAL_QUADRIC_Q_HOLDOUT",
        "basis_primes": train_primes,
        "basis_modulus": str(basis_modulus),
        "basis_modulus_bits": basis_modulus.bit_length(),
        "quadric_primes": SMALL_PRIMES,
        "quadric_modulus": str(quadric_modulus),
        "quadric_modulus_bits": quadric_modulus.bit_length(),
        "holdout_prime": HOLDOUT_PRIME,
        "basis_holdout_mismatch_count": len(basis_mismatches),
        "quadric_holdout_mismatch_count": len(quadric_mismatches),
        "max_basis_numerator_digits": max(len(str(abs(x.numerator))) for x in nonzero_basis),
        "max_basis_denominator_digits": max(len(str(x.denominator)) for x in nonzero_basis),
        "adjoint_slots_j_m_i": modular[HOLDOUT_PRIME]["adjoint_slots_j_m_i"],
        "canonical_vectors": [[encode(x) for x in row] for row in basis],
        "quadric_pairs": [list(pair) for pair in pairs],
        "quadric_coefficients": [encode(x) for x in quadric],
        "twice_symmetric_matrix": [[encode(x) for x in row] for row in twice_symmetric],
        "twice_symmetric_determinant": encode(det),
        "determinant_squarefree_part": squarefree_part,
        "determinant_over_squarefree_part_sqrt": encode(square_root),
        "ruling_field": "Q(sqrt(4873)) = Q(sqrt(11*443))",
        "branch_orientation_field": "Q(sqrt(-23))",
        "fields_equal": False,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({
        key: payload[key]
        for key in (
            "status", "basis_modulus_bits", "quadric_modulus_bits",
            "holdout_prime", "basis_holdout_mismatch_count",
            "quadric_holdout_mismatch_count", "determinant_squarefree_part",
            "ruling_field", "fields_equal",
        )
    }, indent=2))


if __name__ == "__main__":
    main()
