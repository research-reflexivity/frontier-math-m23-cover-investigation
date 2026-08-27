#!/usr/bin/env python3
"""Emit a self-contained Magma certificate for the canonical quadric."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUTPUT = ROOT / "verification" / "verify_canonical_quadric.m"


def rational(entry: dict[str, int]) -> str:
    numerator = int(entry["numerator"])
    denominator = int(entry["denominator"])
    return str(numerator) if denominator == 1 else f"({numerator}/{denominator})"


def magma_rows(name: str, rows: list[list[str]]) -> str:
    body = ",\n".join("    [" + ",".join(row) + "]" for row in rows)
    return f"{name} := [\n{body}\n];\n"


def render() -> str:
    payload = json.loads((DATA / "canonical_quadric_Q.json").read_text())
    if payload.get("status") != "PASS_CANONICAL_QUADRIC_Q_HOLDOUT":
        raise SystemExit("canonical-quadric reconstruction does not carry PASS status")
    zlines = [
        line for line in (DATA / "Fint_coefficients_Z.json").read_text().splitlines()
        if line
    ]
    source = json.loads(zlines[-1])
    if len(source) != 24 or any(len(row) != 9 for row in source):
        raise SystemExit("unexpected source coefficient-table dimensions")
    vectors = payload["canonical_vectors"]
    if len(vectors) != 4 or any(len(row) != 88 for row in vectors):
        raise SystemExit("unexpected canonical basis dimensions")
    slots = payload["adjoint_slots_j_m_i"]
    if len(slots) != 88:
        raise SystemExit("unexpected adjoint slot count")
    coefficients = payload["quadric_coefficients"]
    if len(coefficients) != 10:
        raise SystemExit("unexpected quadric coefficient count")

    header = r'''// Self-contained exact Magma certificate for the canonical quadric.
// Generated deterministically by
// scripts/emit_canonical_quadric_magma_certificate.py.  Do not edit by hand.
//
// It checks the unique quadratic relation among the four reconstructed
// canonical adjoints in Q(T)[V]/(F), its determinant, and the resulting
// ruling field Q(sqrt(4873)) != Q(sqrt(-23)).

'''
    tables = "\n".join([
        magma_rows("F_V_THEN_T", [[str(int(x)) for x in row] for row in source]),
        magma_rows("ADJOINT_SLOTS", [[str(int(x)) for x in row] for row in slots]),
        magma_rows("CANONICAL_VECTORS", [[rational(x) for x in row] for row in vectors]),
        "QUADRIC_COEFFICIENTS := ["
        + ",".join(rational(x) for x in coefficients)
        + "];\n",
    ])
    program = r'''
Q := Rationals();
KT<T> := FunctionField(Q);
RV<V> := PolynomialRing(KT);

function PolynomialInT(coefficients)
    value := KT!0;
    for i in [1..#coefficients] do
        value +:= (KT!coefficients[i]) * T^(i - 1);
    end for;
    return value;
end function;

F := RV![ PolynomialInT(row) : row in F_V_THEN_T ];
Fhat := F / LeadingCoefficient(F);
assert Degree(Fhat) eq 23;
assert IsMonic(Fhat);
D := T^2 + 23;

A := [ RV!0 : a in [1..4] ];
for a in [1..4] do
    for index in [1..88] do
        j := ADJOINT_SLOTS[index][1];
        m := ADJOINT_SLOTS[index][2];
        i := ADJOINT_SLOTS[index][3];
        A[a] +:= (KT!CANONICAL_VECTORS[a][index]) * T^i * D^m * V^j;
    end for;
end for;

// The final four RREF-free columns form the identity matrix.
for a in [1..4] do
    for b in [1..4] do
        assert CANONICAL_VECTORS[a][84 + b] eq (a eq b select 1 else 0);
    end for;
end for;

pairs := [ <i,j> : i,j in [1..4] | i le j ];
products := [ (A[pair[1]] * A[pair[2]]) mod Fhat : pair in pairs ];
relation := RV!0;
for index in [1..10] do
    relation +:= (KT!QUADRIC_COEFFICIENTS[index]) * products[index];
end for;
assert relation eq 0;

M := ZeroMatrix(KT, 23, 10);
for row in [0..22] do
    for column in [1..10] do
        M[row + 1,column] := Coefficient(products[column], row);
    end for;
end for;
assert Rank(M) eq 9;
assert M * Matrix(KT, 10, 1, QUADRIC_COEFFICIENTS) eq ZeroMatrix(KT, 23, 1);
identity_and_uniqueness_ok :=
    relation eq 0 and Rank(M) eq 9 and
    M * Matrix(KT, 10, 1, QUADRIC_COEFFICIENTS) eq ZeroMatrix(KT, 23, 1);

B := ZeroMatrix(Q, 4, 4);
for index in [1..10] do
    i := pairs[index][1];
    j := pairs[index][2];
    coefficient := Q!QUADRIC_COEFFICIENTS[index];
    if i eq j then
        B[i,j] := 2 * coefficient;
    else
        B[i,j] := coefficient;
        B[j,i] := coefficient;
    end if;
end for;
determinant := Determinant(B);
assert determinant eq 644454138716416151027888/970299;
assert determinant eq 4873 * (38141181796/3267)^2;
assert 4873 eq 11 * 443;
assert not IsSquare(Q!(-4873/23));
ruling_field_ok :=
    determinant eq 644454138716416151027888/970299 and
    determinant eq 4873 * (38141181796/3267)^2 and
    4873 eq 11 * 443 and not IsSquare(Q!(-4873/23));
if identity_and_uniqueness_ok and ruling_field_ok then
    print "PASS_MAGMA_CANONICAL_QUADRIC_IDENTITY_AND_UNIQUENESS";
    print "PASS_MAGMA_RULING_FIELD_Q_SQRT_4873_DIFFERS_FROM_Q_SQRT_MINUS_23";
    print "PASS_CANONICAL_QUADRIC_MAGMA_CERTIFICATE";
end if;
'''
    return header + tables + program


def main() -> None:
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    expected = render()
    if args.write:
        OUTPUT.write_text(expected)
        print(f"WROTE {OUTPUT.relative_to(ROOT)} ({len(expected.encode())} bytes)")
        return
    if not OUTPUT.is_file() or OUTPUT.read_text() != expected:
        raise SystemExit(
            "stale generated Magma certificate; run "
            "scripts/emit_canonical_quadric_magma_certificate.py --write"
        )
    print("PASS deterministic canonical-quadric Magma certificate rendering")


if __name__ == "__main__":
    main()
