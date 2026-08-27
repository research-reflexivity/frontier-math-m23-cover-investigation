#!/usr/bin/env python3
"""Render the canonical JSON coefficient table as a deterministic GP input."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data" / "Fint_coefficients_Z.json"
OUTPUT = ROOT / "data" / "Fint_coefficients_Z.gp"


def t_polynomial(row: list[int]) -> str:
    terms: list[str] = []
    for degree, coefficient in enumerate(row):
        if coefficient == 0:
            continue
        if degree == 0:
            terms.append(str(coefficient))
        elif degree == 1:
            terms.append(f"({coefficient})*T")
        else:
            terms.append(f"({coefficient})*T^{degree}")
    return "+".join(terms) if terms else "0"


def render() -> str:
    table = json.loads(SOURCE.read_text(encoding="utf-8"))
    if not isinstance(table, list) or len(table) != 24:
        raise SystemExit("expected 24 rows indexed by V-degree")
    terms: list[str] = []
    for v_degree, raw_row in enumerate(table):
        if not isinstance(raw_row, list) or len(raw_row) != 9:
            raise SystemExit(f"row {v_degree} does not have 9 T-coefficients")
        if any(type(value) is not int for value in raw_row):
            raise SystemExit(f"row {v_degree} contains a non-integer coefficient")
        coefficient = t_polynomial(raw_row)
        if coefficient == "0":
            continue
        if v_degree == 0:
            terms.append(f"  ({coefficient})")
        elif v_degree == 1:
            terms.append(f"  ({coefficient})*V")
        else:
            terms.append(f"  ({coefficient})*V^{v_degree}")
    if not terms or table[23][8] == 0:
        raise SystemExit("coefficient table does not have full bidegree (23,8)")
    header = (
        "\\\\ Deterministically generated from Fint_coefficients_Z.json.\n"
        "\\\\ Canonical indices: outer row = V-degree, inner row = T-degree.\n"
        "\\\\ Repository notation: F is the integral HJLPPZ polynomial; "
        "Fhat = F/(T^2+23)^4.\n"
        "F = ("
    )
    return header + "+".join(term.strip() for term in terms) + ");\nFhat = F/(T^2+23)^4;\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    expected = render()
    if args.write:
        OUTPUT.write_text(expected, encoding="utf-8")
        print(f"WROTE {OUTPUT.relative_to(ROOT)}")
        return
    if not OUTPUT.is_file() or OUTPUT.read_text(encoding="utf-8") != expected:
        raise SystemExit(
            "stale generated GP input; run scripts/render_fint_gp.py --write"
        )
    print("PASS deterministic source-model GP rendering")


if __name__ == "__main__":
    main()
