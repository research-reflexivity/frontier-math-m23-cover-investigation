#!/usr/bin/env python3
"""Recompute the 88-to-4 adjoint kernel and optimal pencil modulo 31."""

from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
import recover_canonical_mod31 as core


def main():
    core.P = 31
    lines = (ROOT / "data" / "node_graph_mod31.inc").read_text().splitlines()
    h = core.parse_singular_univariate(lines[0].split("=", 1)[1])
    s = core.mod_poly(core.parse_singular_univariate(lines[1].split("=", 1)[1]), h)
    assert len(h) == 85 and h[-1] == 1

    slots = core.canonical_slots()
    assert len(slots) == 88
    d = [23, 0, 1]
    s_powers = [[1] + [0] * 83]
    for _ in range(21):
        s_powers.append(core.mul_mod(s_powers[-1], s, h))
    d_powers = [core.mod_poly(core.poly_pow(d, exponent), h) for exponent in range(4)]
    columns = []
    for j, exponent, i in slots:
        factor = core.mul_mod(d_powers[exponent], [0] * i + [1], h)
        columns.append(core.mul_mod(factor, s_powers[j], h))
    evaluation = [[columns[column][row] for column in range(88)] for row in range(84)]
    rank, canonical, pivots = core.rref_nullspace(evaluation)
    assert rank == 84 and len(canonical) == 4 and pivots == list(range(84))

    slot_index = {slot: index for index, slot in enumerate(slots)}
    boundary = []
    for kind in ("odd", "constant"):
        row = []
        for vector in canonical:
            c0 = vector[slot_index[(16, 2, 0)]]
            c1 = vector[slot_index[(16, 2, 1)]]
            c2 = vector[slot_index[(16, 2, 2)]]
            row.append(c1 if kind == "odd" else (c0 - 23 * c2) % 31)
        boundary.append(row)
    pencil_rank, pencil, pencil_pivots = core.rref_nullspace(boundary)
    assert pencil_rank == 2 and len(pencil) == 2 and pencil_pivots == [0, 1]

    result = {
        "status": "PASS_ADJOINT_AND_PENCIL_MOD31",
        "boundary_regular_adjoint_dimension": 88,
        "node_evaluation_rank": 84,
        "canonical_dimension": 4,
        "boundary_evaluation_rank": 2,
        "degree4_pencil_dimension": 2,
        "canonical_rref_pivots": pivots,
        "pencil_rref_pivots": pencil_pivots,
    }
    (ROOT / "verification" / "adjoint_verification_summary.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print(result["status"])
    print("adjoints=88 node_rank=84 canonical=4 boundary_rank=2 pencil=2")


if __name__ == "__main__":
    main()
