#!/usr/bin/env python3
"""Recover the genus-4 canonical system and optimal pencil modulo 31.

The bidegree-(23,8) plane model has 84 affine nodes.  Its toric adjoint
space has dimension 154.  Regularity at the two totally ramified boundary
points cuts this to an explicit 88-dimensional space; vanishing at the 84
nodes should then leave the four canonical differentials.
"""

from __future__ import annotations

import json
import re
from itertools import combinations_with_replacement
from pathlib import Path

P = 31
HERE = Path(__file__).resolve().parent.parent / "data"


def inv(a: int) -> int:
    return pow(a % P, P - 2, P)


def trim(a: list[int]) -> list[int]:
    while len(a) > 1 and a[-1] % P == 0:
        a.pop()
    return [x % P for x in a]


def add(a: list[int], b: list[int]) -> list[int]:
    c = [0] * max(len(a), len(b))
    for i, x in enumerate(a):
        c[i] += x
    for i, x in enumerate(b):
        c[i] += x
    return trim(c)


def mul(a: list[int], b: list[int]) -> list[int]:
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] = (c[i + j] + x * y) % P
    return trim(c)


def mod_poly(a: list[int], h: list[int]) -> list[int]:
    a = trim(a[:])
    n = len(h) - 1
    assert h[-1] == 1
    while len(a) - 1 >= n:
        c = a[-1]
        k = len(a) - 1 - n
        if c:
            for i in range(n + 1):
                a[i + k] = (a[i + k] - c * h[i]) % P
        a = trim(a)
    return a + [0] * (n - len(a))


def mul_mod(a: list[int], b: list[int], h: list[int]) -> list[int]:
    return mod_poly(mul(a, b), h)


def parse_singular_univariate(expr: str) -> list[int]:
    expr = expr.strip().rstrip(";").replace("-", "+-")
    out: dict[int, int] = {}
    for raw in expr.split("+"):
        term = raw.strip()
        if not term:
            continue
        if "T" not in term:
            coeff, exponent = int(term), 0
        else:
            left, right = term.split("T", 1)
            if left in ("", "+"):
                coeff = 1
            elif left == "-":
                coeff = -1
            else:
                coeff = int(left)
            exponent = int(right) if right else 1
        out[exponent] = (out.get(exponent, 0) + coeff) % P
    ans = [0] * (max(out) + 1)
    for exponent, coeff in out.items():
        ans[exponent] = coeff
    return trim(ans)


def rref_nullspace(matrix: list[list[int]]) -> tuple[int, list[list[int]], list[int]]:
    if not matrix:
        return 0, [], []
    a = [[x % P for x in row] for row in matrix]
    rows, cols = len(a), len(a[0])
    pivots: list[int] = []
    r = 0
    for c in range(cols):
        pivot = next((i for i in range(r, rows) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        scale = inv(a[r][c])
        a[r] = [(scale * x) % P for x in a[r]]
        for i in range(rows):
            if i != r and a[i][c]:
                q = a[i][c]
                a[i] = [(x - q * y) % P for x, y in zip(a[i], a[r])]
        pivots.append(c)
        r += 1
        if r == rows:
            break
    free = [c for c in range(cols) if c not in pivots]
    basis: list[list[int]] = []
    for f in free:
        v = [0] * cols
        v[f] = 1
        for i, c in enumerate(pivots):
            v[c] = (-a[i][f]) % P
        basis.append(v)
    return len(pivots), basis, pivots


def poly_pow(a: list[int], n: int) -> list[int]:
    out = [1]
    for _ in range(n):
        out = mul(out, a)
    return out


def canonical_slots() -> list[tuple[int, int, int]]:
    slots: list[tuple[int, int, int]] = []
    for j in range(0, 5):
        slots.extend((j, 0, i) for i in range(7))
    for j in range(5, 11):
        slots.extend((j, 1, i) for i in range(5))
    for j in range(11, 17):
        slots.extend((j, 2, i) for i in range(3))
    for j in range(17, 22):
        slots.append((j, 3, 0))
    assert len(slots) == 88
    return slots


def expand_numerator(v: list[int], slots: list[tuple[int, int, int]]) -> dict[tuple[int, int], int]:
    # Dictionary (V exponent, T exponent) -> coefficient, expanding D^m.
    out: dict[tuple[int, int], int] = {}
    d = [23, 0, 1]
    powers = [poly_pow(d, m) for m in range(4)]
    for c, (j, m, i) in zip(v, slots):
        if not c:
            continue
        for k, q in enumerate(powers[m]):
            if q:
                key = (j, i + k)
                out[key] = (out.get(key, 0) + c * q) % P
    return {k: c for k, c in out.items() if c}


def singular_poly(terms: dict[tuple[int, int], int], vname: str = "V") -> str:
    pieces: list[str] = []
    for (j, i), c in sorted(terms.items(), key=lambda item: (item[0][0], item[0][1]), reverse=True):
        mon = []
        if i:
            mon.append("T" if i == 1 else f"T^{i}")
        if j:
            mon.append(vname if j == 1 else f"{vname}^{j}")
        m = "*".join(mon) if mon else "1"
        pieces.append(f"{c}*{m}")
    return " + ".join(pieces) if pieces else "0"


def main() -> None:
    lines = (HERE / "node_graph_mod31.inc").read_text().splitlines()
    h = parse_singular_univariate(lines[0].split("=", 1)[1])
    s = parse_singular_univariate(lines[1].split("=", 1)[1])
    assert len(h) == 85 and h[-1] == 1
    s = mod_poly(s, h)

    d = [23, 0, 1]
    slots = canonical_slots()
    s_powers = [[1] + [0] * 83]
    for _ in range(21):
        s_powers.append(mul_mod(s_powers[-1], s, h))
    d_powers = [mod_poly(poly_pow(d, m), h) for m in range(4)]

    columns: list[list[int]] = []
    for j, m, i in slots:
        factor = mul_mod(d_powers[m], [0] * i + [1], h)
        columns.append(mul_mod(factor, s_powers[j], h))
    matrix = [[columns[c][r] for c in range(88)] for r in range(84)]
    rank, canonical, canonical_pivots = rref_nullspace(matrix)
    assert rank == 84 and len(canonical) == 4

    # L(K-b-c): the only possible order-zero boundary contribution is the
    # D^2*c(T)*V^16 slot.  Vanishing at T=15 and T=16 imposes two conditions.
    slot_index = {slot: i for i, slot in enumerate(slots)}
    eval_rows: list[list[int]] = []
    # If c(T)=c0+c1*T+c2*T^2 is the V^16/D^2 boundary coefficient,
    # vanishing at both conjugate points is exactly c1=0, c0=23*c2.
    for kind in ("odd", "constant"):
        row = []
        for vec in canonical:
            c0 = vec[slot_index[(16, 2, 0)]]
            c1 = vec[slot_index[(16, 2, 1)]]
            c2 = vec[slot_index[(16, 2, 2)]]
            row.append(c1 if kind == "odd" else (c0 - 23*c2) % P)
        eval_rows.append(row)
    pencil_rank, pencil_combinations, pencil_pivots = rref_nullspace(eval_rows)
    assert pencil_rank == 2 and len(pencil_combinations) == 2
    pencil: list[list[int]] = []
    for comb in pencil_combinations:
        pencil.append([
            sum(comb[k] * canonical[k][i] for k in range(4)) % P
            for i in range(88)
        ])

    canonical_terms = [expand_numerator(v, slots) for v in canonical]
    pencil_terms = [expand_numerator(v, slots) for v in pencil]
    payload = {
        "status": "PASS_CANONICAL_4_AND_OPTIMAL_PENCIL_2_MOD_31",
        "prime": P,
        "node_scheme_degree": 84,
        "boundary_regular_adjoint_dimension": 88,
        "node_evaluation_rank": rank,
        "canonical_dimension": len(canonical),
        "canonical_rref_pivots": canonical_pivots,
        "boundary_evaluation_rank_on_canonical_space": pencil_rank,
        "degree4_pencil_dimension": len(pencil),
        "pencil_rref_pivots": pencil_pivots,
        "slots": [list(x) for x in slots],
        "canonical_vectors": canonical,
        "pencil_vectors": pencil,
        "canonical_numerators": [singular_poly(x) for x in canonical_terms],
        "pencil_numerators": [singular_poly(x) for x in pencil_terms],
        "theorem": (
            "For omega=A(T,V)dT/F_V, toric interior support gives 154 candidates. "
            "At each of the two total-ramification points ord(T-t0)=23, ord(V)=-4, "
            "and ord(dT/F_V)=18. Boundary regularity is therefore equivalent to the "
            "displayed D^m divisibilities, leaving 88 dimensions. The 84 reduced affine "
            "nodes impose rank 84, leaving H^0(K) of dimension 4. Vanishing at both "
            "boundary points has rank 2, producing the two-dimensional degree-4 pencil."
        ),
    }
    (HERE / "canonical_pencil_mod31.json").write_text(json.dumps(payload, indent=2) + "\n")
    (HERE / "canonical_pencil_mod31.inc").write_text(
        "\n".join([
            *(f"poly A{i} = {singular_poly(x)};" for i, x in enumerate(canonical_terms)),
            *(f"poly J{i} = {singular_poly(x)};" for i, x in enumerate(pencil_terms)),
        ]) + "\n"
    )
    print(payload["status"])
    print(f"rank={rank} canonical_dimension={len(canonical)}")
    print(f"pencil_rank={pencil_rank} pencil_dimension={len(pencil)}")


if __name__ == "__main__":
    main()
