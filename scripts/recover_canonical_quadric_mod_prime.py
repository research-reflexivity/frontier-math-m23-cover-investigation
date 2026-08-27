#!/usr/bin/env python3
"""Recover the RREF canonical basis and its quadric at one good prime.

This is the canonical-only part of the modular optimal-model computation.
It deliberately stops before the much more expensive degree-(23,4)
elimination.  The quadric is first found from rational points; the companion
Sage verifier checks the resulting identity in F_p(T)[V]/(F).
"""

from __future__ import annotations

import argparse
import json
import subprocess
import tempfile
from pathlib import Path

import recover_canonical_mod31 as core

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"


def terms_to_singular(table: list[list[int]], p: int) -> str:
    pieces: list[str] = []
    for j in range(len(table) - 1, -1, -1):
        for i in range(len(table[j]) - 1, -1, -1):
            c = table[j][i] % p
            if not c:
                continue
            factors = [str(c)]
            if i:
                factors.append("T" if i == 1 else f"T^{i}")
            if j:
                factors.append("V" if j == 1 else f"V^{j}")
            pieces.append("*".join(factors))
    return " + ".join(pieces) if pieces else "0"


def run_singular(script: Path) -> str:
    result = subprocess.run(
        ["Singular", "-q", script.name],
        cwd=script.parent,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if result.returncode != 0 or "PASS_NODE_GRAPH" not in result.stdout:
        raise RuntimeError(f"Singular failed ({result.returncode}):\n{result.stdout}")
    return result.stdout.strip().splitlines()[-1]


def evaluate_table(table: list[list[int]], t: int, v: int, p: int) -> int:
    answer = 0
    for row in reversed(table):
        coefficient = 0
        for c in reversed(row):
            coefficient = (coefficient * t + c) % p
        answer = (answer * v + coefficient) % p
    return answer


def evaluate_adjoint(
    vector: list[int], slots: list[tuple[int, int, int]], t: int, v: int, p: int
) -> int:
    d = (t * t + 23) % p
    answer = 0
    for c, (j, m, i) in zip(vector, slots):
        answer = (answer + c * pow(v, j, p) * pow(d, m, p) * pow(t, i, p)) % p
    return answer


def determinant(matrix: list[list[int]], p: int) -> int:
    a = [[x % p for x in row] for row in matrix]
    answer = 1
    for c in range(len(a)):
        pivot = next((r for r in range(c, len(a)) if a[r][c]), None)
        if pivot is None:
            return 0
        if pivot != c:
            a[c], a[pivot] = a[pivot], a[c]
            answer = -answer
        q = a[c][c]
        answer = answer * q % p
        qinv = pow(q, p - 2, p)
        for r in range(c + 1, len(a)):
            scale = a[r][c] * qinv % p
            for j in range(c, len(a)):
                a[r][j] = (a[r][j] - scale * a[c][j]) % p
    return answer % p


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("prime", type=int)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--basis-only",
        action="store_true",
        help="recover the canonical RREF basis but skip O(p^2) point enumeration",
    )
    args = parser.parse_args()
    p = args.prime
    if p in (2, 3, 23):
        raise SystemExit("choose a prime away from 2, 3, and 23")
    core.P = p

    zlines = [
        line for line in (DATA / "Fint_coefficients_Z.json").read_text().splitlines()
        if line
    ]
    ztable = json.loads(zlines[-1])
    table = [[int(c) % p for c in row] for row in ztable]

    with tempfile.TemporaryDirectory(prefix=f"m23_quadric_{p}_") as dirname:
        work = Path(dirname)
        (work / "Fp.inc").write_text(f"poly f = {terms_to_singular(table, p)};\n")
        script = work / "nodes.sing"
        script.write_text(f'''option(redSB);
ring r = {p},(V,T),lp;
< "Fp.inc";
ideal J = f,diff(f,T),diff(f,V);
ideal G = std(J);
if (size(G) != 2) {{ ERROR("node ideal is not a two-generator graph"); }}
poly h84 = G[1];
poly graph = G[2];
if (subst(h84,V,0) != h84) {{ h84=G[2]; graph=G[1]; }}
poly s = -subst(graph,V,0);
if (deg(h84,intvec(0,1)) != 84) {{ ERROR("node eliminant degree is not 84"); }}
if (gcd(h84,diff(h84,T)) != 1) {{ ERROR("node eliminant is not squarefree"); }}
poly Hess = diff(diff(f,T),T)*diff(diff(f,V),V)-diff(diff(f,T),V)^2;
poly HessNodes = reduce(Hess,G);
if (gcd(HessNodes,h84) != 1) {{ ERROR("not all singularities are ordinary nodes"); }}
write(":w node_graph.inc","poly h84 = "+string(h84)+";");
write(":a node_graph.inc","poly s = "+string(s)+";");
print("PASS_NODE_GRAPH degree=84 ordinary=1 squarefree=1");
exit;
''')
        node_log = run_singular(script)
        lines = (work / "node_graph.inc").read_text().splitlines()

    h = core.parse_singular_univariate(lines[0].split("=", 1)[1])
    s = core.parse_singular_univariate(lines[1].split("=", 1)[1])
    if len(h) != 85 or h[-1] != 1:
        raise RuntimeError("unexpected node polynomial")
    s = core.mod_poly(s, h)

    slots = core.canonical_slots()
    d = [23 % p, 0, 1]
    s_powers = [[1] + [0] * 83]
    for _ in range(21):
        s_powers.append(core.mul_mod(s_powers[-1], s, h))
    d_powers = [core.mod_poly(core.poly_pow(d, m), h) for m in range(4)]
    columns: list[list[int]] = []
    for j, m, i in slots:
        factor = core.mul_mod(d_powers[m], [0] * i + [1], h)
        columns.append(core.mul_mod(factor, s_powers[j], h))
    evaluation = [[columns[c][r] for c in range(88)] for r in range(84)]
    rank, canonical, pivots = core.rref_nullspace(evaluation)
    if rank != 84 or len(canonical) != 4 or pivots != list(range(84)):
        raise RuntimeError(f"unexpected canonical RREF data: rank={rank}, pivots={pivots}")

    payload = {
        "status": "PASS_CANONICAL_BASIS_MOD_PRIME",
        "prime": p,
        "node_log": node_log,
        "canonical_rank": rank,
        "canonical_rref_pivots": pivots,
        "adjoint_slots_j_m_i": [list(x) for x in slots],
        "canonical_vectors": canonical,
    }
    if not args.basis_only:
        pairs = [(i, j) for i in range(4) for j in range(i, 4)]
        quadric_rows: list[list[int]] = []
        curve_points = 0
        for t in range(p):
            for v in range(p):
                if evaluate_table(table, t, v, p):
                    continue
                curve_points += 1
                values = [evaluate_adjoint(a, slots, t, v, p) for a in canonical]
                if any(values):
                    quadric_rows.append([values[i] * values[j] % p for i, j in pairs])
        quadric_rank, quadrics, quadric_pivots = core.rref_nullspace(quadric_rows)
        if quadric_rank != 9 or len(quadrics) != 1:
            raise RuntimeError(
                f"rational-point evaluation did not isolate one quadric: rank={quadric_rank}"
            )
        q = quadrics[0]
        symmetric_twice = [[0] * 4 for _ in range(4)]
        for c, (i, j) in zip(q, pairs):
            if i == j:
                symmetric_twice[i][j] = 2 * c % p
            else:
                symmetric_twice[i][j] = c
                symmetric_twice[j][i] = c
        det = determinant(symmetric_twice, p)
        payload.update({
            "status": "CANDIDATE_CANONICAL_QUADRIC_MOD_PRIME",
            "affine_curve_points": curve_points,
            "usable_canonical_images": len(quadric_rows),
            "quadric_rank_from_points": quadric_rank,
            "quadric_rref_pivots": quadric_pivots,
            "quadric_pairs": [list(x) for x in pairs],
            "quadric_coefficients": q,
            "twice_symmetric_matrix": symmetric_twice,
            "twice_symmetric_determinant": det,
            "determinant_legendre": pow(det, (p - 1) // 2, p),
            "minus_23_legendre": pow(-23 % p, (p - 1) // 2, p),
        })
    serialized = json.dumps(payload, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    else:
        print(serialized, end="")


if __name__ == "__main__":
    main()
