#!/usr/bin/env python3
"""Run the compact exact certificate for the reconstructed (23,4) model."""

from __future__ import annotations

import json
import math
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"


def run(command: list[str]) -> str:
    result = subprocess.run(
        command,
        cwd=HERE,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if result.returncode != 0 or "FAIL_" in result.stdout:
        raise RuntimeError(f"certificate command failed: {command}\n{result.stdout}")
    return result.stdout


def main() -> None:
    singular = os.environ.get("SINGULAR", "Singular")
    equation = json.loads((DATA / "optimal_23_4_Z.json").read_text())
    table = [[int(entry) for entry in row] for row in equation["coefficients_T_then_W"]]
    content = 0
    for row in table:
        for entry in row:
            content = math.gcd(content, abs(entry))
    if content != 1 or table[4][23] == 0:
        raise RuntimeError("equation is not primitive of full bidegree")

    modular = json.loads((DATA / "modular" / "f4_mod_31.json").read_text())
    scale = pow(table[4][23], -1, 31)
    if [[entry * scale % 31 for entry in row] for row in table] != modular["coefficients_T_then_W"]:
        raise RuntimeError("Z-equation does not reduce to the certified mod-31 equation")

    with ThreadPoolExecutor(max_workers=2) as executor:
        exact_future = executor.submit(
            run, [singular, "--cpus=1", "--threads=1", "-q", "verify_function_identity_Q.sing"]
        )
        finite_future = executor.submit(
            run, [singular, "--cpus=1", "--threads=1", "-q", "verify_f4_mod31.sing"]
        )
        exact_log = exact_future.result()
        finite_log = finite_future.result()
    if "PASS_EXACT_FUNCTION_FIELD_IDENTITY" not in exact_log:
        raise RuntimeError(exact_log)
    if "PASS_F4_MOD31_CERTIFICATE" not in finite_log:
        raise RuntimeError(finite_log)

    payload = {
        "status": "PASS_OPTIMAL_23_4_EXACT_CERTIFICATE",
        "primitive_Z_polynomial": True,
        "degree_W": 23,
        "degree_T": 4,
        "exact_identity_in_QT_mod_F": True,
        "reduction_prime": 31,
        "irreducible_mod_31": True,
        "divides_defining_resultant_mod_31": True,
        "consequence": (
            "F23_4 is irreducible over Q and Q(T,W)=Q(T,V) under "
            "W=J0/J1; hence it defines the same degree-23 regular M23 cover."
        ),
        "exact_log": exact_log.strip().splitlines(),
        "finite_log": finite_log.strip().splitlines(),
    }
    (HERE / "verification_summary.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(payload["status"])
    print(f"degree_W={payload['degree_W']} degree_T={payload['degree_T']}")
    print("exact_identity=1 irreducible_over_Q=1 same_function_field=1")


if __name__ == "__main__":
    main()
