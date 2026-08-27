#!/usr/bin/env python3
"""Check the recorded long exact third-fibre certificate and its inputs."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SUMMARY_PATH = ROOT / "verification" / "hurwitz_degree23_third_fiber_sage_summary.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


summary = json.loads(SUMMARY_PATH.read_text())
assert summary["status"] == "passed_exact_characteristic_zero"
certificate = ROOT / summary["certificate"]
assert sha256(certificate) == summary["certificate_sha256"]
for relative_path, expected in summary["input_sha256"].items():
    assert sha256(ROOT / relative_path) == expected

results = summary["results"]
assert results == {
    "branch_gcd_degree": 14,
    "branch_resultant_degree": 42,
    "common_gcd_degree": 6,
    "extra_gcd_degree": 8,
    "extra_gcd_squarefree": True,
    "generic_gcd_degree": 6,
    "generic_resultant_degree": 42,
    "plane_curve_degree": 6,
    "plane_denominator_degree": 7,
    "plane_numerator_degree": 7,
}
print("PASS recorded exact characteristic-zero third-fibre certificate")
