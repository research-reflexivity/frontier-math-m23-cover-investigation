#!/usr/bin/env python3
"""Bind the recorded Sydney Magma run to the Galois-closure certificate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
CERTIFICATE = HERE / "verify_hurwitz_galois_closure.m"
SUMMARY = HERE / "hurwitz_galois_closure_magma_summary.json"
EXPECTED_MARKERS = [
    "PASS_MAGMA_HURWITZ_TRACE_SEXTIC_GALOIS_GROUP_S6",
    "PASS_MAGMA_HURWITZ_ABSOLUTE_FIELD_IS_TRACE_FIELD_TIMES_Q_SQRT_MINUS_23",
    "PASS_HURWITZ_GALOIS_CLOSURE_MAGMA_CERTIFICATE",
]


def main():
    payload = json.loads(SUMMARY.read_text())
    raw = CERTIFICATE.read_bytes()
    assert payload["engine"] == "Magma"
    assert payload["magma_version"] == "2.29-9"
    assert payload["runner"] == "University of Sydney Magma Calculator"
    assert payload["status"] == "PASS_INDEPENDENT_MAGMA_CERTIFICATE"
    assert payload["markers"] == EXPECTED_MARKERS
    assert payload["input_path"] == "verification/verify_hurwitz_galois_closure.m"
    assert payload["input_bytes"] == len(raw)
    assert payload["input_sha256"] == hashlib.sha256(raw).hexdigest()
    assert payload["total_seconds"] > 0
    assert payload["memory_megabytes"] > 0
    print("PASS recorded Hurwitz Galois-closure Magma run matches the input hash")


if __name__ == "__main__":
    main()
