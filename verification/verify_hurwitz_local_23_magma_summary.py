#!/usr/bin/env python3
"""Bind the recorded Sydney Magma run to the local-23 certificate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
CERTIFICATE = HERE / "verify_hurwitz_local_23.m"
SUMMARY = HERE / "hurwitz_local_23_magma_summary.json"
EXPECTED_MARKERS = [
    "PASS_MAGMA_HURWITZ_LOCAL_PRIME_DATA_1_PLUS_2_PLUS_4",
    "PASS_MAGMA_HURWITZ_LOCAL_DECOMPOSITION_GROUP_V4",
    "PASS_HURWITZ_LOCAL_23_MAGMA_CERTIFICATE",
]


def main():
    payload = json.loads(SUMMARY.read_text())
    raw = CERTIFICATE.read_bytes()
    assert payload["engine"] == "Magma"
    assert payload["magma_version"] == "2.29-9"
    assert payload["runner"] == "University of Sydney Magma Calculator"
    assert payload["status"] == "PASS_INDEPENDENT_MAGMA_CERTIFICATE"
    assert payload["markers"] == EXPECTED_MARKERS
    assert payload["input_path"] == "verification/verify_hurwitz_local_23.m"
    assert payload["input_bytes"] == len(raw)
    assert payload["input_sha256"] == hashlib.sha256(raw).hexdigest()
    assert payload["total_seconds"] > 0
    assert payload["memory_megabytes"] > 0
    print("PASS recorded local-23 Magma run matches the input hash")


if __name__ == "__main__":
    main()
