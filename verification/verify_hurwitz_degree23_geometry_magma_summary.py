#!/usr/bin/env python3
"""Bind the recorded Sydney Magma run to the generated geometry certificate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
CERTIFICATE = HERE / "verify_hurwitz_degree23_geometry.m"
SUMMARY = HERE / "hurwitz_degree23_geometry_magma_summary.json"
EXPECTED_MARKERS = [
    "PASS Magma: critical-fibre length 15 in all 12 embeddings modulo 863153",
    "PASS Magma: seven base points plus eight residual ramification points",
    "PASS Magma: collision-free residue embeddings [ 2, 3, 4, 5, 9, 12 ]",
    "PASS Magma: good reduction excludes characteristic-zero projected-point collisions",
]


def main():
    payload = json.loads(SUMMARY.read_text())
    raw = CERTIFICATE.read_bytes()
    assert payload["engine"] == "Magma"
    assert payload["magma_version"] == "2.29-9"
    assert payload["runner"] == "University of Sydney Magma Calculator"
    assert payload["status"] == "PASS_INDEPENDENT_MAGMA_GEOMETRY_CERTIFICATE"
    assert payload["markers"] == EXPECTED_MARKERS
    assert payload["prime"] == 863153
    assert payload["residue_embedding_count"] == 12
    assert payload["collision_free_residue_embeddings"] == [2, 3, 4, 5, 9, 12]
    assert payload["input_path"] == "verification/verify_hurwitz_degree23_geometry.m"
    assert payload["input_bytes"] == len(raw)
    assert payload["input_sha256"] == hashlib.sha256(raw).hexdigest()
    assert payload["total_seconds"] > 0
    assert payload["memory_megabytes"] > 0
    print("PASS recorded Hurwitz-geometry Magma run matches the input hash")


if __name__ == "__main__":
    main()
