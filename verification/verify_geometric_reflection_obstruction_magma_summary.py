#!/usr/bin/env python3
"""Check the recorded Magma run for the geometric-reflection obstruction."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
SUMMARY = HERE / "geometric_reflection_obstruction_magma_summary.json"
SOURCE = HERE / "certify_geometric_reflection_obstruction.m"
MARKER = "PASS_GEOMETRIC_REFLECTION_OBSTRUCTION_MAGMA"


def main() -> None:
    payload = json.loads(SUMMARY.read_text())
    source = SOURCE.read_bytes()
    assert payload["engine"] == "Magma"
    assert payload["version"] == "V2.29-9"
    assert payload["status"] == "pass"
    assert payload["input_path"] == "verification/certify_geometric_reflection_obstruction.m"
    assert payload["input_bytes"] == len(source)
    assert payload["input_sha256"] == hashlib.sha256(source).hexdigest()
    assert payload["output"][-1] == MARKER
    assert f'print "{MARKER}";' in source.decode()
    print("PASS recorded geometric-reflection Magma run matches the input hash")


if __name__ == "__main__":
    main()
