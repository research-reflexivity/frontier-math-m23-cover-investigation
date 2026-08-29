#!/usr/bin/env python3
"""Bind the recorded Sydney Magma connector runs to their exact inputs."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CASES = {
    "fano_affine_odd_fixed_point_magma_summary.json": {
        "input": "notes/certify_fano_affine_odd_fixed_point_lemma.m",
        "markers": ["PASS_FANO_AFFINE_ODD_FIXED_POINT_LEMMA_MAGMA"],
    },
    "pinched_tag_finite_identities_magma_summary.json": {
        "input": "notes/certify_pinched_tag_finite_identities.m",
        "markers": ["PASS_PINCHED_TAG_FINITE_IDENTITIES"],
    },
    "wild_parameter_orientation_magma_summary.json": {
        "input": "notes/certify_wild_parameter_orientation.m",
        "markers": ["PASS_WILD_PARAMETER_ORIENTATION"],
    },
    "pointed_relative_bockstein_magma_summary.json": {
        "input": "notes/audit_pointed_relative_bockstein.m",
        "markers": ["PASS_POINTED_RELATIVE_BOCKSTEIN_AUDIT"],
    },
    "log_quadratic_orientation_line_magma_summary.json": {
        "input": "notes/audit_log_quadratic_orientation_line.m",
        "markers": [
            "log_quadratic_orientation_line_magma_audit=PASS",
            "normalization_telescope_audit=PASS",
        ],
    },
}


def main() -> None:
    for summary_name, expected in CASES.items():
        summary_path = ROOT / "verification" / summary_name
        payload = json.loads(summary_path.read_text())
        input_path = ROOT / expected["input"]
        raw = input_path.read_bytes()

        assert payload["engine"] == "Magma"
        assert payload["magma_version"] == "2.29-9"
        assert payload["runner"] == "University of Sydney Magma Calculator"
        assert payload["status"] == "PASS_INDEPENDENT_MAGMA_CERTIFICATE"
        assert payload["input_path"] == expected["input"]
        assert payload["input_bytes"] == len(raw)
        assert payload["input_sha256"] == hashlib.sha256(raw).hexdigest()
        assert payload["markers"] == expected["markers"]
        assert payload["total_seconds"] > 0
        assert payload["memory_megabytes"] > 0
        assert payload["seed"] > 0

        output = payload["output"]
        assert output
        assert not any("error" in line.lower() for line in output)
        for marker in expected["markers"]:
            assert marker in output
            assert f'print "{marker}";' in raw.decode("utf-8")

    print("PASS recorded connector Magma runs match the current input hashes")


if __name__ == "__main__":
    main()
