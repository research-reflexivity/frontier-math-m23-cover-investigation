#!/usr/bin/env python3
"""Check that the recorded Magma run applies to the current certificate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SUMMARY_PATH = ROOT / "verification" / "magma_verification_summary.json"
GAUSS_SUMMARY_PATH = (
    ROOT / "verification" / "gauss_prolongation_magma_summary.json"
)
QUADRIC_SUMMARY_PATH = (
    ROOT / "verification" / "canonical_quadric_magma_summary.json"
)
EXPECTED_OUTPUT = [
    "PASS_MAGMA_EXACT_FUNCTION_FIELD_IDENTITY",
    "PASS_MAGMA_IRREDUCIBLE_MOD_31",
    "PASS_OPTIMAL_23_4_MAGMA_CERTIFICATE",
]


def main() -> None:
    summary = json.loads(SUMMARY_PATH.read_text())

    assert summary["status"] == "PASS_INDEPENDENT_MAGMA_CERTIFICATE"
    assert summary["primitive_over_Z"] is True
    assert summary["degree_T"] == 4
    assert summary["degree_W"] == 23
    assert summary["exact_identity_in_QT_mod_Fhat"] is True
    assert summary["irreducible_mod_31"] is True
    assert summary["output"] == EXPECTED_OUTPUT

    input_path = (ROOT / summary["input_file"]).resolve()
    input_path.relative_to(ROOT)
    certificate = input_path.read_bytes()
    assert len(certificate) == summary["input_bytes"]
    assert hashlib.sha256(certificate).hexdigest() == summary["input_sha256"]

    source = certificate.decode("utf-8")
    for marker in EXPECTED_OUTPUT:
        assert f'print "{marker}";' in source

    gauss = json.loads(GAUSS_SUMMARY_PATH.read_text())
    gauss_marker = "PASS_GAUSS_PROLONGATION_OBSTRUCTION_MAGMA_CERTIFICATE"
    assert gauss["status"] == "PASS_INDEPENDENT_MAGMA_CERTIFICATE"
    assert gauss["magma_version"] == "2.29-9"
    assert gauss["fixed_ids"] == [3, 6, 7]
    assert gauss["decomposition_group_order"] == 253
    assert gauss["gauss_prolongation_count"] == 40320
    assert gauss["anchor_gauge_orders"] == [32, 32, 32]
    assert gauss["reflection_pair_count"] == 20160
    assert gauss["pair_gauge_orbit_counts"] == [1260, 1260, 1260]
    assert gauss["xr_anchors_form_one_simultaneous_conjugacy_orbit"] is True
    assert gauss["output"] == [gauss_marker]

    gauss_input_path = (ROOT / gauss["input_file"]).resolve()
    gauss_input_path.relative_to(ROOT)
    gauss_certificate = gauss_input_path.read_bytes()
    assert len(gauss_certificate) == gauss["input_bytes"]
    assert hashlib.sha256(gauss_certificate).hexdigest() == gauss["input_sha256"]
    assert f'print "{gauss_marker}";' in gauss_certificate.decode("utf-8")

    quadric = json.loads(QUADRIC_SUMMARY_PATH.read_text())
    quadric_markers = [
        "PASS_MAGMA_CANONICAL_QUADRIC_IDENTITY_AND_UNIQUENESS",
        "PASS_MAGMA_RULING_FIELD_Q_SQRT_4873_DIFFERS_FROM_Q_SQRT_MINUS_23",
        "PASS_CANONICAL_QUADRIC_MAGMA_CERTIFICATE",
    ]
    assert quadric["status"] == "PASS_INDEPENDENT_MAGMA_CERTIFICATE"
    assert quadric["magma_version"] == "2.29-9"
    assert quadric["canonical_quadric_unique"] is True
    assert quadric["determinant_squarefree_part"] == 4873
    assert quadric["ruling_field"] == "Q(sqrt(4873)) = Q(sqrt(11*443))"
    assert quadric["branch_orientation_field"] == "Q(sqrt(-23))"
    assert quadric["fields_equal"] is False
    assert quadric["output"] == quadric_markers

    quadric_input_path = (ROOT / quadric["input_file"]).resolve()
    quadric_input_path.relative_to(ROOT)
    quadric_certificate = quadric_input_path.read_bytes()
    assert len(quadric_certificate) == quadric["input_bytes"]
    assert hashlib.sha256(quadric_certificate).hexdigest() == quadric["input_sha256"]
    quadric_source = quadric_certificate.decode("utf-8")
    for marker in quadric_markers:
        assert f'print "{marker}";' in quadric_source

    print("PASS recorded Magma runs match the current certificate hashes")


if __name__ == "__main__":
    main()
