#!/usr/bin/env python3
"""Verify the compact recorded Hurwitz tail-bound summary."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / "verification" / "hurwitz_tail_summary.json"


def close(left: float, right: float, relative: float = 2e-12) -> bool:
    return abs(left - right) <= relative * max(abs(left), abs(right), 1e-300)


def main() -> None:
    payload = json.loads(SUMMARY.read_text())
    assert payload["status"] == "PASS_RECORDED_HURWITZ_UNIFORM_TAIL_BOUND"
    assert payload["certified_covering_radius_upper"] < 0.471

    for relative, digest in payload["source_sha256"].items():
        actual = hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
        assert actual == digest, f"tail source hash changed: {relative}"

    records = payload["classes"]
    assert [record["class_id"] for record in records] == list(range(1, 8))
    certified_digits = []
    for record in records:
        assert record["model_terms"] == 700
        assert record["samples"] == 1280
        assert len(record["anchor_positions_patch_mode"]) == 4
        assert len(record["all_mode_residual_uppers"]) == 4
        assert len(record["all_mode_residual_balls"]) == 4
        assert len(record["branch_jet_row_l2_error_uppers_0_through_19"]) == 20
        assert len(record["stability_sha256"]) == 64
        assert len(record["model_sha256"]) == 64

        residuals = record["all_mode_residual_uppers"]
        assert max(residuals) == record["maximum_all_mode_residual_upper"]
        assert max(residuals) < 1e-90
        external = record["anchor_external_forcing_upper"]
        assert external < 1e-140
        factor = record["all_residual_to_low_solution_factor_upper"]
        low_errors = [factor * (residual + external) for residual in residuals]
        combined = math.sqrt(math.fsum(value * value for value in low_errors))
        # The generator applies a deliberate 1+1e-12 outward inflation.
        assert combined <= record["four_column_low_coefficient_error_l2_upper"]
        assert close(
            record["four_column_low_coefficient_error_l2_upper"],
            combined * (1 + 1e-12),
        )

        rho = record["rho"]
        outer = record["outer_radius"]
        leading_error = record["four_column_low_coefficient_error_l2_upper"] / rho**3
        assert close(leading_error, record["model_vs_true_leading_matrix_error_upper"])
        sigma = record["true_branch_leading_sigma_minimum_lower"]
        assert leading_error < sigma
        inverse = 1 / (sigma - leading_error)
        assert close(inverse, record["model_branch_leading_inverse_norm_upper"])

        branch_outer = record["branch_normalized_outer_sup_norm_upper"]
        expected_rows = []
        for mode in range(20):
            raw = record["four_column_low_coefficient_error_l2_upper"] / rho**mode
            true_row = 2 * branch_outer / outer**mode
            expected_rows.append(
                (raw + true_row * leading_error) * inverse * (1 + 1e-12)
            )
        for expected, stored in zip(
            expected_rows,
            record["branch_jet_row_l2_error_uppers_0_through_19"],
        ):
            assert close(expected, stored)
        maximum = max(expected_rows)
        assert close(
            maximum,
            record["maximum_first_20_branch_jet_row_l2_error_upper"],
        )
        assert maximum < 1e-70
        digits = math.floor(-math.log10(maximum))
        assert digits == record["certified_decimal_digits_first_20_branch_jet_rows"]
        certified_digits.append(digits)

    assert min(certified_digits) == payload[
        "minimum_certified_decimal_digits_first_20_branch_jet_rows"
    ]
    print("PASS recorded Arb atlas cover and seven finite-Q stability bounds")
    print("PASS all 1280 DFT modes and first 20 branch-jet error bounds")
    print(f"PASS at least {min(certified_digits)} certified branch-jet digits")


if __name__ == "__main__":
    main()
