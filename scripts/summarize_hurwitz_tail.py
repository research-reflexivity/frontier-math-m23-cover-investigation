#!/usr/bin/env python3
"""Combine the recorded geometry, stability, and all-mode residual bounds."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE_FILES = [
    "scripts/acb_fft.pyx",
    "scripts/analyze_hurwitz_tail_stability.py",
    "scripts/certify_hurwitz_acb.py",
    "scripts/compute_hurwitz_acb_model.py",
    "scripts/hurwitz_high_precision.py",
    "scripts/summarize_hurwitz_tail.py",
    "verification/verify_hurwitz_tail_geometry.py",
]
OUTWARD = 1 + 1e-12
INWARD = 1 - 1e-12


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def load(path: Path) -> tuple[dict[str, object], str]:
    data = path.read_bytes()
    return json.loads(data), hashlib.sha256(data).hexdigest()


def anchor_pairs(indices: list[int], mode_count: int) -> list[list[int]]:
    return [[index // mode_count, index % mode_count] for index in indices]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--geometry", type=Path, required=True)
    parser.add_argument("--stability", type=Path, nargs=7, required=True)
    parser.add_argument("--models", type=Path, nargs=7, required=True)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()

    geometry, geometry_sha = load(arguments.geometry)
    require(
        geometry["status"] == "PASS_HURWITZ_UNIFORM_ATLAS_COVER",
        "geometry certificate did not pass",
    )
    require(
        geometry["certified_covering_radius_upper"] < 0.471,
        "atlas covering radius is too large",
    )

    stabilities = {}
    stability_hashes = {}
    for path in arguments.stability:
        payload, digest = load(path)
        class_id = int(payload["class_id"])
        require(class_id not in stabilities, "duplicate stability class")
        require(
            payload["status"]
            == "PASS_HURWITZ_TAIL_LOW_MODE_STABILITY_CERTIFICATE",
            f"class {class_id} stability certificate did not pass",
        )
        stabilities[class_id] = payload
        stability_hashes[class_id] = digest

    models = {}
    model_hashes = {}
    for path in arguments.models:
        payload, digest = load(path)
        class_id = int(payload["class_id"])
        require(class_id not in models, "duplicate model class")
        require(
            payload["status"] == "PASS_MIXED_PRECISION_NUMERICAL_CANONICAL_QUADRIC",
            f"class {class_id} model did not pass",
        )
        models[class_id] = payload
        model_hashes[class_id] = digest

    require(sorted(stabilities) == list(range(1, 8)), "missing stability class")
    require(sorted(models) == list(range(1, 8)), "missing model class")

    records = []
    for class_id in range(1, 8):
        stability = stabilities[class_id]
        model = models[class_id]
        low = stability["left_inverse_certificate"]
        sup_norm = stability["normalized_sup_norm_certificate"]
        require(stability["low_terms"] == 60, "wrong low cutoff")
        require(stability["samples"] == 1280, "wrong stability sample count")
        require(model["samples"] == 1280, "wrong model sample count")
        require(model["terms"] < model["samples"], "model cutoff must be below Q")
        require(model["rho"] == stability["rho"], "rho changed between runs")
        require(
            anchor_pairs(stability["anchor_indices"], 61)
            == anchor_pairs(model["anchor_indices"], model["terms"] + 1),
            "semantic anchor positions changed between runs",
        )

        residuals = [
            (
                float(item["residual_norm_midpoint"])
                + float(item["residual_norm_radius"])
            )
            * OUTWARD
            for item in model["all_dft_mode_residuals"]
        ]
        require(
            all(item["output_mode_count"] == 1280
                for item in model["all_dft_mode_residuals"]),
            "not all DFT modes were evaluated",
        )
        external = OUTWARD * float(
            sup_norm["q1280_anchor_normalized_external_forcing_l2_upper"]
        )
        low_factor = OUTWARD * float(
            low["certified_all_residual_to_low_solution_factor_upper"]
        )
        low_errors = [low_factor * (residual + external) for residual in residuals]
        column_error_l2 = math.sqrt(math.fsum(value * value for value in low_errors))
        column_error_l2 *= OUTWARD

        rho = float(stability["rho"])
        outer = float(stability["outer_radius"])
        leading_error = column_error_l2 / rho**3
        true_leading_sigma = INWARD * float(
            sup_norm["certified_branch_leading_sigma_minimum_lower"]
        )
        require(leading_error < true_leading_sigma, "model leading jet may be singular")
        model_leading_inverse = 1 / (true_leading_sigma - leading_error)
        branch_outer = OUTWARD * float(
            sup_norm["branch_normalized_outer_sup_norm_upper"]
        )
        jet_errors = []
        for mode in range(20):
            raw_row_error = column_error_l2 / rho**mode
            true_branch_row_norm = 2 * branch_outer / outer**mode
            normalized_error = (
                raw_row_error + true_branch_row_norm * leading_error
            ) * model_leading_inverse
            jet_errors.append(normalized_error * OUTWARD)
        maximum_jet_error = max(jet_errors)
        require(maximum_jet_error < 1e-70, "fewer than 70 certified jet digits")

        records.append(
            {
                "class_id": class_id,
                "model_terms": int(model["terms"]),
                "samples": int(model["samples"]),
                "rho": rho,
                "outer_radius": outer,
                "anchor_positions_patch_mode": anchor_pairs(
                    model["anchor_indices"], model["terms"] + 1
                ),
                "maximum_all_mode_residual_upper": max(residuals),
                "all_mode_residual_uppers": residuals,
                "all_mode_residual_balls": [
                    item["residual_norm_ball"]
                    for item in model["all_dft_mode_residuals"]
                ],
                "anchor_external_forcing_upper": external,
                "all_residual_to_low_solution_factor_upper": low_factor,
                "four_column_low_coefficient_error_l2_upper": column_error_l2,
                "true_branch_leading_sigma_minimum_lower": true_leading_sigma,
                "model_vs_true_leading_matrix_error_upper": leading_error,
                "model_branch_leading_inverse_norm_upper": model_leading_inverse,
                "branch_normalized_outer_sup_norm_upper": branch_outer,
                "branch_jet_row_l2_error_uppers_0_through_19": jet_errors,
                "maximum_first_20_branch_jet_row_l2_error_upper": maximum_jet_error,
                "certified_decimal_digits_first_20_branch_jet_rows": int(
                    math.floor(-math.log10(maximum_jet_error))
                ),
                "stability_sha256": stability_hashes[class_id],
                "model_sha256": model_hashes[class_id],
            }
        )

    result = {
        "status": "PASS_RECORDED_HURWITZ_UNIFORM_TAIL_BOUND",
        "scope": (
            "Arb atlas cover, a posteriori low-block inverse, finite-Q Schur "
            "bound, outer-sup bootstrap, all-mode Acb residuals, and certified "
            "errors for the first 20 branch-normalized Taylor rows"
        ),
        "geometry_sha256": geometry_sha,
        "source_sha256": {
            relative: hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
            for relative in SOURCE_FILES
        },
        "certified_covering_radius_upper": geometry[
            "certified_covering_radius_upper"
        ],
        "classes": records,
        "minimum_certified_decimal_digits_first_20_branch_jet_rows": min(
            record["certified_decimal_digits_first_20_branch_jet_rows"]
            for record in records
        ),
    }
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    print(rendered, end="")
    if arguments.output:
        arguments.output.write_text(rendered)


if __name__ == "__main__":
    main()
