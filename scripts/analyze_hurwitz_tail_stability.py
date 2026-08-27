#!/usr/bin/env sage-python
"""Low-mode stability and outer-sup certificate for the Hurwitz tail proof.

This script builds the dense low-mode block of the sampled multi-centre
Hejhal operator, augments it by four normalization rows, and measures the
singular gap.  It also bounds the coupling from all omitted input modes by
using the routewise target radii.

Without ``--certify-left-inverse`` the output is only a LAPACK diagnostic.
With that flag, the LAPACK inverse is treated merely as a proposed matrix and
is checked a posteriori using explicit IEEE-754 error bounds and Acb-enclosed
route geometry.  The proof and norm conventions are recorded in
``HURWITZ_TAIL_BOUND.md``.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import qr, svdvals
from sage.all import ComplexBallField, I

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import hurwitz_high_precision as hp  # noqa: E402
import certify_hurwitz_acb as acb  # noqa: E402

OUTWARD_FLOAT_INFLATION = 1 + 1e-12


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def dense_low_mode_matrix(
    routes: hp.RouteTable, backend: str = "direct"
) -> np.ndarray:
    """Return ``T-I`` for all modes through ``routes.terms``."""

    patch_count = len(routes.atlas.centers)
    mode_count = routes.terms + 1
    dimension = patch_count * mode_count
    matrix = np.zeros((dimension, dimension), dtype=np.complex128)
    if backend == "fft":
        sample_rows = np.arange(routes.samples)
        for source in range(patch_count):
            values = np.zeros(
                (routes.samples, dimension), dtype=np.complex128
            )
            target = routes.targets[source]
            base = routes.target_base[source]
            factor = routes.factor[source]
            powers = np.ones(routes.samples, dtype=np.complex128)
            for mode in range(mode_count):
                values[sample_rows, target * mode_count + mode] = factor * powers
                powers *= base
            matrix[
                source * mode_count : (source + 1) * mode_count
            ] = np.fft.fft(values, axis=0)[:mode_count] / routes.samples
    elif backend == "direct":
        samples = np.arange(routes.samples, dtype=np.float64)
        modes = np.arange(mode_count, dtype=np.float64)
        fourier = np.exp(
            -2j * np.pi * modes[:, None] * samples[None, :] / routes.samples
        ) / routes.samples
        for source in range(patch_count):
            for target in range(patch_count):
                selected = np.flatnonzero(routes.targets[source] == target)
                if not len(selected):
                    continue
                values = np.empty(
                    (len(selected), mode_count), dtype=np.complex128
                )
                values[:, 0] = routes.factor[source, selected]
                for mode in range(1, mode_count):
                    values[:, mode] = (
                        values[:, mode - 1]
                        * routes.target_base[source, selected]
                    )
                matrix[
                    source * mode_count : (source + 1) * mode_count,
                    target * mode_count : (target + 1) * mode_count,
                ] = fourier[:, selected] @ values
    else:
        raise ValueError(f"unknown dense matrix backend: {backend}")
    matrix -= np.eye(dimension, dtype=np.complex128)
    return matrix


def gamma_bound(operation_count: int) -> float:
    unit_roundoff = 2.0**-53
    product = operation_count * unit_roundoff
    require(product < 1, "floating-point gamma bound overflowed")
    return product / (1 - product)


def induced_two_norm_upper(matrix: np.ndarray) -> float:
    absolute = np.abs(matrix)
    row_sum = float(np.max(np.sum(absolute, axis=1)))
    column_sum = float(np.max(np.sum(absolute, axis=0)))
    inflation = 1 + gamma_bound(max(matrix.shape))
    return math.sqrt(row_sum * column_sum) * inflation


def maximum_fourier_weight_error(
    mode_count: int, sample_count: int, precision: int
) -> float:
    field = ComplexBallField(precision)
    pi = field.pi()
    maximum = 0.0
    for mode in range(mode_count):
        for sample in range(sample_count):
            angle = -2 * pi * mode * sample / sample_count
            exact = (
                angle.cos() + field(I) * angle.sin()
            ) / sample_count
            approximate = np.exp(
                -2j * np.pi * mode * sample / sample_count
            ) / sample_count
            error = (
                exact - acb.decimal_ball(field, approximate)
            ).abs().upper()
            maximum = max(maximum, float(error))
    return maximum


def route_geometry_operator_error(
    routes: hp.RouteTable,
    refined_routes: hp.RouteTable,
    refinement: dict[str, float],
) -> dict[str, float]:
    mode_count = routes.terms + 1
    base_error = refinement[
        "maximum_target_base_error_upper_vs_complex128"
    ]
    factor_error = refinement[
        "maximum_factor_error_upper_vs_complex128"
    ]
    base_bound = max(
        float(np.abs(routes.target_base).max()),
        float(np.abs(refined_routes.target_base).max()),
    ) + base_error
    factor_bound = float(np.abs(routes.factor).max()) + factor_error
    base_bound *= OUTWARD_FLOAT_INFLATION
    factor_bound *= OUTWARD_FLOAT_INFLATION
    coefficient_norm = math.sqrt(
        math.fsum(base_bound ** (2 * mode) for mode in range(mode_count))
    )
    derivative_norm = math.sqrt(
        math.fsum(
            mode * mode * base_bound ** (2 * (mode - 1))
            for mode in range(1, mode_count)
        )
    )
    row_functional_error = (
        factor_error * coefficient_norm
        + factor_bound * base_error * derivative_norm
    )
    target_counts = np.bincount(
        routes.targets.ravel(), minlength=len(routes.atlas.centers)
    )
    replication_factor = math.sqrt(
        float(target_counts.max()) / routes.samples
    )
    return {
        "base_error_upper": base_error,
        "factor_error_upper": factor_error,
        "base_bound": base_bound,
        "factor_bound": factor_bound,
        "row_functional_error": row_functional_error,
        "maximum_target_occurrence_count": int(target_counts.max()),
        "replication_factor": replication_factor,
        "operator_error_bound": (
            replication_factor
            * row_functional_error
            * OUTWARD_FLOAT_INFLATION
        ),
    }


def route_column_norms_with_uniform_errors(
    routes: hp.RouteTable,
    refinement: dict[str, float],
) -> tuple[float, float]:
    """Outward bound for the exact low/high route-column norms."""

    base_error = refinement[
        "maximum_target_base_error_upper_vs_complex128"
    ]
    factor_error = refinement[
        "maximum_factor_error_upper_vs_complex128"
    ]
    radius = (
        np.abs(routes.target_base) + base_error
    ) * OUTWARD_FLOAT_INFLATION
    factor = (
        np.abs(routes.factor) + factor_error
    ) * OUTWARD_FLOAT_INFLATION
    require(float(radius.max()) < 1, "route base error reaches the unit circle")
    denominator = 1 - radius * radius
    high_rows = (
        factor
        * factor
        * radius ** (2 * (routes.terms + 1))
        / denominator
    )
    low_rows = (
        factor
        * factor
        * (1 - radius ** (2 * (routes.terms + 1)))
        / denominator
    )
    patch_count = len(routes.atlas.centers)
    high_by_target = np.zeros(patch_count, dtype=np.float64)
    low_by_target = np.zeros(patch_count, dtype=np.float64)
    for target in range(patch_count):
        mask = routes.targets == target
        # Every summand has already been rounded outward through the uniform
        # Acb error.  Inflate the positive floating sum by gamma_n as well.
        count = int(np.count_nonzero(mask))
        inflation = 1 + gamma_bound(max(1, count))
        # With the 1/Q-normalized DFT, Parseval contributes exactly one
        # factor 1/Q to the squared l2 norm.  The samples already occur in
        # the positive sum, so a second division by Q would be incorrect.
        high_by_target[target] = (
            math.fsum(high_rows[mask].ravel())
            * inflation
            / routes.samples
        )
        low_by_target[target] = (
            math.fsum(low_rows[mask].ravel())
            * inflation
            / routes.samples
        )
    return float(np.sqrt(low_by_target.max())), float(
        np.sqrt(high_by_target.max())
    )


def floating_left_inverse_certificate(
    augmented: np.ndarray,
    routes: hp.RouteTable,
    precision: int,
) -> dict[str, object]:
    from scipy.linalg import svd

    # The SVD is used only to propose a double-precision left inverse.  The
    # certificate below verifies that particular stored matrix a posteriori;
    # it does not trust the singular values returned by LAPACK.
    left_vectors, singular_values, right_vectors = svd(
        augmented,
        full_matrices=False,
        overwrite_a=False,
        check_finite=False,
        lapack_driver="gesdd",
    )
    left_inverse = (
        right_vectors.conjugate().T / singular_values[None, :]
    ) @ left_vectors.conjugate().T
    stored_residual = (
        left_inverse @ augmented
        - np.eye(augmented.shape[1], dtype=np.complex128)
    )
    left_inverse_norm = induced_two_norm_upper(left_inverse)
    stored_residual_norm = induced_two_norm_upper(stored_residual)

    absolute_product = np.abs(left_inverse) @ np.abs(augmented)
    multiplication_roundoff = (
        gamma_bound(8 * augmented.shape[0])
        * induced_two_norm_upper(absolute_product)
    )

    acb.verify_route_membership(routes)
    refined_routes, refinement = acb.reevaluate_route_midpoints(
        precision, routes
    )
    geometry = route_geometry_operator_error(
        routes, refined_routes, refinement
    )
    fourier_weight_error = maximum_fourier_weight_error(
        routes.terms + 1, routes.samples, precision
    )
    dimension = augmented.shape[1]
    route_factor_bound = geometry["factor_bound"]
    # A deliberately coarse Frobenius bound for construction of the direct
    # DFT matrix.  The factor 64 covers complex products, the ordered dot
    # product, and the iterative powers in each block.
    direct_arithmetic_error = (
        dimension
        * gamma_bound(64 * routes.samples)
        * route_factor_bound
    )
    fourier_weight_operator_error = (
        dimension
        * routes.samples
        * fourier_weight_error
        * route_factor_bound
    )
    matrix_error = (
        geometry["operator_error_bound"]
        + direct_arithmetic_error
        + fourier_weight_operator_error
    )
    eta = (
        stored_residual_norm
        + multiplication_roundoff
        + left_inverse_norm * matrix_error
    )
    require(eta < 1, "a posteriori left-inverse residual is not a contraction")
    certified_sigma_lower = (1 - eta) / left_inverse_norm
    certified_low_norm, certified_high_norm = (
        route_column_norms_with_uniform_errors(routes, refinement)
    )
    require(certified_high_norm < 1, "certified high block is not contractive")
    certified_schur_perturbation = (
        certified_high_norm
        * certified_low_norm
        / (1 - certified_high_norm)
    )
    certified_schur_margin = (
        certified_sigma_lower - certified_schur_perturbation
    )
    require(certified_schur_margin > 0, "certified Schur margin is not positive")
    high_to_low_residual_factor = (
        certified_high_norm / (1 - certified_high_norm)
    )
    certified_low_solution_factor = math.hypot(
        1, high_to_low_residual_factor
    ) / certified_schur_margin
    certified_high_solution_factor = (
        1 + certified_low_norm * certified_low_solution_factor
    ) / (1 - certified_high_norm)
    certified_full_inverse_bound = math.hypot(
        certified_low_solution_factor, certified_high_solution_factor
    )
    return {
        "status": "PASS_A_POSTERIORI_LEFT_INVERSE_IEEE754_ACB_GEOMETRY",
        "scope": (
            "Acb route geometry plus an a posteriori IEEE-754 left-inverse "
            "bound for the direct-DFT low matrix"
        ),
        "precision_bits": precision,
        "left_inverse_norm_upper": left_inverse_norm,
        "stored_left_inverse_residual_norm_upper": stored_residual_norm,
        "left_inverse_product_roundoff_bound": multiplication_roundoff,
        "route_refinement": refinement,
        "route_geometry_operator_error": geometry,
        "maximum_scaled_fourier_weight_error_upper": fourier_weight_error,
        "direct_matrix_arithmetic_error_bound": direct_arithmetic_error,
        "fourier_weight_operator_error_bound": fourier_weight_operator_error,
        "total_exact_matrix_vs_stored_matrix_error_bound": matrix_error,
        "left_inverse_contraction_eta": eta,
        "certified_augmented_sigma_minimum_lower": certified_sigma_lower,
        "certified_low_input_operator_norm_upper": certified_low_norm,
        "certified_high_input_operator_norm_upper": certified_high_norm,
        "certified_schur_perturbation_upper": certified_schur_perturbation,
        "certified_schur_margin_lower": certified_schur_margin,
        "certified_high_residual_to_low_solution_factor_upper": (
            high_to_low_residual_factor / certified_schur_margin
        ),
        "certified_all_residual_to_low_solution_factor_upper": (
            certified_low_solution_factor
        ),
        "certified_all_residual_to_high_solution_factor_upper": (
            certified_high_solution_factor
        ),
        "certified_full_inverse_norm_upper": certified_full_inverse_bound,
    }


def ball_matrix_sigma_lower(entries: list[list]) -> dict[str, float]:
    """Certify a singular-value lower bound for a small square ball matrix."""

    dimension = len(entries)
    require(dimension > 0, "empty ball matrix")
    require(
        all(len(row) == dimension for row in entries),
        "ball matrix is not square",
    )
    field = entries[0][0].parent()
    midpoint = np.asarray(
        [
            [
                complex(
                    float(value.real().center()),
                    float(value.imag().center()),
                )
                for value in row
            ]
            for row in entries
        ],
        dtype=np.complex128,
    )
    entry_errors = [
        float(
            (
                value
                - acb.decimal_ball(field, midpoint[row, column])
            ).abs().upper()
        )
        for row, values in enumerate(entries)
        for column, value in enumerate(values)
    ]
    matrix_error = math.sqrt(math.fsum(error * error for error in entry_errors))
    inverse = np.linalg.inv(midpoint)
    inverse_norm = induced_two_norm_upper(inverse)
    residual = inverse @ midpoint - np.eye(dimension, dtype=np.complex128)
    residual_norm = induced_two_norm_upper(residual)
    absolute_product = np.abs(inverse) @ np.abs(midpoint)
    multiplication_roundoff = (
        gamma_bound(8 * dimension)
        * induced_two_norm_upper(absolute_product)
    )
    eta = (
        residual_norm
        + multiplication_roundoff
        + inverse_norm * matrix_error
    )
    require(eta < 1, "small ball matrix inverse is not certified")
    return {
        "stored_inverse_norm_upper": inverse_norm,
        "stored_inverse_residual_norm_upper": residual_norm,
        "inverse_product_roundoff_bound": multiplication_roundoff,
        "ball_matrix_vs_stored_matrix_error_bound": matrix_error,
        "inverse_contraction_eta": eta,
        "sigma_minimum_lower": (1 - eta) / inverse_norm,
    }


def normalized_outer_sup_norm_certificate(
    basis: np.ndarray,
    anchors: np.ndarray,
    routes: hp.RouteTable,
    left_inverse: dict[str, object],
    covering_radius: float,
    outer_radius: float,
    full_terms: int,
    precision: int,
) -> dict[str, object]:
    """Close the sup-norm bootstrap and certify branch-jet normalization."""

    require(
        covering_radius <= 0.471,
        "this certificate expects the separately certified radius 0.471",
    )
    patch_count = len(routes.atlas.centers)
    mode_count = routes.terms + 1
    sigma_lower = float(
        left_inverse["certified_augmented_sigma_minimum_lower"]
    )
    high_norm = float(
        left_inverse["certified_high_input_operator_norm_upper"]
    )

    rho = routes.rho
    coefficient_ratio = rho / outer_radius
    evaluation_ratio = covering_radius / rho
    outer_target_ratio = covering_radius / outer_radius
    omitted_coefficient_norm = (
        math.sqrt(patch_count)
        * coefficient_ratio**mode_count
        / math.sqrt(1 - coefficient_ratio**2)
    )
    source_alias_norm = (
        math.sqrt(patch_count * mode_count)
        * coefficient_ratio**routes.samples
        / (1 - coefficient_ratio**routes.samples)
    )
    low_evaluation_norm = math.sqrt(
        math.fsum(
            evaluation_ratio ** (2 * mode) for mode in range(mode_count)
        )
    )
    inner_tail = (
        outer_target_ratio**mode_count / (1 - outer_target_ratio)
    )

    # The exact Acb value is below 79628 for the present atlas.  Use the
    # round integer 80000 in the inequalities, and verify that assertion.
    ball = acb.build_ball_geometry(precision)
    centers = acb.atlas_centers(ball, routes.atlas)
    y_maximum = max(center.imag().upper() for center in centers)
    y_minimum = min(center.imag().lower() for center in centers)
    real_field = y_maximum.parent()
    outer_ball = real_field(str(outer_radius))
    exact_transition_bound = (
        y_maximum / y_minimum / (1 - outer_ball**2)
    )
    transition_bound = 80000.0
    require(
        exact_transition_bound < transition_bound,
        "round outer transition bound is too small",
    )

    feedback = transition_bound * (
        low_evaluation_norm
        / sigma_lower
        * (high_norm * omitted_coefficient_norm + source_alias_norm)
        + inner_tail
    )
    require(feedback < 1, "outer sup-norm bootstrap is not contractive")
    anchor_outer_sup_norm = (
        transition_bound
        * low_evaluation_norm
        / sigma_lower
        / (1 - feedback)
    )
    anchor_outer_sup_norm *= 1 + 1e-12

    anchor_basis = basis @ np.linalg.inv(basis[anchors])
    refined_routes, _ = acb.reevaluate_route_midpoints(precision, routes)
    evaluator = acb.CachedAcbHejhalEvaluator(precision, refined_routes)
    residual_uppers = []
    anchor_coefficients = []
    for column in range(4):
        coefficients = acb.high_precision_coefficients(
            anchor_basis[:, column],
            [int(value) for value in anchors],
            column,
            precision + 64,
        )
        anchor_coefficients.append(coefficients)
        _, metrics = evaluator.evaluate(coefficients)
        residual_uppers.append(
            metrics["residual_norm_midpoint"]
            + metrics["residual_norm_radius"]
        )

    low_coefficient_errors = [
        (
            residual
            + (
                high_norm * omitted_coefficient_norm
                + source_alias_norm
            )
            * anchor_outer_sup_norm
        )
        / sigma_lower
        for residual in residual_uppers
    ]
    maximum_low_coefficient_error = max(low_coefficient_errors)

    polynomial = branch_normalized_polynomial_bound(
        basis, routes, evaluation_ratio, anchors
    )
    computed_leading_sigma_diagnostic = float(
        polynomial["anchor_normalized_leading_sigma_minimum"]
    )
    branch_patch = patch_count - 1
    branch_start = branch_patch * mode_count
    q_rotation = (
        evaluator.field(I)
        * (evaluator.ball.pi - evaluator.ball.pi / 23)
    ).exp()
    leading_balls = [
        [
            acb.decimal_ball(
                evaluator.field,
                anchor_coefficients[column][branch_start + row],
            )
            / (evaluator.rho * q_rotation) ** row
            for column in range(4)
        ]
        for row in range(4)
    ]
    computed_leading_certificate = ball_matrix_sigma_lower(leading_balls)
    computed_leading_sigma_lower = float(
        computed_leading_certificate["sigma_minimum_lower"]
    )
    # Four columns and the largest q-jet rescaling rho^{-3} give this
    # Frobenius-norm perturbation bound for the 4-by-4 leading matrix.
    leading_matrix_error = (
        2 * rho**-3 * maximum_low_coefficient_error
    )
    certified_leading_sigma = (
        computed_leading_sigma_lower - leading_matrix_error
    )
    require(
        certified_leading_sigma > 0,
        "branch leading matrix is not certified invertible",
    )
    certified_leading_inverse_norm = 1 / certified_leading_sigma
    branch_outer_sup_norm = (
        2 * certified_leading_inverse_norm * anchor_outer_sup_norm
    )
    branch_outer_sup_norm *= 1 + 1e-12

    actual_target_radius_upper = float(
        left_inverse["route_refinement"]["maximum_target_radius_upper"]
    )
    actual_target_ratio = actual_target_radius_upper / outer_radius
    n480_pointwise_target_tail_per_unit = (
        actual_target_ratio ** (full_terms + 1)
        / (1 - actual_target_ratio)
    )
    route_factor_bound = float(
        left_inverse["route_geometry_operator_error"]["factor_bound"]
    )
    n480_sampled_target_l2_per_unit = (
        math.sqrt(patch_count)
        * route_factor_bound
        * n480_pointwise_target_tail_per_unit
    )

    # In the square Q-mode collocation system, coefficients with indices
    # 0,...,Q-1 are unknowns.  The genuine external target tail therefore
    # starts at Q, while the source side contributes the usual DFT aliases
    # c_{k+lQ}.  Parseval supplies sqrt(patch_count) for the target routes;
    # summing the geometric aliases in l2 supplies 1/sqrt(1-r^2).
    q_target_log_upper = (
        0.5 * math.log(patch_count)
        + math.log(route_factor_bound)
        + routes.samples * math.log(actual_target_ratio)
        - math.log1p(-actual_target_ratio)
        + 1e-12
    )
    # The present value is below the smallest positive normal binary64.
    # Retain a deliberately larger positive upper bound instead of emitting
    # zero, which would be mathematically false even though this term is
    # negligible beside the source alias.
    q_target_tail_l2_per_unit = (
        sys.float_info.min
        if q_target_log_upper < math.log(sys.float_info.min)
        else math.exp(q_target_log_upper)
    )
    source_ratio = rho / outer_radius
    q_source_alias_l2_per_unit = (
        math.sqrt(patch_count)
        * source_ratio**routes.samples
        / (
            (1 - source_ratio**routes.samples)
            * math.sqrt(1 - source_ratio**2)
        )
    )
    return {
        "status": "PASS_HURWITZ_NORMALIZED_OUTER_SUP_NORM_CERTIFICATE",
        "scope": (
            "normalized outer sup norm, branch-jet invertibility, N=480 "
            "target tail, and Q=1280 source alias"
        ),
        "certified_covering_radius_used": covering_radius,
        "outer_radius": outer_radius,
        "exact_outer_transition_factor_ball": str(exact_transition_bound),
        "outer_transition_factor_upper_used": transition_bound,
        "omitted_coefficient_l2_per_unit_outer_sup_norm": (
            omitted_coefficient_norm
        ),
        "source_alias_l2_per_unit_outer_sup_norm": source_alias_norm,
        "low_evaluation_l2_norm": low_evaluation_norm,
        "inner_tail_per_unit_outer_sup_norm": inner_tail,
        "sup_norm_bootstrap_feedback": feedback,
        "anchor_normalized_outer_sup_norm_upper": anchor_outer_sup_norm,
        "anchor_basis_finite_residual_norm_uppers": residual_uppers,
        "anchor_basis_low_coefficient_error_uppers": low_coefficient_errors,
        "computed_branch_leading_sigma_minimum_diagnostic": (
            computed_leading_sigma_diagnostic
        ),
        "computed_branch_leading_matrix_certificate": (
            computed_leading_certificate
        ),
        "branch_leading_matrix_error_upper": leading_matrix_error,
        "certified_branch_leading_sigma_minimum_lower": (
            certified_leading_sigma
        ),
        "certified_branch_leading_inverse_norm_upper": (
            certified_leading_inverse_norm
        ),
        "branch_normalized_outer_sup_norm_upper": branch_outer_sup_norm,
        "actual_target_radius_upper": actual_target_radius_upper,
        "full_terms": full_terms,
        "n480_pointwise_target_tail_per_unit_outer_sup_norm": (
            n480_pointwise_target_tail_per_unit
        ),
        "n480_sampled_target_l2_per_unit_outer_sup_norm": (
            n480_sampled_target_l2_per_unit
        ),
        "n480_branch_normalized_pointwise_target_tail_upper": (
            n480_pointwise_target_tail_per_unit * branch_outer_sup_norm
        ),
        "n480_branch_normalized_sampled_target_l2_upper": (
            n480_sampled_target_l2_per_unit * branch_outer_sup_norm
        ),
        "q1280_target_tail_l2_per_unit_outer_sup_norm": (
            q_target_tail_l2_per_unit
        ),
        "q1280_source_alias_l2_per_unit_outer_sup_norm": (
            q_source_alias_l2_per_unit
        ),
        "q1280_anchor_normalized_external_forcing_l2_upper": (
            (q_target_tail_l2_per_unit + q_source_alias_l2_per_unit)
            * anchor_outer_sup_norm
        ),
        "q1280_branch_normalized_external_forcing_l2_upper": (
            (q_target_tail_l2_per_unit + q_source_alias_l2_per_unit)
            * branch_outer_sup_norm
        ),
        "acb_evaluation_backend": evaluator.evaluation_backend,
        "acb_fourier_backend": evaluator.fourier_backend,
    }


def route_column_norms(
    routes: hp.RouteTable, low_terms: int
) -> tuple[float, float]:
    """Bound low- and high-input columns in coefficient ``l2`` norm."""

    patch_count = len(routes.atlas.centers)
    radius = np.abs(routes.target_base)
    factor = np.abs(routes.factor)
    denominator = 1 - radius * radius
    high_rows = (
        factor
        * factor
        * radius ** (2 * (low_terms + 1))
        / denominator
    )
    low_rows = (
        factor
        * factor
        * (1 - radius ** (2 * (low_terms + 1)))
        / denominator
    )
    high_by_target = np.zeros(patch_count, dtype=np.float64)
    low_by_target = np.zeros(patch_count, dtype=np.float64)
    for target in range(patch_count):
        mask = routes.targets == target
        high_by_target[target] = math.fsum(high_rows[mask].ravel()) / routes.samples
        low_by_target[target] = math.fsum(low_rows[mask].ravel()) / routes.samples
    return float(np.sqrt(low_by_target.max())), float(
        np.sqrt(high_by_target.max())
    )


def branch_normalized_polynomial_bound(
    basis: np.ndarray,
    routes: hp.RouteTable,
    evaluation_ratio: float,
    anchors: np.ndarray,
) -> dict[str, object]:
    patch_count = len(routes.atlas.centers)
    mode_count = routes.terms + 1
    branch_patch = patch_count - 1
    require(
        routes.atlas.labels[branch_patch][0] == "b",
        "final chart is not the order-23 branch chart",
    )
    block = basis[
        branch_patch * mode_count : (branch_patch + 1) * mode_count
    ]
    rotation = np.exp(1j * (np.pi - np.pi / 23))
    powers = np.arange(mode_count, dtype=np.float64)
    branch_series = block / (
        (routes.rho * rotation) ** powers[:, None]
    )
    leading = branch_series[:4]
    normalized = basis @ np.linalg.inv(leading)
    coefficient_blocks = normalized.reshape(patch_count, mode_count, 4)
    polynomial_bounds = np.sum(
        np.abs(coefficient_blocks)
        * evaluation_ratio ** powers[None, :, None],
        axis=1,
    )
    anchor_basis = basis @ np.linalg.inv(basis[anchors])
    anchor_block = anchor_basis[
        branch_patch * mode_count : (branch_patch + 1) * mode_count
    ]
    anchor_branch_series = anchor_block / (
        (routes.rho * rotation) ** powers[:, None]
    )
    anchor_leading = anchor_branch_series[:4]
    anchor_coefficient_blocks = anchor_basis.reshape(
        patch_count, mode_count, 4
    )
    anchor_polynomial_bounds = np.sum(
        np.abs(anchor_coefficient_blocks)
        * evaluation_ratio ** powers[None, :, None],
        axis=1,
    )
    return {
        "branch_leading_condition": float(np.linalg.cond(leading)),
        "maximum_scaled_coefficient": float(np.abs(normalized).max()),
        "maximum_weighted_polynomial_by_column": [
            float(value) for value in polynomial_bounds.max(axis=0)
        ],
        "maximum_weighted_polynomial": float(polynomial_bounds.max()),
        "anchor_normalized_leading_inverse_two_norm": float(
            np.linalg.norm(np.linalg.inv(anchor_leading), ord=2)
        ),
        "anchor_normalized_leading_sigma_minimum": float(
            np.linalg.svd(anchor_leading, compute_uv=False)[-1]
        ),
        "anchor_normalized_maximum_weighted_polynomial": float(
            anchor_polynomial_bounds.max()
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--class-id", type=int, choices=range(1, 8), required=True)
    parser.add_argument("--low-terms", type=int, default=60)
    parser.add_argument("--samples", type=int, default=1280)
    parser.add_argument("--rho", type=float, default=float(hp.DEFAULT_RHO))
    parser.add_argument("--covering-radius", type=float, default=0.471)
    parser.add_argument("--outer-radius", type=float, default=0.99)
    parser.add_argument("--full-terms", type=int, default=480)
    parser.add_argument(
        "--matrix-backend", choices=("direct", "fft"), default="direct"
    )
    parser.add_argument("--certify-left-inverse", action="store_true")
    parser.add_argument("--precision", type=int, default=192)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    require(arguments.samples >= arguments.low_terms + 1, "too few samples")
    require(
        0 < arguments.covering_radius < arguments.rho < arguments.outer_radius < 1,
        "need covering radius < rho < outer radius < 1",
    )

    started = time.perf_counter()
    geometry = hp.triangle_geometry()
    atlas = hp.build_atlas(geometry)
    routes = hp.build_route_table(
        arguments.class_id,
        arguments.low_terms,
        arguments.samples,
        arguments.rho,
        geometry,
        atlas,
    )
    operator = hp.MultiCentreHejhalOperator(routes)
    basis, _, solver = hp.solve_fixed_space(operator)
    _, _, pivot = qr(basis.T, mode="economic", pivoting=True)
    anchors = np.asarray(pivot[:4], dtype=np.int64)

    matrix_started = time.perf_counter()
    matrix = dense_low_mode_matrix(routes, backend=arguments.matrix_backend)
    matrix_seconds = time.perf_counter() - matrix_started
    identity_rows = np.eye(matrix.shape[1], dtype=np.complex128)[anchors]
    augmented = np.vstack((matrix, identity_rows))
    singular_started = time.perf_counter()
    singular_values = svdvals(
        augmented,
        overwrite_a=not arguments.certify_left_inverse,
        check_finite=False,
    )
    singular_seconds = time.perf_counter() - singular_started
    sigma_minimum = float(singular_values[-1])
    sigma_next = float(singular_values[-2])

    low_norm, high_norm = route_column_norms(routes, arguments.low_terms)
    require(high_norm < 1, "high-mode target block is not a contraction")
    schur_perturbation = high_norm * low_norm / (1 - high_norm)
    schur_margin = sigma_minimum - schur_perturbation
    require(schur_margin > 0, "numerical Schur margin is not positive")
    low_solution_factor = 1 / (schur_margin * (1 - high_norm))
    high_solution_factor = (
        1 + low_norm * low_solution_factor
    ) / (1 - high_norm)
    full_inverse_bound_diagnostic = math.hypot(
        low_solution_factor, high_solution_factor
    )

    polynomial = branch_normalized_polynomial_bound(
        basis,
        routes,
        arguments.covering_radius / arguments.rho,
        anchors,
    )
    imaginary_parts = [center.imag for center in atlas.centers]
    outer_transition_factor = (
        max(imaginary_parts)
        / min(imaginary_parts)
        / (1 - arguments.outer_radius**2)
    )
    outer_ratio = arguments.covering_radius / arguments.outer_radius
    outer_tail_at_low_cutoff = (
        outer_ratio ** (arguments.low_terms + 1) / (1 - outer_ratio)
    )
    bootstrap_contraction = outer_transition_factor * outer_tail_at_low_cutoff
    require(
        bootstrap_contraction < 1,
        "outer-radius sup-norm bootstrap is not a contraction",
    )
    outer_sup_norm_diagnostic = (
        outer_transition_factor
        * polynomial["maximum_weighted_polynomial"]
        / (1 - bootstrap_contraction)
    )

    actual_target_radius = float(routes.diagnostics["maximum_target_radius"])
    target_ratio_outer = actual_target_radius / arguments.outer_radius
    target_tail = (
        target_ratio_outer ** (arguments.full_terms + 1)
        / (1 - target_ratio_outer)
    )
    source_ratio_outer = arguments.rho / arguments.outer_radius
    source_alias = (
        source_ratio_outer**arguments.samples
        / (1 - source_ratio_outer**arguments.samples)
    )

    left_inverse_certificate = None
    normalized_sup_norm_certificate = None
    if arguments.certify_left_inverse:
        require(
            arguments.matrix_backend == "direct",
            "left-inverse certificate requires the direct matrix backend",
        )
        left_inverse_certificate = floating_left_inverse_certificate(
            augmented, routes, arguments.precision
        )
        normalized_sup_norm_certificate = normalized_outer_sup_norm_certificate(
            basis,
            anchors,
            routes,
            left_inverse_certificate,
            arguments.covering_radius,
            arguments.outer_radius,
            arguments.full_terms,
            arguments.precision,
        )

    result = {
        "status": "PASS_NUMERICAL_HURWITZ_TAIL_STABILITY_DIAGNOSTIC",
        "scope": (
            "double-precision low-block singular gap and Schur-tail "
            "diagnostic; interval singular-value certification still required"
        ),
        "class_id": arguments.class_id,
        "low_terms": arguments.low_terms,
        "full_terms": arguments.full_terms,
        "samples": arguments.samples,
        "rho": arguments.rho,
        "covering_radius": arguments.covering_radius,
        "outer_radius": arguments.outer_radius,
        "patch_count": len(atlas.centers),
        "low_dimension": matrix.shape[1],
        "matrix_backend": arguments.matrix_backend,
        "anchor_indices": [int(value) for value in anchors],
        "low_augmented_sigma_minimum": sigma_minimum,
        "low_augmented_sigma_next": sigma_next,
        "low_augmented_sigma_maximum": float(singular_values[0]),
        "low_augmented_condition": float(singular_values[0] / singular_values[-1]),
        "high_input_operator_norm_bound": high_norm,
        "low_input_operator_norm_bound": low_norm,
        "schur_perturbation_bound": schur_perturbation,
        "schur_margin": schur_margin,
        "full_inverse_bound_diagnostic": full_inverse_bound_diagnostic,
        "actual_maximum_target_radius": actual_target_radius,
        "actual_maximum_target_base": actual_target_radius / arguments.rho,
        "maximum_route_factor": float(np.abs(routes.factor).max()),
        "outer_transition_factor_bound": outer_transition_factor,
        "outer_tail_at_low_cutoff": outer_tail_at_low_cutoff,
        "bootstrap_contraction": bootstrap_contraction,
        "outer_sup_norm_diagnostic": outer_sup_norm_diagnostic,
        "target_tail_at_full_cutoff_per_unit_outer_sup_norm": target_tail,
        "source_dft_alias_per_unit_outer_sup_norm": source_alias,
        "branch_normalization": polynomial,
        "solver": solver,
        "matrix_wall_seconds": matrix_seconds,
        "singular_values_wall_seconds": singular_seconds,
        "total_wall_seconds": time.perf_counter() - started,
    }
    if left_inverse_certificate is not None:
        result["left_inverse_certificate"] = left_inverse_certificate
        result["normalized_sup_norm_certificate"] = (
            normalized_sup_norm_certificate
        )
        result["status"] = "PASS_HURWITZ_TAIL_LOW_MODE_STABILITY_CERTIFICATE"
        result["scope"] = (
            "Acb-certified route geometry, a posteriori IEEE-754 low-mode "
            "left inverse, infinite high-mode Schur bound, normalized outer "
            "sup norm, and N=480/Q=1280 tail constants"
        )
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    print(rendered, end="")
    if arguments.output:
        arguments.output.write_text(rendered)


if __name__ == "__main__":
    main()
