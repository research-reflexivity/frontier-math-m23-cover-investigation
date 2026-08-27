#!/usr/bin/env sage-python
"""Compute one mixed-precision Acb canonical quadric for an M23 cover."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from sage.all import ComplexField, I, matrix, vector

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import certify_hurwitz_acb as acb  # noqa: E402
import hurwitz_high_precision as hp  # noqa: E402


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def refine_vector(
    operator: hp.MultiCentreHejhalOperator,
    evaluator: acb.CachedAcbHejhalEvaluator,
    basis: np.ndarray,
    anchors: list[int],
    column: int,
    precision: int,
    neumann_iterations: int,
    refine_rounds: int,
) -> tuple[list, dict[str, object]]:
    coefficients = acb.high_precision_coefficients(
        basis[:, column], anchors, column, precision + 64
    )
    residual, initial_metrics = evaluator.evaluate(coefficients)
    coefficient_field = ComplexField(precision + 64)
    anchor_set = set(anchors)
    rounds: list[dict[str, object]] = []
    final_metrics = initial_metrics
    for round_index in range(refine_rounds):
        correction, correction_metrics = acb.projected_neumann_correction(
            operator,
            basis,
            anchors,
            acb.residual_midpoints(residual),
            neumann_iterations,
        )
        for index, value in enumerate(correction):
            if index in anchor_set:
                continue
            coefficients[index] += (
                coefficient_field(format(float(value.real), ".17g"))
                + coefficient_field(I)
                * coefficient_field(format(float(value.imag), ".17g"))
            )
        residual, final_metrics = evaluator.evaluate(coefficients)
        rounds.append(
            {
                "round": round_index + 1,
                **correction_metrics,
                "residual_norm_after": final_metrics["residual_norm_midpoint"],
            }
        )
    return coefficients, {
        "column": column,
        "initial_residual_norm": initial_metrics["residual_norm_midpoint"],
        "final_residual_norm": final_metrics["residual_norm_midpoint"],
        "final_residual_norm_ball": final_metrics["residual_norm_ball"],
        "rounds": rounds,
    }


def branch_series(
    refined_basis: list[list],
    routes: hp.RouteTable,
    precision: int,
) -> tuple[list[list], dict[str, object], object]:
    field = ComplexField(precision + 64)
    branch_patch = len(routes.atlas.centers) - 1
    require(routes.atlas.labels[branch_patch][0] == "b", "missing b chart")
    start = branch_patch * (routes.terms + 1)
    rho = field(format(routes.rho, ".17g"))
    theta = field.pi() - field.pi() / 23
    q_rotation = (field(I) * theta).exp()
    raw = matrix(
        field,
        routes.terms + 1,
        4,
        lambda row, column: refined_basis[column][start + row]
        / (rho * q_rotation) ** row,
    )
    leading = raw[:4, :]
    normalized = raw * leading.inverse()
    echelon_error = max(
        abs(normalized[row, column] - (1 if row == column else 0))
        for row in range(4)
        for column in range(4)
    )
    return (
        [
            [normalized[row, column] for column in range(4)]
            for row in range(routes.terms + 1)
        ],
        {
            "branch_patch": branch_patch,
            "echelon_error": str(echelon_error),
            "q_rotation": str(q_rotation),
        },
        leading.inverse(),
    )


def marked_branch_point(
    refined_basis: list[list],
    routes: hp.RouteTable,
    basis_change,
    A,
    precision: int,
    label: str,
    series_row: int,
) -> dict[str, object]:
    """Read a marked point from one optional central atlas chart."""

    field = ComplexField(precision + 64)
    patches = [
        index
        for index, atlas_label in enumerate(routes.atlas.labels)
        if atlas_label[0] == label
    ]
    require(len(patches) == 1, f"the atlas must contain exactly one {label} chart")
    patch = patches[0]
    start = patch * (routes.terms + 1)
    raw_point = matrix(
        field,
        1,
        4,
        lambda _, column: refined_basis[column][start + series_row],
    )
    canonical_point = (raw_point * basis_change).list()
    require(abs(canonical_point[0]) > field("1e-30"), "c point left y_0 chart")
    canonical_point = [value / canonical_point[0] for value in canonical_point]
    scaled_point = [
        value * field(A) ** index for index, value in enumerate(canonical_point)
    ]
    scaled_point = [value / scaled_point[0] for value in scaled_point]
    return {
        "patch": patch,
        "series_row": series_row,
        "canonical_coordinates_x0_equal_1": [str(value) for value in canonical_point],
        "scale_free_coordinates_y0_equal_1": [str(value) for value in scaled_point],
    }


def opposite_branch_point(
    refined_basis: list[list],
    routes: hp.RouteTable,
    basis_change,
    A,
    precision: int,
) -> dict[str, object]:
    """Read the opposite order-23 point from its optional central chart."""

    return marked_branch_point(
        refined_basis, routes, basis_change, A, precision, "c", 0
    )


def order_two_branch_point(
    refined_basis: list[list],
    routes: hp.RouteTable,
    basis_change,
    A,
    precision: int,
) -> dict[str, object]:
    """Read one point over the order-two vertex.

    The identity sheet is fixed by the order-two branch cycle, so pulled-back
    differentials vanish to first order at the elliptic centre.  Their common
    first Taylor coefficient, rather than the constant coefficient, gives the
    canonical point.
    """

    return marked_branch_point(
        refined_basis, routes, basis_change, A, precision, "a", 1
    )


def convolve(left: list, right: list, terms: int) -> list:
    field = left[0].parent()
    result = [field(0) for _ in range(terms)]
    for first, left_value in enumerate(left[:terms]):
        for second, right_value in enumerate(right[: terms - first]):
            result[first + second] += left_value * right_value
    return result


def canonical_quadric(series: list[list], precision: int) -> dict[str, object]:
    field = ComplexField(precision + 64)
    jet_terms = min(len(series), 20)
    monomials = [(i, j) for i in range(4) for j in range(i, 4)]
    columns = [
        convolve(
            [row[first] for row in series],
            [row[second] for row in series],
            jet_terms,
        )
        for first, second in monomials
    ]
    relation_matrix = matrix(
        field,
        jet_terms,
        len(monomials),
        lambda row, column: columns[column][row],
    )
    normalized_index = monomials.index((0, 2))
    unknown = [index for index in range(len(monomials)) if index != normalized_index]
    approximate = np.asarray(
        [
            [complex(relation_matrix[row, column]) for column in unknown]
            for row in range(jet_terms)
        ],
        dtype=np.complex128,
    )
    approximate /= np.maximum(
        np.linalg.norm(approximate, axis=1)[:, None], 1e-300
    )
    from scipy.linalg import qr

    _, _, row_pivots = qr(approximate.T, mode="economic", pivoting=True)
    selected_rows = [int(value) for value in row_pivots[: len(unknown)]]
    solve_scales = {
        row: max(abs(relation_matrix[row, column]) for column in unknown)
        for row in selected_rows
    }
    square = matrix(
        field,
        [
            [
                relation_matrix[row, column] / solve_scales[row]
                for column in unknown
            ]
            for row in selected_rows
        ],
    )
    right_hand_side = vector(
        field,
        [
            -relation_matrix[row, normalized_index] / solve_scales[row]
            for row in selected_rows
        ],
    )
    solved = square.solve_right(right_hand_side)
    coefficients = [field(0) for _ in monomials]
    coefficients[normalized_index] = field(1)
    for index, value in zip(unknown, solved):
        coefficients[index] = value
    residual = relation_matrix * vector(field, coefficients)
    maximum_residual = max(abs(value) for value in residual)
    residual_norm = sum(abs(value) ** 2 for value in residual).sqrt()
    row_norms = [
        sum(abs(relation_matrix[row, column]) ** 2 for column in range(len(monomials))).sqrt()
        for row in range(jet_terms)
    ]
    balanced_residual = [
        abs(residual[row]) / row_norms[row]
        for row in range(jet_terms)
        if row_norms[row] != 0
    ]
    return {
        "monomials": [list(monomial) for monomial in monomials],
        "coefficients": [str(value) for value in coefficients],
        "A": str(coefficients[monomials.index((0, 3))]),
        "selected_jet_rows": selected_rows,
        "jet_terms": jet_terms,
        "maximum_jet_residual": str(maximum_residual),
        "jet_residual_norm": str(residual_norm),
        "maximum_balanced_jet_residual": str(max(balanced_residual)),
        "balanced_jet_residual_norm": str(
            sum(value**2 for value in balanced_residual).sqrt()
        ),
    }


def canonical_cubic(
    series: list[list], quadric: dict[str, object], precision: int
) -> dict[str, object]:
    """Choose a Galois-compatible Petri cubic modulo quadric multiples.

    The four coefficients indexed by ``gauge_monomials`` kill the four
    multiples ``x_i Q``.  The coefficient of ``x_1*x_2*x_3`` is then set to
    one.  Unlike an orthogonal-complement convention, these coordinate
    conditions are algebraic and can be applied identically at every
    embedding of the Hurwitz algebra.
    """

    field = ComplexField(precision + 64)
    # Twenty jets already determine the five-dimensional cubic kernel.  The
    # higher local coefficients are the first to feel the finite Fourier
    # cutoff, so using a longer prefix can hide the fifth relation even when
    # the four exact quadric multiples remain visible.
    jet_terms = min(len(series), 20)
    quadratic_monomials = [
        tuple(monomial) for monomial in quadric["monomials"]
    ]
    quadric_coefficients = [
        field(coefficient) for coefficient in quadric["coefficients"]
    ]
    cubic_monomials = [
        (i, j, k)
        for i in range(4)
        for j in range(i, 4)
        for k in range(j, 4)
    ]
    columns = [
        convolve(
            convolve(
                [row[first] for row in series],
                [row[second] for row in series],
                jet_terms,
            ),
            [row[third] for row in series],
            jet_terms,
        )
        for first, second, third in cubic_monomials
    ]
    relation_matrix = matrix(
        field,
        jet_terms,
        len(cubic_monomials),
        lambda row, column: columns[column][row],
    )

    monomial_index = {
        monomial: index for index, monomial in enumerate(cubic_monomials)
    }
    gauge_monomials = [(0, 0, 2), (0, 1, 2), (0, 2, 2), (0, 2, 3)]
    gauge_indices = [monomial_index[monomial] for monomial in gauge_monomials]
    multiplier_coefficients = [
        [field(0) for _ in range(4)] for _ in cubic_monomials
    ]
    for variable in range(4):
        for coefficient, pair in zip(quadric_coefficients, quadratic_monomials):
            triple = tuple(sorted((*pair, variable)))
            multiplier_coefficients[monomial_index[triple]][variable] += coefficient
    gauge_matrix = matrix(
        field,
        [[multiplier_coefficients[index][variable] for variable in range(4)]
         for index in gauge_indices],
    )
    require(gauge_matrix.det() != 0, "Petri cubic gauge is not transverse")

    # Locate a numerically safe coefficient of the one-dimensional quotient
    # relation.  This diagnostic is used below only to choose the affine
    # chart; the final normalization is still imposed by an exact coordinate
    # equation in the high-precision solve.
    full_approximate = np.asarray(
        [
            [complex(relation_matrix[row, column]) for column in range(20)]
            for row in range(jet_terms)
        ],
        dtype=np.complex128,
    )
    full_approximate /= np.maximum(
        np.linalg.norm(full_approximate, axis=1)[:, None], 1e-300
    )
    _, cubic_singular_values, cubic_vh = np.linalg.svd(
        full_approximate, full_matrices=False
    )
    cubic_kernel = cubic_vh.conjugate().T[:, -5:]
    multiplier_approximate = np.asarray(
        [[complex(value) for value in row] for row in multiplier_coefficients],
        dtype=np.complex128,
    )
    multiplier_orthogonal, _ = np.linalg.qr(multiplier_approximate)
    overlap = multiplier_orthogonal.conjugate().T @ cubic_kernel
    _, _, overlap_vh = np.linalg.svd(overlap, full_matrices=True)
    quotient_relation = cubic_kernel @ overlap_vh.conjugate().T[:, -1]
    quotient_relation -= multiplier_orthogonal @ (
        multiplier_orthogonal.conjugate().T @ quotient_relation
    )
    gauge_correction = np.linalg.solve(
        multiplier_approximate[gauge_indices, :],
        -quotient_relation[gauge_indices],
    )
    gauged_relation = quotient_relation + multiplier_approximate @ gauge_correction
    normalization_candidates = [
        index for index in range(20) if index not in set(gauge_indices)
    ]
    diagnostic_normalized_index = max(
        normalization_candidates, key=lambda index: abs(gauged_relation[index])
    )
    normalized_monomial = (1, 2, 3)
    normalized_index = monomial_index[normalized_monomial]
    require(
        abs(gauged_relation[normalized_index]) > 1e-4,
        "fixed Petri cubic normalization coefficient is numerically zero",
    )
    fixed_indices = set(gauge_indices + [normalized_index])
    unknown = [
        index for index in range(len(cubic_monomials)) if index not in fixed_indices
    ]

    approximate = np.asarray(
        [
            [complex(relation_matrix[row, column]) for column in unknown]
            for row in range(jet_terms)
        ],
        dtype=np.complex128,
    )
    approximate /= np.maximum(
        np.linalg.norm(approximate, axis=1)[:, None], 1e-300
    )
    from scipy.linalg import qr

    _, _, row_pivots = qr(approximate.T, mode="economic", pivoting=True)
    selected_rows = [int(value) for value in row_pivots[: len(unknown)]]
    solve_scales = {
        row: max(abs(relation_matrix[row, column]) for column in unknown)
        for row in selected_rows
    }
    square = matrix(
        field,
        [
            [
                relation_matrix[row, column] / solve_scales[row]
                for column in unknown
            ]
            for row in selected_rows
        ],
    )
    right_hand_side = vector(
        field,
        [
            -relation_matrix[row, normalized_index] / solve_scales[row]
            for row in selected_rows
        ],
    )
    solved = square.solve_right(right_hand_side)
    coefficients = [field(0) for _ in cubic_monomials]
    coefficients[normalized_index] = field(1)
    for index, value in zip(unknown, solved):
        coefficients[index] = value
    residual = relation_matrix * vector(field, coefficients)
    row_norms = [
        sum(
            abs(relation_matrix[row, column]) ** 2
            for column in range(len(cubic_monomials))
        ).sqrt()
        for row in range(jet_terms)
    ]
    balanced_residual = [
        abs(residual[row]) / row_norms[row]
        for row in range(jet_terms)
        if row_norms[row] != 0
    ]

    return {
        "monomials": [list(monomial) for monomial in cubic_monomials],
        "coefficients": [str(value) for value in coefficients],
        "gauge_monomials": [list(monomial) for monomial in gauge_monomials],
        "normalized_monomial": list(normalized_monomial),
        "normalization_coefficient_before_scaling": str(
            gauged_relation[normalized_index]
        ),
        "largest_diagnostic_monomial": list(
            cubic_monomials[diagnostic_normalized_index]
        ),
        "approximate_gauged_coefficients": [
            [float(value.real), float(value.imag)] for value in gauged_relation
        ],
        "approximate_smallest_singular_values": [
            float(value) for value in cubic_singular_values[-7:]
        ],
        "approximate_gauged_residual_norm": float(
            np.linalg.norm(full_approximate @ gauged_relation)
            / np.linalg.norm(gauged_relation)
        ),
        "gauge_determinant": str(gauge_matrix.det()),
        "selected_jet_rows": selected_rows,
        "jet_terms": jet_terms,
        "maximum_jet_residual": str(max(abs(value) for value in residual)),
        "jet_residual_norm": str(
            sum(abs(value) ** 2 for value in residual).sqrt()
        ),
        "maximum_balanced_jet_residual": str(max(balanced_residual)),
        "balanced_jet_residual_norm": str(
            sum(value**2 for value in balanced_residual).sqrt()
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--class-id", type=int, choices=range(1, 8), default=6)
    parser.add_argument("--terms", type=int, default=120)
    parser.add_argument("--samples", type=int, default=320)
    parser.add_argument("--rho", type=float, default=float(hp.DEFAULT_RHO))
    parser.add_argument("--precision", type=int, default=384)
    parser.add_argument("--neumann-iterations", type=int, default=160)
    parser.add_argument("--refine-rounds", type=int, default=1)
    parser.add_argument("--include-c-chart", action="store_true")
    parser.add_argument("--include-a-chart", action="store_true")
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    require(arguments.terms >= 23, "canonical quadric needs at least 24 jets")

    geometry = hp.triangle_geometry()
    atlas = hp.build_atlas(
        geometry,
        include_c=arguments.include_c_chart,
        include_a=arguments.include_a_chart,
    )
    routes = hp.build_route_table(
        arguments.class_id,
        arguments.terms,
        arguments.samples,
        arguments.rho,
        geometry,
        atlas,
    )
    acb.verify_route_membership(routes)
    routes, midpoint_metrics = acb.reevaluate_route_midpoints(
        arguments.precision, routes
    )
    operator = hp.MultiCentreHejhalOperator(routes)
    basis, anchors = acb.normalized_fixed_basis(operator)
    evaluator = acb.CachedAcbHejhalEvaluator(arguments.precision, routes)
    refined_basis: list[list] = []
    refinement: list[dict[str, object]] = []
    for column in range(4):
        coefficients, metrics = refine_vector(
            operator,
            evaluator,
            basis,
            anchors,
            column,
            arguments.precision,
            arguments.neumann_iterations,
            arguments.refine_rounds,
        )
        refined_basis.append(coefficients)
        refinement.append(metrics)
        print(
            f"refined column {column}: residual {metrics['final_residual_norm']:.6e}",
            file=sys.stderr,
            flush=True,
        )
    all_mode_residuals: list[dict[str, object]] = []
    for column, coefficients in enumerate(refined_basis):
        _, metrics = evaluator.evaluate(coefficients, all_modes=True)
        all_mode_residuals.append(
            {
                "column": column,
                "residual_norm_ball": metrics["residual_norm_ball"],
                "residual_norm_midpoint": metrics["residual_norm_midpoint"],
                "residual_norm_radius": metrics["residual_norm_radius"],
                "output_mode_count": metrics["output_mode_count"],
            }
        )
    series, branch_metrics, basis_change = branch_series(
        refined_basis, routes, arguments.precision
    )
    branch_metrics["series_coefficient_order"] = (
        "rows q^0 through q^19; four branch-normalized canonical forms"
    )
    branch_metrics["series_first_20"] = [
        [str(value) for value in row] for row in series[:20]
    ]
    quadric = canonical_quadric(series, arguments.precision)
    cubic = canonical_cubic(series, quadric, arguments.precision)
    c_point = None
    if arguments.include_c_chart:
        c_point = opposite_branch_point(
            refined_basis,
            routes,
            basis_change,
            quadric["A"],
            arguments.precision,
        )
    a_point = None
    if arguments.include_a_chart:
        a_point = order_two_branch_point(
            refined_basis,
            routes,
            basis_change,
            quadric["A"],
            arguments.precision,
        )
    result = {
        "status": "PASS_MIXED_PRECISION_NUMERICAL_CANONICAL_QUADRIC",
        "scope": (
            "high-precision midpoint model with Acb-certified finite basis "
            "residuals; no uniform tail or algebraic-recognition claim"
        ),
        "class_id": arguments.class_id,
        "terms": arguments.terms,
        "samples": arguments.samples,
        "rho": arguments.rho,
        "precision_bits": arguments.precision,
        "patch_count": len(atlas.centers),
        "dimension": operator.dimension,
        "anchor_indices": anchors,
        "route_midpoint_refinement": midpoint_metrics,
        "acb_evaluation_backend": evaluator.evaluation_backend,
        "acb_fourier_backend": evaluator.fourier_backend,
        "basis_refinement": refinement,
        "all_dft_mode_residuals": all_mode_residuals,
        "branch": branch_metrics,
        "opposite_branch_point": c_point,
        "order_two_branch_point": a_point,
        "canonical_quadric": quadric,
        "petri_cubic": cubic,
    }
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    print(rendered, end="")
    if arguments.output:
        arguments.output.write_text(rendered)


if __name__ == "__main__":
    main()
