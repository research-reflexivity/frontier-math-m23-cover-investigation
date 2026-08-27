#!/usr/bin/env sage-python
"""Two-precision Arb/Acb check of the finite multi-centre Hejhal equations.

Run this file with ``sage -python``.  All transcendental geometry is rebuilt
from the exact signature ``(2,23,23)``, the rational atlas parameters, and
the recorded triangle-reduction words.  The finite Fourier residual is then
evaluated in Sage's ``ComplexBallField`` (FLINT/Acb) at two precisions.

This certifies the finite equations for the exported decimal coefficients.
It does not, by itself, bound the omitted Taylor tail or prove that rounded
double-precision coefficients are exact algebraic numbers.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import math
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
from sage.all import ComplexBallField, ComplexField, I, RealBallField

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import hurwitz_high_precision as hp  # noqa: E402
from compute_hurwitz_covers import (  # noqa: E402
    X_REPRESENTATIVES,
    Y,
    compose_right,
    permutation_power,
)

try:  # Optional Cython bridge to FLINT's quasi-linear Acb DFT.
    import acb_fft  # type: ignore[import-not-found]  # noqa: E402
except ImportError:
    acb_fft = None


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def matrix_identity(field):
    return ((field(1), field(0)), (field(0), field(1)))


def matrix_multiply(left, right):
    return tuple(
        tuple(
            sum(left[row][middle] * right[middle][column] for middle in range(2))
            for column in range(2)
        )
        for row in range(2)
    )


def matrix_inverse(matrix):
    determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    return (
        (matrix[1][1] / determinant, -matrix[0][1] / determinant),
        (-matrix[1][0] / determinant, matrix[0][0] / determinant),
    )


def matrix_power(matrix, exponent: int, field):
    if exponent < 0:
        return matrix_power(matrix_inverse(matrix), -exponent, field)
    result = matrix_identity(field)
    base = matrix
    while exponent:
        if exponent & 1:
            result = matrix_multiply(result, base)
        base = matrix_multiply(base, base)
        exponent >>= 1
    return result


def mobius(matrix, z):
    return (matrix[0][0] * z + matrix[0][1]) / (
        matrix[1][0] * z + matrix[1][1]
    )


def disc_coordinate(z, center):
    return (z - center) / (z - center.conjugate())


def upper_half_plane_coordinate(w, center):
    return (center - w * center.conjugate()) / (1 - w)


@dataclass(frozen=True)
class BallGeometry:
    field: object
    pi: object
    mu: object
    delta_a: tuple
    delta_b: tuple
    vertex_b: object
    vertex_c: object


def build_ball_geometry(precision: int) -> BallGeometry:
    field = ComplexBallField(precision)
    pi = field.pi()
    angle = pi / 23
    sine = angle.sin()
    cosine = angle.cos()
    lam = cosine / sine
    mu = lam + (lam * lam - 1).sqrt()
    delta_a = ((field(0), field(1)), (-field(1), field(0)))
    delta_b = (
        (cosine, mu * sine),
        (-sine / mu, cosine),
    )
    real_c = (mu * mu - 1) / (2 * mu * cosine / sine)
    imag_c = (1 - real_c * real_c).sqrt()
    vertex_b = field(I) * mu
    vertex_c = real_c + field(I) * imag_c
    return BallGeometry(
        field, pi, mu, delta_a, delta_b, vertex_b, vertex_c
    )


def poincare_to_klein(w):
    squared_absolute = w.real() * w.real() + w.imag() * w.imag()
    return 2 * w / (1 + squared_absolute)


def klein_to_poincare(k):
    squared_absolute = k.real() * k.real() + k.imag() * k.imag()
    return k / (1 + (1 - squared_absolute).sqrt())


def atlas_centers(ball: BallGeometry, atlas: hp.Atlas) -> list:
    field = ball.field
    vertex_a = field(I)
    disk_b = disc_coordinate(ball.vertex_b, vertex_a)
    disk_c = disc_coordinate(ball.vertex_c, vertex_a)
    klein_b = poincare_to_klein(disk_b)
    klein_c = poincare_to_klein(disk_c)
    bases: list = []
    for parameter in atlas.edge_parameters:
        scalar = field(parameter.numerator) / parameter.denominator
        klein = (1 - scalar) * klein_c + scalar * klein_b
        bases.append(
            upper_half_plane_coordinate(klein_to_poincare(klein), vertex_a)
        )
    centers: list = []
    for kind, orbit, _, _, exponent in atlas.labels:
        if kind == "b":
            centers.append(ball.vertex_b)
        elif kind == "c":
            centers.append(ball.vertex_c)
        elif kind == "a":
            centers.append(vertex_a)
        else:
            centers.append(
                mobius(
                    matrix_power(ball.delta_b, exponent, field), bases[orbit]
                )
            )
    return centers


def route_gamma(ball: BallGeometry, coset_exponent: int, word):
    field = ball.field
    delta = matrix_identity(field)
    for generator, exponent in word:
        matrix = ball.delta_a if generator == "a" else ball.delta_b
        delta = matrix_multiply(matrix_power(matrix, exponent, field), delta)
    return matrix_multiply(
        matrix_power(ball.delta_b, coset_exponent, field), delta
    )


def decimal_ball(field, value):
    if not isinstance(value, (complex, float, np.complexfloating, np.floating)):
        return field(value)
    real = field(format(float(value.real), ".17g"))
    imag = field(format(float(value.imag), ".17g"))
    return real + field(I) * imag


def verify_route_membership(routes: hp.RouteTable) -> None:
    sigma_a = X_REPRESENTATIVES[routes.class_id - 1]
    identity = np.arange(len(Y), dtype=np.int64)
    for index, (coset_exponent, word) in enumerate(
        zip(routes.coset_exponents.ravel(), routes.triangle_words)
    ):
        delta = identity.copy()
        for generator, exponent in word:
            permutation = sigma_a if generator == "a" else Y
            delta = compose_right(permutation_power(permutation, exponent), delta)
        gamma = compose_right(permutation_power(Y, int(coset_exponent)), delta)
        require(int(gamma[0]) == 0, f"route {index} is not in the subgroup")


def reevaluate_route_midpoints(
    precision: int, routes: hp.RouteTable
) -> tuple[hp.RouteTable, dict[str, float]]:
    """Replace rounded reduction geometry by Acb-computed midpoints."""

    ball = build_ball_geometry(precision)
    field = ball.field
    centers = atlas_centers(ball, routes.atlas)
    rho = field(format(routes.rho, ".17g"))
    target_base = np.empty_like(routes.target_base)
    factor = np.empty_like(routes.factor)
    maximum_target_base_change = 0.0
    maximum_factor_change = 0.0
    maximum_target_base_error_upper = 0.0
    maximum_factor_error_upper = 0.0
    maximum_target_radius_upper = 0.0
    row = 0
    for source, center in enumerate(centers):
        for sample in range(routes.samples):
            angle = 2 * ball.pi * sample / routes.samples
            source_w = rho * (angle.cos() + field(I) * angle.sin())
            z = upper_half_plane_coordinate(source_w, center)
            gamma = route_gamma(
                ball,
                int(routes.coset_exponents[source, sample]),
                routes.triangle_words[row],
            )
            reduced = mobius(gamma, z)
            target = int(routes.targets[source, sample])
            target_w = disc_coordinate(reduced, centers[target])
            ball_base = target_w / rho
            automorphy = (gamma[1][0] * z + gamma[1][1]) ** -2
            ball_factor = (
                automorphy * (1 - target_w) ** 2 / (1 - source_w) ** 2
            )
            target_base[source, sample] = complex(
                float(ball_base.real().center()),
                float(ball_base.imag().center()),
            )
            factor[source, sample] = complex(
                float(ball_factor.real().center()),
                float(ball_factor.imag().center()),
            )
            maximum_target_base_change = max(
                maximum_target_base_change,
                abs(target_base[source, sample] - routes.target_base[source, sample]),
            )
            maximum_factor_change = max(
                maximum_factor_change,
                abs(factor[source, sample] - routes.factor[source, sample]),
            )
            maximum_target_base_error_upper = max(
                maximum_target_base_error_upper,
                float(
                    (
                        ball_base
                        - decimal_ball(
                            field, routes.target_base[source, sample]
                        )
                    ).abs().upper()
                ),
            )
            maximum_factor_error_upper = max(
                maximum_factor_error_upper,
                float(
                    (
                        ball_factor
                        - decimal_ball(field, routes.factor[source, sample])
                    ).abs().upper()
                ),
            )
            maximum_target_radius_upper = max(
                maximum_target_radius_upper, float(target_w.abs().upper())
            )
            row += 1
    require(
        maximum_target_radius_upper < routes.rho,
        "Acb-refined route leaves its Cauchy circle",
    )
    refined = dataclasses.replace(routes, target_base=target_base, factor=factor)
    return refined, {
        "precision_bits": precision,
        "maximum_target_base_change": maximum_target_base_change,
        "maximum_factor_change": maximum_factor_change,
        "maximum_target_base_error_upper_vs_complex128": (
            maximum_target_base_error_upper
        ),
        "maximum_factor_error_upper_vs_complex128": maximum_factor_error_upper,
        "maximum_target_radius_upper": maximum_target_radius_upper,
    }


class CachedAcbHejhalEvaluator:
    """Cache ball-valued geometry and Fourier weights across coefficient runs."""

    def __init__(self, precision: int, routes: hp.RouteTable):
        self.precision = precision
        self.routes = routes
        self.ball = build_ball_geometry(precision)
        self.field = self.ball.field
        self.centers = atlas_centers(self.ball, routes.atlas)
        self.rho = self.field(format(routes.rho, ".17g"))
        self.target_base: list = []
        self.factor: list = []
        self.maximum_target_radius = 0.0
        self.maximum_target_base_midpoint_error = 0.0
        self.maximum_factor_midpoint_error = 0.0
        row = 0
        for source, center in enumerate(self.centers):
            for sample in range(routes.samples):
                angle = 2 * self.ball.pi * sample / routes.samples
                source_w = self.rho * (
                    angle.cos() + self.field(I) * angle.sin()
                )
                z = upper_half_plane_coordinate(source_w, center)
                gamma = route_gamma(
                    self.ball,
                    int(routes.coset_exponents[source, sample]),
                    routes.triangle_words[row],
                )
                reduced = mobius(gamma, z)
                target = int(routes.targets[source, sample])
                target_w = disc_coordinate(reduced, self.centers[target])
                base = target_w / self.rho
                automorphy = (gamma[1][0] * z + gamma[1][1]) ** -2
                factor = (
                    automorphy
                    * (1 - target_w) ** 2
                    / (1 - source_w) ** 2
                )
                self.target_base.append(base)
                self.factor.append(factor)
                base_midpoint = complex(
                    float(base.real().center()),
                    float(base.imag().center()),
                )
                factor_midpoint = complex(
                    float(factor.real().center()),
                    float(factor.imag().center()),
                )
                self.maximum_target_base_midpoint_error = max(
                    self.maximum_target_base_midpoint_error,
                    abs(base_midpoint - routes.target_base[source, sample]),
                )
                self.maximum_factor_midpoint_error = max(
                    self.maximum_factor_midpoint_error,
                    abs(factor_midpoint - routes.factor[source, sample]),
                )
                self.maximum_target_radius = max(
                    self.maximum_target_radius, float(target_w.abs().upper())
                )
                row += 1
        require(
            self.maximum_target_radius < routes.rho,
            "cached Acb target disk leaves the Cauchy circle",
        )
        self.evaluation_backend = (
            "compiled_flint_acb_horner" if acb_fft is not None else "python_horner"
        )
        self.fourier_backend = "flint_acb_dft" if acb_fft is not None else "direct"
        self.weights = None
        if acb_fft is None:
            self.weights = [
                [
                    (
                        -2
                        * self.ball.pi
                        * mode
                        * sample
                        / routes.samples
                    ).cos()
                    + self.field(I)
                    * (
                        -2
                        * self.ball.pi
                        * mode
                        * sample
                        / routes.samples
                    ).sin()
                    for sample in range(routes.samples)
                ]
                for mode in range(routes.terms + 1)
            ]

    def evaluate(
        self, coefficients: Sequence, all_modes: bool = False
    ) -> tuple[list, dict[str, object]]:
        started = time.perf_counter()
        routes = self.routes
        coefficient_balls = [
            decimal_ball(self.field, value) for value in coefficients
        ]
        blocks = [
            coefficient_balls[
                patch * (routes.terms + 1) :
                (patch + 1) * (routes.terms + 1)
            ]
            for patch in range(len(self.centers))
        ]
        if acb_fft is not None:
            flat_values = acb_fft.evaluate_routes(
                blocks,
                [int(value) for value in routes.targets.ravel()],
                self.target_base,
                self.factor,
            )
            values = [
                flat_values[
                    source * routes.samples : (source + 1) * routes.samples
                ]
                for source in range(len(self.centers))
            ]
        else:
            values = [
                [self.field(0) for _ in range(routes.samples)]
                for _ in self.centers
            ]
            row = 0
            for source in range(len(self.centers)):
                for sample in range(routes.samples):
                    target = int(routes.targets[source, sample])
                    polynomial = blocks[target][-1]
                    for mode in range(routes.terms - 1, -1, -1):
                        polynomial = (
                            polynomial * self.target_base[row]
                            + blocks[target][mode]
                        )
                    values[source][sample] = self.factor[row] * polynomial
                    row += 1
        residual: list = []
        squared_norm = RealBallField(self.precision)(0)
        maximum_component_radius = 0.0
        output_mode_count = (
            routes.samples if all_modes else routes.terms + 1
        )
        for source in range(len(self.centers)):
            transform = None
            if acb_fft is not None:
                transform = acb_fft.dft(values[source])
            for mode in range(output_mode_count):
                if transform is not None:
                    total = transform[mode] / routes.samples
                else:
                    total = self.field(0)
                    for sample in range(routes.samples):
                        if self.weights is not None and mode < len(self.weights):
                            weight = self.weights[mode][sample]
                        else:
                            angle = (
                                -2
                                * self.ball.pi
                                * mode
                                * sample
                                / routes.samples
                            )
                            weight = angle.cos() + self.field(I) * angle.sin()
                        total += values[source][sample] * weight
                    total /= routes.samples
                if mode <= routes.terms:
                    total -= blocks[source][mode]
                residual.append(total)
                absolute = total.abs()
                squared_norm += absolute * absolute
                maximum_component_radius = max(
                    maximum_component_radius, float(total.rad())
                )
        norm = squared_norm.sqrt()
        return residual, {
            "precision_bits": self.precision,
            "residual_norm_ball": str(norm),
            "residual_norm_midpoint": float(norm.center()),
            "residual_norm_radius": float(norm.rad()),
            "maximum_component_radius": maximum_component_radius,
            "output_mode_count": output_mode_count,
            "maximum_target_radius_upper": self.maximum_target_radius,
            "maximum_target_base_midpoint_error_vs_complex128": (
                self.maximum_target_base_midpoint_error
            ),
            "maximum_factor_midpoint_error_vs_complex128": (
                self.maximum_factor_midpoint_error
            ),
            "evaluation_backend": self.evaluation_backend,
            "fourier_backend": self.fourier_backend,
            "wall_seconds": time.perf_counter() - started,
        }


def evaluate_residual(
    precision: int,
    routes: hp.RouteTable,
    coefficients: Sequence,
) -> tuple[list, dict[str, object]]:
    started = time.perf_counter()
    ball = build_ball_geometry(precision)
    field = ball.field
    centers = atlas_centers(ball, routes.atlas)
    rho = field(format(routes.rho, ".17g"))
    coefficient_balls = [decimal_ball(field, value) for value in coefficients]
    blocks = [
        coefficient_balls[
            patch * (routes.terms + 1) : (patch + 1) * (routes.terms + 1)
        ]
        for patch in range(len(centers))
    ]
    values: list[list] = [
        [field(0) for _ in range(routes.samples)] for _ in centers
    ]
    maximum_target_radius = 0.0
    maximum_target_base_midpoint_error = 0.0
    maximum_factor_midpoint_error = 0.0
    row = 0
    for source, center in enumerate(centers):
        for sample in range(routes.samples):
            angle = 2 * ball.pi * sample / routes.samples
            source_w = rho * (angle.cos() + field(I) * angle.sin())
            z = upper_half_plane_coordinate(source_w, center)
            gamma = route_gamma(
                ball,
                int(routes.coset_exponents[source, sample]),
                routes.triangle_words[row],
            )
            reduced = mobius(gamma, z)
            target = int(routes.targets[source, sample])
            target_w = disc_coordinate(reduced, centers[target])
            target_base = target_w / rho
            automorphy = (gamma[1][0] * z + gamma[1][1]) ** -2
            factor = automorphy * (1 - target_w) ** 2 / (1 - source_w) ** 2
            polynomial = blocks[target][-1]
            for mode in range(routes.terms - 1, -1, -1):
                polynomial = polynomial * target_base + blocks[target][mode]
            values[source][sample] = factor * polynomial
            target_radius = float(target_w.abs().upper())
            maximum_target_radius = max(maximum_target_radius, target_radius)
            target_midpoint = complex(
                float(target_base.real().center()),
                float(target_base.imag().center()),
            )
            factor_midpoint = complex(
                float(factor.real().center()),
                float(factor.imag().center()),
            )
            maximum_target_base_midpoint_error = max(
                maximum_target_base_midpoint_error,
                abs(target_midpoint - routes.target_base[source, sample]),
            )
            maximum_factor_midpoint_error = max(
                maximum_factor_midpoint_error,
                abs(factor_midpoint - routes.factor[source, sample]),
            )
            row += 1

    residual: list = []
    squared_norm = RealBallField(precision)(0)
    maximum_component_radius = 0.0
    for source in range(len(centers)):
        for mode in range(routes.terms + 1):
            total = field(0)
            for sample in range(routes.samples):
                angle = -2 * ball.pi * mode * sample / routes.samples
                weight = angle.cos() + field(I) * angle.sin()
                total += values[source][sample] * weight
            total /= routes.samples
            total -= blocks[source][mode]
            residual.append(total)
            absolute = total.abs()
            squared_norm += absolute * absolute
            maximum_component_radius = max(
                maximum_component_radius, float(total.rad())
            )
    norm = squared_norm.sqrt()
    return residual, {
        "precision_bits": precision,
        "residual_norm_ball": str(norm),
        "residual_norm_midpoint": float(norm.center()),
        "residual_norm_radius": float(norm.rad()),
        "maximum_component_radius": maximum_component_radius,
        "maximum_target_radius_upper": maximum_target_radius,
        "maximum_target_base_midpoint_error_vs_complex128": (
            maximum_target_base_midpoint_error
        ),
        "maximum_factor_midpoint_error_vs_complex128": (
            maximum_factor_midpoint_error
        ),
        "wall_seconds": time.perf_counter() - started,
    }


def overlap_failures(low: Sequence, high: Sequence) -> int:
    require(len(low) == len(high), "cross-precision residual lengths differ")
    failures = 0
    for low_value, high_value in zip(low, high):
        if not low_value.real().overlaps(high_value.real()):
            failures += 1
        if not low_value.imag().overlaps(high_value.imag()):
            failures += 1
    return failures


def normalized_fixed_basis(operator: hp.MultiCentreHejhalOperator) -> tuple[np.ndarray, list[int]]:
    from scipy.linalg import qr

    basis, _, _ = hp.solve_fixed_space(operator)
    _, _, pivot = qr(basis.T, mode="economic", pivoting=True)
    anchors = [int(value) for value in pivot[:4]]
    chart = basis[anchors, :]
    normalized = basis @ np.linalg.inv(chart)
    require(
        np.max(np.abs(normalized[anchors, :] - np.eye(4))) < 5e-12,
        "fixed-space chart normalization failed",
    )
    return normalized, anchors


def high_precision_coefficients(
    coefficients: np.ndarray,
    anchors: Sequence[int],
    vector: int,
    precision: int,
) -> list:
    field = ComplexField(precision)
    result = [
        field(format(float(value.real), ".17g"))
        + field(I) * field(format(float(value.imag), ".17g"))
        for value in coefficients
    ]
    for column, anchor in enumerate(anchors):
        result[anchor] = field(1 if column == vector else 0)
    return result


def residual_midpoints(residual: Sequence) -> np.ndarray:
    return np.asarray(
        [
            complex(
                float(value.real().center()),
                float(value.imag().center()),
            )
            for value in residual
        ],
        dtype=np.complex128,
    )


def projected_neumann_correction(
    operator: hp.MultiCentreHejhalOperator,
    normalized_basis: np.ndarray,
    anchors: Sequence[int],
    residual: np.ndarray,
    iterations: int,
) -> tuple[np.ndarray, dict[str, float | int]]:
    """Approximately invert ``I-T`` on the anchor-zero complement."""

    projection_basis = normalized_basis @ np.linalg.inv(
        normalized_basis[list(anchors), :]
    )

    def project(vector: np.ndarray) -> np.ndarray:
        return vector - projection_basis @ vector[list(anchors)]

    right_hand_side = project(residual)
    correction = np.zeros_like(residual)
    for _ in range(iterations):
        pulled = operator.apply_H(correction) + correction
        correction = project(right_hand_side + pulled)
    predicted = residual + operator.apply_H(correction)
    return correction, {
        "iterations": iterations,
        "correction_norm": float(np.linalg.norm(correction)),
        "predicted_residual_norm": float(np.linalg.norm(predicted)),
        "predicted_projected_residual_norm": float(
            np.linalg.norm(project(predicted))
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--class-id", type=int, choices=range(1, 8), default=6)
    parser.add_argument("--terms", type=int, default=20)
    parser.add_argument("--samples", type=int, default=64)
    parser.add_argument("--rho", type=float, default=float(hp.DEFAULT_RHO))
    parser.add_argument("--vector", type=int, choices=range(4), default=0)
    parser.add_argument("--precision-low", type=int, default=128)
    parser.add_argument("--precision-high", type=int, default=192)
    parser.add_argument("--refine-rounds", type=int, default=0)
    parser.add_argument("--neumann-iterations", type=int, default=120)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    require(
        arguments.precision_high > arguments.precision_low,
        "high precision must exceed low precision",
    )

    geometry = hp.triangle_geometry()
    atlas = hp.build_atlas(geometry)
    routes = hp.build_route_table(
        arguments.class_id,
        arguments.terms,
        arguments.samples,
        arguments.rho,
        geometry,
        atlas,
    )
    verify_route_membership(routes)
    routes, midpoint_refinement = reevaluate_route_midpoints(
        arguments.precision_high, routes
    )
    operator = hp.MultiCentreHejhalOperator(routes)
    basis, anchors = normalized_fixed_basis(operator)
    initial_coefficients = basis[:, arguments.vector]
    double_residual_norm = float(
        np.linalg.norm(operator.apply_H(initial_coefficients))
    )
    coefficients: Sequence = initial_coefficients
    if arguments.refine_rounds:
        coefficients = high_precision_coefficients(
            initial_coefficients,
            anchors,
            arguments.vector,
            arguments.precision_high + 64,
        )
    high_residual, high_metrics = evaluate_residual(
        arguments.precision_high, routes, coefficients
    )
    initial_high_metrics = dict(high_metrics)
    refinement: list[dict[str, object]] = []
    coefficient_field = ComplexField(arguments.precision_high + 64)
    for round_index in range(arguments.refine_rounds):
        correction, correction_metrics = projected_neumann_correction(
            operator,
            basis,
            anchors,
            residual_midpoints(high_residual),
            arguments.neumann_iterations,
        )
        for index, value in enumerate(correction):
            coefficients[index] += (
                coefficient_field(format(float(value.real), ".17g"))
                + coefficient_field(I)
                * coefficient_field(format(float(value.imag), ".17g"))
            )
        high_residual, high_metrics = evaluate_residual(
            arguments.precision_high, routes, coefficients
        )
        refinement.append(
            {
                "round": round_index + 1,
                **correction_metrics,
                "acb_residual_norm_after": high_metrics[
                    "residual_norm_midpoint"
                ],
            }
        )
    low_residual, low_metrics = evaluate_residual(
        arguments.precision_low, routes, coefficients
    )
    failures = overlap_failures(low_residual, high_residual)
    require(failures == 0, "cross-precision Acb residual balls do not overlap")
    require(
        high_metrics["maximum_target_radius_upper"] < arguments.rho,
        "Acb target disk is not contained in the Cauchy circle",
    )
    result = {
        "status": "PASS_ACB_FINITE_MULTICENTRE_HEJHAL_TWO_PRECISIONS",
        "scope": (
            "finite Fourier equations after mixed-precision coefficient "
            "refinement; no omitted-tail or algebraic-recognition claim"
            if arguments.refine_rounds
            else "finite Fourier equations for decimal complex128 coefficients; "
            "no omitted-tail or algebraic-recognition claim"
        ),
        "class_id": arguments.class_id,
        "terms": arguments.terms,
        "samples": arguments.samples,
        "rho": arguments.rho,
        "patch_count": len(atlas.centers),
        "dimension": operator.dimension,
        "vector": arguments.vector,
        "anchor_indices": anchors,
        "route_subgroup_membership_failures": 0,
        "route_midpoint_refinement": midpoint_refinement,
        "complex128_residual_norm": double_residual_norm,
        "initial_high_precision": initial_high_metrics,
        "mixed_precision_refinement": refinement,
        "low_precision": low_metrics,
        "high_precision": high_metrics,
        "cross_precision_component_overlap_failures": failures,
    }
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    print(rendered, end="")
    if arguments.output:
        arguments.output.write_text(rendered)


if __name__ == "__main__":
    main()
