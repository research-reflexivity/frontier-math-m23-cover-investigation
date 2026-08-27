#!/usr/bin/env sage-python
"""Certify the geometric input to the Hurwitz Taylor-tail estimate.

Run with ``sage -python``.  The certificate covers both halves of the
``(2,23,23)`` triangle by geodesic mesh cells in the Klein model.  For each
cell it chooses one mesh vertex and one atlas centre, then uses convexity of
hyperbolic balls and the pseudohyperbolic triangle inequality to bound every
point of the cell.  All transcendental coordinates and all reported upper
bounds are evaluated with Arb balls.

This proves the uniform atlas-covering part of the automorphic tail bound.
It deliberately does not certify the finite Hejhal matrix's normalized
inverse; that separate, finite-dimensional stability check is described in
``HURWITZ_TAIL_BOUND.md``.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from sage.all import ComplexBallField, I, RealBallField

HERE = Path(__file__).resolve().parent
SCRIPTS = HERE.parent / "scripts"
sys.path.insert(0, str(SCRIPTS))

import certify_hurwitz_acb as acb  # noqa: E402
import hurwitz_high_precision as hp  # noqa: E402


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def pseudohyperbolic_double(left: complex, right: complex) -> float:
    return abs((left - right) / (1 - right.conjugate() * left))


def pseudohyperbolic_ball(left, right):
    return ((left - right) / (1 - right.conjugate() * left)).abs()


def klein_to_poincare_ball(value):
    squared_absolute = value.real() ** 2 + value.imag() ** 2
    return value / (1 + (1 - squared_absolute).sqrt())


def mesh_cells(order: int):
    for first in range(order):
        for second in range(order - first):
            yield (
                (first, second),
                (first + 1, second),
                (first, second + 1),
            )
            if first + second <= order - 2:
                yield (
                    (first + 1, second + 1),
                    (first + 1, second),
                    (first, second + 1),
                )


def rational_mesh_point(vertices, first: int, second: int, order: int):
    return (
        vertices[0]
        + (vertices[1] - vertices[0]) * first / order
        + (vertices[2] - vertices[0]) * second / order
    )


def ball_mesh_point(vertices, first: int, second: int, order: int, field):
    return (
        vertices[0]
        + (vertices[1] - vertices[0]) * field(first) / order
        + (vertices[2] - vertices[0]) * field(second) / order
    )


def half_triangle_cover(
    order: int,
    double_vertices: list[complex],
    ball_vertices: list,
    double_centers: list[complex],
    ball_centers: list,
    real_field,
) -> dict[str, object]:
    double_points: dict[tuple[int, int], complex] = {}
    ball_points: dict[tuple[int, int], object] = {}
    selected_centers: dict[tuple[int, int], int] = {}
    selected_radius_balls: dict[tuple[int, int], object] = {}

    for first in range(order + 1):
        for second in range(order + 1 - first):
            key = (first, second)
            klein_double = rational_mesh_point(
                double_vertices, first, second, order
            )
            point_double = hp.klein_to_poincare(klein_double)
            klein_ball = ball_mesh_point(
                ball_vertices, first, second, order, real_field
            )
            point_ball = klein_to_poincare_ball(klein_ball)
            distances = [
                pseudohyperbolic_double(point_double, center)
                for center in double_centers
            ]
            selected = int(np.argmin(distances))
            double_points[key] = point_double
            ball_points[key] = point_ball
            selected_centers[key] = selected
            selected_radius_balls[key] = pseudohyperbolic_ball(
                point_ball, ball_centers[selected]
            )

    maximum_bound = real_field(0)
    worst_cell = None
    cell_count = 0
    for keys in mesh_cells(order):
        cell_count += 1
        # Select the best of the three possible base vertices using only
        # midpoint geometry.  The subsequent ball calculation certifies the
        # selected choice; it need not certify that this choice was optimal.
        choices = []
        for base_key in keys:
            base = double_points[base_key]
            epsilon = max(
                pseudohyperbolic_double(base, double_points[other])
                for other in keys
            )
            radius = pseudohyperbolic_double(
                base, double_centers[selected_centers[base_key]]
            )
            choices.append(((radius + epsilon) / (1 + radius * epsilon), base_key))
        _, base_key = min(choices)

        base_ball = ball_points[base_key]
        epsilon_ball = max(
            pseudohyperbolic_ball(base_ball, ball_points[other])
            for other in keys
        )
        radius_ball = selected_radius_balls[base_key]
        cell_bound = (radius_ball + epsilon_ball) / (
            1 + radius_ball * epsilon_ball
        )
        if cell_bound.upper() > maximum_bound.upper():
            maximum_bound = cell_bound
            worst_cell = {
                "vertices": [list(key) for key in keys],
                "base_vertex": list(base_key),
                "atlas_center": selected_centers[base_key],
                "base_to_center": str(radius_ball),
                "cell_pseudohyperbolic_radius": str(epsilon_ball),
                "combined_bound": str(cell_bound),
            }

    return {
        "mesh_vertex_count": len(ball_points),
        "mesh_cell_count": cell_count,
        "covering_radius_ball": str(maximum_bound),
        "covering_radius_upper": float(maximum_bound.upper()),
        "worst_cell": worst_cell,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh-order", type=int, default=160)
    parser.add_argument("--precision", type=int, default=192)
    parser.add_argument("--required-radius", type=float, default=0.471)
    parser.add_argument("--rho", type=float, default=float(hp.DEFAULT_RHO))
    parser.add_argument("--outer-radius", type=float, default=0.99)
    parser.add_argument("--samples", type=int, default=1280)
    parser.add_argument("--terms", type=int, nargs="+", default=[360, 480])
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    require(arguments.mesh_order >= 2, "mesh order must be at least two")
    require(
        0 < arguments.required_radius < arguments.rho < arguments.outer_radius < 1,
        "need 0 < required radius < rho < outer radius < 1",
    )

    real_field = RealBallField(arguments.precision)
    complex_field = ComplexBallField(arguments.precision)
    geometry = hp.triangle_geometry()
    atlas = hp.build_atlas(geometry)
    double_centers = [hp.disc_coordinate(center) for center in atlas.centers]

    ball = acb.build_ball_geometry(arguments.precision)
    vertex_a_ball = complex_field(I)
    ball_atlas_centers = acb.atlas_centers(ball, atlas)
    ball_centers = [
        acb.disc_coordinate(center, vertex_a_ball)
        for center in ball_atlas_centers
    ]

    vertex_a, vertex_b, vertex_c = hp.triangle_vertices(geometry)
    double_halves = [
        [
            hp.poincare_to_klein(hp.disc_coordinate(vertex))
            for vertex in (vertex_a, vertex_b, third)
        ]
        for third in (vertex_c, -vertex_c.real + 1j * vertex_c.imag)
    ]
    reflected_c_ball = -ball.vertex_c.real() + complex_field(I) * ball.vertex_c.imag()
    ball_halves = [
        [
            acb.poincare_to_klein(acb.disc_coordinate(vertex, vertex_a_ball))
            for vertex in (vertex_a_ball, ball.vertex_b, third)
        ]
        for third in (ball.vertex_c, reflected_c_ball)
    ]

    halves = [
        half_triangle_cover(
            arguments.mesh_order,
            double_vertices,
            ball_vertices,
            double_centers,
            ball_centers,
            real_field,
        )
        for double_vertices, ball_vertices in zip(double_halves, ball_halves)
    ]
    covering_upper = max(half["covering_radius_upper"] for half in halves)
    require(
        covering_upper < arguments.required_radius,
        "certified covering radius does not meet the requested bound",
    )

    # Use the declared rational decimal radius in every downstream bound.
    # The float above is only a human-readable projection of the Arb upper
    # endpoint and must not be rounded back into a purportedly sharp ball.
    radius = real_field(str(arguments.required_radius))
    rho = real_field(str(arguments.rho))
    outer = real_field(str(arguments.outer_radius))
    ratio_rho = radius / rho
    ratio_outer = radius / outer
    source_ratio = rho / outer
    tail_bounds = {}
    for terms in arguments.terms:
        require(terms >= 0, "Taylor cutoff must be nonnegative")
        relative_rho = ratio_rho ** (terms + 1) / (1 - ratio_rho)
        relative_outer = ratio_outer ** (terms + 1) / (1 - ratio_outer)
        tail_bounds[str(terms)] = {
            "target_tail_per_unit_M_rho_ball": str(relative_rho),
            "target_tail_per_unit_M_rho_upper": float(relative_rho.upper()),
            "target_tail_per_unit_M_outer_ball": str(relative_outer),
            "target_tail_per_unit_M_outer_upper": float(relative_outer.upper()),
        }
    alias = source_ratio ** arguments.samples / (
        1 - source_ratio ** arguments.samples
    )

    result = {
        "status": "PASS_HURWITZ_UNIFORM_ATLAS_COVER",
        "scope": (
            "Arb-certified all-points atlas cover and geometric-series tail "
            "constants; finite normalized Hejhal inverse not included"
        ),
        "precision_bits": arguments.precision,
        "mesh_order": arguments.mesh_order,
        "patch_count": len(atlas.centers),
        "required_covering_radius": arguments.required_radius,
        "certified_covering_radius_upper": covering_upper,
        "covering_radius_used_in_tail_bounds": arguments.required_radius,
        "rho": arguments.rho,
        "outer_radius": arguments.outer_radius,
        "covering_to_rho_ratio_ball": str(ratio_rho),
        "covering_to_outer_ratio_ball": str(ratio_outer),
        "dft_sample_count": arguments.samples,
        "source_alias_per_unit_M_outer_ball": str(alias),
        "source_alias_per_unit_M_outer_upper": float(alias.upper()),
        "tail_bounds": tail_bounds,
        "halves": halves,
    }
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    print(rendered, end="")
    if arguments.output:
        arguments.output.write_text(rendered)


if __name__ == "__main__":
    main()
