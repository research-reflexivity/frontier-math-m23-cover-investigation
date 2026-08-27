#!/usr/bin/env python3
"""Multi-centre, matrix-free Hejhal pilot for the seven M23 covers.

The original numerical pilot uses one Taylor centre in each of the 23 long
triangle patches.  That is an excellent double-precision calculation, but
its covering radius is too close to one for a practical 100-digit solve.
This module uses two 23-point orbits of centres on the long ``b--c`` edge,
together with the point above ``b``.  It then applies the Fourier-projected
Hejhal operator without forming its dense matrix.

The default edge parameters are the small rational approximations
``27/500`` and ``677/1000``.  A deterministic Klein-model mesh observes a
worst radius below 0.46, compared with 0.87082 for the original atlas.
The rational parameters are intentional: the same atlas can be reconstructed
with Arb/Acb without treating rounded centre coordinates as exact input.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Sequence

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from compute_hurwitz_covers import (  # noqa: E402
    B_ORDER,
    DEGREE,
    X_REPRESENTATIVES,
    Y,
    TriangleGeometry,
    canonical_relations,
    cyclic_coset_exponents,
    cyclic_coset_matrices,
    disc_coordinate,
    mobius,
    reduce_to_subgroup,
    triangle_geometry,
    upper_half_plane_coordinate,
)


DEFAULT_EDGE_PARAMETERS = (Fraction(27, 500), Fraction(677, 1000))
DEFAULT_RHO = Fraction(18, 25)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def triangle_vertices(geometry: TriangleGeometry) -> tuple[complex, complex, complex]:
    """Return the elliptic vertices ``a,b,c`` in the upper half-plane."""

    cot_a = 1 / math.tan(math.pi / 2)
    cot_b = 1 / math.tan(math.pi / B_ORDER)
    real_c = (geometry.mu**2 - 1) / (2 * (cot_a + geometry.mu * cot_b))
    imag_c = math.sqrt(1 - (real_c - cot_a) ** 2)
    return 1j, geometry.mu * 1j, real_c + imag_c * 1j


def poincare_to_klein(w: complex) -> complex:
    return 2 * w / (1 + abs(w) ** 2)


def klein_to_poincare(k: complex) -> complex:
    require(abs(k) < 1 + 1e-13, "Klein point lies outside the disk")
    return k / (1 + math.sqrt(max(0.0, 1 - abs(k) ** 2)))


def klein_edge_center(
    geometry: TriangleGeometry, parameter: Fraction | float
) -> complex:
    """Interpolate from ``c`` to ``b`` along their hyperbolic geodesic.

    Geodesics are Euclidean segments in the Klein model, so a rational
    interpolation parameter has an unambiguous Acb reconstruction.
    """

    _, b, c = triangle_vertices(geometry)
    klein_b = poincare_to_klein(disc_coordinate(b))
    klein_c = poincare_to_klein(disc_coordinate(c))
    value = float(parameter)
    klein = (1 - value) * klein_c + value * klein_b
    return upper_half_plane_coordinate(klein_to_poincare(klein))


@dataclass(frozen=True)
class Atlas:
    centers: np.ndarray
    labels: tuple[tuple[str, int, int, int, int], ...]
    edge_parameters: tuple[Fraction, ...]


def build_atlas(
    geometry: TriangleGeometry,
    edge_parameters: Sequence[Fraction] = DEFAULT_EDGE_PARAMETERS,
    include_b: bool = True,
    include_c: bool = False,
    include_a: bool = False,
) -> Atlas:
    """Build two edge-centre orbits and optional ramification-point charts."""

    cosets = cyclic_coset_matrices(Y, geometry.delta_b)
    exponents = cyclic_coset_exponents(Y)
    centers: list[complex] = []
    labels: list[tuple[str, int, int, int, int]] = []
    parameters = tuple(edge_parameters)
    for orbit_index, parameter in enumerate(parameters):
        base = klein_edge_center(geometry, parameter)
        for sheet, representative in enumerate(cosets):
            centers.append(mobius(representative, base))
            labels.append(
                (
                    "edge",
                    orbit_index,
                    parameter.numerator,
                    parameter.denominator,
                    exponents[sheet],
                )
            )
    if include_a:
        centers.append(1j)
        labels.append(("a", 0, 1, 1, 0))
    if include_c:
        centers.append(triangle_vertices(geometry)[2])
        labels.append(("c", 0, 1, 1, 0))
    if include_b:
        centers.append(geometry.mu * 1j)
        labels.append(("b", 0, 1, 1, 0))
    array = np.asarray(centers, dtype=np.complex128)
    require(
        len(array)
        == DEGREE * len(parameters) + int(include_b) + int(include_c) + int(include_a),
        "bad atlas size",
    )
    return Atlas(array, tuple(labels), parameters)


def sample_triangle_domain(
    geometry: TriangleGeometry, mesh_order: int
) -> np.ndarray:
    """Sample both halves of the doubled triangle in Klein barycentrics."""

    require(mesh_order >= 2, "mesh order must be at least two")
    a, b, c = triangle_vertices(geometry)
    samples: list[complex] = []
    for third in (c, -c.real + 1j * c.imag):
        vertices = [
            poincare_to_klein(disc_coordinate(vertex))
            for vertex in (a, b, third)
        ]
        for first in range(mesh_order + 1):
            for second in range(mesh_order + 1 - first):
                u = first / mesh_order
                v = second / mesh_order
                klein = (
                    (1 - u - v) * vertices[0]
                    + u * vertices[1]
                    + v * vertices[2]
                )
                samples.append(
                    upper_half_plane_coordinate(klein_to_poincare(klein))
                )
    return np.asarray(samples, dtype=np.complex128)


def atlas_covering_diagnostics(
    geometry: TriangleGeometry, atlas: Atlas, mesh_order: int = 160
) -> dict[str, object]:
    points = sample_triangle_domain(geometry, mesh_order)
    minimum_radii = np.full(points.shape, np.inf, dtype=np.float64)
    for center in atlas.centers:
        minimum_radii = np.minimum(
            minimum_radii, np.abs((points - center) / (points - center.conjugate()))
        )
    worst_index = int(np.argmax(minimum_radii))
    covering_radius = float(minimum_radii[worst_index])
    rho = float(DEFAULT_RHO)
    digits = 100
    estimated_terms = math.ceil(
        digits * math.log(10) / math.log(rho / covering_radius)
    )
    return {
        "patch_count": int(len(atlas.centers)),
        "edge_parameters": [
            [parameter.numerator, parameter.denominator]
            for parameter in atlas.edge_parameters
        ],
        "mesh_order": mesh_order,
        "mesh_point_count": int(len(points)),
        "maximum_radius_on_mesh": covering_radius,
        "worst_mesh_point": [
            float(points[worst_index].real),
            float(points[worst_index].imag),
        ],
        "default_rho": rho,
        "covering_to_cauchy_ratio": covering_radius / rho,
        "geometric_terms_for_100_digits": estimated_terms,
        "original_single_orbit_radius": geometry.domain_radius,
    }


@dataclass(frozen=True)
class RouteTable:
    class_id: int
    terms: int
    samples: int
    rho: float
    weight: int
    atlas: Atlas
    targets: np.ndarray
    target_base: np.ndarray
    factor: np.ndarray
    coset_exponents: np.ndarray
    triangle_words: tuple[tuple[tuple[str, int], ...], ...]
    diagnostics: dict[str, object]


def build_route_table(
    class_id: int,
    terms: int,
    samples: int,
    rho: float,
    geometry: TriangleGeometry,
    atlas: Atlas,
    weight: int = 2,
) -> RouteTable:
    """Reduce all source-circle samples and select their nearest chart."""

    require(1 <= class_id <= 7, "class id must lie between one and seven")
    require(samples >= terms + 1, "Fourier grid is too short for the cutoff")
    require(0 < rho < 1, "Cauchy radius must lie in the unit disk")
    sigma_a = X_REPRESENTATIVES[class_id - 1]
    cosets = cyclic_coset_matrices(Y, geometry.delta_b)
    patch_count = len(atlas.centers)
    targets = np.empty((patch_count, samples), dtype=np.int64)
    target_base = np.empty((patch_count, samples), dtype=np.complex128)
    factor = np.empty((patch_count, samples), dtype=np.complex128)
    coset_exponents = np.empty((patch_count, samples), dtype=np.int64)
    triangle_words: list[tuple[tuple[str, int], ...]] = []
    maximum_target_radius = 0.0
    minimum_target_margin = 1.0
    maximum_iterations = 0
    maximum_roundtrip_error = 0.0
    angles = 2 * math.pi * np.arange(samples, dtype=np.float64) / samples
    source_w = rho * np.exp(1j * angles)

    for source, center in enumerate(atlas.centers):
        for sample, w in enumerate(source_w):
            z = upper_half_plane_coordinate(complex(w), complex(center))
            reduction = reduce_to_subgroup(
                z, geometry, sigma_a, Y, cosets
            )
            radii = np.abs(
                (reduction.reduced_z - atlas.centers)
                / (reduction.reduced_z - np.conjugate(atlas.centers))
            )
            target = int(np.argmin(radii))
            target_w = disc_coordinate(
                reduction.reduced_z, complex(atlas.centers[target])
            )
            c_entry = reduction.gamma[1, 0]
            d_entry = reduction.gamma[1, 1]
            automorphy = (c_entry * z + d_entry) ** -weight
            targets[source, sample] = target
            target_base[source, sample] = target_w / rho
            factor[source, sample] = (
                automorphy * (1 - target_w) ** weight / (1 - w) ** weight
            )
            coset_exponents[source, sample] = reduction.coset_exponent
            triangle_words.append(reduction.triangle_word)
            target_radius = abs(target_w)
            maximum_target_radius = max(maximum_target_radius, target_radius)
            minimum_target_margin = min(minimum_target_margin, rho - target_radius)
            maximum_iterations = max(
                maximum_iterations, reduction.triangle_iterations
            )
            reconstructed = upper_half_plane_coordinate(
                target_w, complex(atlas.centers[target])
            )
            maximum_roundtrip_error = max(
                maximum_roundtrip_error,
                abs(reduction.reduced_z - reconstructed),
            )

    require(minimum_target_margin > 0, "atlas does not contain a reduced sample")
    diagnostics: dict[str, object] = {
        "class_id": class_id,
        "terms": terms,
        "modes_per_patch": terms + 1,
        "samples": samples,
        "rho": rho,
        "weight": weight,
        "patch_count": patch_count,
        "operator_dimension": patch_count * (terms + 1),
        "maximum_target_radius": maximum_target_radius,
        "minimum_target_margin": minimum_target_margin,
        "maximum_triangle_iterations": maximum_iterations,
        "maximum_reduction_roundtrip_error": maximum_roundtrip_error,
    }
    return RouteTable(
        class_id,
        terms,
        samples,
        rho,
        weight,
        atlas,
        targets,
        target_base,
        factor,
        coset_exponents,
        tuple(triangle_words),
        diagnostics,
    )


class MultiCentreHejhalOperator:
    """Square Fourier projection of the multi-centre consistency equations."""

    def __init__(self, routes: RouteTable, rhs_batch_size: int = 4):
        self.routes = routes
        self.patch_count = len(routes.atlas.centers)
        self.mode_count = routes.terms + 1
        self.samples = routes.samples
        self.dimension = self.patch_count * self.mode_count
        self.rhs_batch_size = max(1, rhs_batch_size)
        self.targets = routes.targets.ravel()
        self.target_base = routes.target_base.ravel()
        self.factor = routes.factor.ravel()
        self.target_rows = [
            np.flatnonzero(self.targets == target)
            for target in range(self.patch_count)
        ]

    def _matrix_input(self, vectors: np.ndarray) -> tuple[np.ndarray, bool]:
        array = np.asarray(vectors, dtype=np.complex128)
        squeeze = array.ndim == 1
        if squeeze:
            array = array[:, None]
        require(
            array.ndim == 2 and array.shape[0] == self.dimension,
            "wrong operator input shape",
        )
        return array, squeeze

    @staticmethod
    def _restore(array: np.ndarray, squeeze: bool) -> np.ndarray:
        return array[:, 0] if squeeze else array

    def apply_H(self, vectors: np.ndarray) -> np.ndarray:
        matrix, squeeze = self._matrix_input(vectors)
        result = np.empty_like(matrix)
        for start in range(0, matrix.shape[1], self.rhs_batch_size):
            stop = min(start + self.rhs_batch_size, matrix.shape[1])
            block = matrix[:, start:stop].reshape(
                self.patch_count, self.mode_count, -1
            )
            values = block[self.targets, self.mode_count - 1, :].copy()
            for mode in range(self.mode_count - 2, -1, -1):
                values *= self.target_base[:, None]
                values += block[self.targets, mode, :]
            values *= self.factor[:, None]
            values = values.reshape(self.patch_count, self.samples, -1)
            projected = np.fft.fft(values, axis=1)[
                :, : self.mode_count, :
            ] / self.samples
            result[:, start:stop] = (projected - block).reshape(
                self.dimension, -1
            )
        return self._restore(result, squeeze)

    def apply_H_adjoint(self, vectors: np.ndarray) -> np.ndarray:
        matrix, squeeze = self._matrix_input(vectors)
        result = np.empty_like(matrix)
        for start in range(0, matrix.shape[1], self.rhs_batch_size):
            stop = min(start + self.rhs_batch_size, matrix.shape[1])
            block = matrix[:, start:stop].reshape(
                self.patch_count, self.mode_count, -1
            )
            spectrum = np.zeros(
                (self.patch_count, self.samples, stop - start),
                dtype=np.complex128,
            )
            spectrum[:, : self.mode_count, :] = block
            sample_dual = np.fft.ifft(spectrum, axis=1).reshape(
                self.patch_count * self.samples, -1
            )
            weighted = sample_dual * self.factor.conjugate()[:, None]
            pulled = np.zeros_like(block)
            conjugate_power = np.ones_like(self.target_base)
            for mode in range(self.mode_count):
                contributions = weighted * conjugate_power[:, None]
                for target, rows in enumerate(self.target_rows):
                    if len(rows):
                        pulled[target, mode, :] = np.sum(
                            contributions[rows, :], axis=0
                        )
                conjugate_power *= self.target_base.conjugate()
            result[:, start:stop] = (pulled - block).reshape(
                self.dimension, -1
            )
        return self._restore(result, squeeze)

    def scipy_linear_operator(self):
        from scipy.sparse.linalg import LinearOperator

        return LinearOperator(
            (self.dimension, self.dimension),
            matvec=self.apply_H,
            rmatvec=self.apply_H_adjoint,
            matmat=self.apply_H,
            rmatmat=self.apply_H_adjoint,
            dtype=np.complex128,
        )

    def adjoint_error(self, seed: int = 23) -> float:
        generator = np.random.default_rng(seed)
        left = generator.normal(size=self.dimension) + 1j * generator.normal(
            size=self.dimension
        )
        right = generator.normal(size=self.dimension) + 1j * generator.normal(
            size=self.dimension
        )
        lhs = np.vdot(left, self.apply_H(right))
        rhs = np.vdot(self.apply_H_adjoint(left), right)
        return float(abs(lhs - rhs) / max(1.0, abs(lhs), abs(rhs)))


def solve_fixed_space(
    operator: MultiCentreHejhalOperator,
    block_size: int = 10,
    iterations: int = 100,
) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    """Find the four-plane as the fixed space of the compact pullback ``T``.

    Asking a scalar Krylov method directly for the smallest singular values
    of ``T-I`` is unnecessarily difficult, and a scalar start cannot recover
    a nearly repeated eigenvalue reliably.  The federalist Cauchy operator
    ``T`` is compact, while the weight-two forms have eigenvalue one, so block
    subspace iteration isolates the desired fixed space quickly.  We still
    verify the returned vectors by applying ``T-I`` itself.
    """

    require(block_size > 4, "need at least one eigenvalue beyond genus four")
    started = time.perf_counter()
    generator = np.random.default_rng(23)
    subspace = generator.normal(
        size=(operator.dimension, block_size)
    ) + 1j * generator.normal(
        size=(operator.dimension, block_size)
    )
    subspace, _ = np.linalg.qr(subspace)
    for _ in range(iterations):
        pulled = operator.apply_H(subspace) + subspace
        subspace, _ = np.linalg.qr(pulled)
    pulled = operator.apply_H(subspace) + subspace
    rayleigh = subspace.conjugate().T @ pulled
    eigenvalues, eigenvectors = np.linalg.eig(rayleigh)
    order = np.argsort(np.abs(eigenvalues - 1))
    eigenvalues = eigenvalues[order]
    selected = subspace @ eigenvectors[:, order[:4]]
    basis, _ = np.linalg.qr(selected)
    residual_norms = [
        float(np.linalg.norm(operator.apply_H(basis[:, column])))
        for column in range(4)
    ]
    fixed_distances = np.abs(eigenvalues - 1)
    diagnostics = {
        "dominant_eigenvalues_by_distance_from_one": [
            [float(value.real), float(value.imag)] for value in eigenvalues
        ],
        "distances_from_one": [float(value) for value in fixed_distances],
        "fixed_space_separation_ratio": float(
            fixed_distances[4] / max(fixed_distances[3], 1e-300)
        ),
        "basis_residual_norms": residual_norms,
        "block_size": block_size,
        "subspace_iterations": iterations,
        "solver_wall_seconds": time.perf_counter() - started,
    }
    return basis, eigenvalues, diagnostics


def branch_series_from_atlas_basis(
    basis: np.ndarray,
    routes: RouteTable,
) -> tuple[np.ndarray, dict[str, object]]:
    """Read the order-23 local series directly from the central atlas chart."""

    branch_patch = len(routes.atlas.centers) - 1
    require(
        routes.atlas.labels[branch_patch][0] == "b",
        "the final atlas chart is not the order-23 centre",
    )
    block = basis[
        branch_patch * (routes.terms + 1) :
        (branch_patch + 1) * (routes.terms + 1),
        :,
    ]
    theta = math.pi / 23
    q_rotation = complex(math.cos(math.pi - theta), math.sin(math.pi - theta))
    powers = np.arange(routes.terms + 1, dtype=np.float64)
    coefficients = block / (
        (routes.rho * q_rotation) ** powers[:, None]
    )
    leading = coefficients[:4, :]
    condition = float(np.linalg.cond(leading))
    normalized = coefficients @ np.linalg.inv(leading)
    diagnostics: dict[str, object] = {
        "branch_patch": branch_patch,
        "branch_leading_matrix_condition": condition,
        "branch_echelon_error": float(
            np.max(np.abs(normalized[:4, :] - np.eye(4)))
        ),
        "q_rotation": [q_rotation.real, q_rotation.imag],
    }
    return normalized, diagnostics


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--class-id", type=int, choices=range(1, 8), default=6)
    parser.add_argument("--terms", type=int, default=30)
    parser.add_argument("--samples", type=int, default=96)
    parser.add_argument("--rho", type=float, default=float(DEFAULT_RHO))
    parser.add_argument("--mesh-order", type=int, default=160)
    parser.add_argument("--atlas-only", action="store_true")
    parser.add_argument("--solve", action="store_true")
    parser.add_argument("--canonical", action="store_true")
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    if arguments.canonical:
        arguments.solve = True
        require(arguments.terms >= 19, "canonical extraction needs 20 jets")

    geometry = triangle_geometry()
    atlas = build_atlas(geometry)
    result: dict[str, object] = {
        "status": "PASS_MULTICENTRE_ATLAS_GEOMETRY",
        "signature": [2, 23, 23],
        "atlas": atlas_covering_diagnostics(
            geometry, atlas, mesh_order=arguments.mesh_order
        ),
    }
    if not arguments.atlas_only:
        routes = build_route_table(
            arguments.class_id,
            arguments.terms,
            arguments.samples,
            arguments.rho,
            geometry,
            atlas,
        )
        operator = MultiCentreHejhalOperator(routes)
        adjoint_error = operator.adjoint_error()
        require(adjoint_error < 2e-12, "matrix-free adjoint check failed")
        result["routes"] = routes.diagnostics
        result["operator_adjoint_error"] = adjoint_error
        result["status"] = "PASS_MULTICENTRE_MATRIX_FREE_OPERATOR"
        if arguments.solve:
            basis, _, solver = solve_fixed_space(operator)
            result["solver"] = solver
            result["status"] = "PASS_MULTICENTRE_GENUS_FOUR_SOLVE"
            if arguments.canonical:
                series, branch = branch_series_from_atlas_basis(basis, routes)
                result["branch"] = branch
                result["canonical_relations"] = canonical_relations(series)
                result["status"] = "PASS_MULTICENTRE_CANONICAL_MODEL"
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    print(rendered, end="")
    if arguments.output:
        arguments.output.write_text(rendered)


if __name__ == "__main__":
    main()
