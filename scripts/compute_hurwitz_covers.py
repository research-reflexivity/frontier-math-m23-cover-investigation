#!/usr/bin/env python3
"""Numerical pilot for the seven (2A,23A,23B) M23 covers.

This is a clean-room implementation of the cocompact triangle-group method
of Klug--Musty--Schiavone--Voight.  It deliberately starts with the smallest
useful milestone: recover the four-dimensional space of weight-two forms for
each Nielsen triple.  Algebraic reconstruction is only attempted after that
analytic calibration has passed for the distinguished rational cover.

For the numerical solve we use the signature (2,23,23), with sigma_a=x and
sigma_b=y.  The 23 cosets are represented by powers of delta_b.  Following
the ``federalist'' variant of the algorithm, we expand simultaneously at the
23 centers alpha_i*z_a; every patch then only needs to cover one triangle.

The script needs NumPy.  In the Codex workspace it can be run with the Python
runtime reported by ``codex_app.load_workspace_dependencies``.
"""

from __future__ import annotations

import argparse
import cmath
import json
import math
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np


A_ORDER = 2
B_ORDER = 23
C_ORDER = 23
DEGREE = 23


def permutation_from_cycles(cycles: Iterable[Sequence[int]]) -> np.ndarray:
    """Return a zero-based right-action permutation from one-based cycles."""

    result = np.arange(DEGREE, dtype=np.int64)
    for cycle in cycles:
        for source, target in zip(cycle, cycle[1:] + cycle[:1]):
            result[source - 1] = target - 1
    return result


Y = permutation_from_cycles(
    [
        (
            1,
            2,
            11,
            10,
            16,
            9,
            6,
            3,
            23,
            19,
            20,
            14,
            21,
            17,
            4,
            8,
            22,
            5,
            18,
            15,
            13,
            7,
            12,
        )
    ]
)


X_REPRESENTATIVES = [
    permutation_from_cycles(
        [(4, 16), (5, 10), (6, 21), (8, 19), (9, 18), (11, 23), (13, 14), (17, 22)]
    ),
    permutation_from_cycles(
        [(4, 19), (5, 9), (6, 17), (8, 16), (10, 18), (11, 14), (13, 23), (21, 22)]
    ),
    permutation_from_cycles(
        [(4, 6), (5, 9), (7, 21), (8, 15), (11, 19), (12, 18), (13, 20), (14, 17)]
    ),
    permutation_from_cycles(
        [(3, 17), (5, 11), (7, 18), (8, 16), (9, 21), (12, 19), (14, 22), (15, 23)]
    ),
    permutation_from_cycles(
        [(3, 11), (5, 17), (6, 20), (7, 18), (8, 19), (10, 13), (12, 16), (14, 22)]
    ),
    permutation_from_cycles(
        [(3, 21), (4, 16), (9, 22), (10, 15), (11, 20), (12, 19), (13, 18), (14, 17)]
    ),
    permutation_from_cycles(
        [(3, 8), (5, 12), (6, 15), (9, 10), (11, 16), (13, 21), (17, 19), (20, 23)]
    ),
]


def compose_right(first: np.ndarray, second: np.ndarray) -> np.ndarray:
    """Right-action product: i^(first*second)=(i^first)^second."""

    return second[first]


def inverse_permutation(permutation: np.ndarray) -> np.ndarray:
    result = np.empty_like(permutation)
    result[permutation] = np.arange(len(permutation), dtype=np.int64)
    return result


def permutation_power(permutation: np.ndarray, exponent: int) -> np.ndarray:
    if exponent < 0:
        return permutation_power(inverse_permutation(permutation), -exponent)
    result = np.arange(len(permutation), dtype=np.int64)
    base = permutation.copy()
    while exponent:
        if exponent & 1:
            result = compose_right(result, base)
        base = compose_right(base, base)
        exponent >>= 1
    return result


def matrix_power(matrix: np.ndarray, exponent: int) -> np.ndarray:
    result = np.linalg.matrix_power(matrix, exponent)
    determinant = float(np.linalg.det(result))
    if determinant <= 0:
        raise ArithmeticError(f"nonpositive determinant {determinant}")
    return result / math.sqrt(determinant)


def mobius(matrix: np.ndarray, z: complex) -> complex:
    return (matrix[0, 0] * z + matrix[0, 1]) / (matrix[1, 0] * z + matrix[1, 1])


def disc_coordinate(z: complex, center: complex = 1j) -> complex:
    return (z - center) / (z - center.conjugate())


def upper_half_plane_coordinate(w: complex, center: complex = 1j) -> complex:
    return (center - w * center.conjugate()) / (1 - w)


def argument_0_2pi(z: complex) -> float:
    value = math.atan2(z.imag, z.real)
    return value if value >= 0 else value + 2 * math.pi


@dataclass(frozen=True)
class TriangleGeometry:
    mu: float
    delta_a: np.ndarray
    delta_b: np.ndarray
    disc_A: np.ndarray
    domain_radius: float


def triangle_geometry() -> TriangleGeometry:
    a = A_ORDER
    b = B_ORDER
    c = C_ORDER
    lambda_value = (
        math.cos(math.pi / a) * math.cos(math.pi / b) + math.cos(math.pi / c)
    ) / (math.sin(math.pi / a) * math.sin(math.pi / b))
    mu = lambda_value + math.sqrt(lambda_value * lambda_value - 1)
    delta_a = np.array(
        [
            [math.cos(math.pi / a), math.sin(math.pi / a)],
            [-math.sin(math.pi / a), math.cos(math.pi / a)],
        ],
        dtype=np.float64,
    )
    delta_b = np.array(
        [
            [math.cos(math.pi / b), mu * math.sin(math.pi / b)],
            [-math.sin(math.pi / b) / mu, math.cos(math.pi / b)],
        ],
        dtype=np.float64,
    )
    disc_A = np.array([[mu + 1, 1 - mu], [1 - mu, mu + 1]], dtype=np.float64)

    cot_a = 1 / math.tan(math.pi / a)
    cot_b = 1 / math.tan(math.pi / b)
    real_c = (mu * mu - 1) / (2 * (cot_a + mu * cot_b))
    imag_c = math.sqrt(1 / math.sin(math.pi / a) ** 2 - (real_c - cot_a) ** 2)
    vertices = [1j, mu * 1j, real_c + imag_c * 1j, -real_c + imag_c * 1j]
    domain_radius = max(abs(disc_coordinate(vertex)) for vertex in vertices)
    return TriangleGeometry(mu, delta_a, delta_b, disc_A, domain_radius)


@dataclass(frozen=True)
class Reduction:
    delta: np.ndarray
    delta_permutation: np.ndarray
    reduced_z: complex
    iterations: int
    word: tuple[tuple[str, int], ...]


def reduce_to_triangle(
    z: complex,
    geometry: TriangleGeometry,
    sigma_a: np.ndarray,
    sigma_b: np.ndarray,
    maximum_iterations: int = 1000,
) -> Reduction:
    """Algorithm 3.11, tracking both matrix and permutation word."""

    delta = np.eye(2, dtype=np.float64)
    delta_permutation = np.arange(DEGREE, dtype=np.int64)
    current_z = z
    word: list[tuple[str, int]] = []
    for iteration in range(1, maximum_iterations + 1):
        w = disc_coordinate(current_z)
        alpha = argument_0_2pi(w)
        exponent_a = -math.floor(A_ORDER * alpha / (2 * math.pi) + 0.5)
        if exponent_a % A_ORDER == 0:
            exponent_a = 0
        if exponent_a:
            power_matrix = matrix_power(geometry.delta_a, exponent_a)
            current_z = mobius(power_matrix, current_z)
            delta = power_matrix @ delta
            delta_permutation = compose_right(
                permutation_power(sigma_a, exponent_a), delta_permutation
            )
            word.append(("a", exponent_a))

        w = disc_coordinate(current_z)
        transformed = -(
            (geometry.disc_A[0, 0] * w + geometry.disc_A[0, 1])
            / (geometry.disc_A[1, 0] * w + geometry.disc_A[1, 1])
        )
        beta = argument_0_2pi(transformed)
        exponent_b = -math.floor(B_ORDER * beta / (2 * math.pi) + 0.5)
        if exponent_b % B_ORDER == 0:
            exponent_b = 0
        if exponent_b == 0:
            return Reduction(
                delta,
                delta_permutation,
                current_z,
                iteration,
                tuple(word),
            )
        power_matrix = matrix_power(geometry.delta_b, exponent_b)
        current_z = mobius(power_matrix, current_z)
        delta = power_matrix @ delta
        delta_permutation = compose_right(
            permutation_power(sigma_b, exponent_b), delta_permutation
        )
        word.append(("b", exponent_b))
    raise RuntimeError("triangle reduction did not terminate")


def cyclic_coset_matrices(
    cycle: np.ndarray, generator: np.ndarray
) -> list[np.ndarray]:
    """Return alpha_i with 1^(cycle^k)=i, indexed by i."""

    representatives: list[np.ndarray | None] = [None] * DEGREE
    for exponent in range(DEGREE):
        sheet = int(permutation_power(cycle, exponent)[0])
        representatives[sheet] = matrix_power(generator, exponent)
    if any(matrix is None for matrix in representatives):
        raise ArithmeticError("the chosen permutation did not generate all cosets")
    return [matrix for matrix in representatives if matrix is not None]


@lru_cache(maxsize=None)
def _cyclic_coset_exponents(cycle_key: tuple[int, ...]) -> tuple[int, ...]:
    cycle = np.asarray(cycle_key, dtype=np.int64)
    representatives: list[int | None] = [None] * DEGREE
    for exponent in range(DEGREE):
        sheet = int(permutation_power(cycle, exponent)[0])
        representatives[sheet] = exponent
    if any(exponent is None for exponent in representatives):
        raise ArithmeticError("the chosen permutation did not generate all cosets")
    return tuple(int(exponent) for exponent in representatives)


def cyclic_coset_exponents(cycle: np.ndarray) -> list[int]:
    """Return the power used for each sheet in ``cyclic_coset_matrices``."""

    return list(_cyclic_coset_exponents(tuple(int(value) for value in cycle)))


def side_pairing_generators(
    geometry: TriangleGeometry,
    sigma_a: np.ndarray,
    sigma_b: np.ndarray,
    cosets: Sequence[np.ndarray],
) -> list[np.ndarray]:
    """Build the nontrivial alpha_j*epsilon*alpha_i^-1 side pairings."""

    generators: list[np.ndarray] = []
    moves = [
        (geometry.delta_a, sigma_a),
        (matrix_power(geometry.delta_a, -1), permutation_power(sigma_a, -1)),
        (geometry.delta_b, sigma_b),
        (matrix_power(geometry.delta_b, -1), permutation_power(sigma_b, -1)),
    ]
    for source in range(DEGREE):
        for move_matrix, move_permutation in moves:
            target = int(move_permutation[source])
            pairing = cosets[source] @ move_matrix @ np.linalg.inv(cosets[target])
            determinant = float(np.linalg.det(pairing))
            pairing = pairing / math.sqrt(determinant)
            if np.linalg.norm(pairing - np.eye(2)) < 1e-11 or np.linalg.norm(
                pairing + np.eye(2)
            ) < 1e-11:
                continue
            if not any(
                min(np.linalg.norm(pairing - old), np.linalg.norm(pairing + old)) < 1e-10
                for old in generators
            ):
                generators.append(pairing)
    return generators


def improve_to_dirichlet_domain(
    z: complex,
    gamma: np.ndarray,
    generators: Sequence[np.ndarray],
    maximum_iterations: int = 1000,
) -> tuple[complex, np.ndarray, int]:
    """Greedily move to the side-pairing translate nearest the center i."""

    current_z = mobius(gamma, z)
    current_radius = abs(disc_coordinate(current_z))
    for iteration in range(maximum_iterations):
        best_radius = current_radius
        best_z = current_z
        best_generator: np.ndarray | None = None
        for generator in generators:
            candidate_z = mobius(generator, current_z)
            candidate_radius = abs(disc_coordinate(candidate_z))
            if candidate_radius < best_radius - 5e-13:
                best_radius = candidate_radius
                best_z = candidate_z
                best_generator = generator
        if best_generator is None:
            return current_z, gamma, iteration
        gamma = best_generator @ gamma
        determinant = float(np.linalg.det(gamma))
        gamma = gamma / math.sqrt(determinant)
        current_z = best_z
        current_radius = best_radius
    raise RuntimeError("Dirichlet improvement did not terminate")


@dataclass(frozen=True)
class SubgroupReduction:
    gamma: np.ndarray
    reduced_z: complex
    reduced_w: complex
    triangle_iterations: int
    coset: int
    coset_exponent: int
    triangle_word: tuple[tuple[str, int], ...]


def reduce_to_subgroup(
    z: complex,
    geometry: TriangleGeometry,
    sigma_a: np.ndarray,
    sigma_b: np.ndarray,
    coset_matrices: Sequence[np.ndarray],
    generators: Sequence[np.ndarray] | None = None,
) -> SubgroupReduction:
    triangle = reduce_to_triangle(z, geometry, sigma_a, sigma_b)
    inverse_delta = inverse_permutation(triangle.delta_permutation)
    coset = int(inverse_delta[0])
    coset_exponent = cyclic_coset_exponents(sigma_b)[coset]
    gamma = coset_matrices[coset] @ triangle.delta
    reduced_z = mobius(gamma, z)
    dirichlet_iterations = 0
    if generators:
        reduced_z, gamma, dirichlet_iterations = improve_to_dirichlet_domain(
            z, gamma, generators
        )
    reduced_w = disc_coordinate(reduced_z)
    return SubgroupReduction(
        gamma=gamma,
        reduced_z=reduced_z,
        reduced_w=reduced_w,
        triangle_iterations=triangle.iterations + dirichlet_iterations,
        coset=coset,
        coset_exponent=coset_exponent,
        triangle_word=triangle.word,
    )


def build_hejhal_matrix(
    class_id: int,
    terms: int,
    samples: int,
    rho: float,
    geometry: TriangleGeometry,
    weight: int = 2,
) -> tuple[np.ndarray, dict[str, float | int]]:
    sigma_a = X_REPRESENTATIVES[class_id - 1]
    sigma_b = Y
    cosets = cyclic_coset_matrices(sigma_b, geometry.delta_b)
    centers = [mobius(representative, 1j) for representative in cosets]
    angles = 2 * math.pi * np.arange(samples, dtype=np.float64) / samples
    source_w = rho * np.exp(1j * angles)
    block_size = terms + 1
    dimension = DEGREE * block_size
    matrix = np.zeros((dimension, dimension), dtype=np.complex128)
    powers = np.arange(block_size, dtype=np.int64)
    maximum_reduction_error = 0.0
    maximum_target_radius = 0.0
    maximum_iterations = 0

    for source_coset, center in enumerate(centers):
        row_slice = slice(source_coset * block_size, (source_coset + 1) * block_size)
        for w in source_w:
            z = upper_half_plane_coordinate(complex(w), center)
            reduction = reduce_to_subgroup(
                z, geometry, sigma_a, sigma_b, cosets
            )
            target_coset = reduction.coset
            target_center = centers[target_coset]
            target_w = disc_coordinate(reduction.reduced_z, target_center)
            c_entry = reduction.gamma[1, 0]
            d_entry = reduction.gamma[1, 1]
            automorphy = (c_entry * z + d_entry) ** -weight
            source_vector = automorphy / (
                (w / rho) ** powers * (1 - w) ** weight
            )
            target_vector = (
                (target_w / rho) ** powers * (1 - target_w) ** weight
            )
            column_slice = slice(
                target_coset * block_size, (target_coset + 1) * block_size
            )
            matrix[row_slice, column_slice] += (
                np.outer(source_vector, target_vector) / samples
            )
            maximum_target_radius = max(maximum_target_radius, abs(target_w))
            maximum_iterations = max(
                maximum_iterations, reduction.triangle_iterations
            )
            membership_error = abs(
                mobius(reduction.gamma, z)
                - upper_half_plane_coordinate(target_w, target_center)
            )
            maximum_reduction_error = max(
                maximum_reduction_error, membership_error
            )

    diagnostics: dict[str, float | int] = {
        "class_id": class_id,
        "terms": terms,
        "samples": samples,
        "rho": rho,
        "weight": weight,
        "patch_count": DEGREE,
        "matrix_dimension": dimension,
        "domain_radius": geometry.domain_radius,
        "maximum_target_radius": maximum_target_radius,
        "maximum_triangle_iterations": maximum_iterations,
        "maximum_reduction_roundtrip_error": maximum_reduction_error,
    }
    return matrix, diagnostics


def solve_weight_two(
    class_id: int,
    terms: int,
    samples: int,
    rho: float,
    geometry: TriangleGeometry,
    weight: int = 2,
) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    matrix, diagnostics = build_hejhal_matrix(
        class_id, terms, samples, rho, geometry, weight=weight
    )
    residual = matrix - np.eye(matrix.shape[0], dtype=np.complex128)
    _, singular_values, vh = np.linalg.svd(residual, full_matrices=False)
    order = np.argsort(singular_values)
    basis = vh.conj().T[:, order[:4]]
    ordered_singular_values = singular_values[order]
    diagnostics.update(
        {
            "smallest_singular_values": [
                float(value) for value in ordered_singular_values[:10]
            ],
            "genus_four_gap_ratio": float(
                ordered_singular_values[4] / ordered_singular_values[3]
            ),
            "basis_residual_norms": [
                float(np.linalg.norm(residual @ basis[:, column]))
                for column in range(4)
            ],
        }
    )
    return basis, ordered_singular_values, diagnostics


def order_23_branch_series(
    class_id: int,
    basis: np.ndarray,
    terms: int,
    samples: int,
    rho_main: float,
    rho_branch: float,
    geometry: TriangleGeometry,
) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    """Transfer the federalist basis to the HJLPPZ q-coordinate at b."""

    sigma_a = X_REPRESENTATIVES[class_id - 1]
    sigma_b = Y
    cosets = cyclic_coset_matrices(sigma_b, geometry.delta_b)
    centers = [mobius(representative, 1j) for representative in cosets]
    branch_center = geometry.mu * 1j
    theta = math.pi / 23
    # This sends c to the positive real axis and a to argument -pi/23,
    # exactly the disk normalization in HJLPPZ.
    q_rotation = cmath.exp(1j * (math.pi - theta))
    block_size = terms + 1
    dimension = DEGREE * block_size
    transfer = np.zeros((block_size, dimension), dtype=np.complex128)
    powers = np.arange(block_size, dtype=np.int64)
    maximum_target_radius = 0.0
    maximum_roundtrip_error = 0.0

    for sample in range(samples):
        q = rho_branch * cmath.exp(2j * math.pi * sample / samples)
        branch_w = q / q_rotation
        z = upper_half_plane_coordinate(branch_w, branch_center)
        reduction = reduce_to_subgroup(z, geometry, sigma_a, sigma_b, cosets)
        target_center = centers[reduction.coset]
        target_w = disc_coordinate(reduction.reduced_z, target_center)
        c_entry = reduction.gamma[1, 0]
        d_entry = reduction.gamma[1, 1]
        automorphy = (c_entry * z + d_entry) ** -2
        # The derivative dz/dq contributes a common constant and the factor
        # (1-branch_w)^-2.  The common constant disappears on echelonizing.
        source_vector = automorphy / (
            (q / rho_branch) ** powers * (1 - branch_w) ** 2
        )
        target_vector = (
            (target_w / rho_main) ** powers * (1 - target_w) ** 2
        )
        column_slice = slice(
            reduction.coset * block_size, (reduction.coset + 1) * block_size
        )
        transfer[:, column_slice] += (
            np.outer(source_vector, target_vector) / samples
        )
        maximum_target_radius = max(maximum_target_radius, abs(target_w))
        maximum_roundtrip_error = max(
            maximum_roundtrip_error,
            abs(
                mobius(reduction.gamma, z)
                - upper_half_plane_coordinate(target_w, target_center)
            ),
        )

    balanced_branch = transfer @ basis
    coefficients = balanced_branch / (
        rho_branch ** np.arange(block_size, dtype=np.float64)[:, None]
    )
    leading_matrix = coefficients[:4, :]
    condition = float(np.linalg.cond(leading_matrix))
    coordinate_change = np.linalg.inv(leading_matrix)
    normalized = coefficients @ coordinate_change
    diagnostics: dict[str, object] = {
        "branch_rho": rho_branch,
        "branch_samples": samples,
        "branch_target_maximum_radius": maximum_target_radius,
        "branch_roundtrip_error": maximum_roundtrip_error,
        "branch_leading_matrix_condition": condition,
        "branch_echelon_error": float(
            np.max(np.abs(normalized[:4, :] - np.eye(4)))
        ),
        "q_rotation": [q_rotation.real, q_rotation.imag],
    }
    return normalized, coordinate_change, diagnostics


def normalized_relation(
    vector: np.ndarray, tolerance: float = 1e-8, preferred_index: int | None = None
) -> np.ndarray:
    if preferred_index is not None and abs(vector[preferred_index]) > tolerance * np.max(
        np.abs(vector)
    ):
        return vector / vector[preferred_index]
    for coefficient in vector:
        if abs(coefficient) > tolerance * np.max(np.abs(vector)):
            return vector / coefficient
    raise ArithmeticError("relation is numerically zero")


def canonical_relations(series: np.ndarray) -> dict[str, object]:
    """Recover the canonical quadric and the Petri cubic from local series."""

    # The high Taylor coefficients are the first quantities contaminated by
    # the finite Hejhal cutoff, and their natural size grows rapidly.  The
    # canonical ideal is already determined by 20 jets, so use that stable
    # prefix and balance each coefficient equation to unit row norm.
    precision = min(series.shape[0], 20)
    quadratic_monomials = [(i, j) for i in range(4) for j in range(i, 4)]
    quadratic_series = np.column_stack(
        [
            np.convolve(series[:, i], series[:, j])[:precision]
            for i, j in quadratic_monomials
        ]
    )
    quadratic_balanced = quadratic_series / np.maximum(
        np.linalg.norm(quadratic_series, axis=1)[:, None], 1e-300
    )
    _, quadratic_singular_values, quadratic_vh = np.linalg.svd(
        quadratic_balanced, full_matrices=False
    )
    quadric = normalized_relation(quadratic_vh.conj().T[:, -1], preferred_index=2)

    cubic_monomials = [
        (i, j, k)
        for i in range(4)
        for j in range(i, 4)
        for k in range(j, 4)
    ]
    cubic_series = np.column_stack(
        [
            np.convolve(np.convolve(series[:, i], series[:, j]), series[:, k])[
                :precision
            ]
            for i, j, k in cubic_monomials
        ]
    )
    cubic_balanced = cubic_series / np.maximum(
        np.linalg.norm(cubic_series, axis=1)[:, None], 1e-300
    )
    _, cubic_singular_values, cubic_vh = np.linalg.svd(
        cubic_balanced, full_matrices=False
    )
    cubic_kernel = cubic_vh.conj().T[:, -5:]

    monomial_index = {monomial: index for index, monomial in enumerate(cubic_monomials)}
    quadric_multiples = np.zeros((len(cubic_monomials), 4), dtype=np.complex128)
    for variable in range(4):
        for coefficient, pair in zip(quadric, quadratic_monomials):
            triple = tuple(sorted((*pair, variable)))
            quadric_multiples[monomial_index[triple], variable] += coefficient
    q_orthogonal, _ = np.linalg.qr(quadric_multiples)
    overlap = q_orthogonal.conj().T @ cubic_kernel
    _, _, overlap_vh = np.linalg.svd(overlap, full_matrices=True)
    cubic = cubic_kernel @ overlap_vh.conj().T[:, -1]
    cubic -= q_orthogonal @ (q_orthogonal.conj().T @ cubic)
    cubic = normalized_relation(cubic, preferred_index=2)

    def render(values: np.ndarray) -> list[list[float]]:
        return [[float(value.real), float(value.imag)] for value in values]

    return {
        "quadratic_monomials": [list(monomial) for monomial in quadratic_monomials],
        "quadric_coefficients": render(quadric),
        "canonical_jet_precision": precision,
        "quadric_series_residual": float(
            np.linalg.norm(quadratic_balanced @ quadric)
        ),
        "quadric_smallest_singular_values": [
            float(value) for value in quadratic_singular_values[-3:]
        ],
        "cubic_monomials": [list(monomial) for monomial in cubic_monomials],
        "petri_cubic_coefficients": render(cubic),
        "petri_cubic_series_residual": float(
            np.linalg.norm(cubic_balanced @ cubic)
        ),
        "cubic_smallest_singular_values": [
            float(value) for value in cubic_singular_values[-7:]
        ],
        "cubic_quadric_multiple_overlap": float(
            np.linalg.norm(q_orthogonal.conj().T @ cubic)
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--class-id", type=int, choices=range(1, 8), default=6)
    parser.add_argument("--all-classes", action="store_true")
    parser.add_argument("--terms", type=int, default=30)
    parser.add_argument("--samples", type=int, default=92)
    parser.add_argument("--rho", type=float)
    parser.add_argument("--canonical", action="store_true")
    parser.add_argument("--branch-rho", type=float, default=0.70)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--geometry-only", action="store_true")
    arguments = parser.parse_args()

    geometry = triangle_geometry()
    rho = arguments.rho or min(0.95, geometry.domain_radius + 0.025)
    header = {
        "signature": [A_ORDER, B_ORDER, C_ORDER],
        "mu": geometry.mu,
        "domain_radius": geometry.domain_radius,
        "rho": rho,
    }
    if arguments.geometry_only:
        print(json.dumps(header, indent=2, sort_keys=True))
        return
    if arguments.canonical and arguments.terms < 19:
        raise ValueError("canonical extraction requires at least 20 Taylor coefficients")
    if rho <= geometry.domain_radius:
        raise ValueError("rho must strictly contain each triangular patch")
    class_results: list[dict[str, object]] = []
    class_ids = range(1, 8) if arguments.all_classes else [arguments.class_id]
    for class_id in class_ids:
        basis, _, diagnostics = solve_weight_two(
            class_id, arguments.terms, arguments.samples, rho, geometry
        )
        class_result: dict[str, object] = dict(diagnostics)
        if arguments.canonical:
            branch_series, _, branch_diagnostics = order_23_branch_series(
                class_id,
                basis,
                arguments.terms,
                arguments.samples,
                rho,
                arguments.branch_rho,
                geometry,
            )
            class_result.update(branch_diagnostics)
            class_result["canonical_relations"] = canonical_relations(branch_series)
        class_results.append(class_result)
    result: dict[str, object]
    if arguments.all_classes:
        result = {**header, "classes": class_results}
    else:
        result = {**header, **class_results[0]}
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    print(rendered, end="")
    if arguments.output:
        arguments.output.write_text(rendered)


if __name__ == "__main__":
    main()
