#!/usr/bin/env sage-python
"""Certify branch-cycle continuation for the degree-one, class-6 map.

Each accepted target segment is covered by 23 pairwise disjoint root disks.
A fixed-slope complex interval-Newton estimate proves a contraction on every
disk uniformly for the whole segment.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from sage.all import ComplexBallField, PolynomialRing, QQ, RealBallField


ROOT = Path(__file__).resolve().parents[1]
MODEL_PATH = ROOT / "data" / "optimal_23_4_Z.json"


def qcomplex(real, imag=0):
    return (QQ(real), QQ(imag))


def midpoint(left, right):
    return ((left[0] + right[0]) / 2, (left[1] + right[1]) / 2)


class DegreeOneCover:
    def __init__(self, precision):
        self.precision = precision
        self.real_field = RealBallField(precision)
        self.field = ComplexBallField(precision)
        self.ring = PolynomialRing(self.field, "w")
        self.table = json.loads(MODEL_PATH.read_text())["coefficients_T_then_W"]
        self.sqrt_minus_23 = self.field(23).sqrt() * self.field.gen(0)
        # The integral optimal coordinate has roots of size about 1/32 on
        # the chosen base fibre.  Track z=32*w to improve root conditioning.
        self.root_scale = self.field(32)
        self.beta_coefficients = [
            [self.field.zero() for _ in range(24)] for _ in range(5)
        ]
        for t_degree in range(5):
            for beta_degree in range(5):
                binomial_coefficient = sum(
                    (-1) ** left
                    * math.comb(t_degree, left)
                    * math.comb(4 - t_degree, beta_degree - left)
                    for left in range(
                        max(0, beta_degree - (4 - t_degree)),
                        min(t_degree, beta_degree) + 1,
                    )
                )
                for w_degree in range(24):
                    self.beta_coefficients[beta_degree][w_degree] += (
                        self.field(self.table[t_degree][w_degree])
                        * self.sqrt_minus_23**t_degree
                        * binomial_coefficient
                        / self.root_scale**w_degree
                    )
        self.beta_polynomials = [
            self.ring(coefficients) for coefficients in self.beta_coefficients
        ]

    def point(self, value):
        return self.field(self.real_field(value[0]), self.real_field(value[1]))

    def segment_ball(self, left, right):
        real = self.real_field(left[0]).union(self.real_field(right[0]))
        imag = self.real_field(left[1]).union(self.real_field(right[1]))
        return self.field(real, imag)

    def polynomial(self, beta):
        polynomial = self.beta_polynomials[4]
        for beta_degree in range(3, -1, -1):
            polynomial = polynomial * beta + self.beta_polynomials[beta_degree]
        if polynomial.degree() != 23 or polynomial[23].contains_zero():
            raise ValueError("degree-23 leading coefficient not certified")
        return polynomial

    def roots(self, beta):
        polynomial = self.polynomial(beta)
        roots = (polynomial / polynomial[23]).roots(multiplicities=False)
        if len(roots) != 23:
            raise ValueError(f"isolated {len(roots)} roots, expected 23")
        if any(
            roots[first].overlaps(roots[second])
            for first in range(23)
            for second in range(first)
        ):
            raise ValueError("root balls are not pairwise disjoint")
        return roots

    def centered_coefficients(self, beta, root):
        """Coefficients of P(beta+db, root+dw) in db and dw."""

        beta_shifted = []
        for beta_order in range(5):
            value = self.ring.zero()
            for beta_degree in range(beta_order, 5):
                value += (
                    self.beta_polynomials[beta_degree]
                    * math.comb(beta_degree, beta_order)
                    * beta ** (beta_degree - beta_order)
                )
            beta_shifted.append(value)
        variable = self.ring.gen()
        result = []
        for polynomial in beta_shifted:
            coefficients = list(polynomial(variable + root))
            result.append(coefficients + [self.field.zero()] * (24 - len(coefficients)))
        return result

    def value_from_centered_coefficients(self, delta_beta, coefficients):
        value = self.field.zero()
        for beta_order in range(4, -1, -1):
            value = value * delta_beta + coefficients[beta_order][0]
        return value

    def derivative_from_centered_coefficients(
        self, delta_beta, delta_root, coefficients
    ):
        value = self.field.zero()
        for beta_order in range(4, -1, -1):
            root_series = self.field.zero()
            for root_order in range(22, -1, -1):
                root_series = (
                    root_series * delta_root
                    + (root_order + 1) * coefficients[beta_order][root_order + 1]
                )
            value = value * delta_beta + root_series
        return value


def unique_overlap_map(left, right):
    matches = []
    for value in left:
        candidates = [index for index, other in enumerate(right) if value.overlaps(other)]
        if len(candidates) != 1:
            return None
        matches.append(candidates[0])
    if len(set(matches)) != len(matches):
        return None
    return matches


def ball_center(field, value):
    return field(value.real().mid(), value.imag().mid())


def contained_in_disk(value, center, radius):
    return (value - center).abs().upper() < radius


def certified_tubes(cover, left, right):
    middle = midpoint(left, right)
    middle_value = cover.point(middle)
    middle_roots = cover.roots(middle_value)
    centers = [ball_center(cover.field, value) for value in middle_roots]
    separations = []
    for index, center in enumerate(centers):
        separation = min(
            (center - other).abs().lower()
            for other_index, other in enumerate(centers)
            if other_index != index
        )
        separations.append(separation)

    parameter_ball = cover.segment_ball(left, right)
    delta_parameter = parameter_ball - middle_value
    tubes = []
    for root_index, (center, separation) in enumerate(zip(centers, separations)):
        centered = cover.centered_coefficients(middle_value, center)
        slope = cover.field.one() / centered[0][1]
        eta = (
            slope
            * cover.value_from_centered_coefficients(
                delta_parameter, centered
            )
        ).abs().upper()
        accepted = None
        for exponent in range(2, 31):
            radius = separation / (2**exponent)
            real_box = center.real().add_error(radius)
            imag_box = center.imag().add_error(radius)
            root_box = cover.field(real_box, imag_box)
            derivative_family = cover.derivative_from_centered_coefficients(
                delta_parameter, root_box - center, centered
            )
            contraction = (
                cover.field.one() - slope * derivative_family
            ).abs().upper()
            if (
                contraction < 1
                and eta + contraction * radius < radius
                and contained_in_disk(middle_roots[root_index], center, radius)
            ):
                accepted = (center, radius)
                break
        if accepted is None:
            raise ValueError("no fixed-slope interval Newton radius succeeded")
        tubes.append(accepted)

    for first, (center, radius) in enumerate(tubes):
        for other_center, other_radius in tubes[:first]:
            if not radius + other_radius < (center - other_center).abs().lower():
                raise ValueError("Newton disks are not pairwise disjoint")
    return tubes


def disk_containment_map(roots, tubes):
    matches = []
    for value in roots:
        candidates = [
            index
            for index, (center, radius) in enumerate(tubes)
            if contained_in_disk(value, center, radius)
        ]
        if len(candidates) != 1:
            return None
        matches.append(candidates[0])
    if len(set(matches)) != len(matches):
        return None
    return matches


def continue_segment(cover, left, right, left_roots, depth, statistics):
    statistics["attempted_tubes"] += 1
    try:
        tubes = certified_tubes(cover, left, right)
        right_roots = cover.roots(cover.point(right))
        left_to_tube = disk_containment_map(left_roots, tubes)
        right_to_tube = disk_containment_map(right_roots, tubes)
        if left_to_tube is not None and right_to_tube is not None:
            inverse_left_to_tube = {
                tube: left_index for left_index, tube in enumerate(left_to_tube)
            }
            ordered = [None] * 23
            for right_index, tube_index in enumerate(right_to_tube):
                ordered[inverse_left_to_tube[tube_index]] = right_roots[right_index]
            statistics["certified_tubes"] += 1
            statistics["maximum_depth"] = max(statistics["maximum_depth"], depth)
            return ordered
    except (ArithmeticError, ValueError):
        pass
    if depth >= 24:
        raise RuntimeError(f"failed to certify segment {left} -> {right}")
    middle = midpoint(left, right)
    middle_roots = continue_segment(
        cover, left, middle, left_roots, depth + 1, statistics
    )
    return continue_segment(
        cover, middle, right, middle_roots, depth + 1, statistics
    )


def continue_path(cover, vertices):
    initial = cover.roots(cover.point(vertices[0]))
    current = initial
    statistics = {"attempted_tubes": 0, "certified_tubes": 0, "maximum_depth": 0}
    for segment_index, (left, right) in enumerate(zip(vertices, vertices[1:]), 1):
        current = continue_segment(cover, left, right, current, 0, statistics)
        print(
            "path segment",
            segment_index,
            "of",
            len(vertices) - 1,
            statistics,
            flush=True,
        )
    permutation = unique_overlap_map(current, initial)
    if permutation is None:
        raise RuntimeError("closed path did not return to the isolated base fibre")
    return [value + 1 for value in permutation], statistics


def cycle_lengths(permutation):
    seen = set()
    lengths = []
    for start in range(1, len(permutation) + 1):
        if start in seen:
            continue
        length = 0
        value = start
        while value not in seen:
            seen.add(value)
            length += 1
            value = permutation[value - 1]
        lengths.append(length)
    return sorted(lengths)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--precision", type=int, default=256)
    arguments = parser.parse_args()
    cover = DegreeOneCover(arguments.precision)
    base = qcomplex(QQ(1) / 2, 2)
    top_zero = qcomplex(0, QQ(1) / 4)
    top_one = qcomplex(1, QQ(1) / 4)
    loops = {
        "zero": [
            base,
            top_zero,
            qcomplex(QQ(-1) / 4, 0),
            qcomplex(0, QQ(-1) / 4),
            qcomplex(QQ(1) / 4, 0),
            top_zero,
            base,
        ],
        "one": [
            base,
            top_one,
            qcomplex(QQ(3) / 4, 0),
            qcomplex(1, QQ(-1) / 4),
            qcomplex(QQ(5) / 4, 0),
            top_one,
            base,
        ],
        "infinity": [
            base,
            qcomplex(QQ(5) / 2, 0),
            qcomplex(QQ(1) / 2, -2),
            qcomplex(QQ(-3) / 2, 0),
            base,
        ],
    }
    for label, vertices in loops.items():
        permutation, statistics = continue_path(cover, vertices)
        print(label, permutation)
        print(label, "cycle_lengths", cycle_lengths(permutation), statistics)


if __name__ == "__main__":
    main()
