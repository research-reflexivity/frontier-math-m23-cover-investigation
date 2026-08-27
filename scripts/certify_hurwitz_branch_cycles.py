#!/usr/bin/env sage-python
"""Certify branch cycles for one embedding of the sextic Hurwitz component."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

from sage.all import ComplexBallField, ComplexField, PolynomialRing, QQ, RealBallField


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import certify_degree_one_branch_cycles as continuation  # noqa: E402
from reconstruct_hurwitz_degree23_maps import build_fields  # noqa: E402


ELIMINANT_PATH = ROOT / "data" / "hurwitz_monodromy_eliminant_candidate.json"
ALGEBRA_PATH = ROOT / "data" / "hurwitz_algebra_candidate.json"
CLASS_IDS = [7, 4, 1, 5, 3, 2]


def rational(record):
    return QQ(record["numerator"]) / QQ(record["denominator"])


class SexticCover(continuation.DegreeOneCover):
    def __init__(self, precision, class_id):
        self.precision = precision
        self.real_field = RealBallField(precision)
        self.field = ComplexBallField(precision)
        self.ring = PolynomialRing(self.field, "x")
        payload = json.loads(ELIMINANT_PATH.read_text())
        algebra = json.loads(ALGEBRA_PATH.read_text())
        fields = build_fields(algebra)
        exact_field = fields["absolute"]
        approximate_field = ComplexField(precision)
        approximate_embeddings = exact_field.embeddings(approximate_field)
        ball_ring = PolynomialRing(self.field, "a_ball")
        ball_roots = ball_ring(exact_field.polynomial()).roots(multiplicities=False)
        unused_roots = list(ball_roots)
        embeddings = []
        for approximate_embedding in approximate_embeddings:
            approximate_root = approximate_embedding(exact_field.gen())
            root = min(
                unused_roots,
                key=lambda value: abs(complex(value.mid()) - complex(approximate_root)),
            )
            unused_roots.remove(root)
            embeddings.append(exact_field.hom([root], self.field, check=False))
        embeddings = [
            embedding
            for embedding in embeddings
            if embedding(fields["sqrt_minus_23_absolute"]).imag() > 0
        ]
        if len(embeddings) != 6:
            raise AssertionError("expected six embeddings over the chosen K0 place")
        embedding = embeddings[CLASS_IDS.index(class_id)]

        def decode(records):
            coefficients = [rational(value) for value in records]
            return exact_field(coefficients)

        exact_coefficients = [
            [decode(payload["coefficients_target_then_x"][i][j]) for j in range(24)]
            for i in range(7)
        ]

        embedded_coefficients = [
            [embedding(exact_coefficients[i][j]) for j in range(24)]
            for i in range(7)
        ]
        base = continuation.qcomplex(QQ(1) / 2, 2)
        preliminary_roots = None
        preliminary_exponent = None
        for scale_exponent in [0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6]:
            scale_rational = (
                QQ(2) ** scale_exponent
                if scale_exponent >= 0
                else QQ(1) / (QQ(2) ** (-scale_exponent))
            )
            self.root_scale = self.field(scale_rational)
            self.beta_coefficients = [
                [
                    embedded_coefficients[i][j] / self.root_scale**j
                    for j in range(24)
                ]
                for i in range(7)
            ]
            self.beta_polynomials = [
                self.ring(coefficients) for coefficients in self.beta_coefficients
            ]
            try:
                preliminary_roots = self.roots(self.point(base))
                preliminary_exponent = scale_exponent
                break
            except ValueError:
                continue
        if preliminary_roots is None:
            raise ValueError("could not condition the base-fibre root isolation")
        median = sorted(float(root.abs().upper()) for root in preliminary_roots)[11]
        scale_exponent = preliminary_exponent + int(round(-math.log2(median)))
        scale_rational = (
            QQ(2) ** scale_exponent
            if scale_exponent >= 0
            else QQ(1) / (QQ(2) ** (-scale_exponent))
        )
        self.root_scale = self.field(scale_rational)
        self.beta_coefficients = [
            [
                embedding(exact_coefficients[i][j]) / self.root_scale**j
                for j in range(24)
            ]
            for i in range(7)
        ]
        self.beta_polynomials = [
            self.ring(coefficients) for coefficients in self.beta_coefficients
        ]
        self.class_id = class_id
        self.embedding_index = CLASS_IDS.index(class_id) + 1
        self.scale_rational = scale_rational

    def polynomial(self, beta):
        polynomial = self.beta_polynomials[6]
        for beta_degree in range(5, -1, -1):
            polynomial = polynomial * beta + self.beta_polynomials[beta_degree]
        if polynomial.degree() != 23 or polynomial[23].contains_zero():
            raise ValueError("degree-23 leading coefficient not certified")
        return polynomial

    def centered_coefficients(self, beta, root):
        beta_shifted = []
        for beta_order in range(7):
            value = self.ring.zero()
            for beta_degree in range(beta_order, 7):
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
        for beta_order in range(6, -1, -1):
            value = value * delta_beta + coefficients[beta_order][0]
        return value

    def derivative_from_centered_coefficients(
        self, delta_beta, delta_root, coefficients
    ):
        value = self.field.zero()
        for beta_order in range(6, -1, -1):
            root_series = self.field.zero()
            for root_order in range(22, -1, -1):
                root_series = (
                    root_series * delta_root
                    + (root_order + 1) * coefficients[beta_order][root_order + 1]
                )
            value = value * delta_beta + root_series
        return value


def loops(radius=QQ(2) / 5, infinity_detour=False):
    base = continuation.qcomplex(QQ(1) / 2, 2)
    top_zero = continuation.qcomplex(0, radius)
    top_one = continuation.qcomplex(1, radius)
    result = {
        "zero": [
            base,
            top_zero,
            continuation.qcomplex(-radius, 0),
            continuation.qcomplex(0, -radius),
            continuation.qcomplex(radius, 0),
            top_zero,
            base,
        ],
        "one": [
            base,
            top_one,
            continuation.qcomplex(1 - radius, 0),
            continuation.qcomplex(1, -radius),
            continuation.qcomplex(1 + radius, 0),
            top_one,
            base,
        ],
        "infinity": [
            base,
            continuation.qcomplex(QQ(5) / 2, 0),
            continuation.qcomplex(QQ(1) / 2, -2),
            continuation.qcomplex(QQ(-3) / 2, 0),
            base,
        ],
    }
    if infinity_detour:
        result["infinity"] = [
            base,
            continuation.qcomplex(QQ(5) / 2, 0),
            continuation.qcomplex(QQ(1) / 2, -2),
            continuation.qcomplex(QQ(-3) / 2, 0),
            continuation.qcomplex(QQ(-3) / 2, 2),
            base,
        ]
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--class-id", type=int, choices=CLASS_IDS, required=True)
    parser.add_argument("--precision", type=int, default=256)
    parser.add_argument("--loop", choices=("zero", "one", "infinity"), action="append")
    parser.add_argument("--radius", default="2/5")
    parser.add_argument("--infinity-detour", action="store_true")
    arguments = parser.parse_args()
    cover = SexticCover(arguments.precision, arguments.class_id)
    print(
        "class",
        cover.class_id,
        "embedding",
        cover.embedding_index,
        "root_scale",
        cover.scale_rational,
        flush=True,
    )
    selected_loops = loops(QQ(arguments.radius), arguments.infinity_detour)
    if arguments.loop:
        selected_loops = {
            label: selected_loops[label] for label in arguments.loop
        }
    for label, vertices in selected_loops.items():
        permutation, statistics = continuation.continue_path(cover, vertices)
        print(label, permutation)
        print(
            label,
            "cycle_lengths",
            continuation.cycle_lengths(permutation),
            statistics,
            flush=True,
        )


if __name__ == "__main__":
    main()
