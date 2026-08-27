#!/usr/bin/env sage-python
"""Recover the exact degree-23 pencil on the reconstructed canonical curves.

For a marked pair ``b,c`` on a genus-four canonical complete intersection,
the two four-dimensional spaces

    H^0(5K-23b),  H^0(5K-23c)

are related by multiplication by the desired rational function.  The unique
multiplier is found by linear algebra in the degree-ten canonical ring.  The
degree-12 component encodes the six previously missing conjugate covers.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from sage.all import PolynomialRing, PowerSeriesRing, QQ, matrix, prod, vector


ROOT = Path(__file__).resolve().parents[1]
ALGEBRA_PATH = ROOT / "data" / "hurwitz_algebra_candidate.json"
CANONICAL_PATH = ROOT / "data" / "hurwitz_canonical_models_candidate.json"
MARKED_PATH = ROOT / "data" / "hurwitz_marked_points_candidate.json"
BRANCH_PATH = ROOT / "data" / "hurwitz_degree23_branch_candidate.json"
DEFAULT_OUTPUT = ROOT / "data" / "hurwitz_degree23_maps_candidate.json"


def rational(record):
    return QQ(record["numerator"]) / QQ(record["denominator"])


def encode_rational(value):
    return {"numerator": int(value.numerator()), "denominator": int(value.denominator())}


def exponent_tuples(total, slots=4):
    if slots == 1:
        return [(total,)]
    return [
        (first, *rest)
        for first in range(total + 1)
        for rest in exponent_tuples(total - first, slots - 1)
    ]


def build_fields(algebra):
    base_ring = PolynomialRing(QQ, "s")
    s = base_ring.gen()
    base_field = QQ.extension(s**2 + 23, "sqrt_minus_23")
    sqrt_minus_23 = base_field.gen()

    def decode_quadratic(record):
        return base_field(rational(record["rational_part"])) + base_field(
            rational(record["sqrt_minus_23_part"])
        ) * sqrt_minus_23

    hurwitz_ring = PolynomialRing(base_field, "J")
    factors = [
        hurwitz_ring([
            decode_quadratic(value)
            for value in factor["coefficients_ascending"]
        ])
        for factor in algebra["exact_factors"]
    ]
    sextic = next(factor for factor in factors if factor.degree() == 6)
    relative_field = base_field.extension(sextic, "j")
    absolute_field = relative_field.absolute_field("a")
    _, to_absolute = absolute_field.structure()
    sqrt_minus_23_absolute = to_absolute(relative_field(sqrt_minus_23))
    return {
        "base": base_field,
        "sqrt_minus_23": sqrt_minus_23,
        "absolute": absolute_field,
        "sqrt_minus_23_absolute": sqrt_minus_23_absolute,
        "integral_basis": absolute_field.integral_basis(),
    }


def reconstruct(component, algebra, canonical, marked, branch, check_basepoints=False):
    fields = build_fields(algebra)
    base_field = fields["base"]
    sqrt_minus_23 = fields["sqrt_minus_23"]
    absolute_field = fields["absolute"]
    integral_basis = fields["integral_basis"]
    field = base_field if component == "degree_one" else absolute_field

    def degree_one_curve(record):
        value = record["degree_one_component"]
        return base_field(rational(value["rational_part"])) + base_field(
            rational(value["sqrt_minus_23_part"])
        ) * sqrt_minus_23

    def degree_one_point(record):
        return base_field(rational(record["rational_part"])) + base_field(
            rational(record["sqrt_minus_23_part"])
        ) * sqrt_minus_23

    def sextic_curve(record):
        return sum(
            rational(value) * basis
            for value, basis in zip(
                record["sextic_component"]["integral_basis_coordinates"],
                integral_basis,
            )
        )

    def sextic_point(record):
        return sum(
            rational(value) * basis
            for value, basis in zip(record["integral_basis_coordinates"], integral_basis)
        )

    quadric_monomials = [tuple(value) for value in canonical["quadric"]["monomials"]]
    quadric_coefficients = [
        field(degree_one_curve(value)) if component == "degree_one" else sextic_curve(value)
        for value in canonical["quadric"]["coefficients"]
    ]
    cubic_monomials = [tuple(value) for value in canonical["petri_cubic"]["monomials"]]
    cubic_coefficients = [
        field(degree_one_curve(value)) if component == "degree_one" else sextic_curve(value)
        for value in canonical["petri_cubic"]["coefficients"]
    ]
    b_point = [field(1), field(0), field(0), field(0)]
    c_point = (
        [
            field(degree_one_point(value))
            for value in marked["c_coordinates"]["degree_one_component"]
        ]
        if component == "degree_one"
        else [
            sextic_point(value)
            for value in marked["c_coordinates"]["sextic_component"]
        ]
    )

    def evaluate_form(point, monomials, coefficients):
        return sum(
            coefficient * prod(point[index] for index in monomial)
            for monomial, coefficient in zip(monomials, coefficients)
        )

    def derivative(coefficients, monomials, variable, point):
        total = field(0)
        for monomial, coefficient in zip(monomials, coefficients):
            multiplicity = monomial.count(variable)
            if multiplicity:
                factors = list(monomial)
                factors.remove(variable)
                total += coefficient * multiplicity * prod(
                    point[index] for index in factors
                )
        return total

    def local_series(point, precision=31):
        for parameter_index in (1, 2, 3):
            solve_indices = [
                index for index in (1, 2, 3) if index != parameter_index
            ]
            jacobian = matrix(field, [
                [
                    derivative(
                        quadric_coefficients,
                        quadric_monomials,
                        index,
                        point,
                    )
                    for index in solve_indices
                ],
                [
                    derivative(
                        cubic_coefficients,
                        cubic_monomials,
                        index,
                        point,
                    )
                    for index in solve_indices
                ],
            ])
            if not jacobian.det():
                continue
            series_ring = PowerSeriesRing(field, "t", default_prec=precision)
            t = series_ring.gen()
            values = [series_ring(value).add_bigoh(precision) for value in point]
            values[parameter_index] += t
            for order in range(1, precision):
                residual = vector(field, [
                    evaluate_form(values, quadric_monomials, quadric_coefficients)[order],
                    evaluate_form(values, cubic_monomials, cubic_coefficients)[order],
                ])
                correction = jacobian.solve_right(-residual)
                for index, value in zip(solve_indices, correction):
                    values[index] += series_ring(value) * t**order
            if (
                evaluate_form(values, quadric_monomials, quadric_coefficients) == 0
                and evaluate_form(values, cubic_monomials, cubic_coefficients) == 0
            ):
                return values
        raise ValueError("no nonsingular local coordinate choice")

    def relation_matrix(degree):
        monomials = exponent_tuples(degree)
        index = {monomial: position for position, monomial in enumerate(monomials)}
        rows = []
        for tail in exponent_tuples(degree - 2):
            row = [field(0)] * len(monomials)
            for (left, right), coefficient in zip(
                quadric_monomials, quadric_coefficients
            ):
                exponent = list(tail)
                exponent[left] += 1
                exponent[right] += 1
                row[index[tuple(exponent)]] += coefficient
            rows.append(row)
        for tail in exponent_tuples(degree - 3):
            row = [field(0)] * len(monomials)
            for (first, second, third), coefficient in zip(
                cubic_monomials, cubic_coefficients
            ):
                exponent = list(tail)
                exponent[first] += 1
                exponent[second] += 1
                exponent[third] += 1
                row[index[tuple(exponent)]] += coefficient
            rows.append(row)
        return monomials, index, matrix(field, rows)

    def quotient_basis(degree):
        monomials, index, relations = relation_matrix(degree)
        echelon = relations.echelon_form()
        pivots = echelon.pivots()
        free = [position for position in range(len(monomials)) if position not in pivots]
        if len(free) != 6 * degree - 3:
            raise AssertionError("canonical-ring Hilbert function changed")
        return monomials, index, echelon, pivots, free

    monomials_5, _, relations_5, pivots_5, free_5 = quotient_basis(5)

    def section_space(point):
        values = local_series(point)
        series = []
        for position in free_5:
            value = values[0].parent()(1)
            for coordinate, exponent in zip(values, monomials_5[position]):
                value *= coordinate**exponent
            series.append(value)
        jets = matrix(
            field,
            23,
            len(series),
            lambda row, column: series[column][row],
        )
        kernel = jets.right_kernel()
        if jets.rank() != 23 or kernel.dimension() != 4:
            raise AssertionError("unexpected 23-jet rank")
        sections = []
        for coefficients in kernel.basis():
            section = vector(field, len(monomials_5))
            for position, value in zip(free_5, coefficients):
                section[position] = value
            sections.append(section)
        return sections, values

    b_sections, b_series = section_space(b_point)
    c_sections, c_series = section_space(c_point)
    monomials_10, index_10, relations_10, pivots_10, free_10 = quotient_basis(10)

    def multiply_sections(left, right):
        result = vector(field, len(monomials_10))
        for left_index, left_value in enumerate(left):
            if not left_value:
                continue
            for right_index, right_value in enumerate(right):
                if right_value:
                    exponent = tuple(
                        a + b
                        for a, b in zip(
                            monomials_5[left_index], monomials_5[right_index]
                        )
                    )
                    result[index_10[exponent]] += left_value * right_value
        return result

    def quotient_coordinates(value):
        result = vector(field, value)
        for row, pivot in zip(relations_10.rows(), pivots_10):
            if result[pivot]:
                result -= result[pivot] * row
        if any(result[pivot] for pivot in pivots_10):
            raise AssertionError("degree-ten reduction failed")
        return vector(field, [result[position] for position in free_10])

    products = [
        [
            quotient_coordinates(multiply_sections(b_sections[row], c_sections[column]))
            for column in range(4)
        ]
        for row in range(4)
    ]
    equations = []
    for left in range(4):
        for right in range(left + 1, 4):
            for coordinate in range(len(free_10)):
                equation = [field(0)] * 16
                for column in range(4):
                    equation[4 * right + column] += products[left][column][coordinate]
                    equation[4 * left + column] -= products[right][column][coordinate]
                equations.append(equation)
    system = matrix(field, equations)
    kernel = system.right_kernel()
    if system.rank() != 15 or kernel.dimension() != 1:
        raise AssertionError("multiplier is not unique")
    multiplier = matrix(field, 4, 4, kernel.basis()[0])
    if not multiplier.det():
        raise AssertionError("multiplier is singular")
    denominator_sections = [
        sum(multiplier[row, column] * c_sections[column] for column in range(4))
        for row in range(4)
    ]
    for left in range(4):
        for right in range(left + 1, 4):
            relation = multiply_sections(
                b_sections[left], denominator_sections[right]
            ) - multiply_sections(b_sections[right], denominator_sections[left])
            if not quotient_coordinates(relation).is_zero():
                raise AssertionError("cross-multiplication identity failed")

    def evaluate_section(section, point):
        return sum(
            value
            * prod(point[index]**exponent for index, exponent in enumerate(monomial))
            for monomial, value in zip(monomials_5, section)
        )

    def section_series(section, values):
        total = values[0].parent()(0)
        for monomial, coefficient in zip(monomials_5, section):
            term = values[0].parent()(coefficient)
            for value, exponent in zip(values, monomial):
                term *= value**exponent
            total += term
        return total

    chosen_row = next(
        row
        for row in range(4)
        if evaluate_section(b_sections[row], c_point)
        and evaluate_section(denominator_sections[row], b_point)
    )
    numerator = b_sections[chosen_row]
    denominator = denominator_sections[chosen_row]
    if section_series(numerator, b_series).valuation() != 23:
        raise AssertionError("numerator does not have order 23 at b")
    if section_series(denominator, c_series).valuation() != 23:
        raise AssertionError("denominator does not have order 23 at c")

    basepoint_diagnostics = None
    if check_basepoints:
        projective_ring = PolynomialRing(
            field, 4, names=("X0", "X1", "X2", "X3"), order="degrevlex"
        )
        coordinates = projective_ring.gens()

        def projective_form(monomials, coefficients):
            return sum(
                projective_ring(value)
                * prod(
                    coordinates[index] ** exponent
                    for index, exponent in enumerate(monomial)
                )
                for monomial, value in zip(monomials, coefficients)
            )

        quadric_form = sum(
            projective_ring(value) * coordinates[left] * coordinates[right]
            for (left, right), value in zip(
                quadric_monomials, quadric_coefficients
            )
        )
        cubic_form = sum(
            projective_ring(value)
            * coordinates[first]
            * coordinates[second]
            * coordinates[third]
            for (first, second, third), value in zip(
                cubic_monomials, cubic_coefficients
            )
        )
        point_ideals = {
            "b": projective_ring.ideal(coordinates[1:]),
            "c": projective_ring.ideal([
                c_point[0] * coordinates[index]
                - c_point[index] * coordinates[0]
                for index in range(1, 4)
            ]),
        }
        section_bases = {
            "b": b_sections,
            "c": denominator_sections,
        }
        basepoint_diagnostics = {}
        for label in ("b", "c"):
            started = time.time()
            forms = [
                projective_form(monomials_5, section)
                for section in section_bases[label]
            ]
            base_ideal = projective_ring.ideal(quadric_form, cubic_form, *forms)
            residual = base_ideal.saturation(point_ideals[label])[0]
            if not residual.is_one():
                raise AssertionError(f"residual base point away from {label}")
            basepoint_diagnostics[label] = {
                "residual_ideal_is_unit": True,
                "seconds": time.time() - started,
            }

    first_nonzero = next(value for value in numerator if value)
    numerator = numerator / first_nonzero
    denominator = denominator / first_nonzero
    lambda_value = None
    if component == "sextic":
        lambda_value = field([
            QQ(value) / QQ(branch["lambda"]["denominator"])
            for value in branch["lambda"]["power_basis_numerators"]
        ])

    field_degree = int(field.degree())

    def encode_element(value):
        coefficients = list(field(value)) + [QQ(0)] * field_degree
        return [encode_rational(coefficient) for coefficient in coefficients[:field_degree]]

    def encode_section(section):
        return [
            {
                "monomial": list(monomial),
                "coefficient_power_basis": encode_element(value),
            }
            for monomial, value in zip(monomials_5, section)
            if value
        ]

    result = {
        "component": component,
        "field_degree": field_degree,
        "chosen_row": chosen_row,
        "degree_5_quotient_dimension": len(free_5),
        "degree_10_quotient_dimension": len(free_10),
        "jet_ranks": {"b": 23, "c": 23},
        "section_space_dimensions": {"b": 4, "c": 4},
        "multiplier_system_rank": int(system.rank()),
        "multiplier_kernel_dimension": int(kernel.dimension()),
        "multiplier_determinant_nonzero": True,
        "orders": {"numerator_at_b": 23, "denominator_at_c": 23},
        "numerator": encode_section(numerator),
        "denominator": encode_section(denominator),
    }
    if lambda_value is not None:
        result["third_branch_value_power_basis"] = encode_element(lambda_value)
        result["normalized_denominator"] = encode_section(lambda_value * denominator)
        result["normalized_map"] = "beta=numerator/normalized_denominator"
    if basepoint_diagnostics is not None:
        result["residual_basepoint_check"] = {
            label: {"residual_ideal_is_unit": True}
            for label in ("b", "c")
        }
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--component",
        choices=("degree_one", "sextic", "all"),
        default="sextic",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check-basepoints", action="store_true")
    arguments = parser.parse_args()
    algebra = json.loads(ALGEBRA_PATH.read_text())
    canonical = json.loads(CANONICAL_PATH.read_text())
    marked = json.loads(MARKED_PATH.read_text())
    branch = json.loads(BRANCH_PATH.read_text())
    components = (
        ("degree_one", "sextic")
        if arguments.component == "all"
        else (arguments.component,)
    )
    results = [
        reconstruct(
            component,
            algebra,
            canonical,
            marked,
            branch,
            check_basepoints=arguments.check_basepoints,
        )
        for component in components
    ]
    payload = {
        "schema": "m23.cover-investigation.degree23-maps-candidate.v1",
        "status": "exact_map_sections_and_exact_third_branch_normalization",
        "construction": "H0(5K-23b) to H0(5K-23c) multiplier in the canonical ring",
        "components": results,
        "third_fiber_certificate": {
            "status": "verified_exact",
            "certificate": "verification/verify_hurwitz_degree23_third_fiber.sage",
            "extra_gcd_degree": 8,
            "extra_gcd_squarefree": True,
            "magma_collision_free_residue_embeddings": [2, 3, 4, 5, 9, 12],
            "magma_geometry_certificate": "verification/verify_hurwitz_degree23_geometry.m",
            "projected_point_collisions_excluded": True,
        },
    }
    arguments.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    for result in results:
        print(
            "PASS",
            result["component"],
            "rank-15 multiplier and exact order-(23,23) map sections",
        )
    print(f"wrote {arguments.output}")


if __name__ == "__main__":
    main()
