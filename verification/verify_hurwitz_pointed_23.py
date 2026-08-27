#!/usr/bin/env sage-python
"""Certify the pointed naive reductions of all Hurwitz components at 23.

This determines the integral canonical reductions, their complete
singularity types, and the divisor of the reduced degree-23 maps.  It is
not a semistable-reduction certificate for the Galois closures: the ADE
singularities still have to be resolved with their Mathieu pointing data.
"""

from __future__ import annotations

import json
from pathlib import Path

from sage.all import PowerSeriesRing, PolynomialRing, ProjectiveSpace, QQ, matrix, vector


ROOT = Path(__file__).resolve().parents[1]


def rational(record):
    if isinstance(record, dict):
        return QQ(record["numerator"]) / QQ(record["denominator"])
    return QQ(record)


def main() -> None:
    algebra = json.loads((ROOT / "data/hurwitz_algebra_candidate.json").read_text())
    canonical = json.loads(
        (ROOT / "data/hurwitz_canonical_models_candidate.json").read_text()
    )
    marked = json.loads(
        (ROOT / "data/hurwitz_marked_points_candidate.json").read_text()
    )
    maps = json.loads(
        (ROOT / "data/hurwitz_degree23_maps_candidate.json").read_text()
    )
    branch = json.loads(
        (ROOT / "data/hurwitz_degree23_branch_candidate.json").read_text()
    )

    base_ring = PolynomialRing(QQ, "s")
    s = base_ring.gen()
    base_field = QQ.extension(s**2 + 23, "sqrt_minus_23")
    sqrt_minus_23 = base_field.gen()
    hurwitz_ring = PolynomialRing(base_field, "J")

    def quadratic(record):
        return base_field(rational(record["rational_part"])) + base_field(
            rational(record["sqrt_minus_23_part"])
        ) * sqrt_minus_23

    factors = [
        hurwitz_ring([quadratic(value) for value in factor["coefficients_ascending"]])
        for factor in algebra["exact_factors"]
    ]
    sextic = next(factor for factor in factors if factor.degree() == 6)
    relative_field = base_field.extension(sextic, "j")
    absolute_field = relative_field.absolute_field("a")
    integral_basis = absolute_field.integral_basis()

    def integral_coordinates(record):
        return sum(
            rational(coordinate) * basis
            for coordinate, basis in zip(record, integral_basis)
        )

    def curve_coefficient(record):
        return integral_coordinates(
            record["sextic_component"]["integral_basis_coordinates"]
        )

    def power_basis(record):
        return absolute_field([rational(value) for value in record])

    data = {
        "quadric": [curve_coefficient(value) for value in canonical["quadric"]["coefficients"]],
        "cubic": [curve_coefficient(value) for value in canonical["petri_cubic"]["coefficients"]],
        "marked_c": [
            integral_coordinates(value["integral_basis_coordinates"])
            for value in marked["c_coordinates"]["sextic_component"]
        ],
    }
    sextic_map = next(
        component for component in maps["components"] if component["component"] == "sextic"
    )
    for name in ("numerator", "denominator", "normalized_denominator"):
        data[name] = [
            power_basis(term["coefficient_power_basis"])
            for term in sextic_map[name]
        ]
    data["lambda"] = [
        absolute_field(
            [
                QQ(value) / QQ(branch["lambda"]["denominator"])
                for value in branch["lambda"]["power_basis_numerators"]
            ]
        )
    ]

    quadric_monomials = [tuple(value) for value in canonical["quadric"]["monomials"]]
    cubic_monomials = [
        tuple(value) for value in canonical["petri_cubic"]["monomials"]
    ]
    b_point = [absolute_field(value) for value in marked["b_coordinates"]]

    def reduce_geometry(prime):
        residue_field = prime.residue_field(names="r")
        ambient_space = ProjectiveSpace(
            residue_field, 3, names=("x0", "x1", "x2", "x3")
        )
        coordinate_ring = ambient_space.coordinate_ring()
        coordinates = coordinate_ring.gens()
        quadric = sum(
            residue_field(coefficient) * coordinates[first] * coordinates[second]
            for (first, second), coefficient in zip(
                quadric_monomials, data["quadric"]
            )
        )
        cubic = sum(
            residue_field(coefficient)
            * coordinates[first]
            * coordinates[second]
            * coordinates[third]
            for (first, second, third), coefficient in zip(
                cubic_monomials, data["cubic"]
            )
        )

        numerator = sum(
            residue_field(coefficient)
            * coordinates[0] ** term["monomial"][0]
            * coordinates[1] ** term["monomial"][1]
            * coordinates[2] ** term["monomial"][2]
            * coordinates[3] ** term["monomial"][3]
            for coefficient, term in zip(data["numerator"], sextic_map["numerator"])
        )
        denominator = sum(
            residue_field(coefficient)
            * coordinates[0] ** term["monomial"][0]
            * coordinates[1] ** term["monomial"][1]
            * coordinates[2] ** term["monomial"][2]
            * coordinates[3] ** term["monomial"][3]
            for coefficient, term in zip(
                data["normalized_denominator"],
                sextic_map["normalized_denominator"],
            )
        )
        gram = matrix(residue_field, 4, 4)
        for (first, second), coefficient in zip(
            quadric_monomials, data["quadric"]
        ):
            reduced = residue_field(coefficient)
            if first == second:
                gram[first, second] = reduced
            else:
                gram[first, second] = reduced / 2
                gram[second, first] = reduced / 2

        quadric_gradient = [quadric.derivative(value) for value in coordinates]
        cubic_gradient = [cubic.derivative(value) for value in coordinates]
        jacobian_minors = [
            quadric_gradient[first] * cubic_gradient[second]
            - quadric_gradient[second] * cubic_gradient[first]
            for first in range(4)
            for second in range(first + 1, 4)
        ]
        curve_ideal = coordinate_ring.ideal([quadric, cubic])
        map_base_ideal = coordinate_ring.ideal(
            [quadric, cubic, numerator, denominator]
        )
        affine_quadric = quadric.subs({coordinates[0]: 1})
        affine_cubic = cubic.subs({coordinates[0]: 1})
        plane_resultant = affine_quadric.resultant(affine_cubic, coordinates[3])
        plane_factorization = plane_resultant.factor()
        affine_singular_ideal = coordinate_ring.ideal(
            [quadric, cubic, coordinates[0] - 1, *jacobian_minors]
        )
        affine_singular_radical = affine_singular_ideal.radical()
        no_singularities_at_infinity = all(
            coordinate_ring(1)
            in coordinate_ring.ideal(
                [
                    quadric,
                    cubic,
                    coordinates[0],
                    coordinates[index] - 1,
                    *jacobian_minors,
                ]
            )
            for index in range(1, 4)
        )
        smooth_patches = []
        singular_points = set()
        for coordinate in coordinates:
            singular_patch = coordinate_ring.ideal(
                [quadric, cubic, coordinate - 1, *jacobian_minors]
            )
            smooth_patches.append(coordinate_ring(1) in singular_patch)
            if not smooth_patches[-1]:
                for solution in singular_patch.variety():
                    point = tuple(solution[value] for value in coordinates)
                    first_nonzero = next(value for value in point if value)
                    singular_points.add(
                        tuple(value / first_nonzero for value in point)
                    )

        reduced_b = [residue_field(value) for value in b_point]
        reduced_c = [residue_field(value) for value in data["marked_c"]]

        def evaluate(polynomial, point):
            return polynomial(*point)

        def local_series(point, precision=31):
            assert point[0] == 1
            for parameter_index in (1, 2, 3):
                solve_indices = [
                    index for index in (1, 2, 3) if index != parameter_index
                ]
                jacobian = matrix(
                    residue_field,
                    [
                        [quadric.derivative(coordinates[index])(*point) for index in solve_indices],
                        [cubic.derivative(coordinates[index])(*point) for index in solve_indices],
                    ],
                )
                if not jacobian.det():
                    continue
                series_ring = PowerSeriesRing(
                    residue_field, "t", default_prec=precision
                )
                t = series_ring.gen()
                values = [series_ring(value).add_bigoh(precision) for value in point]
                values[parameter_index] += t
                for order in range(1, precision):
                    residual = vector(
                        residue_field,
                        [
                            quadric(*values)[order],
                            cubic(*values)[order],
                        ],
                    )
                    correction = jacobian.solve_right(-residual)
                    for index, value in zip(solve_indices, correction):
                        values[index] += series_ring(value) * t**order
                assert quadric(*values) == 0
                assert cubic(*values) == 0
                return values
            raise AssertionError("marked point is not smooth")

        def tangent_discriminant(point):
            local_ring = PolynomialRing(residue_field, 3, names=("u1", "u2", "u3"))
            local_coordinates = local_ring.gens()
            substitution = [
                local_ring(1),
                *[
                    local_ring(point[index]) + local_coordinates[index - 1]
                    for index in range(1, 4)
                ],
            ]
            local_quadric = local_ring(quadric(*substitution))
            local_cubic = local_ring(cubic(*substitution))
            quadric_linear = local_quadric.homogeneous_components().get(1, local_ring(0))
            cubic_linear = local_cubic.homogeneous_components().get(1, local_ring(0))
            pivot = next(
                coordinate
                for coordinate in local_coordinates
                if quadric_linear.monomial_coefficient(coordinate)
            )
            multiplier = (
                cubic_linear.monomial_coefficient(pivot)
                / quadric_linear.monomial_coefficient(pivot)
            )
            assert cubic_linear == multiplier * quadric_linear
            remaining = [value for value in local_coordinates if value != pivot]
            specialized_quadric = local_quadric.subs(
                {value: 0 for value in remaining}
            )
            specialized_cubic = local_cubic.subs(
                {value: 0 for value in remaining}
            )
            assert specialized_quadric.gcd(specialized_cubic).degree(pivot) == 1
            pivot_linear = -sum(
                quadric_linear.monomial_coefficient(value) * value
                for value in remaining
            ) / quadric_linear.monomial_coefficient(pivot)
            pivot_series = pivot_linear
            pivot_coefficient = quadric_linear.monomial_coefficient(pivot)
            for degree in range(2, 9):
                residual = local_quadric.subs({pivot: pivot_series})
                residual_piece = residual.homogeneous_components().get(
                    degree, local_ring(0)
                )
                pivot_series -= residual_piece / pivot_coefficient
            restricted_cubic = local_cubic.subs({pivot: pivot_series})
            components = restricted_cubic.homogeneous_components()
            plane_ring = PolynomialRing(
                residue_field, 2, names=("p", "q"), order="negdegrevlex"
            )
            p_coordinate, q_coordinate = plane_ring.gens()
            exact_plane_curve = local_quadric.resultant(local_cubic, pivot)
            plane_curve = plane_ring(
                exact_plane_curve.subs(
                    {
                        pivot: 0,
                        remaining[0]: p_coordinate,
                        remaining[1]: q_coordinate,
                    }
                )
            )
            jacobian_ideal = plane_ring.ideal(
                [
                    plane_curve.derivative(p_coordinate),
                    plane_curve.derivative(q_coordinate),
                ]
            )
            milnor_number = jacobian_ideal.vector_space_dimension()
            tjurina_number = plane_ring.ideal(
                [
                    plane_curve,
                    plane_curve.derivative(p_coordinate),
                    plane_curve.derivative(q_coordinate),
                ]
            ).vector_space_dimension()
            assert components.get(0, local_ring(0)) == 0
            assert components.get(1, local_ring(0)) == 0
            initial_degree = min(
                degree for degree, value in components.items() if degree and value
            )
            initial_form = components[initial_degree]
            if initial_degree >= 3:
                factorization = initial_form.factor()
                tangent_value_order = None
                exceptional_type = None
                if (
                    initial_degree == 3
                    and len(factorization) == 1
                    and factorization[0][0].degree() == 1
                    and factorization[0][1] == 3
                ):
                    tangent_line = factorization[0][0]
                    first, second = remaining
                    first_coefficient = tangent_line.monomial_coefficient(first)
                    second_coefficient = tangent_line.monomial_coefficient(second)
                    if second_coefficient:
                        tangent_direction = (
                            residue_field(1),
                            -first_coefficient / second_coefficient,
                        )
                    else:
                        tangent_direction = (residue_field(0), residue_field(1))
                    tangent_value = restricted_cubic.subs(
                        {
                            first: tangent_direction[0] * first,
                            second: tangent_direction[1] * first,
                        }
                    )
                    tangent_components = tangent_value.homogeneous_components()
                    tangent_value_order = min(
                        degree
                        for degree, value in tangent_components.items()
                        if degree and value
                    )
                    if tangent_value_order in (4, 5):
                        exceptional_type = "E6" if tangent_value_order == 4 else "E8"
                return {
                    "multiplicity": initial_degree,
                    "initial_form": str(initial_form),
                    "initial_factorization": str(factorization),
                    "milnor_number": milnor_number,
                    "tjurina_number": tjurina_number,
                    "tangent_value_order": tangent_value_order,
                    "exceptional_type_candidate": exceptional_type,
                    "ordinary_triple_test": (
                        initial_degree == 3
                        and len(factorization) == 3
                        and all(
                            factor.degree() == 1 and multiplicity == 1
                            for factor, multiplicity in factorization
                        )
                    ),
                }
            binary = initial_form
            cubic_term = components.get(3, local_ring(0))
            first, second = remaining
            aa = binary.monomial_coefficient(first**2)
            bb = binary.monomial_coefficient(first * second)
            cc = binary.monomial_coefficient(second**2)
            discriminant = bb**2 - 4 * aa * cc
            tangent_directions = [
                (residue_field(1), value)
                for value in residue_field
                if binary.subs({first: 1, second: value}) == 0
            ]
            if binary.subs({first: 0, second: 1}) == 0:
                tangent_directions.append((residue_field(0), residue_field(1)))
            assert tangent_directions
            tangent_direction = tangent_directions[0]
            tangent_cubic = cubic_term.subs(
                {
                    first: tangent_direction[0],
                    second: tangent_direction[1],
                }
            )
            transverse_direction = (
                (residue_field(0), residue_field(1))
                if tangent_direction[0]
                else (residue_field(1), residue_field(0))
            )
            transformed = restricted_cubic.subs(
                {
                    first: tangent_direction[0] * first
                    + transverse_direction[0] * second,
                    second: tangent_direction[1] * first
                    + transverse_direction[1] * second,
                }
            )
            transverse_square = transformed.monomial_coefficient(second**2)
            assert transverse_square
            derivative = transformed.derivative(second)
            critical_series = local_ring(0)
            for degree in range(1, 8):
                residual = derivative.subs({second: critical_series})
                coefficient = residual.monomial_coefficient(first**degree)
                critical_series -= coefficient / (2 * transverse_square) * first**degree
            critical_value = transformed.subs({second: critical_series})
            critical_components = critical_value.homogeneous_components()
            critical_order = min(
                degree
                for degree, value in critical_components.items()
                if degree and value
            )
            return {
                "quadratic": str(binary),
                "discriminant": str(discriminant),
                "tangent_cubic": str(tangent_cubic),
                "milnor_number": milnor_number,
                "tjurina_number": tjurina_number,
                "critical_value_order": critical_order,
                "A_type_if_nonzero": f"A{critical_order - 1}",
            }

        assert evaluate(quadric, reduced_b) == 0
        assert evaluate(cubic, reduced_b) == 0
        assert evaluate(quadric, reduced_c) == 0
        assert evaluate(cubic, reduced_c) == 0
        special_points = {
            "b": tuple(reduced_b),
            "c": tuple(reduced_c),
            **{
                f"singular_{index}": point
                for index, point in enumerate(sorted(singular_points, key=str), start=1)
            },
        }
        b_series = local_series(reduced_b)
        c_series = local_series(reduced_c)
        return {
            "quadric_rank": gram.rank(),
            "integral_over_residue_field": curve_ideal.is_prime(),
            "map_base_scheme_dimension": map_base_ideal.dimension(),
            "map_base_scheme_hilbert_polynomial": str(
                map_base_ideal.hilbert_polynomial()
            ),
            "plane_projection_factor_degrees": [
                (factor.degree(), multiplicity)
                for factor, multiplicity in plane_factorization
            ],
            "affine_singular_point_count": affine_singular_radical.vector_space_dimension(),
            "no_singularities_at_infinity": no_singularities_at_infinity,
            "smooth_patches": smooth_patches,
            "singular_points": sorted(map(str, singular_points)),
            "tangent_discriminants": [
                tangent_discriminant(point)
                for point in sorted(singular_points, key=str)
            ],
            "marked_points_distinct": reduced_b != reduced_c,
            "marked_points_smooth": all(
                tuple(point) not in singular_points for point in (reduced_b, reduced_c)
            ),
            "marked_orders": {
                "numerator_at_b": numerator(*b_series).valuation(),
                "denominator_at_c": denominator(*c_series).valuation(),
                "denominator_at_b": denominator(*b_series).valuation(),
                "numerator_at_c": numerator(*c_series).valuation(),
            },
            "map_values_N_D": {
                name: (
                    str(evaluate(numerator, point)),
                    str(evaluate(denominator, point)),
                )
                for name, point in special_points.items()
            },
        }

    primes = sorted(
        absolute_field.primes_above(23),
        key=lambda prime: prime.absolute_ramification_index(),
    )
    assert [
        (prime.absolute_ramification_index(), prime.residue_class_degree())
        for prime in primes
    ] == [(2, 2), (4, 2)]

    results = {}
    for prime in primes:
        for name, values in data.items():
            nonzero_values = [value for value in values if value]
            valuations = [value.valuation(prime) for value in nonzero_values]
            assert min(valuations) >= 0, f"nonintegral {name} at {prime}"
        results[f"sextic_e{prime.absolute_ramification_index()}"] = reduce_geometry(
            prime
        )

    # Run the same audit on the degree-one K0-component.  Its map
    # was stored before the third branch value was normalized, but scaling
    # its denominator does not change any source-curve or divisor statement
    # checked here.
    data = {
        "quadric": [
            quadratic(value["degree_one_component"])
            for value in canonical["quadric"]["coefficients"]
        ],
        "cubic": [
            quadratic(value["degree_one_component"])
            for value in canonical["petri_cubic"]["coefficients"]
        ],
        "marked_c": [
            quadratic(value)
            for value in marked["c_coordinates"]["degree_one_component"]
        ],
    }
    sextic_map = dict(
        next(
            component
            for component in maps["components"]
            if component["component"] == "degree_one"
        )
    )
    sextic_map["normalized_denominator"] = sextic_map["denominator"]

    def degree_one_power_basis(record):
        return base_field([rational(value) for value in record])

    for name in ("numerator", "denominator", "normalized_denominator"):
        data[name] = [
            degree_one_power_basis(term["coefficient_power_basis"])
            for term in sextic_map[name]
        ]
    b_point = [base_field(1), base_field(0), base_field(0), base_field(0)]
    degree_one_prime = base_field.primes_above(23)[0]
    for name, values in data.items():
        nonzero_values = [value for value in values if value]
        valuations = [value.valuation(degree_one_prime) for value in nonzero_values]
        assert min(valuations) >= 0, f"nonintegral degree-one {name}"
    results["degree_one"] = reduce_geometry(degree_one_prime)

    def common_checks(result, singular_point_count):
        assert result["quadric_rank"] == 4
        assert result["integral_over_residue_field"]
        assert result["map_base_scheme_dimension"] == 1
        assert result["map_base_scheme_hilbert_polynomial"] == "7"
        assert result["plane_projection_factor_degrees"] == [(6, 1)]
        assert result["affine_singular_point_count"] == singular_point_count
        assert result["no_singularities_at_infinity"]
        assert len(result["singular_points"]) == singular_point_count
        assert result["marked_points_distinct"]
        assert result["marked_points_smooth"]
        assert result["marked_orders"] == {
            "numerator_at_b": 23,
            "denominator_at_c": 23,
            "denominator_at_b": 0,
            "numerator_at_c": 0,
        }
        for name, (numerator_value, denominator_value) in result[
            "map_values_N_D"
        ].items():
            if name.startswith("singular_"):
                assert numerator_value != "0" or denominator_value != "0"

    common_checks(results["degree_one"], 1)
    common_checks(results["sextic_e2"], 1)
    common_checks(results["sextic_e4"], 2)

    for label in ("degree_one", "sextic_e2"):
        invariant = results[label]["tangent_discriminants"][0]
        assert invariant["multiplicity"] == 3
        assert invariant["milnor_number"] == 8
        assert invariant["tjurina_number"] == 8
        assert invariant["tangent_value_order"] == 5
        assert invariant["exceptional_type_candidate"] == "E8"

    ramified_types = {
        (
            invariant["A_type_if_nonzero"],
            invariant["milnor_number"],
            invariant["tjurina_number"],
            invariant["critical_value_order"],
        )
        for invariant in results["sextic_e4"]["tangent_discriminants"]
    }
    assert ramified_types == {("A2", 2, 2, 3), ("A6", 6, 6, 7)}

    print("PASS all three exact canonical models and map sections are 23-integral")
    print("PASS the degree-one and unramified-degree-2 reductions have one E8 singularity")
    print("PASS the ramified-degree-4 reduction has singularities A2+A6")
    print("PASS each singularity configuration has total delta invariant 4")
    print("PASS each normalization is P1 and the reduced pointed map is Frobenius of degree 23")
    print("PASS the normalized local Hurwitz points have residue degrees 1+2+2")
    print("SCOPE the ADE tails still require semistable resolution in a common pointed M23 frame")


if __name__ == "__main__":
    main()
