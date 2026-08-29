#!/usr/bin/env sage-python
"""Explore the mixed-characteristic germs behind the certified ADE reductions.

This is a research script rather than a certificate.  It reuses the setup in
verification/verify_hurwitz_pointed_23.py and captures the local variables of
its three calls to reduce_geometry.  Keeping this separate prevents the
exploratory local-field calculations from weakening the exact certificate.
"""

from __future__ import annotations

import runpy
import sys
import os
from fractions import Fraction
from pathlib import Path

from sage.all import GF, PolynomialRing, PowerSeriesRing, matrix, vector


ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT / "verification" / "verify_hurwitz_pointed_23.py"


captured = []
_residue_lift_cache = {}


def trace_reduce_geometry(frame, event, arg):
    if frame.f_code.co_name == "reduce_geometry" and event == "return":
        captured.append(
            {
                name: frame.f_locals[name]
                for name in (
                    "prime",
                    "residue_field",
                    "quadric",
                    "cubic",
                    "numerator",
                    "denominator",
                    "ordered_singular_points",
                    "data",
                    "sextic_map",
                    "quadric_monomials",
                    "cubic_monomials",
                )
            }
        )
    return trace_reduce_geometry


sys.settrace(trace_reduce_geometry)
try:
    runpy.run_path(str(TARGET), run_name="__main__")
finally:
    sys.settrace(None)

assert len(captured) == 3


def residue_vector(value, residue_field):
    value = residue_field(value)
    coordinates = list(value.polynomial())
    coordinates += [GF(23)(0)] * (residue_field.degree() - len(coordinates))
    coordinates = coordinates[: residue_field.degree()]
    return vector(GF(23), coordinates)


def residue_lift(value, prime, residue_field):
    """Choose an integral number-field lift of one residue-field element."""

    cache_key = (id(prime), str(residue_field(value)))
    if cache_key in _residue_lift_cache:
        return _residue_lift_cache[cache_key]

    number_field = prime.number_field()
    integral_basis = number_field.integral_basis()
    image_matrix = matrix(
        GF(23),
        [
            residue_vector(residue_field(basis), residue_field)
            for basis in integral_basis
        ],
    ).transpose()
    pivots = image_matrix.pivots()
    assert len(pivots) == residue_field.degree()
    square = image_matrix.matrix_from_columns(pivots)
    coefficients = square.solve_right(residue_vector(value, residue_field))
    answer = sum(
        number_field(int(coefficient)) * integral_basis[pivot]
        for coefficient, pivot in zip(coefficients, pivots)
    )
    assert residue_field(answer) == residue_field(value)
    assert answer.valuation(prime) >= 0
    _residue_lift_cache[cache_key] = answer
    return answer


def exact_plane_germ(local, point):
    prime = local["prime"]
    number_field = prime.number_field()
    residue_field = local["residue_field"]
    ring = PolynomialRing(number_field, 4, names=("x0", "x1", "x2", "x3"))
    x0, x1, x2, x3 = ring.gens()
    coordinates = ring.gens()
    quadric = sum(
        number_field(coefficient) * coordinates[first] * coordinates[second]
        for (first, second), coefficient in zip(
            local["quadric_monomials"], local["data"]["quadric"]
        )
    )
    cubic = sum(
        number_field(coefficient)
        * coordinates[first]
        * coordinates[second]
        * coordinates[third]
        for (first, second, third), coefficient in zip(
            local["cubic_monomials"], local["data"]["cubic"]
        )
    )
    affine_quadric = quadric.subs({x0: 1})
    affine_cubic = cubic.subs({x0: 1})
    plane = affine_quadric.resultant(affine_cubic, x3)
    lift_x1 = residue_lift(point[1], prime, residue_field)
    lift_x2 = residue_lift(point[2], prime, residue_field)
    germ = plane.subs({x1: x1 + lift_x1, x2: x2 + lift_x2})
    return germ


def normalized_terms(germ, prime, residue_field):
    uniformizer = prime.gens_reduced()[0]
    terms = []
    for exponent, coefficient in germ.dict().items():
        if not coefficient:
            continue
        valuation = coefficient.valuation(prime)
        unit = coefficient / uniformizer**valuation
        residue = residue_field(unit)
        assert residue
        terms.append(
            (
                valuation,
                sum(exponent),
                exponent[1],
                exponent[2],
                str(residue),
            )
        )
    return sorted(terms)


def tangent_normalize(germ, prime, residue_field):
    number_field = prime.number_field()
    residue_ring = PolynomialRing(residue_field, 2, names=("P", "Q"))
    P, Q = residue_ring.gens()
    special = residue_ring(0)
    for exponent, coefficient in germ.dict().items():
        if coefficient and coefficient.valuation(prime) == 0:
            special += residue_field(coefficient) * (
                P**exponent[1] * Q**exponent[2]
            )
    minimum_degree = min(
        sum(exponent)
        for exponent, coefficient in special.dict().items()
        if coefficient
    )
    tangent_cone = special.homogeneous_components()[minimum_degree]
    tangent_factorization = tangent_cone.factor()
    tangent_line = next(
        factor
        for factor, multiplicity in tangent_factorization
        if factor.degree() == 1 and multiplicity >= 2
    )
    aa = tangent_line.monomial_coefficient(P)
    bb = tangent_line.monomial_coefficient(Q)
    lift_aa = residue_lift(aa, prime, residue_field)
    lift_bb = residue_lift(bb, prime, residue_field)
    ring = germ.parent()
    _, X, Y, _ = ring.gens()
    if aa:
        # P=bX+a^(-1)Y, Q=-aX makes aP+bQ=Y.
        coordinate_substitution = {
            X: lift_bb * X + (number_field(1) / lift_aa) * Y,
            Y: -lift_aa * X,
        }
    else:
        coordinate_substitution = {X: X, Y: Y / lift_bb}
    transformed = germ.subs(coordinate_substitution)
    transformed_special = residue_ring(0)
    for exponent, coefficient in transformed.dict().items():
        if coefficient and coefficient.valuation(prime) == 0:
            transformed_special += residue_field(coefficient) * (
                P**exponent[1] * Q**exponent[2]
            )
    return {
        "minimum_degree": minimum_degree,
        "tangent_cone": tangent_cone,
        "tangent_factorization": tangent_factorization,
        "transformed": transformed,
        "transformed_special": transformed_special,
        "coordinate_substitution": coordinate_substitution,
    }


def exact_map_numerator_at_singularity(local, point, normalization):
    """Pull beta-beta(tau) back to the exact plane projection.

    A linear subresultant expresses x3 rationally in x1,x2 on the canonical
    complete intersection.  Clearing its denominator in numerator and
    denominator does not affect the local map because that denominator is a
    unit at the singular point.
    """

    prime = local["prime"]
    number_field = prime.number_field()
    residue_field = local["residue_field"]
    ring = PolynomialRing(number_field, 4, names=("x0", "x1", "x2", "x3"))
    x0, x1, x2, x3 = ring.gens()
    coordinates = ring.gens()
    quadric = sum(
        number_field(coefficient) * coordinates[first] * coordinates[second]
        for (first, second), coefficient in zip(
            local["quadric_monomials"], local["data"]["quadric"]
        )
    )
    cubic = sum(
        number_field(coefficient)
        * coordinates[first]
        * coordinates[second]
        * coordinates[third]
        for (first, second, third), coefficient in zip(
            local["cubic_monomials"], local["data"]["cubic"]
        )
    )
    map_terms = local["sextic_map"]

    def map_section(name):
        return sum(
            number_field(coefficient)
            * x0 ** term["monomial"][0]
            * x1 ** term["monomial"][1]
            * x2 ** term["monomial"][2]
            * x3 ** term["monomial"][3]
            for coefficient, term in zip(local["data"][name], map_terms[name])
        )

    numerator = map_section("numerator").subs({x0: 1})
    denominator = map_section("normalized_denominator").subs({x0: 1})
    affine_quadric = quadric.subs({x0: 1})
    affine_cubic = cubic.subs({x0: 1})
    subresultants = affine_quadric.subresultants(affine_cubic, x3)
    linear = next(
        value for value in reversed(subresultants) if value.degree(x3) == 1
    )
    aa = linear.derivative(x3)
    bb = linear.subs({x3: 0})
    assert residue_field(aa.subs({x1: residue_lift(point[1], prime, residue_field), x2: residue_lift(point[2], prime, residue_field)}))
    x3_expression = -bb / aa
    clearing_degree = max(numerator.degree(x3), denominator.degree(x3))
    plane_numerator = ring(aa**clearing_degree * numerator.subs({x3: x3_expression}))
    plane_denominator = ring(aa**clearing_degree * denominator.subs({x3: x3_expression}))

    lift_x1 = residue_lift(point[1], prime, residue_field)
    lift_x2 = residue_lift(point[2], prime, residue_field)
    beta_residue = residue_field(
        numerator.subs(
            {
                x1: lift_x1,
                x2: lift_x2,
                x3: residue_lift(point[3], prime, residue_field),
            }
        )
        / denominator.subs(
            {
                x1: lift_x1,
                x2: lift_x2,
                x3: residue_lift(point[3], prime, residue_field),
            }
        )
    )
    beta_lift = residue_lift(beta_residue, prime, residue_field)
    map_numerator = plane_numerator - beta_lift * plane_denominator
    translated = map_numerator.subs({x1: x1 + lift_x1, x2: x2 + lift_x2})
    transformed = translated.subs(normalization["coordinate_substitution"])
    denominator_translated = plane_denominator.subs(
        {x1: x1 + lift_x1, x2: x2 + lift_x2}
    )
    denominator_transformed = denominator_translated.subs(
        normalization["coordinate_substitution"]
    )
    return transformed, denominator_transformed, beta_lift, beta_residue


def weighted_initial_form(polynomial, prime, residue_field, weight_p, weight_q):
    """Return the terms on the lowest mixed-characteristic weight face."""

    uniformizer = prime.gens_reduced()[0]
    records = []
    for exponent, coefficient in polynomial.dict().items():
        if not coefficient:
            continue
        valuation = int(coefficient.valuation(prime))
        weight = Fraction(valuation) + weight_p * exponent[1] + weight_q * exponent[2]
        residue = residue_field(coefficient / uniformizer**valuation)
        records.append((weight, valuation, exponent[1], exponent[2], str(residue)))
    minimum = min(record[0] for record in records)
    return minimum, [record for record in records if record[0] == minimum]


def centered_map_initial_form(
    map_numerator,
    map_denominator,
    curve,
    prime,
    residue_field,
    weight_p,
    weight_q,
    stop_weight=None,
    curve_base_exponents=((5, 0), (0, 3)),
    max_steps=20,
):
    """Remove curve-equation multiples and constant target displacements."""

    _, P, Q, _ = map_numerator.parent().gens()
    current = map_numerator
    steps = []
    for _ in range(max_steps):
        weight, initial = weighted_initial_form(
            current, prime, residue_field, weight_p, weight_q
        )
        if stop_weight is not None and weight >= stop_weight:
            return current, weight, initial, steps
        if all(record[2] == record[3] == 0 for record in initial):
            correction = (
                current.constant_coefficient()
                / map_denominator.constant_coefficient()
            )
            current -= correction * map_denominator
            steps.append(("target", weight, int(correction.valuation(prime))))
            continue

        reduced = False
        for _, _, exponent_p, exponent_q, _ in initial:
            offsets = []
            for base_p, base_q in curve_base_exponents:
                if exponent_p >= base_p and exponent_q >= base_q:
                    offsets.append(
                        (
                            exponent_p - base_p,
                            exponent_q - base_q,
                            P**base_p * Q**base_q,
                        )
                    )
            for offset_p, offset_q, base_monomial in offsets:
                monomial = P**offset_p * Q**offset_q
                ratio = (
                    current.monomial_coefficient(monomial * base_monomial)
                    / curve.monomial_coefficient(base_monomial)
                )
                candidate = current - ratio * monomial * curve
                candidate_weight, candidate_initial = weighted_initial_form(
                    candidate, prime, residue_field, weight_p, weight_q
                )
                old_support = {(record[2], record[3]) for record in initial}
                new_support = {
                    (record[2], record[3]) for record in candidate_initial
                }
                if candidate_weight > weight or (
                    candidate_weight == weight
                    and len(new_support & old_support) < len(old_support)
                ):
                    current = candidate
                    steps.append(
                        (
                            "curve",
                            weight,
                            offset_p,
                            offset_q,
                            int(ratio.valuation(prime)),
                        )
                    )
                    reduced = True
                    break
            if reduced:
                break
        if not reduced:
            return current, weight, initial, steps
    raise AssertionError(f"map-face centering did not terminate: {steps}")


def graded_initial_form(polynomial, prime, residue_field, weight_p, weight_q):
    """Return one mixed-characteristic initial form with a formal pi."""

    graded_ring = PolynomialRing(residue_field, 3, names=("S", "P", "Q"))
    S, P, Q = graded_ring.gens()
    uniformizer = prime.gens_reduced()[0]
    records = []
    for exponent, coefficient in polynomial.dict().items():
        if not coefficient:
            continue
        valuation = int(coefficient.valuation(prime))
        weight = Fraction(valuation) + weight_p * exponent[1] + weight_q * exponent[2]
        residue = residue_field(coefficient / uniformizer**valuation)
        records.append((weight, valuation, exponent[1], exponent[2], residue))
    minimum = min(record[0] for record in records)
    initial = sum(
        residue * S**valuation * P**exponent_p * Q**exponent_q
        for weight, valuation, exponent_p, exponent_q, residue in records
        if weight == minimum
    )
    return minimum, initial


def graded_centered_map_initial_form(
    map_numerator,
    map_denominator,
    curve,
    prime,
    residue_field,
    weight_p,
    weight_q,
    stop_weight,
    max_steps=30,
):
    """Center the map by exact lifts of whole associated-graded quotients."""

    number_field = prime.number_field()
    uniformizer = prime.gens_reduced()[0]
    _, P_exact, Q_exact, _ = map_numerator.parent().gens()
    curve_weight, curve_initial = graded_initial_form(
        curve, prime, residue_field, weight_p, weight_q
    )
    S, P, Q = curve_initial.parent().gens()
    current = map_numerator
    steps = []
    for _ in range(max_steps):
        weight, initial = graded_initial_form(
            current, prime, residue_field, weight_p, weight_q
        )
        if weight >= stop_weight:
            return current, weight, initial, steps, curve_initial
        quotient, remainder = initial.quo_rem(curve_initial)
        if quotient:
            quotient_lift = number_field(0)
            for exponent, coefficient in quotient.dict().items():
                quotient_lift += (
                    residue_lift(coefficient, prime, residue_field)
                    * uniformizer**exponent[0]
                    * P_exact**exponent[1]
                    * Q_exact**exponent[2]
                )
            current -= quotient_lift * curve
            steps.append(("curve", weight, len(quotient.dict())))
            continue
        if remainder and all(
            exponent[1] == exponent[2] == 0
            for exponent, coefficient in remainder.dict().items()
            if coefficient
        ):
            correction = (
                current.constant_coefficient()
                / map_denominator.constant_coefficient()
            )
            current -= correction * map_denominator
            steps.append(("target", weight, int(correction.valuation(prime))))
            continue
        if remainder:
            return current, weight, initial, steps, curve_initial
        raise AssertionError("zero graded initial form")
    raise AssertionError(f"graded map-face centering did not terminate: {steps}")


def truncate_weight(polynomial, weights, cutoff):
    """Discard monomials of weighted degree at least cutoff."""

    return polynomial.parent()(
        {
            exponent: coefficient
            for exponent, coefficient in polynomial.dict().items()
            if coefficient
            and sum(weight * power for weight, power in zip(weights, exponent))
            < cutoff
        }
    )


def residue_expansion(polynomial, prime, residue_field, weights, cutoff):
    """Expand exact coefficients pi-adically in k[S,P,Q] below a weight."""

    graded_ring = PolynomialRing(residue_field, 3, names=("S", "P", "Q"))
    S, P, Q = graded_ring.gens()
    uniformizer = prime.gens_reduced()[0]
    maximum_valuation = (cutoff - 1) // weights[0] + 1
    answer = graded_ring(0)
    for exponent, coefficient in polynomial.dict().items():
        if not coefficient:
            continue
        assert exponent[0] == exponent[3] == 0
        remainder = coefficient
        while remainder:
            valuation = int(remainder.valuation(prime))
            if valuation >= maximum_valuation:
                break
            monomial_weight = (
                weights[0] * valuation
                + weights[1] * exponent[1]
                + weights[2] * exponent[2]
            )
            if monomial_weight >= cutoff:
                break
            digit = residue_field(remainder / uniformizer**valuation)
            answer += digit * S**valuation * P**exponent[1] * Q**exponent[2]
            remainder -= residue_lift(digit, prime, residue_field) * uniformizer**valuation
            assert not remainder or remainder.valuation(prime) > valuation
    return truncate_weight(answer, weights, cutoff)


def weight_initial(polynomial, weights):
    records = [
        (sum(weight * power for weight, power in zip(weights, exponent)), exponent, coefficient)
        for exponent, coefficient in polynomial.dict().items()
        if coefficient
    ]
    minimum = min(record[0] for record in records)
    return minimum, polynomial.parent()(
        {exponent: coefficient for weight, exponent, coefficient in records if weight == minimum}
    )


def critical_coordinate_model(exact_curve, prime, residue_field, weights, cutoff):
    """Put a double-tangent germ in critical Q-coordinate below cutoff."""

    expanded = residue_expansion(
        exact_curve, prime, residue_field, weights, cutoff
    )
    S, P, Q = expanded.parent().gens()
    derivative = expanded.derivative(Q)
    transverse_unit = derivative.derivative(Q).constant_coefficient()
    assert transverse_unit
    critical_section = expanded.parent()(0)
    for _ in range(cutoff):
        residual = truncate_weight(
            derivative.subs({Q: critical_section}), weights, cutoff
        )
        if not residual:
            break
        critical_section = truncate_weight(
            critical_section - residual / transverse_unit, weights, cutoff
        )
    else:
        raise AssertionError("critical-coordinate iteration did not converge")
    assert not truncate_weight(
        derivative.subs({Q: critical_section}), weights, cutoff
    )
    shifted = truncate_weight(
        expanded.subs({Q: Q + critical_section}), weights, cutoff
    )
    return shifted, critical_section


def centered_residual_map(
    exact_numerator,
    exact_denominator,
    critical_section,
    critical_curve,
    prime,
    residue_field,
    weights,
    stop_weight,
    cutoff,
):
    """Center a residual map in the critical coordinate entirely over k."""

    numerator = residue_expansion(
        exact_numerator, prime, residue_field, weights, cutoff
    )
    denominator = residue_expansion(
        exact_denominator, prime, residue_field, weights, cutoff
    )
    S, P, Q = numerator.parent().gens()
    numerator = truncate_weight(
        numerator.subs({Q: Q + critical_section}), weights, cutoff
    )
    denominator = truncate_weight(
        denominator.subs({Q: Q + critical_section}), weights, cutoff
    )
    curve_weight, curve_initial = weight_initial(critical_curve, weights)
    denominator_unit = denominator.constant_coefficient()
    assert denominator_unit
    steps = []
    for _ in range(cutoff):
        weight, initial = weight_initial(numerator, weights)
        if weight >= stop_weight:
            return numerator, weight, initial, steps, curve_initial
        quotient, remainder = initial.quo_rem(curve_initial)
        if quotient:
            numerator = truncate_weight(
                numerator - quotient * critical_curve, weights, cutoff
            )
            steps.append(("curve", weight, len(quotient.dict())))
            continue
        if remainder and all(
            exponent[1] == exponent[2] == 0
            for exponent, coefficient in remainder.dict().items()
            if coefficient
        ):
            correction = remainder / denominator_unit
            numerator = truncate_weight(
                numerator - correction * denominator, weights, cutoff
            )
            steps.append(("target", weight, len(correction.dict())))
            continue
        return numerator, weight, initial, steps, curve_initial
    raise AssertionError(f"residual map centering did not terminate: {steps}")


def critical_discriminant_terms(germ, prime, residue_field):
    """Return the exact Q-discriminant terms and their lower Newton hull.

    At a double tangent the derivative with respect to Q has one simple root
    reducing to Q=0.  All its other roots contribute a unit, so the lower
    Newton polygon of Res_Q(F,F_Q) agrees with that of the critical value
    F(P,h(P)), up to multiplication by a unit power series.
    """

    ring = germ.parent()
    _, P, Q, _ = ring.gens()
    first_derivative = germ.derivative(Q)
    discriminant = germ.resultant(first_derivative, Q)
    uniformizer = prime.gens_reduced()[0]
    terms = []
    for exponent_tuple, coefficient in discriminant.dict().items():
        if not coefficient:
            continue
        assert exponent_tuple[0] == exponent_tuple[2] == exponent_tuple[3] == 0
        exponent = exponent_tuple[1]
        valuation = coefficient.valuation(prime)
        residue = residue_field(coefficient / uniformizer**valuation)
        terms.append((exponent, valuation, str(residue)))
    terms.sort()

    # Monotone-chain lower convex hull.  Collinear intermediate points are
    # retained separately below by reporting all terms lying on each face.
    hull = []
    for exponent, valuation, _ in terms:
        point = (exponent, valuation)
        while len(hull) >= 2:
            x0, y0 = hull[-2]
            x1, y1 = hull[-1]
            x2, y2 = point
            cross = (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0)
            if cross > 0:
                break
            hull.pop()
        hull.append(point)
    return discriminant, terms, hull


for index, local in enumerate(captured, start=1):
    if os.environ.get("M23_ADE_ONLY_RATIONAL") and index != 3:
        continue
    if os.environ.get("M23_ADE_ONLY_A2") and index != 2:
        continue
    if os.environ.get("M23_ADE_ONLY_UNRAMIFIED") and index != 1:
        continue
    if os.environ.get("M23_ADE_ONLY_A6") and index != 2:
        continue
    prime = local["prime"]
    residue_field = local["residue_field"]
    number_field = prime.number_field()
    print(
        "LOCAL",
        index,
        "absolute_e=",
        prime.absolute_ramification_index(),
        "residue_degree=",
        prime.residue_class_degree(),
        "number_field_degree=",
        number_field.degree(),
    )
    print("  residue_modulus=", residue_field.modulus())
    print("  singular_points=", local["ordered_singular_points"])
    print(
        "  prime_generators=",
        [(str(value), value.valuation(prime)) for value in prime.gens_reduced()],
        "two_generator=",
        [(str(value), value.valuation(prime)) for value in prime.gens_two()],
    )
    print(
        "  prime_methods=",
        [
            name
            for name in (
                "reduce",
                "residue_field",
                "uniformizer",
                "valuation",
                "gens_reduced",
                "gens_two",
            )
            if hasattr(prime, name)
        ],
    )
    print(
        "  residue_methods=",
        [
            name
            for name in ("lift", "teichmuller", "modulus")
            if hasattr(residue_field, name)
        ],
    )
    for point_index, point in enumerate(local["ordered_singular_points"], start=1):
        if os.environ.get("M23_ADE_ONLY_A2") and point_index != 1:
            continue
        if os.environ.get("M23_ADE_ONLY_A6") and point_index != 2:
            continue
        germ = exact_plane_germ(local, point)
        normalization = tangent_normalize(germ, prime, residue_field)
        terms = normalized_terms(normalization["transformed"], prime, residue_field)
        print(
            "  GERM",
            point_index,
            "multiplicity=",
            normalization["minimum_degree"],
            "tangent_cone=",
            normalization["tangent_cone"],
            "factorization=",
            normalization["tangent_factorization"],
        )
        print("    transformed_special=", normalization["transformed_special"])
        print("  GERM", point_index, "lowest_valuation_terms=")
        for term in terms[:48]:
            print("   ", term)
        map_tail_case = None
        if index == 3 and point_index == 1:
            map_tail_case = {
                "weights": (Fraction(1, 5), Fraction(1, 3)),
                "stop_weight": Fraction(23, 15),
                "curve_bases": ((5, 0), (0, 3)),
            }
        elif index == 1 and point_index == 1:
            map_tail_case = {
                "weights": (Fraction(1, 5), Fraction(1, 3)),
                "stop_weight": Fraction(23, 15),
                "curve_bases": ((5, 0), (0, 3)),
            }
        elif index == 2 and point_index == 1:
            map_tail_case = {
                "weights": (Fraction(1, 2), Fraction(3, 4)),
                "stop_weight": Fraction(23, 4),
                "curve_bases": ((3, 0), (0, 2)),
            }
        elif index == 2 and point_index == 2:
            map_tail_case = {
                "critical_weights": (7, 2, 7),
                "stop_weight": 23,
                "cutoff": 31,
            }
        if map_tail_case is not None:
            (
                map_numerator,
                map_denominator,
                beta_lift,
                beta_residue,
            ) = exact_map_numerator_at_singularity(local, point, normalization)
            print("    beta_residue=", beta_residue)
            if "critical_weights" in map_tail_case:
                critical_curve, critical_section = critical_coordinate_model(
                    normalization["transformed"],
                    prime,
                    residue_field,
                    map_tail_case["critical_weights"],
                    map_tail_case["cutoff"],
                )
                (
                    centered_numerator,
                    centered_weight,
                    centered_initial,
                    centering_steps,
                    curve_initial,
                ) = centered_residual_map(
                    map_numerator,
                    map_denominator,
                    critical_section,
                    critical_curve,
                    prime,
                    residue_field,
                    map_tail_case["critical_weights"],
                    map_tail_case["stop_weight"],
                    map_tail_case["cutoff"],
                )
                print("    critical_section=", critical_section)
            else:
                map_weight, map_initial = weighted_initial_form(
                    map_numerator,
                    prime,
                    residue_field,
                    *map_tail_case["weights"],
                )
                print("    map_weight=", map_weight)
                print("    map_initial_terms=", map_initial)
                (
                    centered_numerator,
                    centered_weight,
                    centered_initial,
                    centering_steps,
                    curve_initial,
                ) = graded_centered_map_initial_form(
                    map_numerator,
                    map_denominator,
                    normalization["transformed"],
                    prime,
                    residue_field,
                    *map_tail_case["weights"],
                    map_tail_case["stop_weight"],
                )
            print("    centering_steps=", centering_steps)
            print("    curve_initial=", curve_initial)
            print("    centered_map_weight=", centered_weight)
            print("    centered_map_initial_terms=", centered_initial)
            S, P, Q = curve_initial.parent().gens()
            tail_curve = curve_initial.subs({S: 1})
            tail_map = centered_initial.subs({S: 1})
            _, tail_map_reduced = tail_map.quo_rem(tail_curve)
            print("    tail_curve=", tail_curve)
            print("    tail_map_reduced=", tail_map_reduced)
            if "critical_weights" in map_tail_case:
                special_ring = PolynomialRing(
                    residue_field, 2, names=("P", "Q")
                )
                PP, QQ = special_ring.gens()

                def special_reduction(exact_polynomial):
                    answer = special_ring(0)
                    for exponent, coefficient in exact_polynomial.dict().items():
                        if coefficient and coefficient.valuation(prime) == 0:
                            answer += residue_field(coefficient) * PP**exponent[1] * QQ**exponent[2]
                    return answer

                special_curve_exact = special_reduction(
                    normalization["transformed"]
                )
                special_map_exact = special_reduction(map_numerator)
                special_resultant = special_curve_exact.resultant(
                    special_map_exact, QQ
                )
                special_p_order = min(
                    exponent[0]
                    for exponent, coefficient in special_resultant.dict().items()
                    if coefficient
                )
                print("    special_map_intersection_P_order=", special_p_order)
                outer_weights = (15, 2, 7)
                outer_curve, outer_critical_section = critical_coordinate_model(
                    normalization["transformed"],
                    prime,
                    residue_field,
                    outer_weights,
                    32,
                )
                _, outer_curve_initial = weight_initial(
                    outer_curve, outer_weights
                )
                print("    outer_curve_initial=", outer_curve_initial)

                h0_graded = outer_critical_section.subs(
                    {outer_critical_section.parent().gens()[0]: 0}
                )
                h0 = special_ring(0)
                for exponent, coefficient in h0_graded.dict().items():
                    if coefficient:
                        assert exponent[0] == exponent[2] == 0
                        h0 += residue_field(coefficient) * PP**exponent[1]
                shifted_special_curve = special_curve_exact.subs({QQ: QQ + h0})
                shifted_special_map = special_map_exact.subs({QQ: QQ + h0})
                q2_coefficient = shifted_special_curve.monomial_coefficient(QQ**2)
                p7_coefficient = shifted_special_curve.monomial_coefficient(PP**7)
                leading_square = -p7_coefficient / q2_coefficient
                extension_ring = PolynomialRing(residue_field, "z")
                z = extension_ring.gen()
                parameter_field = residue_field.extension(
                    z**2 - leading_square, "w"
                )
                w = parameter_field.gen()
                series_ring = PowerSeriesRing(
                    parameter_field, "s", default_prec=32
                )
                s = series_ring.gen()

                def evaluate_series(polynomial, p_value, q_value):
                    return sum(
                        parameter_field(coefficient)
                        * p_value**exponent[0]
                        * q_value**exponent[1]
                        for exponent, coefficient in polynomial.dict().items()
                        if coefficient
                    )

                v = series_ring(w)
                shifted_derivative = shifted_special_curve.derivative(QQ)
                for _ in range(6):
                    curve_value = (
                        evaluate_series(
                            shifted_special_curve, s**2, s**7 * v
                        )
                        / s**14
                    )
                    derivative_value = (
                        evaluate_series(
                            shifted_derivative, s**2, s**7 * v
                        )
                        / s**7
                    )
                    v = (v - curve_value / derivative_value).add_bigoh(32)
                parameterized_map = evaluate_series(
                    shifted_special_map, s**2, s**7 * v
                )
                print("    special_parameterized_map_order=", parameterized_map.valuation())
                print("    special_parameterized_map_lead=", parameterized_map[23])
                outer_x23 = parameterized_map[23]
                inner_x8 = parameter_field(
                    tail_map_reduced.monomial_coefficient(P**4)
                )
                normalized_x8 = inner_x8 / outer_x23
                assert outer_x23 and inner_x8 and normalized_x8
                print(
                    "    outer_residual_map=",
                    f"X^23 + ({normalized_x8})*X^8",
                )
                print("    outer_zero_fibre=8+15 with fifteen simple nonzero points")
        if normalization["minimum_degree"] == 2:
            discriminant, discriminant_terms, lower_hull = critical_discriminant_terms(
                normalization["transformed"], prime, residue_field
            )
            print("    discriminant_terms=", discriminant_terms)
            print("    discriminant_lower_hull=", lower_hull)
