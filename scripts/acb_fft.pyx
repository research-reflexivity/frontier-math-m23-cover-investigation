# cython: language_level=3
"""Thin Sage/FLINT bridge for fast discrete Fourier transforms of Acb balls."""

from sage.libs.flint.acb cimport (
    _acb_vec_clear,
    _acb_vec_init,
    acb_add,
    acb_mul,
    acb_set,
)
from sage.libs.flint.acb_dft cimport acb_dft, acb_dft_inverse
from sage.libs.flint.acb_poly cimport _acb_poly_evaluate_vec_fast
from sage.libs.flint.types cimport acb_ptr, slong
from sage.rings.complex_arb cimport ComplexBall


def dft(values, bint inverse=False):
    """Return the FLINT Acb DFT of a nonempty list of Sage complex balls."""
    cdef slong length = len(values)
    cdef slong index
    cdef slong precision
    cdef acb_ptr source
    cdef acb_ptr target
    cdef ComplexBall template
    cdef ComplexBall value
    cdef ComplexBall output
    cdef list result

    if length == 0:
        return []
    template = values[0]
    precision = template.parent().precision()
    source = _acb_vec_init(length)
    target = _acb_vec_init(length)
    try:
        for index in range(length):
            value = values[index]
            if value.parent() is not template.parent():
                raise ValueError("all Acb values must have the same parent")
            acb_set(source + index, value.value)
        if inverse:
            acb_dft_inverse(target, source, length, precision)
        else:
            acb_dft(target, source, length, precision)
        result = []
        for index in range(length):
            output = template._new()
            acb_set(output.value, target + index)
            result.append(output)
        return result
    finally:
        _acb_vec_clear(source, length)
        _acb_vec_clear(target, length)


def evaluate_many(coefficients, points):
    """Evaluate one Acb polynomial at many Acb points using FLINT's product tree."""
    cdef slong coefficient_count = len(coefficients)
    cdef slong point_count = len(points)
    cdef slong index
    cdef slong precision
    cdef acb_ptr polynomial
    cdef acb_ptr inputs
    cdef acb_ptr outputs
    cdef ComplexBall template
    cdef ComplexBall value
    cdef ComplexBall output
    cdef list result

    if point_count == 0:
        return []
    if coefficient_count == 0:
        raise ValueError("the coefficient list must be nonempty")
    template = coefficients[0]
    precision = template.parent().precision()
    polynomial = _acb_vec_init(coefficient_count)
    inputs = _acb_vec_init(point_count)
    outputs = _acb_vec_init(point_count)
    try:
        for index in range(coefficient_count):
            value = coefficients[index]
            if value.parent() is not template.parent():
                raise ValueError("all Acb values must have the same parent")
            acb_set(polynomial + index, value.value)
        for index in range(point_count):
            value = points[index]
            if value.parent() is not template.parent():
                raise ValueError("all Acb values must have the same parent")
            acb_set(inputs + index, value.value)
        _acb_poly_evaluate_vec_fast(
            outputs, polynomial, coefficient_count, inputs, point_count, precision
        )
        result = []
        for index in range(point_count):
            output = template._new()
            acb_set(output.value, outputs + index)
            result.append(output)
        return result
    finally:
        _acb_vec_clear(polynomial, coefficient_count)
        _acb_vec_clear(inputs, point_count)
        _acb_vec_clear(outputs, point_count)


def evaluate_routes(blocks, targets, bases, factors):
    """Evaluate target-selected Acb polynomials routewise by compiled Horner."""
    cdef slong patch_count = len(blocks)
    cdef slong coefficient_count
    cdef slong route_count = len(targets)
    cdef slong patch
    cdef slong route
    cdef slong mode
    cdef slong target_index
    cdef slong precision
    cdef acb_ptr polynomial
    cdef acb_ptr outputs
    cdef ComplexBall template
    cdef ComplexBall value
    cdef ComplexBall base
    cdef ComplexBall factor
    cdef ComplexBall output
    cdef list result

    if route_count == 0:
        return []
    if patch_count == 0:
        raise ValueError("the block list must be nonempty")
    coefficient_count = len(blocks[0])
    if coefficient_count == 0:
        raise ValueError("the coefficient blocks must be nonempty")
    if len(bases) != route_count or len(factors) != route_count:
        raise ValueError("route arrays have different lengths")
    template = blocks[0][0]
    precision = template.parent().precision()
    polynomial = _acb_vec_init(patch_count * coefficient_count)
    outputs = _acb_vec_init(route_count)
    try:
        for patch in range(patch_count):
            if len(blocks[patch]) != coefficient_count:
                raise ValueError("coefficient blocks have different lengths")
            for mode in range(coefficient_count):
                value = blocks[patch][mode]
                acb_set(polynomial + patch * coefficient_count + mode, value.value)
        for route in range(route_count):
            target_index = targets[route]
            if target_index < 0 or target_index >= patch_count:
                raise ValueError("route target is out of range")
            base = bases[route]
            factor = factors[route]
            acb_set(
                outputs + route,
                polynomial + (target_index + 1) * coefficient_count - 1,
            )
            for mode in range(coefficient_count - 2, -1, -1):
                acb_mul(outputs + route, outputs + route, base.value, precision)
                acb_add(
                    outputs + route,
                    outputs + route,
                    polynomial + target_index * coefficient_count + mode,
                    precision,
                )
            acb_mul(outputs + route, outputs + route, factor.value, precision)
        result = []
        for route in range(route_count):
            output = template._new()
            acb_set(output.value, outputs + route)
            result.append(output)
        return result
    finally:
        _acb_vec_clear(polynomial, patch_count * coefficient_count)
        _acb_vec_clear(outputs, route_count)
