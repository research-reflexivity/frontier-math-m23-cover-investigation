#!/usr/bin/env sage-python
"""Stress-test the distinction between beta(0) and the third branch value.

The general valuation argument is written in the adversarial audit note.
This test reads existing finite Hensel residues as a probe and checks them
again at the stated finite precision. It is not independent reconstruction,
an all-orders incidence proof, or a new Galois-group computation.
"""
import json
from pathlib import Path

from sage.all import QQ
from harmonic_finite_equations import fixed_family, mixed


def order(a):
    coefficients = list(mixed.K(a))
    return min(2 * c.valuation(23) + j
               for j, c in enumerate(coefficients) if c)


record = json.loads(Path(__file__).with_name(
    "harmonic_local_recognition_checkpoint.json").read_text())
precision = 9
assert record["precision"] > precision
point = [mixed.K([QQ(c) for c in a]) for a in record["point"]]
g = fixed_family(point, precision)
assert not any(g.multiplier_candidate())
print("PASS_FINITE_GLOBAL_FUNCTION_RELATIONS_MOD_PI", precision, flush=True)
R = g.ramification_factor()
assert order(R[0]) >= 2
assert g.numerator(1, 0, 0, 0) == g.denominator(1, 0, 0, 0)


def scalar_and_remainders(factor):
    h, u, cut, _ = g.curve_in_quotient(factor**2)
    N = cut(g.numerator(1, u, h, g.z))
    D = cut(g.denominator(1, u, h, g.z))
    plain = cut(N-D)
    correction = g.A.coefficient(plain[0]/D[0])
    adjusted = cut(plain-correction*D)
    assert adjusted[0] == 0
    return correction, plain, adjusted


correction, plain, adjusted = scalar_and_remainders(R)
assert correction and order(correction) == 5
assert plain[0] and order(plain[0]) == 5
print("THIRD_BRANCH_VALUE_MINUS_BETA_AT_ORIGIN", correction,
      "order", order(correction),
      "leading_pi_digit", g.A.digit(correction, 5), flush=True)
print("UNADJUSTED_CONSTANT_REMAINDER", plain[0],
      "order", order(plain[0]), flush=True)

# The last two digits of R are not known from the divided differential.
# They must not be silently treated as zero when extracting this scalar.
perturbation = (mixed.pi**(precision-2)*(1+g.z**3+g.z**7)
                + mixed.pi**(precision-1)*(2+g.z**4))
perturbed, _, _ = scalar_and_remainders(R+perturbation)
assert g.A.coefficient(perturbed-correction) == 0
print("PASS_UNKNOWN_RAMIFICATION_DIGITS_DO_NOT_CHANGE_SCALAR_MOD_PI",
      precision, flush=True)
print("PASS_BRANCH_VALUE_NORMALIZATION_AUDIT", flush=True)
print("SCOPE: finite-residue probe, not all-orders construction; beta(0)=1 does not make the third branch value exactly 1", flush=True)
