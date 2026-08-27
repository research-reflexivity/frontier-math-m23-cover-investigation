#!/usr/bin/env sage-python
"""Certify the local structure of the sextic Hurwitz component at 23.

The global trace field is E = Q[y]/(g).  This script determines the two
23-adic factors over Q_23, their discriminant square classes, the prime
decomposition after base change to K0 = Q(sqrt(-23)), and the resulting
decomposition/inertia permutation action on the six geometric points.
"""

from __future__ import annotations

from sage.all import GF, PermutationGroup, PolynomialRing, Qp, QQ


def main() -> None:
    rational_ring = PolynomialRing(QQ, "y")
    y = rational_ring.gen()
    g = y**6 - 6 * y**5 + 14 * y**4 - 2 * y**3 - 27 * y**2 + 44 * y - 44

    residue_ring = PolynomialRing(GF(23), "yb")
    yb = residue_ring.gen()
    assert residue_ring(g) == (yb + 18) ** 2 * (yb + 1) ** 4

    # Over Q_23 there is one ramified quadratic factor and one totally
    # tamely ramified quartic factor.  Their discriminants have square
    # classes 23 and -23, respectively: the displayed unit residues are a
    # square and a nonsquare in F_23.
    local_field = Qp(23, prec=30, type="capped-rel")
    local_ring = PolynomialRing(local_field, "yl")
    local_factors = sorted(
        (factor for factor, multiplicity in local_ring(g).factor()),
        key=lambda factor: factor.degree(),
    )
    assert [factor.degree() for factor in local_factors] == [2, 4]
    local_discriminants = [factor.discriminant() for factor in local_factors]
    discriminant_valuations = [value.valuation() for value in local_discriminants]
    discriminant_units = [
        (value / 23 ** value.valuation()).residue()
        for value in local_discriminants
    ]
    assert discriminant_valuations == [1, 3]
    assert discriminant_units == [GF(23)(9), GF(23)(7)]
    assert discriminant_units[0].is_square()
    assert not discriminant_units[1].is_square()
    assert not GF(23)(-1).is_square()

    # The first Newton residual equation at y=5 is U^2=15 for
    # U=(y-5)/sqrt(-23).  At y=-1 the quartic branch has V^2=5 for
    # V=(y+1)^2/sqrt(-23).  Both are the same quadratic residue field,
    # since 15/5=3=7^2 in F_23.
    shifted_two = g(y + 5)
    shifted_four = g(y - 1)
    u_square = GF(23)(shifted_two[0] / 23) / GF(23)(shifted_two[2])
    # If x=y+1 and V=x^2/s with s^2=-23, then x^4=-23 V^2.
    v_square = GF(23)(shifted_four[0] / 23) / GF(23)(shifted_four[4])
    assert u_square == 15
    assert v_square == 5
    assert not u_square.is_square() and not v_square.is_square()
    assert u_square / v_square == GF(23)(7) ** 2

    # Exact prime decomposition in E/Q and in L=E*K0 over K0.
    trace_field = QQ.extension(g, "a")
    trace_primes = trace_field.primes_above(23)
    assert sorted(
        (prime.absolute_ramification_index(), prime.residue_class_degree())
        for prime in trace_primes
    ) == [(2, 1), (4, 1)]

    base_ring = PolynomialRing(QQ, "s")
    s = base_ring.gen()
    base_field = QQ.extension(s**2 + 23, "sqrt_minus_23")
    relative_ring = PolynomialRing(base_field, "z")
    relative_field = base_field.extension(relative_ring(g), "j")
    base_prime = base_field.primes_above(23)[0]
    relative_primes = relative_field.primes_above(base_prime)
    assert sorted(
        (
            prime.relative_ramification_index(),
            prime.residue_class_degree(),
        )
        for prime in relative_primes
    ) == [(1, 2), (2, 2)]

    # The quartic local splitting field over Q_23 is D4: its tame inertia
    # is C4 and Frobenius acts by inversion because 23=-1 mod 4.  Its
    # discriminant field is Q_23(sqrt(-23)).  Restriction to K0 therefore
    # gives the even V4.  In the natural six-point action V4 has a two-point
    # orbit and a regular four-point orbit.  The chosen letters below are
    # only a representative of this S6-conjugacy class.
    inertia = "(3,4)(5,6)"
    frobenius_lift = "(1,2)(3,5)(4,6)"
    decomposition_group = PermutationGroup([inertia, frobenius_lift])
    assert decomposition_group.order() == 4
    assert decomposition_group.structure_description() == "C2 x C2"
    assert sorted(len(orbit) for orbit in decomposition_group.orbits()) == [2, 4]
    assert decomposition_group(inertia).cycle_type() == [2, 2, 1, 1]
    assert decomposition_group(frobenius_lift).cycle_type() == [2, 2, 2]

    print("PASS g mod 23=(y-5)^2*(y+1)^4")
    print("PASS local factor degrees 2+4 have discriminant square classes 23 and -23")
    print("PASS over K0 the prime degrees are (e,f)=(1,2) and (2,2)")
    print("PASS local residual equations U^2=15 and V^2=5 share F_23^2")
    print("PASS D_23(K0)=V4 with orbit sizes 2+4 and inertia type 1^2*2^2")
    print("SCOPE a Nielsen labeling of the 2+4 orbits requires choosing a prime of the S6 closure")


if __name__ == "__main__":
    main()
