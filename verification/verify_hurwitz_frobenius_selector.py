#!/usr/bin/env sage-python
"""Check the finite-field application of the local component criterion.

The theorem's lifting and finite-etale assertions are proved in main.tex.
This script checks the residue polynomial, Frobenius equalizer, quadratic
splitting, and an explicit counterexample to local-to-global descent.
It does NOT construct relative incidence correspondences, validate a
quadratic determinant line, or independently match Nielsen classes.
Run make verify-hurwitz-frobenius-selector to also recompute the local
decomposition and pointed reductions on which this application depends.
"""

from sage.all import GF, PolynomialRing, QQ


def main():
    ring = PolynomialRing(GF(23), "U")
    U = ring.gen()
    factors = [U - 16, U**2 + 1, U**2 + U + 1]
    assert all(f.is_irreducible() for f in factors)
    R23 = factors[0] * factors[1] * factors[2]
    assert R23.is_squarefree()
    B = ring.quotient(R23, "u")
    u = B.gen()
    e = 2*u**4 + 2*u**3 + 4*u**2 + 2*u + 3
    assert e**2 == e
    assert e == -(u**23 - u)**22
    assert (U**23 - U).gcd(R23).monic() == U - 16
    assert (U**(23**2) - U) % R23 == 0
    assert e == 1 - 21*(u**2 + 1)*(u**2 + u + 1)

    # In each residue factor, verify idempotent value and splitting of
    # z^2+1 directly, not merely the parity of the recorded residue degree.
    residue_degrees = []
    splitting = []
    for polynomial, expected in zip(factors, [0, 1, 1]):
        field = GF(23) if polynomial.degree() == 1 else GF(
            23**2, name="a", modulus=polynomial
        )
        value = field(16) if polynomial.degree() == 1 else field.gen()
        assert e.lift()(value) == expected
        zring = PolynomialRing(field, "z")
        z = zring.gen()
        quadratic = z**2 + 1
        linear_factors = sum(
            int(g.degree() == 1) * int(m) for g, m in quadratic.factor()
        )
        assert linear_factors == (0 if expected == 0 else 2)
        residue_degrees.append(int(field.degree()))
        splitting.append(int(linear_factors == 2))
    assert residue_degrees == [1, 2, 2]
    assert splitting == [0, 1, 1]

    # This list is supplied by verify_hurwitz_local_23.py; it is not a
    # fresh local-field factorization inside this finite-field script.
    ef = [(1, 1), (1, 2), (2, 2)]
    assert sum(ei*fi for ei, fi in ef) == 7
    assert sum(ei for ei, fi in ef if fi == 1) == 1
    assert sum(ei*fi for ei, fi in ef if fi % 2 == 0) == 6

    # A locally unique rank-one factor does not imply a global Q-point.
    qring = PolynomialRing(QQ, "X")
    X = qring.gen()
    assert (X**3 - 2).is_irreducible()
    mod5 = PolynomialRing(GF(5), "x")
    x = mod5.gen()
    assert x**3 - 2 == (x - 3)*(x**2 + 3*x + 4)
    assert (x**2 + 3*x + 4).is_irreducible()
    assert (x**3 - 2).is_squarefree()
    print("PASS: Frobenius fixed locus = linear residue component")
    print("PASS: -(u^23-u)^22 = e_sing; quadratic split locus = [0,1,1]")
    print("PASS: component ranks 1+6; local-to-global counterexample")
    print("Scope: finite-field application, not a specialization proof")


if __name__ == "__main__":
    main()
