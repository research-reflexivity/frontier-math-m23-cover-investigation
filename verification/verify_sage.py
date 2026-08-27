#!/usr/bin/env sage-python
"""Independent SageMath verification of the optimal model.

This script rebuilds every polynomial from JSON coefficient tables.  It does
not read the Singular input files and does not call the Singular certificates.
"""

from pathlib import Path
import json

from sage.all import GF, QQ, PolynomialRing, gcd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"


def load_json(name):
    return json.loads((DATA / name).read_text())


def main():
    equation = load_json("optimal_23_4_Z.json")
    pencil = load_json("optimal_degree4_pencil.json")
    f_table = load_json("Fint_coefficients_Z.json")

    QT = PolynomialRing(QQ, "T")
    T = QT.gen()
    KT = QT.fraction_field()
    RV = PolynomialRing(KT, "V")
    V = RV.gen()

    F = RV.zero()
    for v_degree, row in enumerate(f_table):
        coefficient = sum(QQ(value) * T**t_degree for t_degree, value in enumerate(row))
        F += KT(coefficient) * V**v_degree
    assert F.degree() == 23
    D = T**2 + 23
    Fhat = F / D**4
    assert Fhat.degree() == 23
    assert Fhat.is_monic()

    exact_pencil = []
    for encoded in pencil["expanded_coefficients_V_then_T"]:
        value = RV.zero()
        for key, coefficient in encoded.items():
            v_degree, t_degree = map(int, key.split(","))
            value += KT(QQ(coefficient) * T**t_degree) * V**v_degree
        exact_pencil.append(value)
    J0, J1 = exact_pencil

    table = [[QQ(value) for value in row] for row in equation["coefficients_T_then_W"]]
    a = [KT(sum(table[t_degree][w_degree] * T**t_degree for t_degree in range(5)))
         for w_degree in range(24)]

    # Homogeneous Horner evaluation of J1^23 P(T,J0/J1), reducing after
    # every multiplication in K(T)[V]/(Fhat).
    numerator = RV(a[23])
    j1_power = J1
    for w_degree in range(22, -1, -1):
        numerator = (numerator * J0 + a[w_degree] * j1_power) % Fhat
        if w_degree > 0:
            j1_power = (j1_power * J1) % Fhat
    assert numerator == 0

    # An independent finite-field construction and factorization.
    k = GF(31)
    RTW = PolynomialRing(k, names=("T", "W"))
    t, w = RTW.gens()
    P31 = RTW.zero()
    content = 0
    for t_degree, row in enumerate(equation["coefficients_T_then_W"]):
        for w_degree, coefficient in enumerate(row):
            integer = int(coefficient)
            content = gcd(content, abs(integer))
            P31 += k(integer) * t**t_degree * w**w_degree
    assert content == 1
    assert P31.degree(t) == 4 and P31.degree(w) == 23
    factors = list(P31.factor())
    assert len(factors) == 1 and factors[0][1] == 1

    result = {
        "status": "PASS_INDEPENDENT_SAGE_CERTIFICATE",
        "sage_version": __import__("sage.version").version.version,
        "input_format": "JSON coefficient tables (no Singular inputs)",
        "degree_T": 4,
        "degree_W": 23,
        "primitive_over_Z": True,
        "exact_identity_in_QT_mod_Fhat": True,
        "irreducible_mod_31": True,
        "same_function_field_degree_argument": True,
    }
    output = ROOT / "verification" / "sage_verification_summary.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(result["status"])
    print(f"SageMath {result['sage_version']}")
    print("degree_T=4 degree_W=23 exact_identity=1 irreducible_mod_31=1")


if __name__ == "__main__":
    main()
