#!/usr/bin/env python3
"""Export a standalone independent-arithmetic Magma check.

The supplied sections are witnesses, not trusted conclusions. Magma checks
their jets, base divisors, multiplier identities and critical algebra.
Generating this file does NOT constitute a Magma execution certificate.
"""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def scalar(pair):
    a, b = pair
    return f"(({a})+({b})*pi)"


def form(rows):
    return " +\n".join(scalar(row["coefficient"]) + "*" + "*".join(
        f"x{i}^{e}" for i, e in enumerate(row["monomial"]) if e)
        for row in rows) or "P!0"


CHECKS = r'''
I := ideal<P | q,c>;
xx := [x0,x1,x2,x3];
minors := [];
for i in [1..3] do
    for j in [i+1..4] do
        Append(~minors,Derivative(q,i)*Derivative(c,j)-Derivative(q,j)*Derivative(c,i));
    end for;
end for;
for j in [1..4] do
    assert 1 in ideal<P | [q,c,xx[j]-1] cat minors>;
end for;
assert Evaluate(c,[P!0,x1,x2,P!0]) eq
    Evaluate(q,[P!0,x1,x2,P!0]) *
    (-(5+pi)/458*x1-(771+111*pi)/209764*x2);
assert 1-4*Rationals()!12192/52441 eq Rationals()!3673/52441;
print "PASS_MAGMA_HARMONIC_SMOOTH_CURVE_AND_TWO_POINT_INCIDENCE";

function MarkJets(which,q,c)
    K := BaseRing(Parent(q)); S<t> := PowerSeriesRing(K,24);
    if which eq 1 then
        v := [S!1,t,S!0,S!0]; unknown := [3,4];
    else
        v := [S!0,S!0,t,S!1]; unknown := [1,2];
    end if;
    p := [Coefficient(f,0):f in v];
    J := Matrix(K,[[Evaluate(Derivative(f,j),p):j in unknown]:f in [q,c]]);
    assert Determinant(J) ne 0;
    for n in [1..23] do
        err := Vector(K,[Coefficient(Evaluate(f,v),n):f in [q,c]]);
        delta := -err*Transpose(J^-1);
        for j in [1..2] do v[unknown[j]] +:= delta[j]*t^n; end for;
    end for;
    assert &and[Valuation(Evaluate(f,v)) ge 24:f in [q,c]];
    assert [Valuation(f):f in v] eq (which eq 1 select [0,1,2,3] else [3,2,1,0]);
    return v;
end function;

quintics := Monomials((x0+x1+x2+x3)^5);
for which in [1..2] do
    ss := which eq 1 select SA else SB;
    vv := MarkJets(which,q,c);
    assert Rank(Matrix(K,[[MonomialCoefficient(NormalForm(f,I),m):m in quintics]:f in ss])) eq 4;
    assert &and[Valuation(Evaluate(f,vv)) ge 23:f in ss];
    assert &or[Coefficient(Evaluate(f,vv),23) ne 0:f in ss];
    mark := which eq 1 select ideal<P|x1,x2,x3> else ideal<P|x0,x1,x2>;
    assert 1 in Saturation(ideal<P|[q,c] cat ss>,mark);
end for;

products := [[NormalForm(a*b,I):b in SB]:a in SA];
decics := Monomials((x0+x1+x2+x3)^10);
columns := [];
for row in [1..4] do
    for col in [1..4] do
        column := [];
        for j in [2..4] do
            f := P!0;
            if row eq j then f +:= products[1][col]; end if;
            if row eq 1 then f -:= products[j][col]; end if;
            column cat:= [MonomialCoefficient(f,m):m in decics];
        end for;
        Append(~columns,column);
    end for;
end for;
system := Matrix(K,columns);
ker := Nullspace(system);
assert Dimension(ker) eq 1;
witness := Basis(ker)[1]; assert IsZero(witness*system);
M := Matrix(K,4,4,Eltseq(witness)); assert Determinant(M) ne 0;
DD := [&+[M[i,j]*SB[j]:j in [1..4]]:i in [1..4]];
for i in [1..3] do
    for j in [i+1..4] do assert SA[i]*DD[j]-SA[j]*DD[i] in I; end for;
end for;
left := NormalForm(N*DD[1],I); right := NormalForm(D*SA[1],I);
assert left ne 0 and right ne 0;
mon := Monomials(right)[1];
scale := MonomialCoefficient(left,mon)/MonomialCoefficient(right,mon);
assert scale ne 0 and left eq scale*right;
print "PASS_MAGMA_HARMONIC_DEGREE23_FUNCTION";

R<a,b,d> := PolynomialRing(K,3,"grevlex");
aff := [a,R!1,b,d];
Q := Evaluate(q,aff); C := Evaluate(c,aff);
NN := Evaluate(N,aff); Den := Evaluate(D,aff);
critical := Determinant(Matrix(R,3,3,
    [Derivative(Q,j):j in [1..3]] cat [Derivative(C,j):j in [1..3]] cat
    [Den*Derivative(NN,j)-NN*Derivative(Den,j):j in [1..3]]));
crit := Saturation(ideal<R|Q,C,critical>,ideal<R|NN*Den>);
assert QuotientDimension(crit) eq 8;
assert IsRadical(crit);
assert NN-u*Den in crit;
assert u ne 0;
print "PASS_MAGMA_HARMONIC_EIGHT_DISTINCT_EQUAL_CRITICAL_VALUES";
// Genus four and two order-23 fibres leave exactly eight in Riemann--Hurwitz.
print "PASS_MAGMA_HARMONIC_PASSPORT_23_23_2POWER8_1POWER7";
print "SCOPE: exact curve, incidence, function and passport; no local Hensel audit or new Arb run";
'''


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default=str(ROOT / "verification/verify_harmonic_model_generated.m"))
    parser.add_argument("--check", action="store_true", help="Check deterministic output without rewriting it")
    args = parser.parse_args()
    path = ROOT / "notes/harmonic_independent_exact_model_certificate.json"
    raw = path.read_bytes()
    data = json.loads(raw)
    pieces = ["// Generated witnesses; independent Magma execution is still required.\n",
              f"// Input SHA256: {hashlib.sha256(raw).hexdigest()}\n",
              "Qp<p> := PolynomialRing(Rationals());\nK<pi> := NumberField(p^2+23);\n",
              'P<x0,x1,x2,x3> := PolynomialRing(K,4,"lex");\n']
    for key, name in [("quadric", "q"), ("cubic", "c"), ("numerator", "N"), ("denominator", "D")]:
        pieces.append(f"{name} := {form(data[key])};\n")
    for key, name in [("sections_A", "SA"), ("sections_B", "SB")]:
        pieces.append(name + " := [\n" + ",\n".join(form(f) for f in data[key]) + "];\n")
    pieces.append("u := " + scalar(data["third_branch_value"]) + ";\n")
    pieces.append(CHECKS)
    output = Path(args.output)
    if args.check:
        if output.read_text() != "".join(pieces):
            raise SystemExit("Generated harmonic Magma input is stale")
        print("PASS harmonic Magma input matches the current exact witnesses; no Magma execution performed")
        return
    output.write_text("".join(pieces))
    print(f"Prepared {output} ({output.stat().st_size} bytes); GENERATION_ONLY_NOT_AN_EXECUTION")


if __name__ == "__main__":
    main()
