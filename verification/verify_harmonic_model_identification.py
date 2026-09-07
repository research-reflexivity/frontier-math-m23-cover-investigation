#!/usr/bin/env sage-python
"""Identify the independently reconstructed harmonic map exactly.

This is a SECOND, comparison certificate: unlike the construction it
deliberately reads the earlier certified model. It gives a projective
isomorphism and checks the rational functions, transferring the existing
M23 monodromy certificate. It does not replay Arb continuation.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path
from sage.all import PolynomialRing, PowerSeriesRing, QQ, matrix, vector, prod

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"notes"))
import certify_harmonic_recognized_model as new
K=new.K;P=new.P;xs=list(P.gens())


def rat(v):return QQ(v["numerator"])/QQ(v["denominator"])


def decode(v,sign=1):return K(rat(v["rational_part"]))+sign*new.pi*rat(v["sqrt_minus_23_part"])


def load_old(sign):
    can=json.loads((ROOT/"data/hurwitz_canonical_models_candidate.json").read_text())
    marks=json.loads((ROOT/"data/hurwitz_marked_points_candidate.json").read_text())
    equations=[]
    for name in ["quadric","petri_cubic"]:
        r=can[name]
        equations.append(sum(decode(a["degree_one_component"],sign)*prod(xs[i] for i in m)
                             for a,m in zip(r["coefficients"],r["monomials"])))
    points=[list(map(K,marks["b_coordinates"])),
            [decode(a,sign) for a in marks["c_coordinates"]["degree_one_component"]]]
    rec=json.loads((ROOT/"data/hurwitz_degree23_maps_candidate.json").read_text())
    rec=next(r for r in rec["components"] if r["component"]=="degree_one")
    sections=[]
    for name in ["numerator","denominator"]:
        sections.append(sum((rat(r["coefficient_power_basis"][0])+sign*new.pi*rat(r["coefficient_power_basis"][1]))
                            *prod(x**e for x,e in zip(xs,r["monomial"])) for r in rec[name]))
    return equations,points,sections


def jets(equations,point):
    for parameter in [1,2,3]:
        solve=[i for i in [1,2,3] if i!=parameter]
        jac=matrix(K,[[f.derivative(xs[i])(*point) for i in solve] for f in equations])
        if not jac.det():continue
        PS=PowerSeriesRing(K,"t",default_prec=4);t=PS.gen()
        values=[PS(a).add_bigoh(4) for a in point];values[parameter]+=t
        for n in range(1,4):
            err=vector(K,[f(*values)[n] for f in equations])
            delta=jac.solve_right(-err)
            for i,a in zip(solve,delta):values[i]+=a*t**n
        assert all(f(*values).valuation()>=4 for f in equations)
        return values
    raise AssertionError("Singular marked point")


def identify(sign):
    equations,points,sections=load_old(sign);q,c=equations
    jj=[jets(equations,p) for p in points]
    forms=[]
    for i in range(4):
        rows=[[f[n] for f in jj[0]] for n in range(i)]
        rows += [[f[n] for f in jj[1]] for n in range(3-i)]
        ker=matrix(K,rows).right_kernel().basis_matrix();assert ker.nrows()==1
        forms.append(ker.row(0))
    L=matrix(K,forms);assert L.det()
    old=list(L.inverse()*vector(P,xs))
    fq=P(q(*old));fc=P(c(*old))
    x0,x1,x2,x3=xs
    a,b,d,e=[fq.monomial_coefficient(m) for m in [x0*x3,x0*x2,x1*x2,x1*x1]]
    d1=b/d;kappa=1/(e*d1*d1);d2=1/(kappa*b);d3=1/(kappa*a)
    diagonal=matrix.diagonal(K,[1,d1,d2,d3])
    change=L.inverse()*diagonal
    sub=list(change*vector(P,xs))
    fq=P(kappa*q(*sub));fc=P(c(*sub)).reduce([fq])
    fc/=fc.monomial_coefficient(x0*x0*x3)
    if fq!=new.q or fc!=new.c:return None
    assert P(c(*sub)).reduce([new.q,new.c])==0
    for j,point in zip([0,3],points):
        mapped=change.inverse()*vector(K,point)
        assert mapped[j] and all(mapped[k]==0 for k in range(4) if k!=j)
    return change,[P(f(*sub)) for f in sections]


def main(args):
    result=None
    for sign in [1,-1]:
        result=identify(sign)
        if result is not None:break
    assert result is not None,"Neither coefficient embedding identifies the marked curves"
    change,(oldN,oldD)=result
    stored=json.loads((ROOT/"notes/harmonic_independent_exact_model_certificate.json").read_text())
    def form(rows):return P({tuple(r["monomial"]):new.decode(r["coefficient"]) for r in rows})
    N,D=form(stored["numerator"]),form(stored["denominator"])
    ideal=P.ideal(new.q,new.c)
    left=ideal.reduce(oldN*D);right=ideal.reduce(oldD*N);assert left and right
    mon=right.monomials()[0]
    scalar=left.monomial_coefficient(mon)/right.monomial_coefficient(mon)
    assert scalar and left==scalar*right
    print("PASS_EXACT_HARMONIC_MARKED_CURVE_ISOMORPHISM",flush=True)
    print("COEFFICIENT_EMBEDDING_SIGN",sign,flush=True)
    print("OLD_COORDINATES_FROM_INTRINSIC_COORDINATES",change,flush=True)
    print("PASS_EXACT_DEGREE23_MAP_IDENTITY_UP_TO_TARGET_SCALING",flush=True)
    # Validate the group and provenance of the source monodromy record.
    sys.path.insert(0,str(ROOT/"verification"))
    import verify_hurwitz_branch_cycle_summary as branch
    branch.main()
    print("PASS_HARMONIC_MONODROMY_M23_BY_EXACT_ISOMORPHISM",flush=True)
    record={"status":"PASS_EXACT_ISOMORPHISM_TRANSFERS_CERTIFIED_M23_MONODROMY",
            "coefficient_embedding_sign":sign,
            "old_coordinates_from_new":[[new.encode(a) for a in row] for row in change.rows()],
            "old_map_equals_scalar_times_new":new.encode(scalar),
            "new_arb_continuation_run":False,
            "source_sha256":{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()
                for path in [ROOT/"data/hurwitz_canonical_models_candidate.json",
                             ROOT/"data/hurwitz_marked_points_candidate.json",
                             ROOT/"data/hurwitz_degree23_maps_candidate.json",
                             ROOT/"notes/harmonic_independent_exact_model_certificate.json",
                             ROOT/"verification/hurwitz_branch_cycle_summary.json"]}}
    if args.output:Path(args.output).write_text(json.dumps(record,indent=2)+"\n")


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--output")
    main(p.parse_args())
