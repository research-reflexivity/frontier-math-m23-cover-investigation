#!/usr/bin/env sage-python
"""Recognize a canonical model in an intrinsic two-flag frame.

Each coordinate spans H0(omega(-i A-(3-i)B)), for i=0,...,3.
All input comes from the independently computed local Hensel residues.
Recognition remains conjectural until an exact verification succeeds.
"""
import argparse
import json
from pathlib import Path
from sage.all import GF, PolynomialRing, matrix, vector
from harmonic_finite_equations import fixed_family, mixed
from harmonic_truncated_geometry import Arithmetic
from investigate_harmonic_exact_recognition import decode, encode, recognize


def frame(point,precision):
    g=fixed_family(point,precision);A=g.A;k=GF(23)
    jets=[]
    for mark in g.marks:
        hh,uu,_=g.local_jets(mark,count=3)
        fs=[hh,uu,g.z+mark,g.S.one()]
        jets.append(matrix(g.K,[[f[d] for f in fs] for d in range(3)]))
    forms=[]
    for i in range(4):
        rows=list(jets[0].rows()[:i])+list(jets[1].rows()[:3-i])
        J=matrix(g.K,rows)
        special=matrix(k,[[A.digit(a,0) for a in row] for row in rows])
        pivots=special.pivots();assert len(pivots)==3
        free=next(j for j in range(4) if j not in pivots)
        inv=special.matrix_from_columns(pivots).inverse()
        v=vector(g.K,4);v[free]=1
        for n in range(precision):
            delta=-inv*vector(k,[A.digit(a,n) for a in J*v])
            for j,a in zip(pivots,delta):v[j]=A.coefficient(v[j]+g.pi**n*int(a))
        assert all(A.coefficient(a)==0 for a in J*v)
        forms.append(v)
    L=matrix(g.K,forms)
    assert matrix(k,[[A.digit(a,0) for a in row] for row in L.rows()]).det()
    change=L.inverse().apply_map(A.coefficient)
    P=PolynomialRing(g.K,["x0","x1","x2","x3"],order="lex")
    xs=list(P.gens());x0,x1,x2,x3=xs
    def cut(f):return P({m:A.coefficient(a) for m,a in P(f).dict().items()})
    old=[cut(sum(a*x for a,x in zip(row,xs))) for row in change.rows()]
    H,U,V,W=old
    q=cut(g.q(W,U,H,V));c=cut(g.c(W,U,H,V))
    assert all(q.monomial_coefficient(m)==0 for m in [x0*x0,x0*x1,x2*x3,x3*x3])
    # Diagonal substitution x_i=d_i*y_i, with d0=1, and an independent
    # equation scalar kappa. These four normalizations have determinant
    # one on the character lattice, so require no square roots.
    a,b,d,e=[q.monomial_coefficient(m) for m in [x0*x3,x0*x2,x1*x2,x1*x1]]
    print("FRAME_NORMALIZING_SPECIAL_COEFFICIENTS",[A.digit(v,0) for v in [a,b,d,e]],flush=True)
    assert all(A.digit(v,0) for v in [a,b,d,e])
    d1=A.coefficient(b/d)
    kappa=A.coefficient(1/(e*d1*d1))
    d2=A.coefficient(1/(kappa*b));d3=A.coefficient(1/(kappa*a))
    diag=[g.K(1),d1,d2,d3]
    q=cut(kappa*q(*[a*x for a,x in zip(diag,xs)]))
    c=cut(c(*[a*x for a,x in zip(diag,xs)]))
    assert all(q.monomial_coefficient(m)==1 for m in [x0*x3,x0*x2,x1*x2,x1*x1])
    # Lex leader x0*x2 is monic. Reduce the cubic modulo q and fix its
    # scalar using its first unit coefficient in this fixed monomial order.
    assert q.leading_coefficient()==1
    c=cut(c.reduce([q]))
    unit_mon=next(m for m in c.monomials() if A.digit(c.monomial_coefficient(m),0))
    c=cut(c/c.monomial_coefficient(unit_mon))
    print("FRAME_CUBIC_UNIT_MONOMIAL",unit_mon,flush=True)
    return g,P,q,c,L,diag


def main(args):
    data=json.loads(Path(args.checkpoint).read_text())
    assert data["last_correction_verified"]
    point=list(map(decode,data["point"]))
    # Every input coefficient and mark is known through this precision.
    precision=data["precision"]+3
    g,P,q,c,L,diag=frame(point,precision)
    print("INTRINSIC_FRAME_PRECISION",precision,flush=True)
    for label,f in [("quadric",q),("cubic",c)]:
        for mon in f.monomials():
            a=f.monomial_coefficient(mon)
            print(label,str(mon),encode(a),recognize(a,precision),flush=True)
    if args.output:
        out={"scope":"finite intrinsic-frame residues; not an exact model",
             "precision":precision,
             "quadric":[{"monomial":list(m),"coefficient":encode(a)} for m,a in q.dict().items()],
             "cubic":[{"monomial":list(m),"coefficient":encode(a)} for m,a in c.dict().items()]}
        Path(args.output).write_text(json.dumps(out,indent=2)+"\n")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--checkpoint",default="/private/tmp/m23-harmonic-local-recognition.json")
    p.add_argument("--output")
    main(p.parse_args())
