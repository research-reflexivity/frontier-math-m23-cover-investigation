#!/usr/bin/env sage-python
"""Exact tests of the model recognized from harmonic local jets.

This model was independently recovered in the intrinsic two-flag frame.
No data/ characteristic-zero model or stored map is an input. Rational
recognition proposes these coefficients; only the exact tests below
certify the assertions printed by the script.
"""
import argparse
import json
import time
from pathlib import Path
from sage.all import QQ, NumberField, PolynomialRing, PowerSeriesRing, matrix, vector
from harmonic_truncated_geometry import Arithmetic

K=NumberField(PolynomialRing(QQ,"p").gen()**2+23,"pi");pi=K.gen()
P=PolynomialRing(K,["x0","x1","x2","x3"],order="lex")
x0,x1,x2,x3=P.gens();xs=list(P.gens())
q=x0*x2+x0*x3+x1*x1+x1*x2-QQ(127)/229*x1*x3+QQ(12192)/52441*x2*x2
c=(x0*x0*x3+(69+41*pi)/916*x0*x1*x3
   -(288+3552*pi)/52441*x0*x3*x3
   -(5+pi)/458*x1**3-(3061+569*pi)/209764*x1*x1*x2
   -(8649+1515*pi)/104882*x1*x1*x3
   -(298479+49803*pi)/48035956*x1*x2*x2
   -(1041970+181274*pi)/12008989*x1*x2*x3
   +(2318419+389759*pi)/48035956*x1*x3*x3
   -(2350008+338328*pi)/2750058481*x2**3
   -(55642056+9354216*pi)/2750058481*x2*x2*x3)


def log(*args):print(*args,flush=True)


def encode(a):return [str(QQ(b)) for b in list(K(a))]


def decode(a):return K([QQ(b) for b in a])


def certificate_curve():
    gram=matrix(K,4,4,lambda i,j:q.derivative(xs[i]).derivative(xs[j]))
    assert gram.det()
    log("EXACT_QUADRIC_DETERMINANT",gram.det())
    minors=[q.derivative(xs[i])*c.derivative(xs[j])-q.derivative(xs[j])*c.derivative(xs[i])
            for i in range(4) for j in range(i+1,4)]
    for j in range(4):
        ideal=P.ideal([q,c,xs[j]-1,*minors])
        assert ideal.groebner_basis()==[P.one()]
        log("EXACT_SMOOTH_CHART",j)
    assert q(1,0,0,0)==c(1,0,0,0)==0
    assert q(0,0,0,1)==c(0,0,0,1)==0
    qline=q(0,x1,x2,0);cline=c(0,x1,x2,0)
    quotient=-(5+pi)/458*x1-(771+111*pi)/209764*x2
    assert cline==qline*quotient
    assert QQ(1)-4*QQ(12192)/52441==QQ(3673)/52441
    log("EXACT_INCIDENCE_FACTORIZATION",quotient)
    log("EXACT_INCIDENCE_BINARY_DISCRIMINANT",QQ(3673)/52441)


def local_jets(which,count=24):
    S=PowerSeriesRing(K,"t",default_prec=count);t=S.gen()
    if which=="A":
        point=[S(1),t,S(0),S(0)];unknown=[2,3]
    else:
        point=[S(0),S(0),t,S(1)];unknown=[0,1]
    accuracy=1
    while accuracy<count:
        accuracy=min(count,2*accuracy)
        values=[S(f(*point)).add_bigoh(accuracy) for f in [q,c]]
        J=matrix(S,[[S(f.derivative(xs[j])(*point)).add_bigoh(accuracy)
                     for j in unknown] for f in [q,c]])
        det=J.det();assert det[0]
        delta=[(J[1,1]*values[0]-J[0,1]*values[1])/det,
               (-J[1,0]*values[0]+J[0,0]*values[1])/det]
        for j,a in zip(unknown,delta):
            # Restore an exact polynomial representative before the next
            # Newton step; unknown higher coefficients are set to zero.
            point[j]=S((point[j]-a).truncate())
    assert all(S(f(*point)).truncate(count)==0 for f in [q,c])
    point=[f.add_bigoh(count) for f in point]
    assert [f.valuation() for f in point]==([0,1,2,3] if which=="A" else [3,2,1,0])
    log("EXACT_MARKED_VANISHING_SEQUENCE",which,[f.valuation() for f in point])
    return point


def standard_monomials(gb,degree):
    leaders=[tuple(f.exponents()[0]) for f in gb]
    mons=(sum(xs)**degree).monomials()
    return [m for m in mons if not any(all(a>=b for a,b in zip(m.exponents()[0],e)) for e in leaders)]


def exact_map():
    gb=list(P.ideal(q,c).groebner_basis())
    quintics=standard_monomials(gb,5);decics=standard_monomials(gb,10)
    assert len(quintics)==27 and len(decics)==57
    sections=[]
    for label in ["A","B"]:
        points=local_jets(label)
        values=[f(*points) for f in quintics]
        jets=matrix(K,[[v[j] for v in values] for j in range(23)])
        kernel=jets.right_kernel().basis_matrix();assert kernel.nrows()==4
        ss=[sum(a*m for a,m in zip(row,quintics)) for row in kernel.rows()]
        assert all(jets*row==0 for row in kernel.rows())
        assert any(f(*points)[23] for f in ss)
        mark_ideal=P.ideal([x1,x2,x3] if label=="A" else [x0,x1,x2])
        common=P.ideal([q,c,*ss])
        assert common.saturation(mark_ideal)[0].groebner_basis()==[P.one()]
        log("EXACT_RESIDUAL_BASEPOINT_FREE",label)
        sections.append(ss);log("EXACT_QUINTIC_KERNEL",label,4)
    SA,SB=sections
    products=[[P(a*b).reduce(gb) for b in SB] for a in SA]
    coordinates=[[vector(K,[f.monomial_coefficient(m) for m in decics]) for f in row] for row in products]
    pairs=[(0,1),(0,2),(0,3)];columns=[]
    for row in range(4):
        for col in range(4):
            column=[]
            for a,b in pairs:
                v=vector(K,57)
                if row==b:v+=coordinates[a][col]
                if row==a:v-=coordinates[b][col]
                column.extend(v)
            columns.append(column)
    system=matrix(K,columns).transpose()
    kernel=system.right_kernel().basis_matrix()
    assert kernel.nrows()==1
    M=matrix(K,4,4,kernel.row(0));assert M.det()
    DD=[sum(a*f for a,f in zip(row,SB)) for row in M.rows()]
    for a in range(4):
        for b in range(a+1,4):assert P(SA[a]*DD[b]-SA[b]*DD[a]).reduce(gb)==0
    log("EXACT_GLOBAL_MULTIPLIER_RANK",15)
    log("EXACT_ALL_SIX_CROSS_PRODUCTS",True)
    return gb,SA[0],DD[0],SA,SB


def exact_ramification(N,D):
    """Critical algebra away from the two total-ramification fibres.

    Eight distinct critical points in this chart exhaust the remaining
    Riemann--Hurwitz degree. Thus no claim about an omitted chart is
    needed, and equality of their critical values proves three-pointness.
    """
    from itertools import product
    R=PolynomialRing(K,["a","b","d"],order="degrevlex")
    a,b,d=R.gens();aff=[a,R.one(),b,d]
    Q,C=map(R,[q(*aff),c(*aff)]);NN,DD=map(R,[N(*aff),D(*aff)])
    rows=[[f.derivative(z) for z in R.gens()] for f in [Q,C]]
    rows.append([DD*NN.derivative(z)-NN*DD.derivative(z) for z in R.gens()])
    critical=matrix(R,rows).det()
    log("EXACT_CRITICAL_SATURATION_START")
    I=R.ideal(Q,C,critical).saturation(R.ideal(NN*DD))[0]
    assert I.vector_space_dimension()==8
    log("EXACT_CRITICAL_ALGEBRA_DIMENSION",8)
    nr,dr=I.reduce(NN),I.reduce(DD);assert dr
    mon=dr.monomials()[0]
    value=nr.monomial_coefficient(mon)/dr.monomial_coefficient(mon)
    assert value and nr==value*dr
    log("EXACT_SINGLE_THIRD_BRANCH_VALUE",value)
    # A squarefree characteristic polynomial of a multiplication
    # operator on this eight-dimensional algebra certifies etaleness.
    gb=list(I.groebner_basis());leaders=[tuple(f.exponents()[0]) for f in gb]
    exponents=[e for e in product(range(8),repeat=3)
               if not any(all(x>=y for x,y in zip(e,lm)) for lm in leaders)]
    assert len(exponents)==8
    basis=[R.monomial(*e) for e in exponents]
    for element in [b,a+b+d,a+2*b+3*d]:
        columns=[I.reduce(element*f) for f in basis]
        multiplication=matrix(K,[[f.monomial_coefficient(m) for f in columns] for m in basis])
        polynomial=multiplication.charpoly()
        if polynomial.gcd(polynomial.derivative()).degree()==0:break
    else:raise AssertionError("No certified separating coordinate")
    log("EXACT_CRITICAL_SEPARATING_COORDINATE",element)
    log("EXACT_CRITICAL_CHARACTERISTIC_POLYNOMIAL",polynomial)
    log("PASS_EXACT_THREE_POINT_PASSPORT",["23","23","2^8 1^7"])
    return value,polynomial


def check_recognition(path):
    data=json.loads(Path(path).read_text());A=Arithmetic(K,data["precision"])
    for key,f in [("quadric",q),("cubic",c)]:
        observed=P({tuple(row["monomial"]):decode(row["coefficient"]) for row in data[key]})
        assert all(A.coefficient(a)==0 for a in (f-observed).coefficients())
    log("EXACT_MODEL_MATCHES_LOCAL_FRAME_MOD_PI",data["precision"])


def check_local_slice(path):
    from investigate_harmonic_intrinsic_frame import frame
    data=json.loads(Path(path).read_text());assert data["last_correction_verified"]
    precision=data["precision"]+3
    point=[decode(a) for a in data["point"]]
    g,_,fq,fc,L,diag=frame(point,precision)
    change=matrix.diagonal(K,[1/a for a in diag])*L
    sub=list(change*vector(g.P,[g.H,g.U,g.V,g.W]))
    assert g.normal_form(g.P(q(*sub)))==0
    assert g.normal_form(g.P(c(*sub)))==0
    for label,mark,axis in [("A",g.marks[0],0),("B",g.marks[1],3)]:
        hh,uu,_=g.local_jets(mark,count=3)
        mapped=change*vector(K,[hh[0],uu[0],mark,1])
        assert g.A.digit(mapped[axis],0)
        assert all(g.A.coefficient(mapped[j])==0 for j in range(4) if j!=axis)
        log("EXACT_MARKED_FRAME_MATCH_MOD_PI",label,precision)
    log("EXACT_MODEL_PULLBACK_MATCHES_HARMONIC_SLICE_MOD_PI",precision)


def conjugate(a):
    coefficients=list(K(a))
    return K([coefficients[0],-coefficients[1]])


def conjugate_form(f):
    return P({m:conjugate(a) for m,a in f.dict().items()})


def exact_descent(gb,N,D,value):
    """Explicit descent to Q with a conjugate pair of wild branch values.

    This proves arithmetic descent of this exact map, not its M23
    monodromy. The rational target coordinate does not keep both wild
    branch values individually rational.
    """
    a=-QQ(127)/229;b=QQ(12192)/52441
    T=matrix(K,[[0,0,0,a*b],[0,0,b,0],[0,1,0,0],[1/a,0,0,0]])
    transformed=list(T*vector(P,xs))
    assert T*T==b*matrix.identity(K,4)
    assert P(q(*transformed))==b*q
    cc=P(conjugate_form(c)(*transformed)).reduce([q])
    mu=cc.monomial_coefficient(x0*x0*x3)
    assert P(cc-mu*c).reduce([q])==0 and mu.norm()==b**3
    NN=P(conjugate_form(N)(*transformed));DD=P(conjugate_form(D)(*transformed))
    assert P(N*NN-value*conjugate(value)*D*DD).reduce(gb)==0
    log("EXACT_BRANCH_EXCHANGE_CUBIC_SCALAR",mu)
    log("EXACT_NORMALIZED_MAP_INVERSION",True)
    R=(b/mu)*T
    assert R*R.apply_map(conjugate)==matrix.identity(K,4)
    columns=[]
    for j in range(4):
        for scalar in [K.one(),pi]:
            v=vector(K,4);v[j]=scalar
            diff=R*v.apply_map(conjugate)-v
            columns.append([a for coefficient in diff for a in list(coefficient)])
    kernel=matrix(QQ,columns).transpose().right_kernel().basis_matrix()
    assert kernel.nrows()==4
    B=matrix(K,[list(vector(K,[row[2*j]+pi*row[2*j+1] for j in range(4)]))
                for row in kernel.rows()]).transpose()
    assert B.det() and R*B.apply_map(conjugate)==B
    sub=list(B*vector(P,xs))
    rq=P(q(*sub));rq=rq/rq.leading_coefficient()
    assert all(a==conjugate(a) for a in rq.coefficients())
    rc=P(c(*sub)).reduce([rq]);rc=rc/rc.leading_coefficient()
    assert all(a==conjugate(a) for a in rc.coefficients())
    rational_gb=list(P.ideal(rq,rc).groebner_basis())
    numerator=P((N-value*D)(*sub))
    denominator=P((pi*(N+value*D))(*sub))
    rn=P(numerator*conjugate_form(denominator)).reduce(rational_gb)
    rd=P(denominator*conjugate_form(denominator)).reduce(rational_gb)
    assert rn and rd
    assert all(a==conjugate(a) for a in [*rn.coefficients(),*rd.coefficients()])
    log("EXACT_RATIONAL_DESCENT",True)
    log("RATIONAL_TARGET_BRANCH_LOCUS","t=0 and 23*t^2+1=0")
    return {"quadric":serialize(rq),"cubic":serialize(rc),
            "numerator":serialize(rn),"denominator":serialize(rd),
            "source_change_matrix_rows":[[encode(a) for a in row] for row in B.rows()]}


def serialize(f):
    return [{"monomial":list(m),"coefficient":encode(a)} for m,a in f.dict().items()]


def main(args):
    start=time.monotonic();certificate_curve()
    if args.frame:check_recognition(args.frame)
    if args.checkpoint:check_local_slice(args.checkpoint)
    output={"scope":"exact smooth curve and osculating incidence; map not yet checked",
            "quadric":serialize(q),"cubic":serialize(c)}
    if args.map or args.ramification:
        gb,N,D,SA,SB=exact_map()
        output.update({"scope":"exact smooth curve, incidence, and global cross-product identities",
                       "numerator":serialize(N),"denominator":serialize(D),
                       "sections_A":[serialize(f) for f in SA],"sections_B":[serialize(f) for f in SB]})
        if args.ramification:
            value,poly=exact_ramification(N,D)
            output.update({"scope":"exact smooth three-point map with passport (23),(23),(2^8,1^7), and exact osculating incidence; M23 monodromy not certified here",
                           "third_branch_value":encode(value),
                           "critical_polynomial":[encode(a) for a in poly.list()]})
            if args.descent:
                output["rational_descent"]=exact_descent(gb,N,D,value)
    if args.output:Path(args.output).write_text(json.dumps(output,indent=2)+"\n")
    log("ELAPSED_SECONDS",round(time.monotonic()-start,2))


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--map",action="store_true");p.add_argument("--output")
    p.add_argument("--ramification",action="store_true");p.add_argument("--frame")
    p.add_argument("--descent",action="store_true")
    p.add_argument("--checkpoint")
    args=p.parse_args()
    if args.descent:args.ramification=True
    main(args)
