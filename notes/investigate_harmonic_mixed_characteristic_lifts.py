#!/usr/bin/env sage-python
"""Exact truncated mixed-characteristic degree-23-function lifting tests.

Arithmetic is in O/(pi^N), pi^2=-23, using rational number-field elements
reduced after operations. This is not equal-characteristic k[pi]/pi^N
when N>2. No reconstructed cover data are read.
"""
import argparse
import importlib.util
from pathlib import Path
from sage.all import GF, QQ, NumberField, PolynomialRing, matrix, vector

spec=importlib.util.spec_from_file_location("first_order",Path(__file__).with_name(
    "investigate_harmonic_first_order_hurwitz_constraints.py"))
first=importlib.util.module_from_spec(spec);spec.loader.exec_module(first)
k=GF(23)
K=NumberField(PolynomialRing(QQ,"p").gen()**2+23,"pi")
pi=K.gen()
S=PolynomialRing(K,"z");z=S.gen()
P=PolynomialRing(K,["H","U","V","W"])
H,U,V,W=P.gens()
basis=[P({mon:K(int(a)) for mon,a in f.dict().items()}) for f in first.basis]
q0=H*W-U*U-U*V-V*V
c0=U*W*W-V**3


def run(q,c,marks,precision=3,label="test",verbose=True,diagnostic_candidate=False,
        series_count=26,analytic_candidate=False):
    """Solve the global multiplier equations, or a named off-locus candidate.

    ``analytic_candidate`` solves only the fixed fifteen invertible rows
    at each digit, even if other rows failed earlier. Its returned ratios
    are NOT asserted to glue to a global function. The raw relations are
    returned so callers can impose the missing equations themselves.
    Default and final-digit diagnostic behaviour remain unchanged.
    """
    assert precision>=2
    def rational_residue(a,modulus):
        if modulus==1:return 0
        a=QQ(a)
        assert a.denominator()%23
        value=(int(a.numerator())*pow(int(a.denominator()),-1,modulus))%modulus
        return value if value<=modulus//2 else value-modulus
    def cutcoef(a):
        coordinates=list(K(a))+[QQ.zero(),QQ.zero()]
        return K(rational_residue(coordinates[0],23**((precision+1)//2)))+pi*K(
            rational_residue(coordinates[1],23**(precision//2)))
    def cut(f):return S({j:cutcoef(a) for j,a in S(f).dict().items()})
    assert cutcoef(pi*pi)==cutcoef(-23)
    if precision>=3:assert cutcoef(23)!=0
    def digit(a,n):
        coordinates=list(K(a)/pi**n)+[QQ.zero(),QQ.zero()]
        assert all(not x or x.valuation(23)>=0 for x in coordinates[:2])
        return k(rational_residue(coordinates[0],23))
    def digits(f,n,width):return [digit(f[j],n) for j in range(width)]
    def powers(f,d):
        answer=[S.one()]
        for _ in range(d):answer.append(cut(answer[-1]*f))
        return answer
    def evaluate(form,hh,uu):
        hp,up=powers(hh,5),powers(uu,5)
        total=S.zero()
        for (a,b,d,e),coefficient in form.dict().items():
            total=cut(total+cut(coefficient*cut(hp[a]*up[b])*z**d))
        return total
    qtail,ctail=q-q0,c-c0
    hh,uu=z**6+z**4+z*z,z**3
    for _ in range(2*precision+1):
        uu=cut(z**3-evaluate(ctail,hh,uu))
        hh=cut(uu*uu+uu*z+z*z-evaluate(qtail,hh,uu))
    assert evaluate(q,hh,uu)==evaluate(c,hh,uu)==0
    values=[evaluate(f,hh,uu) for f in basis]
    sections=[]
    section_forms=[]
    for mark,data in zip(marks,first.jet_data):
        _,_,kernel,pivots,inverse,_=data
        shifted=[cut(f(z+mark)) for f in values]
        jets=matrix(K,[[f[j] for f in shifted] for j in range(23)])
        current=[]
        current_forms=[]
        for row in kernel:
            coefficients=vector(K,[K(int(a)) for a in row])
            for n in range(1,precision):
                residual=jets*coefficients
                corrections=inverse*(-vector(k,[digit(a,n) for a in residual]))
                for j,a in zip(pivots,corrections):coefficients[j]+=pi**n*int(a)
                coefficients=vector(K,[cutcoef(a) for a in coefficients])
            assert all(cutcoef(a)==0 for a in jets*coefficients)
            current.append(cut(sum(a*f for a,f in zip(coefficients,values))))
            homogeneous=sum(a*f for a,f in zip(coefficients,basis))
            current_forms.append(P({mon:cutcoef(a) for mon,a in homogeneous.dict().items()}))
        sections.append(current)
        section_forms.append(current_forms)
    SA,SB=sections
    width=60+15*(precision-1)+1
    pairs=first.pairs
    def flatten0(polys):
        assert all(f.degree()<width for f in polys)
        return vector(k,[k(f[j]) for f in polys for j in range(width)])
    columns=[]
    for row in range(4):
        for column in range(4):
            polys=[(first.SA0[a]*first.SB0[column] if b==row else first.R.zero())
                   -(first.SA0[b]*first.SB0[column] if a==row else first.R.zero())
                   for a,b in pairs]
            columns.append(flatten0(polys))
    system=matrix(k,columns).transpose()
    assert system.rank()==15
    pivot_columns=system.pivots()
    injective=system.matrix_from_columns(pivot_columns)
    pivot_rows=injective.transpose().pivots()
    inverse=injective.matrix_from_rows(pivot_rows).inverse()
    multiplier=matrix(K,4,4,[K(int(a)) for a in first.M0.list()])
    def denominators():
        return [cut(sum(a*f for a,f in zip(row,SB))) for row in multiplier.rows()]
    def relations():
        DD=denominators()
        return [cut(SA[a]*DD[b]-SA[b]*DD[a]) for a,b in pairs]
    failure=None
    selected_labels=[divmod(j,width) for j in pivot_rows]
    for n in range(1,precision):
        rr=relations()
        assert all(f.degree()<width for f in rr)
        if analytic_candidate:
            # Only these rows are known to vanish at every preceding
            # order. Dividing the other rows by pi^n would be invalid.
            rhs_selected=-vector(k,[digit(rr[pair][degree],n)
                                    for pair,degree in selected_labels])
            correction=inverse*rhs_selected
            for j,a in zip(pivot_columns,correction):
                row,column=divmod(j,4)
                multiplier[row,column]=cutcoef(multiplier[row,column]+pi**n*int(a))
            continue
        rhs=-vector(k,[a for f in rr for a in digits(f,n,width)])
        correction=inverse*vector(k,[rhs[j] for j in pivot_rows])
        obstruction=rhs-injective*correction
        if obstruction:
            signature=[(j,int(a)) for j,a in enumerate(obstruction) if a]
            if verbose:
                print("MIXED_LIFT",label,"FAIL_AT_PI_POWER",n,
                      "nonzero_obstruction_coordinates",len(signature),
                      "first_coordinates",signature[:8],flush=True)
            failure={"pass":False,"order":n,"obstruction":obstruction}
            if not diagnostic_candidate:return failure
            assert n==precision-1 and n>=2
        for j,a in zip(pivot_columns,correction):
            row,column=divmod(j,4)
            multiplier[row,column]=cutcoef(multiplier[row,column]+pi**n*int(a))
    raw_relations=relations()
    if analytic_candidate:
        assert all(cutcoef(raw_relations[pair][degree])==0
                   for pair,degree in selected_labels)
    elif failure is None:assert all(f==0 for f in raw_relations)
    else:assert any(f!=0 for f in raw_relations)
    DD=denominators()
    chosen=next(j for j,f in enumerate(DD) if digit(f[0],0))
    NN,denominator=SA[chosen],DD[chosen]
    # A local power series at e=0, retaining the z^23 term whose
    # derivative has coefficient 23=-pi^2.
    count=series_count
    coefficients=[]
    for j in range(count):
        coefficients.append(cutcoef((NN[j]-sum(denominator[r]*coefficients[j-r]
                             for r in range(1,j+1)))/denominator[0]))
    unit=coefficients[0]
    coefficients=[cutcoef(a/unit) for a in coefficients]
    beta_series=S(coefficients)
    eta=[];eta1=[];eta_lift=S.zero()
    if precision>=3:
        # Divide beta'/beta by pi^2 coefficientwise, then by -1 because
        # pi^2=-23. The derivative vanishes modulo pi^2 on actual lifts.
        derivative=[cutcoef((j+1)*coefficients[j+1]) for j in range(count-1)]
        quotient=[]
        for j in range(count-1):
            quotient.append(cutcoef(derivative[j]-sum(coefficients[r]*quotient[j-r]
                                                     for r in range(1,j+1))))
        assert all(digit(a,0)==digit(a,1)==0 for a in quotient)
        eta=[-digit(a,2) for a in quotient]
        eta_lift=S([-a/pi**2 for a in quotient])
        if precision>=4:
            # Remove the chosen integral lift of the leading coefficient
            # before extracting the next digit of the divided differential.
            eta1=[-digit(a-pi**2*int(digit(a,2)),3) for a in quotient]
    if verbose and failure is None:
        print("MIXED_LIFT",label,"ANALYTIC_CANDIDATE_MOD_PI" if analytic_candidate else "PASS_MOD_PI",precision,
              "ETA0",first.R(eta) if eta else None,flush=True)
    if analytic_candidate:
        hN=section_forms[0][chosen]/unit
        hD=sum(a*f for a,f in zip(multiplier.row(chosen),section_forms[1]))
        hN=P({mon:cutcoef(a) for mon,a in P(hN).dict().items()})
        hD=P({mon:cutcoef(a) for mon,a in P(hD).dict().items()})
        return {"global_relations_vanish":all(f==0 for f in raw_relations),
                "beta":beta_series,"eta_lift":eta_lift,
                "affine_pair":(cut(NN/unit),denominator),
                "homogeneous_pair":(hN,hD),
                "curve_affine":(hh,uu),"multiplier":multiplier,
                "raw_relations":raw_relations,
                "selected_relation_labels":selected_labels}
    if failure is not None:
        # These coefficients are only a diagnostic polynomial extension of
        # the actual invariant off its zero locus. They are NOT the
        # differential of a globally lifted function when the test fails.
        failure.update({"candidate_eta0":eta,"candidate_eta1":eta1,
                        "candidate_eta_lift":eta_lift,
                        "candidate_affine_pair":(cut(NN/unit),denominator),
                        "curve_affine":(hh,uu)})
        return failure
    homogeneous_denominators=[]
    for row in multiplier.rows():
        f=sum(a*g for a,g in zip(row,section_forms[1]))
        homogeneous_denominators.append(P({mon:cutcoef(a) for mon,a in f.dict().items()}))
    return {"pass":True,"beta":beta_series,"eta0":eta,"eta1":eta1,
            "eta_lift":eta_lift,"affine_pair":(cut(NN/unit),denominator),
            "curve_affine":(hh,uu),"multiplier":multiplier,
            "obstruction":vector(k,len(pairs)*width),
            "homogeneous_pairs":list(zip(section_forms[0],homogeneous_denominators))}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--precision",type=int,default=3)
    parser.add_argument("--lambda-values",type=int,nargs="+",default=[0,1,-1,2])
    args=parser.parse_args()
    results={}
    for value in args.lambda_values:
        results[value]=run(q0,c0+pi*value*H**3,[K(1),K(-1)],args.precision,
                          "c_plus_pi_%s_H3"%value)
    if args.precision==3 and all(v in results for v in [0,1,-1,2]):
        assert all(not r["pass"] and r["order"]==2 for r in results.values())
        constant=results[0]["obstruction"]
        quadratic=results[1]["obstruction"]-constant
        assert all(r["obstruction"]==constant+k(v)**2*quadratic
                   for v,r in results.items())
        assert (constant[8],quadratic[8])==(7,k(-3))
        assert (constant[10],quadratic[10])==(12,5)
        L=PolynomialRing(k,"lam");lam=L.gen()
        ob1,ob2=7-3*lam**2,12+5*lam**2
        assert 5*ob1+3*ob2==2
        print("PURE_H3_SECOND_ORDER_OBSTRUCTIONS",ob1,ob2,flush=True)
        print("OBSTRUCTION_UNIT_IDENTITY",5*ob1+3*ob2,flush=True)
        print("PASS_NO_PURE_H3_SECOND_ORDER_TORSION_LIFT",flush=True)
    print("SCOPE: chosen deformation family; higher curve coefficients not solved")


if __name__=="__main__":main()
