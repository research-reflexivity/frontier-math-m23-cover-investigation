#!/usr/bin/env sage-python
"""Audit the finite global-function and full effective-fibre equations.

The selected functions below use canonical-ring coefficients, not the
central-series coefficients of the preceding checker. Agreement on the
zero locus is tested independently. The accompanying geometric argument
is needed to justify eliminating redundant global-function conditions.
"""
import argparse
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from sage.all import GF,matrix,vector
from harmonic_finite_equations import fixed_family,mixed
from harmonic_truncated_geometry import Arithmetic
from certify_harmonic_universal_canonical_ring import standard_monomials

k=GF(23)


def certify_special_linear_systems():
    g=fixed_family([0]*11,6)
    mons=standard_monomials(g.P,5)
    normal=[g.normal_form(f) for f in g.basis]
    change=matrix(k,[[g.A.digit(f.monomial_coefficient(m),0) for f in normal] for m in mons])
    assert change.nrows()==change.ncols()==27 and change.det()
    for label,polys in [("A",mixed.first.gA),("B",mixed.first.gB)]:
        common=polys[0]
        for f in polys[1:]:common=common.gcd(f)
        assert common.degree()==0 and max(f.degree() for f in polys)==7
        print("SPECIAL_RESIDUAL_SYSTEM",label,"rank",4,"basepoint_free",True,flush=True)
    print("SPECIAL_QUINTIC_BASIS_DETERMINANT",change.det(),flush=True)


def evaluate(point,precision=6):
    g=fixed_family([mixed.K(a) for a in point],precision)
    residual=g.multiplier_candidate();R=g.ramification_factor();E=g.effective_fibre(R)
    g.certify_fibre_elimination(R,E)
    small=Arithmetic(mixed.K,precision-5)
    labels=[("fibre",3,15)]+[("fibre",4,j) for j in range(9,15)]
    labels += [("function",pair,str(mon)) for pair in range(3) for mon in g.decics]
    values=[]
    for _,order,j in labels[:7]:
        g.A.digit(E[j],order);values.append(small.coefficient(E[j]/mixed.pi**order))
    for a in residual:
        g.A.digit(a,5);values.append(small.coefficient(-a/mixed.pi**5))
    incidence=g.incidence()
    inc=[small.coefficient(incidence[j]/mixed.pi**3) for j in range(2)]
    return labels,values,inc,g


def sample(point):
    _,values,inc,_=evaluate(point)
    A=Arithmetic(mixed.K,1)
    return [int(A.digit(a,0)) for a in values],[int(A.digit(a,0)) for a in inc]


def samples(points,workers):
    if workers==1:return list(map(sample,points))
    with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context("spawn")) as pool:
        return list(pool.map(sample,points))


def run(workers,audit,lift_digits):
    certify_special_linear_systems()
    points=[[0]*11]
    for j in range(11):
        e=[0]*11;e[j]=1;points.append(e)
    values=samples(points,workers)
    base=vector(k,values[0][0]);inc0=vector(k,values[0][1])
    J=matrix(k,[vector(k,v[0])-base for v in values[1:]]).transpose()
    B=matrix(k,[vector(k,v[1])-inc0 for v in values[1:]]).transpose()
    labels,_,_,_=evaluate([0]*11)
    chosen=J.transpose().pivots();square=J.matrix_from_rows(chosen)
    constants=vector(k,[base[j] for j in chosen])
    function_jacobian=J[7:,:]
    print("FINITE_GLOBAL_FUNCTION_JACOBIAN_RANK",function_jacobian.rank(),flush=True)
    assert function_jacobian.rank()==4
    print("FINITE_COMPLETE_SELECTED_EQUATIONS",[labels[j] for j in chosen],flush=True)
    print("FINITE_COMPLETE_JACOBIAN",[list(row) for row in square.rows()],flush=True)
    print("FINITE_COMPLETE_CONSTANT",list(constants),flush=True)
    print("FINITE_COMPLETE_DETERMINANT",square.det(),flush=True)
    assert len(chosen)==11 and list(chosen[:7])==list(range(7)) and square.det()
    assert J.matrix_from_rows(chosen[7:]).rank()==4
    assert square.det()==11
    solution=square.solve_right(-constants)
    assert J*solution+base==0 and B*solution+inc0==0
    assert list(solution)==list(map(k,[0,0,16,0,22,21,0,10,0,16,16]))
    print("FINITE_COMPLETE_RESIDUAL_ROOT",list(solution),flush=True)
    if audit:
        points=[]
        for i in range(11):
            e=[0]*11;e[i]=-1;points.append(e)
            for j in range(i+1,11):
                e=[0]*11;e[i]=e[j]=1;points.append(e)
        for point,(val,inc) in zip(points,samples(points,workers)):
            assert vector(k,val)==base+J*vector(k,point)
            assert vector(k,inc)==inc0+B*vector(k,point)
        print("FINITE_COMPLETE_FULL_QUADRATIC_AUDIT",len(points),"PASS",flush=True)
    def centered(a):
        a=int(k(a));return a if a<=11 else a-23
    point=[mixed.K(centered(a)) for a in solution]
    for n in range(1,lift_digits):
        _,v,inc,_=evaluate(point,n+6);A=Arithmetic(mixed.K,n+1)
        correction=-square.inverse()*vector(k,[A.digit(v[j],n) for j in chosen])
        point=[a+mixed.pi**n*centered(b) for a,b in zip(point,correction)]
        _,v,inc,g=evaluate(point,n+6)
        assert all(v[j]==0 for j in chosen)
        assert all(a==0 for a in g.global_residual)
        print("FINITE_COMPLETE_LIFT_DIGIT",n,list(correction),"incidence_div_pi3",inc,flush=True)
    # Same genuine global function on this finite zero locus, computed
    # by the independent older central-series method.
    from certify_harmonic_fixed_hensel_system import evaluate as old_evaluate
    _,old_values,old_inc,old=old_evaluate(point,max(6,lift_digits+5))
    _,values,inc,new=evaluate(point,max(6,lift_digits+5))
    assert all(a==0 for a in values) and all(a==0 for a in old_values)
    assert old["global_relations_vanish"] and not any(new.global_residual)
    assert old["multiplier"]==new.multiplier and old_inc==inc
    print("FINITE_COMPLETE_OLD_METHOD_COMPARISON","PASS",flush=True)
    print("SCOPE: finite complete-map equations; no exact osculating incidence claim",flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--workers",type=int,default=4)
    parser.add_argument("--audit-quadratic",action="store_true")
    parser.add_argument("--lift-digits",type=int,default=1)
    args=parser.parse_args();run(args.workers,args.audit_quadratic,args.lift_digits)
