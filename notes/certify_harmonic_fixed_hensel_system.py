#!/usr/bin/env sage-python
"""A fixed eleven-variable Hensel system, not an incidence certificate.

Only the harmonic special curve and independently determined first jets
are inputs. No reconstructed characteristic-zero cover is read. This
tests a fixed analytic candidate off the global-function locus; all its
unimposed cross-product relations remain visible.
"""
import argparse
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from sage.all import GF,PolynomialRing,matrix,vector
from harmonic_truncated_geometry import Arithmetic,osculating_remainder
from investigate_harmonic_second_incidence_order import mixed,complement

k=GF(23);pi=mixed.pi
ROWS=[
    [1,0,-5,0,0,-7,0,11,0],
    [0,6,0,5,3,0,4,0,3],
    [-1,0,-10,0,0,6,0,2,0],
]


def certify_cubic_slice():
    """The nine cubic directions really complement projective gauge."""
    P=PolynomialRing(k,["H","U","V","W"])
    xs=list(P.gens());H,U,V,W=xs
    q=H*W-U**2-U*V-V**2;c=U*W**2-V**3
    quadratic_monomials=sum(xs)**2
    quadratic_monomials=quadratic_monomials.monomials()
    cubic_monomials=(sum(xs)**3).monomials()
    infinitesimal_quadric=[q.derivative(xs[i])*xs[j] for i in range(4) for j in range(4)]
    quadric_map=matrix(k,[[f.monomial_coefficient(mon) for f in infinitesimal_quadric]
                          for mon in quadratic_monomials])
    orthogonal_lie=quadric_map.right_kernel().basis()
    assert len(orthogonal_lie)==6
    orbit=[]
    for v in orthogonal_lie:
        orbit.append(sum(v[4*i+j]*c.derivative(xs[i])*xs[j] for i in range(4) for j in range(4)))
    gauge=[*orbit,c,*[q*x for x in xs]]
    transverse=[H**3,H**2*U,H*U**2,U**3,H**2*V,H*U*V,U**2*V,H*V**2,U*V**2]
    def columns(polys):
        return matrix(k,[[f.monomial_coefficient(mon) for f in polys] for mon in cubic_monomials])
    assert columns(gauge).rank()==11
    complete=columns(gauge+transverse)
    assert complete.rank()==20
    print("CANONICAL_CUBIC_SLICE",{"orthogonal_dimension":6,"gauge_rank":11,
                                    "full_rank":20,"determinant":int(complete.det())},flush=True)


def evaluate(point,precision=6):
    """Return actual O-valued, normalized fixed equations modulo pi^(N-5).

    Input coordinates are integral number-field elements, not just
    residue-field digits. In particular this supports independent tests
    of the Jacobian on lifts that fail at an earlier precision.
    """
    assert len(point)==11 and precision>=6
    A=Arithmetic(mixed.K,precision)
    c=mixed.c0
    for order,row in enumerate(ROWS,1):
        c+=pi**order*sum(int(a)*f for a,f in zip(row,complement))
    c+=pi**3*point[0]*complement[8]
    c+=pi**4*sum(a*f for a,f in zip(point[1:9],complement[:8]))
    marks=[1+pi+9*pi**2+pi**3*point[9],-1+pi-9*pi**2+pi**3*point[10]]
    result=mixed.run(mixed.q0,c,marks,precision=precision,verbose=False,
                     analytic_candidate=True,series_count=8*(precision-2)+2)
    Rarith=Arithmetic(mixed.K,precision-2)
    R=Rarith.weierstrass(result["eta_lift"],8)
    NN,DD=result["affine_pair"]
    P=A.quo_rem(NN-DD,R**2)[1];Q=A.quo_rem(DD,R**2)[1]
    scalar=A.coefficient(P[0]/Q[0])
    E=A.cut(P-scalar*Q)
    # Independent finite-algebra implementation of the fibre: solve
    # directly in O[z]/(R^2), not by reducing the long central series.
    def mod_square(f):return A.quo_rem(f,R**2)[1]
    ubar=mod_square(A.z**3)
    ctail=c-mixed.c0
    for _ in range(precision+1):
        hbar=mod_square(ubar**2+ubar*A.z+A.z**2)
        ubar=mod_square(A.z**3-ctail(hbar,ubar,A.z,1))
    hbar=mod_square(ubar**2+ubar*A.z+A.z**2)
    assert mod_square(c(hbar,ubar,A.z,1))==0
    nform,dform=result["homogeneous_pair"]
    nbar=mod_square(nform(hbar,ubar,A.z,1))
    dbar=mod_square(dform(hbar,ubar,A.z,1))
    assert nbar==mod_square(NN) and dbar==mod_square(DD)
    assert mod_square(nbar-(1+scalar)*dbar)==E
    # The lowest retained orders of these coefficients differ. Keeping
    # the full ramification factor to N-2 determines every displayed
    # normalized value modulo pi^(N-5), by the z^15/z^8 precision lemma.
    labels=[("fibre",3,15)]
    labels += [("fibre",4,j) for j in range(8,15)]
    labels += [("fibre",5,j) for j in range(1,8)]
    # A fixed z bound rather than the growing truncation width: all
    # possible linear correction rows at the initial stage lie here.
    labels += [("function",pair,j) for pair in range(6) for j in range(136)]
    equations=[]
    small=Arithmetic(mixed.K,precision-5)
    for typ,a,j in labels:
        raw,order=(E[j],a) if typ=="fibre" else (-result["raw_relations"][a][j],5)
        # Asserting integrality checks all lower coefficients, not just
        # the particular digit being sampled.
        A.digit(raw,order)
        equations.append(small.coefficient(raw/pi**order))
    inc,_=osculating_remainder(mixed.q0,c,marks,result["curve_affine"],precision-2)
    incidence=[small.coefficient(inc[j]/pi**3) for j in range(2)]
    # Perturb R's two unavailable digits: the fixed normalized values
    # must not depend on setting those coefficients artificially to zero.
    perturbed=R+pi**(precision-2)*(1+A.z**3+A.z**7)+pi**(precision-1)*(2+A.z**4)
    p=A.quo_rem(NN-DD,perturbed**2)[1];q=A.quo_rem(DD,perturbed**2)[1]
    e=A.cut(p-A.coefficient(p[0]/q[0])*q)
    for typ,order,j in labels[:15]:
        assert small.coefficient((e[j]-E[j])/pi**order)==0
    return labels,equations,incidence,result


def sample(point):
    labels,equations,incidence,_=evaluate([mixed.K(a) for a in point])
    A=Arithmetic(mixed.K,1)
    return [int(A.digit(a,0)) for a in equations],[int(A.digit(a,0)) for a in incidence]


def samples(points,workers):
    if workers==1:return list(map(sample,points))
    with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context("spawn")) as pool:
        return list(pool.map(sample,points))


def run(workers,audit,check_off_locus,lift_digits):
    certify_cubic_slice()
    points=[[0]*11]
    for j in range(11):
        e=[0]*11;e[j]=1;points.append(e)
    values=samples(points,workers)
    base=vector(k,values[0][0]);inc0=vector(k,values[0][1])
    J=matrix(k,[vector(k,v[0])-base for v in values[1:]]).transpose()
    B=matrix(k,[vector(k,v[1])-inc0 for v in values[1:]]).transpose()
    labels,_,_,_=evaluate([mixed.K.zero()]*11)
    chosen=J.transpose().pivots();square=J.matrix_from_rows(chosen)
    constant=vector(k,[base[j] for j in chosen])
    print("FIXED_SYSTEM_SELECTED_EQUATIONS",[labels[j] for j in chosen],flush=True)
    print("FIXED_SYSTEM_JACOBIAN",[list(row) for row in square.rows()],flush=True)
    print("FIXED_SYSTEM_CONSTANT",list(constant),flush=True)
    print("FIXED_SYSTEM_DETERMINANT",square.det(),flush=True)
    assert square.nrows()==square.ncols()==11 and square.det()==16
    assert list(constant)==list(map(k,[0,3,0,5,0,11,0,21,0,20,0]))
    solution=square.solve_right(-constant)
    assert list(solution)==list(map(k,[0,0,16,0,22,21,0,10,0,16,16]))
    assert J*solution+base==0 and B*solution+inc0==0
    print("FIXED_SYSTEM_UNIQUE_RESIDUE_SOLUTION",list(solution),flush=True)
    print("FINITE_FIBRE_ALGEBRA_COMPARISON",len(values),"PASS",flush=True)
    if audit:
        points=[]
        for i in range(11):
            e=[0]*11;e[i]=-1;points.append(e)
            for j in range(i+1,11):
                e=[0]*11;e[i]=e[j]=1;points.append(e)
        for point,(val,inc) in zip(points,samples(points,workers)):
            assert vector(k,val)==base+J*vector(k,point)
            assert vector(k,inc)==inc0+B*vector(k,point)
        print("FIXED_SYSTEM_FULL_QUADRATIC_AUDIT",len(points),"PASS",flush=True)
    # This identity is ONLY the reduction of a prospective analytic
    # syzygy. Higher-order defects are explicitly tested below.
    witnesses=matrix(k,[square.transpose().solve_right(B.row(j)) for j in range(2)])
    assert witnesses*square==B and witnesses*constant==inc0
    print("FIXED_SYSTEM_RESIDUAL_INCIDENCE_WITNESSES",[list(row) for row in witnesses],flush=True)
    if check_off_locus:
        point=[mixed.K(a) for a in solution]
        point[0]+=pi
        labels,lifts,inc,result=evaluate(point,precision=8)
        A=Arithmetic(mixed.K,3)
        selected=vector(mixed.K,[lifts[j] for j in chosen])
        w=matrix(mixed.K,[[int(a) for a in row] for row in witnesses.rows()])
        defect=vector(mixed.K,inc)-w*selected
        defect=[A.coefficient(a) for a in defect]
        print("OFF_LOCUS_GLOBAL_RELATIONS_VANISH",result["global_relations_vanish"],flush=True)
        print("OFF_LOCUS_CONSTANT_WITNESS_DEFECT_MOD_PI3",defect,flush=True)
        assert all(A.digit(a,0)==0 for a in defect)
        assert defect==[mixed.K(-115),mixed.K.zero()]
    if lift_digits>1:
        def centered(a):
            a=int(k(a));return a if a<=11 else a-23
        point=[mixed.K(centered(a)) for a in solution]
        inverse=square.inverse()
        for n in range(1,lift_digits):
            _,lifts,inc,_=evaluate(point,precision=n+6)
            A=Arithmetic(mixed.K,n+1)
            residual=vector(k,[A.digit(lifts[j],n) for j in chosen])
            correction=-inverse*residual
            point=[a+pi**n*centered(b) for a,b in zip(point,correction)]
            _,next_lifts,next_inc,result=evaluate(point,precision=n+6)
            assert all(A.coefficient(next_lifts[j])==0 for j in chosen)
            print("FIXED_HENSEL_DIGIT",n,list(correction),flush=True)
            print("FIXED_HENSEL_UNIMPOSED_CHECK",n,
                  {"all_listed_equations_zero":all(a==0 for a in next_lifts),
                   "global_function_relations_zero":result["global_relations_vanish"],
                   "incidence_div_pi3":next_inc},flush=True)
        print("FIXED_HENSEL_POINT",point,flush=True)
    print("SCOPE: fixed necessary equations and an invertible Jacobian; no exact incidence conclusion",flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--workers",type=int,default=4)
    parser.add_argument("--audit-quadratic",action="store_true")
    parser.add_argument("--check-off-locus",action="store_true")
    parser.add_argument("--lift-digits",type=int,default=1)
    args=parser.parse_args();run(args.workers,args.audit_quadratic,args.check_off_locus,args.lift_digits)
