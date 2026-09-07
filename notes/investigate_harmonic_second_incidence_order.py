#!/usr/bin/env sage-python
"""Next necessary effective-fibre equations, without generic cover input.

Only selected coefficients of beta-1 modulo R^2 are used: they are
independent of the unknown higher Weierstrass-factor digits. This is
not an all-orders existence theorem.
"""
import argparse
import importlib.util
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from sage.all import GF, PolynomialRing, matrix, vector
from harmonic_truncated_geometry import Arithmetic,osculating_remainder

spec=importlib.util.spec_from_file_location("mixed",Path(__file__).with_name(
    "investigate_harmonic_mixed_characteristic_lifts.py"))
mixed=importlib.util.module_from_spec(spec);spec.loader.exec_module(mixed)
k=GF(23);H,U,V,W=mixed.P.gens();pi=mixed.pi
complement=[H**3,H**2*U,H*U**2,U**3,H**2*V,H*U*V,U**2*V,H*V**2,U*V**2]
C1=H**3-5*H*U**2-7*H*U*V+11*H*V**2
C20=6*H**2*U+5*U**3+3*H**2*V+4*U**2*V
A5=Arithmetic(mixed.K,5);A3=Arithmetic(mixed.K,3)


def evaluate(t,point,verbose=False):
    C2=C20+int(t)*U*V**2
    C3=sum(int(a)*f for a,f in zip(point[:9],complement))
    c=mixed.c0+pi*C1+pi**2*C2+pi**3*C3
    marks=[1+pi+pi**2*int(point[9]),-1+pi+pi**2*int(point[10])]
    result=mixed.run(mixed.q0,c,marks,precision=5,verbose=verbose,diagnostic_candidate=True)
    eta=result["eta_lift"] if result["pass"] else result["candidate_eta_lift"]
    R=A3.weierstrass(eta,8)
    NN,DD=result["affine_pair"] if result["pass"] else result["candidate_affine_pair"]
    assert A5.coefficient(NN[0]-DD[0])==0
    remainder=A5.quo_rem(NN-DD,R**2)[1]
    fibre=[]
    # Unknown r3 first affects pi^3*z^15 and pi^4*z^8.
    # Unknown r4 first affects pi^4*z^15. All selected coefficients
    # are therefore determined by R modulo pi^3.
    for order,count in [(2,16),(3,15),(4,8)]:
        fibre.extend(A5.digit(remainder[j],order) for j in range(count))
    residual=vector(k,[*result["obstruction"],*fibre])
    incidence,_=osculating_remainder(mixed.q0,c,marks,result["curve_affine"],3)
    inc=vector(k,[A3.digit(incidence[j],2) for j in range(2)])
    return residual,inc,result,R


def solve(t):
    zero=vector(k,11);base,inc0,_,_=evaluate(t,zero)
    print("NEXT_ORDER_BASE",int(t),"nonzero",bool(base),"incidence",list(inc0),flush=True)
    columns=[];inc_columns=[]
    for j in range(11):
        e=vector(k,11);e[j]=1
        value,inc,_,_=evaluate(t,e)
        columns.append(value-base);inc_columns.append(inc-inc0)
        print("NEXT_ORDER_COLUMN",int(t),j,flush=True)
    system=matrix(k,columns).transpose();incmap=matrix(k,inc_columns).transpose()
    test=vector(k,[2,1,-1,0,1,0,0,1,2,-1,2])
    actual,inc,_,_=evaluate(t,test)
    difference=actual-base-system*test
    if difference:
        print("NONLINEAR_RESIDUAL_COORDINATES",[(j,int(a)) for j,a in enumerate(difference) if a],flush=True)
    assert actual[:-1]==(base+system*test)[:-1]
    assert inc==inc0+incmap*test
    base=base[:-1];system=system[:-1,:]
    rank=system.rank();augmented=system.augment(matrix(k,len(base),1,base)).rank()
    print("NEXT_ORDER_SYSTEM",int(t),"rank",rank,"augmented_rank",augmented,flush=True)
    if rank!=augmented:return None
    solution=system.solve_right(-base)
    residual,inc,result,R=evaluate(t,solution,True)
    assert residual[:-1]==0 and result["pass"]
    kernel=system.right_kernel().basis_matrix()
    variation=incmap*kernel.transpose()
    print("NEXT_ORDER_SOLUTION",int(t),list(solution),"kernel",[list(v) for v in kernel.rows()],flush=True)
    print("REMAINING_LAST_FIBRE_EQUATION_AT_SOLUTION",int(residual[-1]),flush=True)
    print("SECOND_INCIDENCE_OBSTRUCTION",list(inc),"variation_rank",variation.rank(),flush=True)
    print("WEIERSTRASS_RAMIFICATION",R,flush=True)
    return solution,kernel,inc,variation


def parameter_sample(job):
    t,point=job
    residual,inc,_,_=evaluate(k(t),vector(k,point))
    return [int(a) for a in residual[:-1]]


def solve_parameter(workers):
    # Coarse degree bound: beta is truncated modulo pi^5, R is prepared
    # modulo pi^3, and division by R^2 uses at most four positive-order
    # factors. The selected coefficient equations are affine in the
    # third curve jet and second mark motions. Interpolate their t
    # coefficients through degree four, then test extra mixed points.
    jobs=[]
    for t in range(5):
        jobs.append((t,[0]*11))
        for j in range(11):
            point=[0]*11;point[j]=1;jobs.append((t,point))
    if workers>1:
        with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context("spawn")) as pool:
            results=list(pool.map(parameter_sample,jobs))
    else:results=[parameter_sample(job) for job in jobs]
    print("PARAMETER_INTERPOLATION_SAMPLES",len(results),flush=True)
    values=[vector(k,row) for row in results]
    inverse=matrix(k,[[k(t)**j for j in range(5)] for t in range(5)]).inverse()
    bases=matrix(k,[values[12*t] for t in range(5)])
    base_coefficients=inverse*bases
    operator_coefficients=[]
    for j in range(11):
        table=matrix(k,[values[12*t+j+1]-values[12*t] for t in range(5)])
        operator_coefficients.append(inverse*table)
    P=PolynomialRing(k,["t",*["v%d"%j for j in range(11)]],order="degrevlex")
    t,*vs=P.gens()
    coeffs=[*base_coefficients.rows()]
    monomials=[t**j for j in range(5)]
    for j,table in enumerate(operator_coefficients):
        coeffs.extend(table.rows());monomials.extend(t**d*vs[j] for d in range(5))
    coefficient_matrix=matrix(k,coeffs).transpose()
    for value,point in [(5,[2,1,-1,0,1,0,0,1,2,-1,2]),(9,[1]*11),(-1,[0,1,0,0,2,0,0,0,1,3,-2])]:
        arguments=[k(value),*map(k,point)]
        expected=coefficient_matrix*vector(k,[m(*arguments) for m in monomials])
        assert expected==vector(k,parameter_sample((value,point)))
    basis=coefficient_matrix.row_space().basis_matrix()
    equations=[sum(a*m for a,m in zip(row,monomials)) for row in basis]
    ideal=P.ideal(equations)
    print("PARAMETER_EQUATION_COUNT",len(equations),flush=True)
    print("PARAMETER_GROEBNER_BASIS",list(ideal.groebner_basis()),flush=True)
    elimination=ideal.elimination_ideal(vs)
    print("PARAMETER_ELIMINATION",list(elimination.groebner_basis()),flush=True)
    T=PolynomialRing(k,"t")
    polynomials=[T(f(t,*([P.zero()]*11))) for f in elimination.groebner_basis()]
    common=T.zero()
    for f in polynomials:common=common.gcd(f)
    assert common
    roots=common.roots(multiplicities=False)
    print("ADMISSIBLE_SECOND_JET_PARAMETERS",roots,flush=True)
    for value in roots:solve(value)


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--t",type=int,nargs="+",default=[0])
    parser.add_argument("--solve-parameter",action="store_true")
    parser.add_argument("--workers",type=int,default=4)
    args=parser.parse_args()
    if args.solve_parameter:solve_parameter(args.workers)
    else:
        for t in args.t:solve(k(t))
    print("SCOPE: next necessary finite-order effective-fibre equations, not a full Hurwitz lift")
