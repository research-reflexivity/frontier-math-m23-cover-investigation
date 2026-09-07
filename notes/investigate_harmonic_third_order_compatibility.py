#!/usr/bin/env sage-python
"""Test both necessary first jets one mixed-characteristic order further.

The two first jets are computed without reconstructed covers by
investigate_harmonic_combined_jet_ideal.py. Here we impose the actual
divided logarithmic differential, not just existence of a meromorphic
differential with the requested local jets. Passing is still finite-order
evidence only, not a Belyi map or a full M23 monodromy certificate.
"""
import importlib.util
from pathlib import Path
from sage.all import GF, matrix, vector

spec=importlib.util.spec_from_file_location("mixed",Path(__file__).with_name(
    "investigate_harmonic_mixed_characteristic_lifts.py"))
mixed=importlib.util.module_from_spec(spec);spec.loader.exec_module(mixed)
k=GF(23);H,U,V,W=mixed.P.gens();pi=mixed.pi
complement=[H**3,H**2*U,H*U**2,U**3,H**2*V,H*U*V,U**2*V,H*V**2,U*V**2]
first_jet=H**3-5*H*U**2-7*H*U*V+11*H*V**2
desired=vector(k,[0,1,0,15,0,12,0])

for sign in [1,-1]:
    def evaluate(point,verbose=False):
        second=sum(int(a)*f for a,f in zip(point,complement))
        result=mixed.run(mixed.q0,mixed.c0+sign*pi*first_jet+pi**2*second,
                         [1+pi,-1+pi],precision=4,verbose=verbose,
                         diagnostic_candidate=True,label="sign=%s"%sign)
        eta=result["eta1"] if result["pass"] else result["candidate_eta1"]
        return vector(k,[*result["obstruction"],
                         *list(vector(k,eta[:7])-desired)]),result
    zero=vector(k,9);base,_=evaluate(zero)
    print("THIRD_ORDER_SIGN",sign,"BASE_RESIDUAL_NONZERO",bool(base),flush=True)
    columns=[]
    for j in range(9):
        e=vector(k,9);e[j]=1
        value,_=evaluate(e)
        columns.append(value-base)
        print("THIRD_ORDER_COLUMN",sign,j,flush=True)
    system=matrix(k,columns).transpose()
    test=vector(k,[1,2,-1,0,1,0,2,1,-2])
    assert evaluate(test)[0]==base+system*test
    rank=system.rank();augmented=system.augment(matrix(k,len(base),1,base)).rank()
    print("THIRD_ORDER_LINEAR_SYSTEM",sign,"rank",rank,"augmented_rank",augmented,flush=True)
    if rank==augmented:
        solution=system.solve_right(-base)
        assert evaluate(solution,verbose=True)[0]==0
        print("THIRD_ORDER_SOLUTION",sign,list(solution),"free_dimension",9-rank,flush=True)
        print("THIRD_ORDER_KERNEL",sign,[list(v) for v in system.right_kernel().basis()],flush=True)
        print("SECOND_CURVE_JET",sum(int(a)*f for a,f in zip(solution,complement)),flush=True)
    else:
        witness=system.transpose().right_kernel().basis_matrix()*base
        assert witness
        print("THIRD_ORDER_EXCLUDES_FIRST_JET",sign,flush=True)
print("SCOPE: mod pi^4 function and critical differential jets, not all-orders ramification")
