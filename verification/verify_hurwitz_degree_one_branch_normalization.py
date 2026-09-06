#!/usr/bin/env sage-python
"""Compute the third-branch scaling on the degree-one component exactly.

The stored degree-one map is unnormalized. A uniform geometric
cross-ratio comparison must account for that scaling before using 16.
"""
import json
from pathlib import Path
from sage.all import QQ, PolynomialRing, matrix, prod

ROOT=Path(__file__).resolve().parents[1]
can=json.loads((ROOT/"data/hurwitz_canonical_models_candidate.json").read_text())
maps=json.loads((ROOT/"data/hurwitz_degree23_maps_candidate.json").read_text())
U=PolynomialRing(QQ,"s0")
s0=U.gen()
K=QQ.extension(s0*s0+23,"s")
s=K.gen()
def rat(x):
    if isinstance(x,dict):return QQ(x["numerator"])/QQ(x["denominator"])
    return QQ(x)
def coeff(x):
    x=x["degree_one_component"]
    return rat(x["rational_part"])+rat(x["sqrt_minus_23_part"])*s
R=PolynomialRing(K,["x1","x2","x3"])
xs=[R.one(),*R.gens()]
equations=[]
for name in ["quadric","petri_cubic"]:
    table=can[name]
    equations.append(sum(coeff(c)*prod(xs[i] for i in mon)
                         for c,mon in zip(table["coefficients"],table["monomials"])))
rec=next(c for c in maps["components"] if c["component"]=="degree_one")
sections=[]
for name in ["numerator","denominator"]:
    sections.append(sum(K([rat(c) for c in term["coefficient_power_basis"]])*
                        prod(x**e for x,e in zip(xs,term["monomial"]))
                        for term in rec[name]))
N,D=sections
rows=[[f.derivative(x) for x in R.gens()] for f in equations]
rows.append([D*N.derivative(x)-N*D.derivative(x) for x in R.gens()])
critical=matrix(R,rows).det()
I=R.ideal([*equations,critical])
print("saturating the exact critical-point ideal away from N*D=0",flush=True)
I=I.saturation(R.ideal(N*D))[0]
dimension=I.vector_space_dimension()
assert 0<dimension<=8
Nr,Dr=I.reduce(N),I.reduce(D)
assert Dr
mon=next(iter(Dr.dict()))
lam=Nr.monomial_coefficient(R.monomial(*mon))/Dr.monomial_coefficient(R.monomial(*mon))
assert lam and Nr==lam*Dr
prime=K.prime_above(23)
k=prime.residue_field()
assert lam.valuation(prime)==0
assert k(lam)==7
assert k(16)/k(lam)==-1
print("finite critical algebra dimension",dimension,flush=True)
print("exact third branch value lambda:",lam,flush=True)
print("lambda modulo 23:",k(lam),flush=True)
print("normalized distinguished attachment:",k(16)/k(lam),flush=True)
print("PASS_DEGREE_ONE_BRANCH_NORMALIZATION")
