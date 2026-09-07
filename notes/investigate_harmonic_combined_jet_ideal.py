#!/usr/bin/env sage-python
"""Combine necessary differential jets and mixed-characteristic function jets.

No reconstructed covers are input. This is a finite-order investigation,
not a full Hurwitz construction. Off the function-lifting locus, the
candidate differential is used ONLY as a polynomial extension for forming
equations; it is never asserted to be a differential of a global map.
"""
import contextlib
import importlib.util
import io
from pathlib import Path
from sage.all import GF, PolynomialRing, matrix, vector


def load(name,filename):
    spec=importlib.util.spec_from_file_location(name,Path(__file__).with_name(filename))
    module=importlib.util.module_from_spec(spec)
    with contextlib.redirect_stdout(io.StringIO()):spec.loader.exec_module(module)
    return module


critical=load("critical","investigate_harmonic_critical_fibre_incidence.py")
mixed=load("mixed","investigate_harmonic_mixed_characteristic_lifts.py")
k=GF(23);P=critical.P;H,U,V,W=P.gens()
q=H*W-U*U-U*V-V*V;c=U*W*W-V**3
qmons=((H+U+V+W)**2).monomials()
cmons=((H+U+V+W)**3).monomials()
infinitesimal=[]
for r in range(4):
    for s in range(4):infinitesimal.append(q.derivative(P.gen(r))*P.gen(s))
qsystem=matrix(k,[[f.monomial_coefficient(m) for f in infinitesimal]+[-q.monomial_coefficient(m)]
                  for m in qmons])
assert qsystem.right_kernel().dimension()==7
gauge=[]
for row in qsystem.right_kernel().basis():
    mat=matrix(k,4,4,row[:16])
    move=mat*vector(P,P.gens())
    gauge.append(sum(c.derivative(x)*v for x,v in zip(P.gens(),move)))
gauge+=[c,*[q*x for x in P.gens()]]
def coefficient_vector(f):return vector(k,[f.monomial_coefficient(m) for m in cmons])
gauge_matrix=matrix(k,[coefficient_vector(f) for f in gauge])
assert gauge_matrix.rank()==11
span=gauge_matrix.row_space();complement=[]
for mon in cmons:
    v=coefficient_vector(mon)
    if v not in span:
        complement.append(mon)
        span=span+matrix(k,[v]).row_space()
assert len(complement)==9
print("CANONICAL_DEFORMATION_COMPLEMENT",complement,flush=True)

columns=[critical.conditions(0,f,0,0,0) for f in complement]
columns+=[critical.conditions(0,0,1,0,0),critical.conditions(0,0,0,1,0)]
columns+=[critical.conditions(0,0,0,0,f) for f in critical.numerator_basis]
system=matrix(k,columns).transpose()
assert system.ncols()==26 and system.rank()==18
origin=system.solve_right(critical.rhs)
directions=system.right_kernel().basis_matrix()
assert directions.nrows()==8
inc_cols=[critical.incidence.obstruction(0,f)[0] for f in complement]
inc_cols+=[critical.incidence.obstruction(0,0,1,0)[0],critical.incidence.obstruction(0,0,0,1)[0]]
inc_cols+=[vector(k,[0,0]) for _ in critical.numerator_basis]
inc_matrix=matrix(k,inc_cols).transpose()
Q=PolynomialRing(k,["t%d"%j for j in range(8)],order="degrevlex")
ts=Q.gens()
affine=[Q(origin[j])+sum(Q(directions[i,j])*ts[i] for i in range(8))
        for j in range(26)]
incidence_equations=list(inc_matrix*vector(Q,affine))
print("INCIDENCE_LINEAR_EQUATIONS",incidence_equations,flush=True)

cache={}
def evaluate(point):
    key=tuple(int(a) for a in point)
    if key in cache:return cache[key]
    values=origin+directions.transpose()*vector(k,point)
    c1=sum(a*f for a,f in zip(values[:9],complement))
    lifted=mixed.P({mon:mixed.K(int(a)) for mon,a in c1.dict().items()})
    result=mixed.run(mixed.q0,mixed.c0+mixed.pi*lifted,
                     [mixed.K(1),mixed.K(-1)],precision=3,
                     label=str(key),verbose=False,diagnostic_candidate=True)
    eta=result["eta0"] if result["pass"] else result["candidate_eta0"]
    residual=vector(k,[*result["obstruction"],*[eta[j] for j in [0,1,3,6]]])
    cache[key]=residual
    print("QUADRATIC_SAMPLE",len(cache),"point",key,"function_lifts",result["pass"],
          "residual_zero",not bool(residual),flush=True)
    return residual

zero=vector(k,8);constant=evaluate(zero)
linear=[];square=[]
for i in range(8):
    e=vector(k,8);e[i]=1
    plus,minus=evaluate(e),evaluate(-e)
    linear.append((plus-minus)/k(2))
    square.append((plus+minus)/k(2)-constant)
cross=[];pairs=[]
for i in range(8):
    for j in range(i+1,8):
        e=vector(k,8);e[i]=e[j]=1
        cross.append(evaluate(e)-constant-linear[i]-linear[j]-square[i]-square[j])
        pairs.append((i,j))
coefficient_matrix=matrix(k,[constant,*linear,*square,*cross]).transpose()
monomials=[Q.one(),*ts,*[t*t for t in ts],*[ts[i]*ts[j] for i,j in pairs]]
rows=coefficient_matrix.row_space().basis_matrix()
equations=[sum(a*m for a,m in zip(row,monomials)) for row in rows]
print("NUMBER_OF_INDEPENDENT_QUADRATIC_EQUATIONS",len(equations),flush=True)
for eq in equations:print("QUADRATIC_EQUATION",eq,flush=True)
# Extra values independently test the quadratic interpolation. The degree
# bound comes from truncation at pi^3: at most two first-order coefficient
# factors can occur, and every matrix inversion uses a fixed residual unit.
for point in [vector(k,[2,1,0,-1,0,1,0,2]),vector(k,[1]*8)]:
    predicted=coefficient_matrix*vector(k,[m(*point) for m in monomials])
    assert predicted==evaluate(point)
ideal=Q.ideal(equations)
print("COMBINED_IDEAL_GROEBNER_START",flush=True)
gb=ideal.groebner_basis()
print("COMBINED_IDEAL_DIMENSION",ideal.dimension(),flush=True)
print("COMBINED_IDEAL_GROEBNER_BASIS",list(gb),flush=True)
print("INCIDENCE_REMAINDERS",[ideal.reduce(f) for f in incidence_equations],flush=True)
if all(ideal.reduce(f)==0 for f in incidence_equations):
    print("PASS_COMBINED_JETS_FORCE_FIRST_ORDER_INCIDENCE",flush=True)
else:
    print("INCIDENCE_NOT_IN_COMBINED_IDEAL",flush=True)


def extract_tail(map_pair,c1):
    # pi=rho^15, (U,V,W)=(rho^3 X,rho^5 Y,rho^6 Z), H=1.
    # The weight bound 24 is below v_rho(23)=30, so k[pi] arithmetic
    # is legitimate here, unlike in the mixed function-lifting test.
    T=PolynomialRing(k,["p","u","v"]);p,u,v=T.gens()
    def cut(f):return T({mon:a for mon,a in T(f).dict().items()
                         if 15*mon[0]+3*mon[1]+5*mon[2]<24})
    w=u*u+u*v+v*v
    def convert(form):
        total=T.zero()
        for (a,b,d,e),coefficient in form.dict().items():
            cs=list(mixed.K(coefficient))+[mixed.QQ.zero(),mixed.QQ.zero()]
            val=k(cs[0])+p*k(cs[1])
            total=cut(total+cut(val*u**b*v**d*w**e))
        return total
    cubic=cut(u*w*w-v**3+p*c1(1,u,v,w))
    assert cubic.monomial_coefficient(v**3)==-1
    rewrite=cut(cubic+v**3)
    def reduce(f):
        f=cut(f)
        while True:
            eligible=[mon for mon in f.dict() if mon[2]>=3]
            if not eligible:return f
            mon=min(eligible,key=lambda x:(15*x[0]+3*x[1]+5*x[2],-x[2]))
            coefficient=f.dict()[mon]
            old=coefficient*p**mon[0]*u**mon[1]*v**mon[2]
            f=cut(f-old+coefficient*p**mon[0]*u**mon[1]*v**(mon[2]-3)*rewrite)
    N,D=[convert(f) for f in map_pair]
    assert D.constant_coefficient()
    denominator=D.constant_coefficient()
    step=cut(1-D/denominator)
    term=T.one();inverse=T.one()
    for _ in range(24):term=cut(term*step);inverse=cut(inverse+term)
    inverse=cut(inverse/denominator)
    assert cut(D*inverse)==1
    value=reduce(cut(-N*inverse)) # normalized third branch value is one
    centre=value( p,0,0 )
    moving=reduce(value-centre)
    assert moving
    weights={15*a+3*b+5*d for a,b,d in moving.dict()}
    print("NEW_TAIL_NONCONSTANT_WEIGHTS",sorted(weights),flush=True)
    assert weights=={23}
    B=PolynomialRing(k,["X","Y"]);X,Y=B.gens()
    tail=sum(a*X**b*Y**d for (j,b,d),a in moving.dict().items())
    lam=c1.monomial_coefficient(H**3)
    leading=tail.monomial_coefficient(X**6*Y)
    assert lam and leading
    assert tail==leading*Y*(X**6+15*lam*X)
    print("NEW_TAIL_CURVE",Y**3-X**5-lam,"MAP",tail,"CENTRE",centre,flush=True)
    return lam,leading


points=ideal.variety()
assert len(points)==2
for point in points:
    parameters=vector(k,[point[t] for t in ts])
    values=origin+directions.transpose()*parameters
    c1=sum(a*f for a,f in zip(values[:9],complement))
    marks=list(values[9:11])
    print("COMBINED_SOLUTION",list(parameters),"CUBIC_FIRST_JET",c1,
          "MARK_MOTIONS",marks,"INCIDENCE",[f(*parameters) for f in incidence_equations],flush=True)
    assert evaluate(parameters)==0
    lifted=mixed.P({mon:mixed.K(int(a)) for mon,a in c1.dict().items()})
    lifted_map=mixed.run(mixed.q0,mixed.c0+mixed.pi*lifted,
                        [mixed.K(1)+mixed.pi*int(marks[0]),
                         mixed.K(-1)+mixed.pi*int(marks[1])],precision=2,verbose=False)
    assert lifted_map["pass"]
    pair=next(pair for pair in lifted_map["homogeneous_pairs"]
              if pair[1](1,0,0,0) and k(list(pair[1](1,0,0,0))[0]))
    extract_tail(pair,c1)
print("SCOPE: necessary finite-order conditions, not an all-orders existence theorem")
