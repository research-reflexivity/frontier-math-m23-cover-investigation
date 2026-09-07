#!/usr/bin/env sage-python
"""Linearized degree-23 function equations at the harmonic E8 curve.

No reconstructed model or map coefficients are read. This tests only the
specified first-order conditions; it is NOT an all-orders Hurwitz lift.
"""
from sage.all import GF, PolynomialRing, matrix, vector, prod

k=GF(23)
P=PolynomialRing(k,["H","U","V","W"])
H,U,V,W=P.gens()
R=PolynomialRing(k,"z");z=R.gen()
F=z**6+z**4+z*z
monomials=(H+U+V+W)**5
all_monomials=monomials.monomials()
all_values=[f(F,z**3,z,1) for f in all_monomials]
images=matrix(k,[[f[j] for f in all_values] for j in range(31)])
assert images.rank()==27
basis=[all_monomials[j] for j in images.pivots()]
values=[f(F,z**3,z,1) for f in basis]
assert len(basis)==27


def combine(coefficients, polys):
    return sum((a*f for a,f in zip(coefficients,polys)),R.zero())


jet_data=[]
for mark in [k(1),k(-1)]:
    shifted=[f(z+mark) for f in values]
    jets=matrix(k,[[f[j] for f in shifted] for j in range(23)])
    assert jets.rank()==23
    kernel=jets.right_kernel().basis()
    assert len(kernel)==4
    pivots=jets.pivots()
    inverse=jets.matrix_from_columns(pivots).inverse()
    sections=[combine(row,values) for row in kernel]
    assert all(f % (z-mark)**23==0 for f in sections)
    jet_data.append((mark,jets,kernel,pivots,inverse,sections))
SA0,SB0=jet_data[0][-1],jet_data[1][-1]
gA=[f//(z-1)**23 for f in SA0]
gB=[f//(z+1)**23 for f in SB0]
common=gA[0]
for f in gA[1:]:common=common.gcd(f)
assert common.degree()==0 and max(f.degree() for f in SA0)==30
matrixB=matrix(k,[[f[j] for f in gB] for j in range(8)])
M0=matrix(k,[matrixB.solve_right(vector(k,[f[j] for j in range(8)])) for f in gA])
assert M0.det()
D0=[combine(row,SB0) for row in M0.rows()]
beta0=R.fraction_field()((z-1)**23/(z+1)**23)
assert all(R.fraction_field()(f/g)==beta0 for f,g in zip(SA0,D0))

# A first-order product of two quintic sections has z-degree at most 75
# for any quadratic/cubic coefficient perturbation considered below.
WIDTH=76
pairs=[(a,b) for a in range(4) for b in range(a+1,4)]
def flatten(polys):
    assert len(polys)==6 and all(f.degree()<WIDTH for f in polys)
    return vector(k,[f[j] for f in polys for j in range(WIDTH)])

columns=[]
for row in range(4):
    for col in range(4):
        polys=[]
        for a,b in pairs:
            polys.append((SA0[a]*SB0[col] if b==row else R.zero())
                         -(SA0[b]*SB0[col] if a==row else R.zero()))
        columns.append(flatten(polys))
system=matrix(k,columns).transpose()
assert system.rank()==15
assert system*vector(k,list(M0.list()))==0
pivot_columns=system.pivots()
injective=system.matrix_from_columns(pivot_columns)
pivot_rows=injective.transpose().pivots()
solve_square=injective.matrix_from_rows(pivot_rows).inverse()


def direction(q1,c1,alpha=0,beta=0):
    q1,c1=P(q1),P(c1)
    u1=-c1(F,z**3,z,1)
    h1=(2*z**3+z)*u1-q1(F,z**3,z,1)
    values1=[f.derivative(H)(F,z**3,z,1)*h1
             +f.derivative(U)(F,z**3,z,1)*u1 for f in basis]
    sections1=[]
    for motion,data in zip([k(alpha),k(beta)],jet_data):
        mark,jets,kernel,pivots,inverse,sections=data
        shifted=[(f1+motion*f0.derivative())(z+mark)
                 for f0,f1 in zip(values,values1)]
        jets1=matrix(k,[[f[j] for f in shifted] for j in range(23)])
        current=[]
        for row in kernel:
            correction=vector(k,27)
            for j,a in zip(pivots,inverse*(-jets1*row)):
                correction[j]=a
            assert jets*correction+jets1*row==0
            current.append(combine(row,values1)+combine(correction,values))
        sections1.append(current)
    SA1,SB1=sections1
    known_D1=[combine(row,SB1) for row in M0.rows()]
    rhs=-flatten([SA1[a]*D0[b]+SA0[a]*known_D1[b]
                  -SA1[b]*D0[a]-SA0[b]*known_D1[a] for a,b in pairs])
    candidate=solve_square*vector(k,[rhs[j] for j in pivot_rows])
    residual=rhs-injective*candidate
    if residual:
        return residual,None
    correction=vector(k,16)
    for j,a in zip(pivot_columns,candidate):correction[j]=a
    M1=matrix(k,4,4,correction)
    D1=[known+combine(row,SB0) for known,row in zip(known_D1,M1.rows())]
    for a,b in pairs:
        assert SA1[a]*D0[b]+SA0[a]*D1[b]-SA1[b]*D0[a]-SA0[b]*D1[a]==0
    betas1=[R.fraction_field()((n1*d0-n0*d1)/(d0*d0))
            for n0,n1,d0,d1 in zip(SA0,SA1,D0,D1)]
    assert all(f==betas1[0] for f in betas1)
    return residual,betas1[0]


def main():
    results=[]
    labels=[]
    for degree,which in [(2,"q"),(3,"c")]:
        for mon in ((H+U+V+W)**degree).monomials():
            residual,beta1=direction(mon if which=="q" else 0,mon if which=="c" else 0)
            labels.append(which+":"+str(mon))
            results.append((residual,beta1))
            assert not residual and beta1==0
            print("DIRECTION",labels[-1],"torsion_function_lifts",not bool(residual),
                  "dualizing_differential_zero_mod_epsilon2",beta1.derivative()==0,
                  "beta1",beta1,flush=True)
    for a,b,label in [(1,0,"mark_A"),(0,1,"mark_B")]:
        residual,beta1=direction(0,0,a,b)
        labels.append(label);results.append((residual,beta1))
        assert not residual and beta1==0
        print("DIRECTION",label,"torsion_function_lifts",not bool(residual),
              "dualizing_differential_zero_mod_epsilon2",beta1.derivative()==0,
              "beta1",beta1,flush=True)
    obstructions=matrix(k,[residual for residual,beta1 in results]).transpose()
    assert obstructions.rank()==0
    print("TORSION_FUNCTION_LINEAR_OBSTRUCTION_RANK",obstructions.rank(),flush=True)
    print("INPUT=harmonic special canonical curve only; no reconstructed coefficients")
    print("SCOPE=first order, not the full three-point cover conditions")
    print("PASS_FIRST_ORDER_HARMONIC_FUNCTION_TEST")


if __name__=="__main__":main()
