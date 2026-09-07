#!/usr/bin/env sage-python
"""Test critical-fibre jet conditions against first-order osculation.

The meromorphic differential is constructed as a cubic numerator divided
by the quadratic section through the two marks. Regularity is imposed on
the residual length-ten divisor, not inferred from numerical periods.
No reconstructed curve or map is read.
"""
import contextlib
import importlib.util
import io
from pathlib import Path
from sage.all import GF, PolynomialRing, matrix, vector

with contextlib.redirect_stdout(io.StringIO()):
    spec=importlib.util.spec_from_file_location("incidence",Path(__file__).with_name(
        "certify_e8_osculating_deformation_obstructions.py"))
    incidence=importlib.util.module_from_spec(spec);spec.loader.exec_module(incidence)
k=GF(23)
P=PolynomialRing(k,["H","U","V","W"]);H,U,V,W=P.gens()
R=PolynomialRing(k,"z");z=R.gen();F=z**6+z**4+z*z
N0=2*H*V*V-2*U*U*W-2*U*V*W
assert N0(F,z**3,z,1)==2*z**8
# Necessary critical-fibre pattern from beta-1=R8^2*F7. See the note:
# beta_2=2(z^9/9+z^11/11+z^13/13) modulo z^15 and
# beta_2=6*z^7*r1^2 modulo z^15, where r1 is divisible by z.
square_root=11+16*z*z+13*z**4
assert square_root**2 % z**6 == 6+7*z*z+13*z**4
beta2=2*(z**9/k(9)+z**11/k(11)+z**13/k(13))
assert (beta2-6*z**9*square_root**2) % z**15 == 0
beta3=11*z*z+2*z**4+21*z**6
assert ((1-z*z)*beta3.derivative()/k(2)-z*square_root) % z**7 == 0
assert -beta3.derivative() % z**7 == z+15*z**3+12*z**5
assert k(5)**3+k(6)**2 == 0
all_mons=((H+U+V+W)**3).monomials()
images=matrix(k,[[f(F,z**3,z,1)[j] for f in all_mons] for j in range(19)])
numerator_basis=[all_mons[j] for j in images.pivots()]
assert len(numerator_basis)==15

# H=1 near the cusp. Localizing at 1-3w removes the two smooth marked
# points from the degree-twelve divisor V^2-W^2.
T=PolynomialRing(k,["u","v","w","t"]);u,v,w,t=T.gens()
q=w-u*u-u*v-v*v;c=u*w*w-v**3;ell=v*v-w*w
localizer=t*(1-3*w)-1
ideal=T.ideal(q,c,ell,localizer)
assert ideal.vector_space_dimension()==10
n0=N0(1,u,v,w)
lift=n0.lift(ideal)
assert sum(a*f for a,f in zip(lift,ideal.gens()))==n0
normal_basis=ideal.normal_basis()
assert len(normal_basis)==10
def coordinates(f):
    reduced=ideal.reduce(f)
    coords=vector(k,[reduced.monomial_coefficient(m) for m in normal_basis])
    assert sum(a*m for a,m in zip(coords,normal_basis))==reduced
    return coords

def conditions(q1,c1,alpha,beta,n1):
    q1,c1,n1=P(q1),P(c1),P(n1)
    alpha,beta=k(alpha),k(beta)
    s,d=alpha+beta,alpha-beta
    ell1=-s*v*w-d*w*w
    regularity=coordinates(n1(1,u,v,w)-lift[0]*q1(1,u,v,w)
                           -lift[1]*c1(1,u,v,w)-lift[2]*ell1)
    uu1=-c1(F,z**3,z,1)
    hh1=(2*z**3+z)*uu1-q1(F,z**3,z,1)
    numerator1=n1(F,z**3,z,1)+N0.derivative(H)(F,z**3,z,1)*hh1
    numerator1+=N0.derivative(U)(F,z**3,z,1)*uu1
    determinant1=(q1.derivative(H)+c1.derivative(U)
                  +(2*U+V)*c1.derivative(H))(F,z**3,z,1)
    residue_A=numerator1(1)+16*alpha-d-2*determinant1(1)
    residue_B=numerator1(-1)-16*beta-d-2*determinant1(-1)
    # All terms involving N0 itself have order >=8 at e. To order six,
    # the variation of eta=N/(ell*Jacobian) is numerator1/(z^2-1).
    inverse=-(1+z*z+z**4+z**6)
    eta1=(numerator1*inverse) % z**7
    return vector(k,[*regularity,residue_A,residue_B,*[eta1[j] for j in range(7)]])

columns=[];labels=[];inc_columns=[]
for degree,which in [(2,"q"),(3,"c")]:
    for mon in ((H+U+V+W)**degree).monomials():
        q1,c1=(mon,0) if which=="q" else (0,mon)
        columns.append(conditions(q1,c1,0,0,0))
        labels.append(which+":"+str(mon))
        inc_columns.append(incidence.obstruction(q1,c1)[0])
for a,b,label in [(1,0,"mark_A"),(0,1,"mark_B")]:
    columns.append(conditions(0,0,a,b,0));labels.append(label)
    inc_columns.append(incidence.obstruction(0,0,a,b)[0])
for mon in numerator_basis:
    columns.append(conditions(0,0,0,0,mon));labels.append("N:"+str(mon))
    inc_columns.append(vector(k,[0,0]))
system=matrix(k,columns).transpose()
inc_matrix=matrix(k,inc_columns).transpose()
rhs=vector(k,[0]*12+[0,1,0,15,0,12,0])
print("CRITICAL_DIFFERENTIAL_LINEAR_SYSTEM",system.nrows(),system.ncols(),
      "rank",system.rank(),"augmented_rank",system.augment(matrix(k,len(rhs),1,rhs)).rank(),flush=True)
solution=system.solve_right(rhs)
assert system*solution==rhs
kernel=system.right_kernel().basis_matrix()
inc_image=inc_matrix*kernel.transpose()
print("INCIDENCE_VARIATION_ON_CRITICAL_KERNEL_RANK",inc_image.rank(),flush=True)
print("INCIDENCE_AT_ONE_CRITICAL_SOLUTION",list(inc_matrix*solution),flush=True)
if inc_image.rank()==0 and inc_matrix*solution==0:
    row_witnesses=[]
    for row in inc_matrix.rows():
        witness=system.transpose().solve_right(row)
        assert witness*system==row
        assert witness*rhs==0
        row_witnesses.append(witness)
    print("INCIDENCE_ROW_CERTIFICATES",[list(v) for v in row_witnesses],flush=True)
    print("PASS_CRITICAL_FIBRE_FORCES_FIRST_ORDER_INCIDENCE",flush=True)
else:
    print("NOT_FORCED_BY_THESE_DIFFERENTIAL_JETS",flush=True)
print("ONE_SOLUTION",[(label,int(a)) for label,a in zip(labels,solution) if a],flush=True)
print("SCOPE: first-order source/differential conditions, not an all-orders Hurwitz certificate")
