#!/usr/bin/env sage-python
"""Identify cyclic Witt designs in explicit characteristic-23 tail germs.

The test uses the degree-seven heptad polynomial. On the correct M23
design it is a globally invariant function, integral over the finite
target line, with pole order less than one at infinity: it must be
constant. A nonconstant Laurent coefficient therefore excludes a design.
The input that each tail has geometric group M23 remains the separate
exact function-field Galois certificate. No generic reconstructed cover
is read and no conclusion about all-orders osculating incidence is made.
"""
from itertools import combinations
from sage.all import GF, PolynomialRing, PowerSeriesRing, binomial

k=GF(23)
B=PolynomialRing(GF(2),"x");x=B.gen()
generator=x**11+x**9+x**7+x**6+x**5+x+1
assert (x**23-1)%generator==0
blocks=[]
for mask in range(1<<12):
    word=generator*sum(B.gen()**j for j in range(12) if (mask>>j)&1)
    support=tuple(j for j in range(23) if word[j])
    if len(support)==7:blocks.append(support)
assert len(blocks)==253
four_sets=[s for block in blocks for s in combinations(block,4)]
assert len(set(four_sets))==len(four_sets)==8855
other=[tuple(sorted((-j)%23 for j in block)) for block in blocks]
assert set(blocks).isdisjoint(other)
for family in [blocks,other]:
    as_set=set(family)
    assert {tuple(sorted((j+1)%23 for j in b)) for b in family}==as_set
    assert {tuple(sorted((2*j)%23 for j in b)) for b in family}==as_set
print("TWO_CYCLIC_WITT_DESIGNS",len(blocks),"each; disjoint; S(4,7,23) checked",flush=True)


def test_tail(name,pole,jump,lam=1):
    # Both the first possible nonconstant invariant term T^-1 and a
    # second such term are retained. T has pole 253 in this tame
    # uniformizer of the degree-253 local normal closure.
    max_power=3
    # The primitive tail has a mu_7 target-scaling symmetry. For these
    # degree-7r expressions its first possible negative target power is
    # T^-7, so stopping at T^-2 would miss the relevant coefficient.
    precision=7*max_power*pole+(8 if name=="primitive7" else 2)*253+1
    S=PowerSeriesRing(k,"t",default_prec=precision);t=S.gen()
    if name=="primitive7":
        unit=S.one();terms=[(9,k(5),0),(2,k(6),77)]
        # Common central label: delta z = rho*z^(4/11)*(1-z^2)*i,
        # rho^11=pi. With z=s*x, s^7=pi, the primitive leading
        # displacement is delta x = i*x^(4/11).
        label_scale=k(1)
    else:
        lam=k(lam);unit=S.one()
        for _ in range(10):
            unit-=(unit**5*(unit-1)-lam*t**165)/(unit**4*(6*unit-5))
        assert unit**5*(unit-1)==lam*t**165
        terms=[(13,k(7)*lam**2,0),(8,k(2)*lam**3,165),
               (3,k(17)*lam**4,330)]
        # At the new component z~r^-1/xi, r^15=pi, xi=X^2/Y.
        # The same central label gives delta xi=i*xi^(-4/11),
        # hence delta X=3*i*t^-18 when xi=t^-11.
        label_scale=k(3)
    # Residual equation for X_i=t^-pole*(unit+t^jump*v_i).
    # Its reduction is v^23-v, so all 23 roots are simple and lift.
    coefficients={23:S.one()}
    for degree,coefficient,offset in terms:
        for r in range(1,degree+1):
            value=coefficient*binomial(degree,r)*unit**(degree-r)*t**(offset+jump*(r-1))
            coefficients[r]=coefficients.get(r,S.zero())+value
    assert coefficients[1][0]==-1
    def equation(v):return sum(a*v**r for r,a in coefficients.items())
    def derivative(v):return sum(r*a*v**(r-1) for r,a in coefficients.items() if r%23)
    roots=[]
    for index in range(23):
        v=S(label_scale*index)
        for _ in range(11):v-=equation(v)/derivative(v)
        assert equation(v)==0
        assert v[0]==label_scale*index
        roots.append(unit+t**jump*v)
    print("LOCAL_ROOTS",name,"lambda",lam,"precision",precision,"all 23 verified",flush=True)
    target_unit=(1+5*t**154+6*t**231 if name=="primitive7"
                 else 2*(unit**8+15*k(lam)*t**165*unit**3))
    outputs=[]
    for number,family in enumerate([blocks,other]):
        values=[S.zero() for _ in range(max_power)]
        for block in family:
            term=S.one()
            for index in block:term*=roots[index]
            for power in range(1,max_power+1):values[power-1]+=term**power
        any_residual=False
        for power,value in enumerate(values,start=1):
            bound=7*power*pole
            # A globally invariant sum is integral over k[T] and has
            # target degree at most floor(bound/253). Remove that unique
            # possible polynomial before inspecting the Laurent tail.
            polynomial=[]
            for degree in range(bound//253,-1,-1):
                coefficient=value[bound-253*degree]/target_unit[0]**degree
                polynomial.append((degree,int(coefficient)))
                value-=coefficient*t**(bound-253*degree)*target_unit**degree
            residual=[(n-bound,int(value[n])) for n in range(precision) if value[n]]
            if power==1 and residual:
                assert residual[0]==((1771,15) if name=="primitive7" else (759,5))
            print("HEPTAD_INVARIANT",name,"lambda",lam,"design",number,"power",power,
                  "polynomial",polynomial,"residual_count",len(residual),"terms",residual[:16],flush=True)
            any_residual|=bool(residual)
        outputs.append(any_residual)
    # No all-orders constancy is inferred from a finite zero expansion.
    # If exactly one candidate is excluded, the other is identified using
    # the separately certified M23 group and the two-design group fact.
    print("EXCLUDED_BY_NONCONSTANT_TERM",name,"lambda",lam,
          [j for j,terms in enumerate(outputs) if terms],flush=True)
    return outputs


primitive=test_tail("primitive7",11,7)
positive=test_tail("new15",33,15,1)
negative=test_tail("new15",33,15,-1)
assert primitive==positive==[False,True]
assert negative==[True,False]
print("PASS_COMMON_CENTRAL_LABEL_WITT_ALIGNMENT: primitive and lambda=+1 match; lambda=-1 differs")
print("SCOPE: first-order selection under the stated integral/annular model hypotheses; not all-orders incidence")
