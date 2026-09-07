#!/usr/bin/env sage-python
"""A universal monic Groebner basis for the canonical deformation slice.

All identities are verified over Z[a0,...,a8], not at sampled covers.
The result supplies finite, exact global function equations; it does
not prove the osculating incidence relation.
"""
from itertools import combinations
from sage.all import ZZ,PolynomialRing


def monomial(P,exponents):
    ans=P.one()
    for x,n in zip(P.gens(),exponents):ans*=x**n
    return ans


def monic_remainder(f,basis,trace=False):
    P=f.parent();remainder=P.zero();quotients=[P.zero() for _ in basis]
    leaders=[tuple(g.exponents()[0]) for g in basis]
    assert all(g.lc()==1 for g in basis)
    while f:
        exp=tuple(f.exponents()[0]);coeff=f.lc()
        for j,(g,lead) in enumerate(zip(basis,leaders)):
            if all(a>=b for a,b in zip(exp,lead)):
                term=coeff*monomial(P,[a-b for a,b in zip(exp,lead)])
                quotients[j]+=term;f-=term*g
                break
        else:
            term=coeff*monomial(P,exp)
            remainder+=term;f-=term
    return (remainder,quotients) if trace else remainder


def standard_monomials(P,degree):
    # Variable ordering is W,U,H,V. Leading monomials are HW, UW^2,
    # U^3W, and U^5. This explicit description also works in degrees
    # zero through three, where the usual Hilbert polynomial does not.
    W,U,H,V=P.gens();answer=[]
    for w in range(degree+1):
        for u in range(degree-w+1):
            for h in range(degree-w-u+1):
                v=degree-w-u-h
                if (h and w) or (u and w>=2) or (u>=3 and w) or u>=5:continue
                answer.append(W**w*U**u*H**h*V**v)
    return answer


def universal_ring():
    A=PolynomialRing(ZZ,["a%d"%j for j in range(9)])
    P=PolynomialRing(A,["W","U","H","V"],order="lex")
    W,U,H,V=P.gens()
    mons=[H**3,H**2*U,H*U**2,U**3,H**2*V,H*U*V,U**2*V,H*V**2,U*V**2]
    C=sum(a*m for a,m in zip(A.gens(),mons))
    S=U**2+U*V+V**2;T=C-V**3
    q=H*W-S;c=U*W**2+T
    g4=U*W*S+H*T;g5=U*S**2+H**2*T
    return A,P,[q,c,g4,g5]


def main():
    A,P,basis=universal_ring();W,U,H,V=P.gens();q,c,g4,g5=basis
    S=U**2+U*V+V**2
    assert g4==H*c-U*W*q
    assert g5==H*g4-U*S*q
    assert [g.lm() for g in basis]==[H*W,U*W**2,U**3*W,U**5]
    for i,j in combinations(range(4),2):
        expi,expj=map(lambda f:tuple(f.exponents()[0]),[basis[i],basis[j]])
        lcm=tuple(max(a,b) for a,b in zip(expi,expj))
        spoly=monomial(P,[a-b for a,b in zip(lcm,expi)])*basis[i]
        spoly-=monomial(P,[a-b for a,b in zip(lcm,expj)])*basis[j]
        remainder,quotients=monic_remainder(spoly,basis,trace=True)
        assert not remainder
        assert spoly==sum(a*b for a,b in zip(quotients,basis))
        print("UNIVERSAL_BUCHBERGER_PAIR",i,j,"PASS",
              "quotient_term_counts",[len(a.dict()) for a in quotients],flush=True)
    for degree in range(13):
        mons=standard_monomials(P,degree)
        expected=[1,4,9][degree] if degree<3 else 6*degree-3
        assert len(mons)==expected
        assert all(monic_remainder(m,basis)==m for m in mons)
        print("CANONICAL_HILBERT_FUNCTION",degree,len(mons),flush=True)
    # A normal form and an exact ideal-membership witness for a generic
    # product of quintic monomials, with no denominator introduced.
    product=(H**5+U**4*V+W**5)*(H**3*U**2+V**5+U*V**3*W)
    remainder,quotients=monic_remainder(product,basis,trace=True)
    assert product==remainder+sum(a*b for a,b in zip(quotients,basis))
    assert all(m in standard_monomials(P,10) for m in remainder.monomials())
    print("UNIVERSAL_DECIC_NORMAL_FORM",len(remainder.dict()),"terms",
          "coefficient_degree",max(a.total_degree() for a in remainder.coefficients()),flush=True)
    print("PASS_UNIVERSAL_CANONICAL_RING_OVER_Z_PARAMETERS",flush=True)
    print("SCOPE: finite global function equations; no incidence ideal membership asserted",flush=True)


if __name__=="__main__":main()
