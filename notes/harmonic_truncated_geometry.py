"""Truncated mixed-characteristic arithmetic and canonical incidence.

Helpers for research scripts. Coefficients lie in Q(pi), pi^2=-23;
arithmetic retains 23, and only p-unit denominators are allowed.
"""
from sage.all import GF, QQ, PolynomialRing, matrix, vector


class Arithmetic:
    def __init__(self,K,precision):
        self.K=K;self.pi=K.gen();self.N=precision;self.k=GF(23)
        self.S=PolynomialRing(K,"z");self.z=self.S.gen()

    def residue(self,a,modulus):
        if modulus==1:return 0
        a=QQ(a);assert a.denominator()%23
        value=int(a.numerator())*pow(int(a.denominator()),-1,modulus)%modulus
        return value if value<=modulus//2 else value-modulus

    def coefficient(self,a):
        cs=list(self.K(a))+[QQ.zero(),QQ.zero()]
        return self.K(self.residue(cs[0],23**((self.N+1)//2)))+self.pi*self.residue(cs[1],23**(self.N//2))

    def cut(self,f):
        return self.S({j:self.coefficient(a) for j,a in self.S(f).dict().items()})

    def digit(self,a,n):
        cs=list(self.K(a)/self.pi**n)+[QQ.zero(),QQ.zero()]
        assert all(not b or b.valuation(23)>=0 for b in cs[:2])
        return self.k(self.residue(cs[0],23))

    def digit_polynomial(self,f,n):
        T=PolynomialRing(self.k,"z")
        return T({j:self.digit(a,n) for j,a in self.S(f).dict().items()})

    def quo_rem(self,f,g):
        f=self.cut(f);g=self.cut(g);quotient=self.S.zero()
        lead=g.leading_coefficient();assert self.digit(lead,0)
        while f and f.degree()>=g.degree():
            term=self.coefficient(f.leading_coefficient()/lead)*self.z**(f.degree()-g.degree())
            quotient=self.cut(quotient+term)
            f=self.cut(f-term*g)
        return quotient,f

    def weierstrass(self,f,degree):
        """Distinguished factor modulo pi^N, from a sufficiently long series.

        For the input truncated at z^count, count>=degree*N ensures that
        omitted higher terms cannot change the factor's coefficients.
        """
        f=self.cut(f);R=self.z**degree
        special=self.digit_polynomial(f,0)
        assert special.valuation()==degree
        T=special.parent();z0=T.gen()
        unit0=special//z0**degree
        inverse=unit0.inverse_mod(z0**degree)
        for n in range(1,self.N):
            _,rem=self.quo_rem(f,R)
            correction=(self.digit_polynomial(rem,n)*inverse)%(z0**degree)
            R=self.cut(R+self.pi**n*self.S([int(a) for a in correction.list()]))
        assert self.quo_rem(f,R)[1]==0
        return R


def osculating_remainder(q,c,marks,affine,precision):
    """Remainder c|line modulo q|line, using actual order-two canonical jets."""
    A=Arithmetic(q.base_ring(),precision);K=A.K;z=A.z;pi=A.pi
    hh,uu=map(A.S,affine)
    functions=[hh,uu,z,A.S.one()]
    planes=[]
    for mark in marks:
        jets=matrix(K,[[A.coefficient(f.derivative(d)(mark)) for f in functions] for d in range(3)])
        square=jets.matrix_from_columns([1,2,3])
        inverse=matrix(A.k,[[A.digit(a,0) for a in row] for row in square.rows()]).inverse()
        solution=vector(K,4);solution[0]=1
        for n in range(precision):
            residual=jets*solution
            correction=-inverse*vector(A.k,[A.digit(a,n) for a in residual])
            for j,a in enumerate(correction,start=1):solution[j]=A.coefficient(solution[j]+pi**n*int(a))
        assert all(A.coefficient(a)==0 for a in jets*solution)
        planes.append(solution)
    a,b=planes
    ul=A.cut(-((a[2]-b[2])*z+a[3]-b[3])/(a[1]-b[1]))
    hl=A.cut(-a[1]*ul-a[2]*z-a[3])
    ql=A.cut(q(hl,ul,z,1));cl=A.cut(c(hl,ul,z,1))
    _,remainder=A.quo_rem(cl,ql)
    assert remainder.degree()<=1
    return remainder,(hl,ul,ql,cl)
