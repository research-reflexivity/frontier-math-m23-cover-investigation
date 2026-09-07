"""Finite canonical-ring, marked-jet, and effective-fibre calculations.

No central series of growing z-degree is used. Global equality is tested
in the 57-dimensional degree-ten canonical ring. The ramification factor
and its square are handled in finite free polynomial quotient algebras.
These routines define a candidate off the global-function locus; that
candidate is never asserted to be a genuine cover merely by construction.
"""
from sage.all import GF,PolynomialRing,matrix,vector
from harmonic_truncated_geometry import Arithmetic
from certify_harmonic_universal_canonical_ring import standard_monomials,monomial
from investigate_harmonic_second_incidence_order import mixed,complement

k=GF(23)


class FiniteGeometry:
    def __init__(self,c,marks,precision):
        self.N=precision;self.K=mixed.K;self.pi=mixed.pi
        self.A=Arithmetic(self.K,precision);self.S=self.A.S;self.z=self.A.z
        self.P=PolynomialRing(self.K,["W","U","H","V"],order="lex")
        W,U,H,V=self.P.gens();self.W,self.U,self.H,self.V=W,U,H,V
        self.c=self.P(c(H,U,V,W));self.C=self.c-U*W**2+V**3
        self.q=H*W-U**2-U*V-V**2
        ss=U**2+U*V+V**2;tt=self.C-V**3
        self.gb=[self.q,self.c,U*W*ss+H*tt,U*ss**2+H**2*tt]
        self.leaders=[tuple(g.exponents()[0]) for g in self.gb]
        self.basis=[self.P(f(H,U,V,W)) for f in mixed.basis]
        self.decics=standard_monomials(self.P,10)
        assert len(self.decics)==57
        self.marks=marks

    def cut_form(self,f):
        return self.P({m:self.A.coefficient(a) for m,a in self.P(f).dict().items()})

    def normal_form(self,f):
        data={tuple(e):self.A.coefficient(a) for e,a in f.dict().items() if self.A.coefficient(a)}
        out={};reducers=[g.dict() for g in self.gb]
        while data:
            exp=max(data);coefficient=data.pop(exp)
            for g,lead in zip(reducers,self.leaders):
                if all(a>=b for a,b in zip(exp,lead)):
                    delta=tuple(a-b for a,b in zip(exp,lead))
                    for power,a in g.items():
                        power=tuple(power)
                        if power==lead:continue
                        target=tuple(a+b for a,b in zip(delta,power))
                        value=self.A.coefficient(data.get(target,0)-coefficient*a)
                        if value:data[target]=value
                        else:data.pop(target,None)
                    break
            else:
                out[exp]=coefficient
        return self.P(out)

    def local_jets(self,mark,count=23):
        """Implicit curve directly in O[t]/t^count at a marked section."""
        z=self.z;A=self.A
        def cut(f):return A.cut(self.S(f).truncate(count))
        zz=z+mark;uu=cut(zz**3)
        for _ in range(self.N+1):
            hh=cut(uu**2+uu*zz+zz**2)
            uu=cut(zz**3-self.C(1,uu,hh,zz))
        hh=cut(uu**2+uu*zz+zz**2)
        assert cut(self.c(1,uu,hh,zz))==0
        values=[cut(f(1,uu,hh,zz)) for f in self.basis]
        return hh,uu,values

    def kernel_sections(self):
        sections=[];self.mark_jets=[]
        for mark,data in zip(self.marks,mixed.first.jet_data):
            hh,uu,values=self.local_jets(mark)
            self.mark_jets.append((hh,uu))
            _,_,kernel,pivots,inverse,_=data
            jets=matrix(self.K,[[f[j] for f in values] for j in range(23)])
            current=[]
            for row in kernel:
                coefficients=vector(self.K,[int(a) for a in row])
                for n in range(1,self.N):
                    residual=jets*coefficients
                    correction=-inverse*vector(k,[self.A.digit(a,n) for a in residual])
                    for j,a in zip(pivots,correction):coefficients[j]+=self.pi**n*int(a)
                    coefficients=vector(self.K,[self.A.coefficient(a) for a in coefficients])
                assert all(self.A.coefficient(a)==0 for a in jets*coefficients)
                current.append(self.cut_form(sum(a*f for a,f in zip(coefficients,self.basis))))
            sections.append(current)
        self.sections=sections
        return sections

    def multiplier_candidate(self):
        SA,SB=self.kernel_sections()
        products=[[self.normal_form(a*b) for b in SB] for a in SA]
        coordinates=[[vector(self.K,[f.monomial_coefficient(m) for m in self.decics])
                      for f in row] for row in products]
        columns=[]
        # Three pairs suffice: S_A,0 is a non-zero-divisor in the
        # canonical ring, so their vanishing implies the other three.
        self.pairs=[(0,1),(0,2),(0,3)]
        for row in range(4):
            for col in range(4):
                column=[]
                for a,b in self.pairs:
                    value=vector(self.K,57)
                    if b==row:value+=coordinates[a][col]
                    if a==row:value-=coordinates[b][col]
                    column.extend(value)
                columns.append(column)
        system=matrix(self.K,columns).transpose()
        special=matrix(k,[[self.A.digit(a,0) for a in row] for row in system.rows()])
        assert special.rank()==15
        pivots=special.pivots();assert tuple(pivots)==tuple(range(15))
        injective=special.matrix_from_columns(pivots)
        rows=injective.transpose().pivots()
        inverse=injective.matrix_from_rows(rows).inverse()
        coefficients=vector(self.K,matrix.identity(self.K,4).list())
        for n in range(1,self.N):
            residual=system*coefficients
            rhs=-vector(k,[self.A.digit(residual[j],n) for j in rows])
            correction=inverse*rhs
            for j,a in zip(pivots,correction):coefficients[j]+=self.pi**n*int(a)
            coefficients=vector(self.K,[self.A.coefficient(a) for a in coefficients])
        residual=[self.A.coefficient(a) for a in system*coefficients]
        assert all(residual[j]==0 for j in rows)
        M=matrix(self.K,4,4,coefficients)
        DD=[self.cut_form(sum(a*f for a,f in zip(row,SB))) for row in M.rows()]
        chosen=next(j for j,f in enumerate(DD) if self.A.digit(f(1,0,0,0),0))
        unit=self.A.coefficient(SA[chosen](1,0,0,0)/DD[chosen](1,0,0,0))
        self.numerator=self.cut_form(SA[chosen]/unit);self.denominator=DD[chosen]
        assert self.numerator(1,0,0,0)==self.denominator(1,0,0,0)
        self.multiplier=M;self.global_residual=residual
        self.global_matrix=system
        self.selected_multiplier_rows=rows
        return residual

    def quotient(self,R):
        """Operations in the finite free algebra O[z]/R."""
        def cut(f):return self.A.quo_rem(f,R)[1]
        def inverse(f):
            special=self.A.digit_polynomial(cut(f),0)
            z0=special.parent().gen()
            inv=special.inverse_mod(z0**R.degree())
            out=self.S([int(a) for a in inv.list()])
            precision=1
            while precision<self.N:
                out=cut(out*cut(2-cut(f*out)));precision*=2
            assert cut(f*out)==1
            return out
        return cut,inverse

    def curve_in_quotient(self,R):
        cut,inverse=self.quotient(R);z=self.z
        uu=cut(z**3)
        for _ in range(self.N+1):
            hh=cut(uu**2+uu*z+z**2)
            uu=cut(z**3-self.C(1,uu,hh,z))
        hh=cut(uu**2+uu*z+z**2)
        assert cut(self.c(1,uu,hh,z))==0
        return hh,uu,cut,inverse

    def eta_in_quotient(self,R):
        hh,uu,cut,inverse=self.curve_in_quotient(R);z=self.z
        def ev(f):return cut(f(1,uu,hh,z))
        fU=cut(ev(self.c.derivative(self.U))+(2*uu+z)*ev(self.c.derivative(self.H)))
        fz=cut(ev(self.c.derivative(self.V))+(uu+2*z)*ev(self.c.derivative(self.H)))
        du=cut(-fz*inverse(fU));dh=cut((2*uu+z)*du+uu+2*z)
        def dlog(f):
            derivative=cut(ev(f.derivative(self.H))*dh+ev(f.derivative(self.U))*du
                           +ev(f.derivative(self.V)))
            return cut(derivative*inverse(ev(f)))
        dlog_beta=cut(dlog(self.numerator)-dlog(self.denominator))
        assert all(self.A.digit(a,0)==self.A.digit(a,1)==0 for a in dlog_beta.list())
        short=Arithmetic(self.K,self.N-2)
        return short.cut(-dlog_beta/self.pi**2)

    def ramification_factor(self):
        z=self.z;short=Arithmetic(self.K,self.N-2)
        # This verifies the complete special unit needed below, not
        # just the multiplicity of its zero.
        special=self.A.digit_polynomial(self.eta_in_quotient(z**16),0)
        z0=special.parent().gen()
        assert special==(-2*sum(z0**j for j in range(8,16,2)))
        inverse_unit=(z0**2-1)/k(2)
        R=z**8
        for n in range(1,self.N-2):
            remainder=self.eta_in_quotient(R)
            correction=(self.A.digit_polynomial(remainder,n)*inverse_unit)%(z0**8)
            R=short.cut(R+self.pi**n*self.S([int(a) for a in correction.list()]))
        assert self.eta_in_quotient(R)==0
        return R

    def effective_fibre(self,R):
        hh,uu,cut,_=self.curve_in_quotient(R**2)
        NN=cut(self.numerator(1,uu,hh,self.z));DD=cut(self.denominator(1,uu,hh,self.z))
        scalar=self.A.coefficient((NN[0]-DD[0])/DD[0])
        E=cut(NN-(1+scalar)*DD)
        self.fibre_denominator=DD
        return E

    def certify_fibre_elimination(self,R,E):
        """The exact derivative identity eliminates coefficients E1,...,E8."""
        DD=self.fibre_denominator;z=self.z
        cut,_=self.quotient(R)
        assert cut(DD*E.derivative()-DD.derivative()*E)==0
        columns=[cut(DD*j*z**(j-1)-DD.derivative()*z**j) for j in range(1,16)]
        L=matrix(self.K,[[f[j] for f in columns] for j in range(8)])
        low=L[:,:8];high=L[:,8:]
        special=matrix(k,[[self.A.digit(a,0) for a in row] for row in low.rows()])
        assert special.det()==k(40320)*self.A.digit(DD[0],0)**8
        assert special.det()
        rhs=-high*vector(self.K,[E[j] for j in range(9,16)])
        answer=vector(self.K,8);inverse=special.inverse()
        for n in range(self.N):
            residual=low*answer-rhs
            correction=-inverse*vector(k,[self.A.digit(a,n) for a in residual])
            answer=vector(self.K,[self.A.coefficient(a+self.pi**n*int(b)) for a,b in zip(answer,correction)])
        assert list(answer)==[E[j] for j in range(1,9)]
        return int(special.det())

    def incidence(self):
        planes=[]
        for mark,(hh,uu) in zip(self.marks,self.mark_jets):
            functions=[hh,uu,self.z+mark,self.S.one()]
            jets=matrix(self.K,[[f[j] for f in functions] for j in range(3)])
            square=jets.matrix_from_columns([1,2,3])
            inverse=matrix(k,[[self.A.digit(a,0) for a in row] for row in square.rows()]).inverse()
            solution=vector(self.K,[1,0,0,0])
            for n in range(self.N):
                residual=jets*solution
                correction=-inverse*vector(k,[self.A.digit(a,n) for a in residual])
                for j,a in enumerate(correction,1):solution[j]=self.A.coefficient(solution[j]+self.pi**n*int(a))
            assert all(self.A.coefficient(a)==0 for a in jets*solution)
            planes.append(solution)
        a,b=planes;z=self.z
        ul=self.A.cut(-((a[2]-b[2])*z+a[3]-b[3])/(a[1]-b[1]))
        hl=self.A.cut(-a[1]*ul-a[2]*z-a[3])
        ql=self.A.cut(self.q(1,ul,hl,z));cl=self.A.cut(self.c(1,ul,hl,z))
        rem=self.A.quo_rem(cl,ql)[1]
        assert rem.degree()<=1
        return rem


def fixed_family(point,precision):
    from certify_harmonic_fixed_hensel_system import ROWS
    pi=mixed.pi;c=mixed.c0
    for order,row in enumerate(ROWS,1):
        c+=pi**order*sum(int(a)*f for a,f in zip(row,complement))
    c+=pi**3*point[0]*complement[8]
    c+=pi**4*sum(a*f for a,f in zip(point[1:9],complement[:8]))
    marks=[1+pi+9*pi**2+pi**3*point[9],-1+pi-9*pi**2+pi**3*point[10]]
    return FiniteGeometry(c,marks,precision)
