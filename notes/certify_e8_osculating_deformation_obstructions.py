#!/usr/bin/env sage-python
"""Special E8/A2+A6 incidence and the harmonic lifting obstruction.

This reads NO reconstructed cover coefficients. It does not impose the
characteristic-zero three-point-map equations and is not a Hurwitz lifting
certificate. All deformation computations are over the dual numbers in
characteristic 23, also valid modulo pi^2 when pi^2 = -23.
"""
from sage.all import GF, PolynomialRing, matrix, vector

k = GF(23)
E = GF(23**2, "i", modulus=PolynomialRing(k, "v").gen()**2+1)
i = E.gen()


def special_incidence(b):
    """a=1, e=0, cusp=infinity; b runs through mu_4 minus {1}."""
    R = PolynomialRing(E, "z")
    z = R.gen()
    a = E.one()
    assert b**4 == 1 and b != a
    quotient, remainder = (z**8).quo_rem((z-a)*(z-b))
    # Opposite residues, as required by (d beta/beta)/23.
    assert a**8 == b**8
    assert remainder.degree() == 0
    u = quotient[5]/2
    f3 = z**3+u*z**2
    f6 = quotient-quotient[3]*f3-quotient[1]*z-quotient[0]
    basis = [R.one(), z, f3, f6]
    assert f6[5] == 2*u
    products = [(r,s) for r in range(4) for s in range(r,4)]
    mat = matrix(E, [[(basis[r]*basis[s])[d] for r,s in products]
                     for d in range(13)])
    assert mat.right_kernel().dimension() == 1
    P = PolynomialRing(E, ["x0", "x1", "x2", "x3"])
    x0,x1,x2,x3 = P.gens()
    quadric = sum(c*P.gen(r)*P.gen(s) for c,(r,s) in
                  zip(mat.right_kernel().basis()[0], products))
    cubic = x0*x0*x2-x1**3-u*x0*x1*x1
    assert quadric(*basis) == cubic(*basis) == 0
    planes = []
    for mark in [a,b]:
        jets = matrix(E, [[f.derivative(d)(mark) for f in basis]
                           for d in range(3)])
        assert jets.rank() == 3
        plane = jets.right_kernel().basis()[0]
        section = sum(c*f for c,f in zip(plane,basis))
        assert section.derivative(3)(mark) != 0
        assert section(b if mark == a else a) != 0
        planes.append(plane)
    line = matrix(E,planes).right_kernel().basis()
    assert len(line) == 2
    S = PolynomialRing(E,["r","s"])
    r,s = S.gens()
    coords = [r*line[0][j]+s*line[1][j] for j in range(4)]
    qline, cline = quadric(*coords), cubic(*coords)
    common = qline.gcd(cline)
    length = common.total_degree()
    expected = 2 if b == -1 else 0
    assert length == expected, (b, length)
    # Direct normalization check excludes an intersection lost at the cusp.
    sections = [sum(c*f for c,f in zip(plane,basis)) for plane in planes]
    assert sections[0].degree() == sections[1].degree() == 6
    assert sections[0].gcd(sections[1]).degree() == length
    resultant = sections[0].monic().resultant(sections[1].monic())
    assert bool(resultant) == (length == 0)
    assert resultant == (E.zero() if b == -1 else E(16))
    print("SPECIAL_E8", "b", b, "basis", basis,
          "common_degree", length, "monic_section_resultant",resultant,flush=True)
    return length


assert [special_incidence(b) for b in [-E.one(),i,-i]] == [2,0,0]


def a2_a6_special_incidence(b):
    """Semigroups <2,3>, <2,7> at e=0, q=infinity; b^3=1."""
    R = PolynomialRing(E,"z")
    z = R.gen()
    assert b**3 == 1 and b != 1
    quotient, remainder = (z**6).quo_rem((z-1)*(z-b))
    assert remainder.degree() == 0
    u = quotient[3]/2
    f2 = z*z+u*z
    f4 = quotient-quotient[2]*f2-quotient[0]
    assert f4[3] == 2*u
    # Multiply the differential basis z^-2,1,f2,f4 by z^2.
    basis = [R.one(),z*z,z*z*f2,z*z*f4]
    P = PolynomialRing(E,["x0","x1","x2","x3"])
    def relations(degree):
        mons = (sum(P.gens())**degree).monomials()
        mat = matrix(E,[[m(*basis)[j] for m in mons]
                        for j in range(6*degree+1)])
        return [sum(c*m for c,m in zip(row,mons))
                for row in mat.right_kernel().basis()]
    quadrics, cubics = relations(2), relations(3)
    assert len(quadrics) == 1 and len(cubics) == 5
    q = quadrics[0]
    c = next(f for f in cubics if f.reduce([q]) != 0)
    assert q(*basis) == c(*basis) == 0
    planes, sections = [], []
    for mark in [E.one(),b]:
        jets = matrix(E,[[f.derivative(j)(mark) for f in basis]
                         for j in range(3)])
        assert jets.rank() == 3
        plane = jets.right_kernel().basis()[0]
        section = sum(a*f for a,f in zip(plane,basis))
        assert section.derivative(3)(mark) != 0
        assert section(b if mark == 1 else E.one()) != 0
        planes.append(plane)
        sections.append(section)
    line = matrix(E,planes).right_kernel().basis()
    S = PolynomialRing(E,["r","s"])
    r,s = S.gens()
    coords = [r*line[0][j]+s*line[1][j] for j in range(4)]
    common = q(*coords).gcd(c(*coords))
    length = common.total_degree()
    print("SPECIAL_A2_A6", "b",b,"basis",basis,
          "common_degree",length,"normalization_gcd",sections[0].gcd(sections[1]),
          "section_degrees",[f.degree() for f in sections],flush=True)
    assert length == 0
    assert sections[0].gcd(sections[1]).degree() == 0
    assert max(f.degree() for f in sections) == 6
    resultant = sections[0].monic().resultant(sections[1].monic())
    assert resultant == E(21)
    print("A2_A6_MONIC_SECTION_RESULTANT",resultant,flush=True)
    return length


T = PolynomialRing(E,"t")
assert [a2_a6_special_incidence(root) for root,mult in
        (T.gen()**2+T.gen()+1).roots()] == [0,0]

# An arbitrary first-order canonical deformation at the harmonic fibre.
# Affine chart W=1: q=h-u^2-uz-z^2, c=u-z^3. Marks have
# z(A)=1+epsilon*alpha and z(B)=-1+epsilon*beta.
R = PolynomialRing(k,"z")
z = R.gen()
F = z**6+z**4+z*z
P = PolynomialRing(k,["H","U","V","W"])
H,U,V,W = P.gens()


def obstruction(q1, c1, alpha=0, beta=0):
    q1, c1 = P(q1), P(c1)
    uu1 = -c1(F,z**3,z,1)
    hh1 = (2*z**3+z)*uu1-q1(F,z**3,z,1)
    base = [F,z**3,z,R.one()]
    perturb = [hh1,uu1,R.zero(),R.zero()]
    plane_data = []
    for mark,motion in [(k(1),k(alpha)),(-k.one(),k(beta))]:
        jets = matrix(k,[[f.derivative(d)(mark) for f in base]
                        for d in range(3)])
        moving = matrix(k,[[g.derivative(d)(mark)
                           +motion*f.derivative(d+1)(mark)
                           for f,g in zip(base,perturb)] for d in range(3)])
        square = jets.matrix_from_columns([1,2,3])
        p0 = vector(k,[1,*square.solve_right(-jets.column(0))])
        p1 = vector(k,[0,*square.solve_right(-moving*p0)])
        assert jets*p0 == 0 and jets*p1+moving*p0 == 0
        plane_data.append((p0,p1))
    a0,a1 = plane_data[0]
    b0,b1 = plane_data[1]
    assert a0 == vector(k,[1,8,10,2])
    assert b0 == vector(k,[1,-8,-10,2])
    ul0,hl0 = 16*z,R(-2)
    ul1 = -((a1[1]-b1[1])*ul0+(a1[2]-b1[2])*z
            +(a1[3]-b1[3]))/(a0[1]-b0[1])
    hl1 = -a0[1]*ul1-a1[1]*ul0-a1[2]*z-a1[3]
    ql0 = hl0-ul0**2-ul0*z-z*z
    cl0 = ul0-z**3
    assert ql0 == 3*(z*z-16) and cl0 == -z*(z*z-16)
    ql1 = hl1-(2*ul0+z)*ul1+q1(hl0,ul0,z,1)
    cl1 = ul1+c1(hl0,ul0,z,1)
    residual = (cl1+z/k(3)*ql1) % (z*z-16)
    return vector(k,[residual[0],residual[1]]), (hl1,ul1,ql1,cl1)


directions = []
labels = []
for degree,which in [(2,"q"),(3,"c")]:
    for mon in (sum(P.gens())**degree).monomials():
        directions.append(obstruction(mon if which == "q" else 0,
                                      mon if which == "c" else 0)[0])
        labels.append(which+":"+str(mon))
for alpha,beta,label in [(1,0,"mark_A"),(0,1,"mark_B")]:
    directions.append(obstruction(0,0,alpha,beta)[0])
    labels.append(label)
linear_map = matrix(k,directions).transpose()
assert len(labels) == 32
assert linear_map.rank() == 2
mark_jacobian = linear_map.matrix_from_columns([30,31])
assert mark_jacobian == matrix(k,[[3,3],[7,-7]])
assert mark_jacobian.det() == 4
print("MOVING_MARK_JACOBIAN",list(mark_jacobian.rows()),
      "determinant",mark_jacobian.det(),flush=True)
print("DEFORMATION_OBSTRUCTION_RANK",linear_map.rank(),"of",len(labels), flush=True)
for label,direction in zip(labels,directions):
    print("OBSTRUCTION_DIRECTION",label,list(direction),flush=True)

control, detail = obstruction(0,H**3)
assert control == vector(k,[2,0])
assert detail == (z,R(15),12*z,R(7))
mark_correction = mark_jacobian.solve_right(-control)
assert mark_correction == vector(k,[15,15])
assert obstruction(0,H**3,*mark_correction)[0] == 0
print("CONTROL_c_plus_pi_H3",list(control),"line_and_equation_variations",detail,flush=True)
print("CONTROL_VALUES_AT_SPECIAL_COMMON_POINTS",[control[0]+a*control[1]
                                               for a in [k(4),k(-4)]],flush=True)
print("UNIQUE_FIRST_ORDER_MARK_CORRECTION",list(mark_correction),flush=True)

# Coordinate/equation changes must not produce a geometric obstruction.
# New coordinates are (I+epsilon*M) times the old ones, so the new
# equations have first variation -df(Mx); the marks must move with them.
q0 = H*W-U*U-U*V-V*V
c0 = U*W*W-V**3
special_marks = [vector(k,[3,1,1,1]),vector(k,[3,-1,-1,1])]
for row in range(4):
    for column in range(4):
        mat = matrix(k,4,4)
        mat[row,column] = 1
        displacement = mat*vector(P,P.gens())
        q1 = -sum(q0.derivative(x)*dx for x,dx in zip(P.gens(),displacement))
        c1 = -sum(c0.derivative(x)*dx for x,dx in zip(P.gens(),displacement))
        motions = []
        for point in special_marks:
            change = mat*point
            motions.append(change[2]-point[2]*change[3])
        assert obstruction(q1,c1,*motions)[0] == 0, (row,column)
for ell in P.gens():
    assert obstruction(0,ell*q0)[0] == 0
assert obstruction(q0,0)[0] == obstruction(0,c0)[0] == 0
print("INFINITESIMAL_PROJECTIVE_AND_EQUATION_CHANGE_INVARIANCE=PASS",flush=True)
print("SCOPE: these are incidence equations, not the uncomputed Hurwitz deformation ideal")
print("PASS_E8_OSCULATING_DEFORMATION_OBSTRUCTIONS")
