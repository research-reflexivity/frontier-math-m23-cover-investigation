#!/usr/bin/env sage-python
"""Exact global separator on the seven certified canonical M23 curves.

Uses ONLY canonical equations and the two marked points. No reduction at
23, Galois orbit labels, or numerical approximations enter the incidence
test. The model data and seven-cover identification have their existing
independent certificates. See HURWITZ_OSCULATING.md and Section 2 of main.tex.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
from sage.all import PolynomialRing, PowerSeriesRing, matrix, prod, vector

parser = argparse.ArgumentParser()
parser.add_argument("--component", choices=["both","degree_one","sextic"], default="both")
parser.add_argument("--residual-points", action="store_true")
parser.add_argument("--change-coordinates", action="store_true")
parser.add_argument("--quartic-monodromy", action="store_true")
parser.add_argument("--verbose", action="store_true")
args = parser.parse_args()
ROOT = Path(__file__).resolve().parents[1]
records = {}
for key in ["algebra","canonical_models","marked_points"]:
    path = ROOT/"data"/("hurwitz_"+key+"_candidate.json")
    records[key] = json.loads(path.read_text())
    print("INPUT",key,hashlib.sha256(path.read_bytes()).hexdigest(),flush=True)
spec = importlib.util.spec_from_file_location("helpers",ROOT/"scripts/reconstruct_hurwitz_degree23_maps.py")
helpers = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helpers)
fields = helpers.build_fields(records["algebra"])
rat = helpers.rational
canonical = records["canonical_models"]
marked = records["marked_points"]

for component in (["degree_one","sextic"] if args.component == "both" else [args.component]):
    K = fields["base"] if component == "degree_one" else fields["absolute"]
    def dec(rec):
        if component == "degree_one":
            return K(rat(rec["rational_part"]))+rat(rec["sqrt_minus_23_part"])*K.gen()
        return sum(rat(a)*b for a,b in zip(rec["integral_basis_coordinates"],fields["integral_basis"]))
    R = PolynomialRing(K,["x0","x1","x2","x3"])
    xs = R.gens()
    equations = []
    for name in ["quadric","petri_cubic"]:
        table = canonical[name]
        equations.append(sum(dec(c[component+"_component"])*prod(xs[i] for i in m)
            for c,m in zip(table["coefficients"],table["monomials"])))
    q,c = equations
    A = [K(v) for v in marked["b_coordinates"]]
    B = [dec(v) for v in marked["c_coordinates"][component+"_component"]]
    if args.change_coordinates:
        change = matrix(K,[[1,0,0,0],[0,1,1,0],[0,0,1,1],[0,1,0,1]])
        assert change.det() != 0
        substitution = change*vector(R,xs)
        equations = [R(f(*substitution)) for f in equations]
        q,c = equations
        A = list(change.inverse()*vector(K,A))
        B = list(change.inverse()*vector(K,B))
    assert A[0] == B[0] == 1 and A != B
    assert all(f(*point) == 0 for f in equations for point in [A,B])
    qmat = matrix(K,4,4)
    for i in range(4):
        for j in range(4): qmat[i,j] = q.derivative(xs[i]).derivative(xs[j])/2
    assert qmat.det() != 0
    def jets(point):
        for parameter in [1,2,3]:
            solve = [i for i in [1,2,3] if i != parameter]
            jac = matrix(K,[[f.derivative(xs[i])(*point) for i in solve] for f in equations])
            if not jac.det(): continue
            PS = PowerSeriesRing(K,"t",default_prec=4)
            t = PS.gen()
            values = [PS(a).add_bigoh(4) for a in point]
            values[parameter] += t
            for n in range(1,4):
                err = vector(K,[f(*values)[n] for f in equations])
                correction = jac.solve_right(-err)
                for i,a in zip(solve,correction): values[i] += a*t**n
            assert all(f(*values).valuation() >= 4 for f in equations)
            return values
        raise AssertionError("marked point not smooth")
    ja,jb = jets(A),jets(B)
    forms = []
    for point_series in [ja,jb]:
        jet_matrix = matrix(K,[[f[j] for f in point_series] for j in range(3)])
        ker = jet_matrix.right_kernel().basis()
        assert len(ker) == 1
        assert matrix(K,[[f[j] for f in point_series] for j in range(4)]).rank() == 4
        form = sum(a*x for a,x in zip(ker[0],xs))
        assert form(*point_series).valuation() == 3
        forms.append(form)
    la,lb = forms
    assert la(*B) != 0 and lb(*A) != 0
    line_matrix = matrix(K,[[form.monomial_coefficient(x) for x in xs] for form in forms])
    assert line_matrix.rank() == 2
    v,w = line_matrix.right_kernel().basis()
    T = PolynomialRing(K,["u","v"])
    u,z = T.gens()
    line = [v[i]*u+w[i]*z for i in range(4)]
    qq,cc = q(*line),c(*line)
    assert qq != 0 and qq.total_degree() == 2
    gg = qq.gcd(cc)
    common_degree = gg.total_degree()
    expected = 2 if component == "degree_one" else 0
    assert common_degree == expected
    if common_degree == 2:
        assert cc % qq == 0
        discr = qq.monomial_coefficient(u*z)**2-4*qq.monomial_coefficient(u**2)*qq.monomial_coefficient(z**2)
        assert discr != 0
    else:
        assert gg.is_constant()
    print("PASS",component,"exact common divisor degree",common_degree,
          "canonical pencil degree",6-common_degree,flush=True)
    if args.verbose or component == "degree_one":
        print("OSCULATING_FORMS",forms,flush=True)
        print("LINE_QUADRIC",qq,flush=True)
        print("LINE_CUBIC_REMAINDER",cc % qq,flush=True)
    if component == "degree_one" and not args.change_coordinates:
        # Verify the compact homogeneous incidence equations printed in the paper.
        s = K.gen()
        printed_a = xs[3]
        printed_b = (xs[0] + (387*s-51)/1508*xs[1]
                     + (4851*s-9547)/19604*xs[2]
                     + (30411308*s-56608155)/53582633*xs[3])
        printed_quadric = (xs[1]**2 + (387*s+1457)/1508*xs[1]*xs[2]
                           - (217389*s+189591)/568516*xs[2]**2)
        assert la == printed_a and lb == printed_b
        printed_on_line = printed_quadric(*line)
        assert printed_on_line != 0 and qq % printed_on_line == 0
        print("PASS compact osculating-plane and secant equations in main.tex", flush=True)
    if args.residual_points and common_degree == 2:
        residual = []
        for form,opposite,point in [(la,lb,A),(lb,la,B)]:
            found = None
            for j in [1,2,3]:
                remove = opposite*(xs[j]-point[j]*xs[0])
                reduced = R.ideal(equations+[form]).saturation(R.ideal(remove))[0]
                if reduced.hilbert_polynomial() != 1: continue
                affine = R.ideal(list(reduced.gens())+[xs[0]-1])
                if affine.vector_space_dimension() != 1: continue
                candidate = [K(affine.reduce(x)) for x in xs]
                assert all(f(*candidate) == 0 for f in equations+[form])
                assert opposite(*candidate) != 0
                assert form(*jets(candidate)).valuation() == 1
                found = candidate
                break
            assert found is not None
            residual.append(found)
        DD,CC = residual
        assert len({tuple(p) for p in [A,B,CC,DD]}) == 4
        print("RESIDUAL_C",CC,flush=True)
        print("RESIDUAL_D",DD,flush=True)
        print("PASS divisor(la/lb)=3A+D-3B-C; hence [C-D]=3[A-B]",flush=True)
    if args.quartic_monodromy and common_degree == 2:
        assert not args.change_coordinates
        assert la == xs[3] and lb.monomial_coefficient(xs[0]) == 1
        P3 = PolynomialRing(K,["U","V","W"])
        U,V,W = P3.gens()
        chart = [W-lb.monomial_coefficient(xs[1])*U-lb.monomial_coefficient(xs[2])*V
                 -lb.monomial_coefficient(xs[3]),U,V,P3.one()]
        fq,fc = [f(*chart) for f in equations]
        elim = fq.resultant(fc,U)
        assert elim.degree(U) == 0 and elim.degree(V) == 4
        base = PolynomialRing(K,"w")
        target = base.gen()
        PR = PolynomialRing(base,"v")
        fibre = PR.gen()
        ff = sum(a*target**powers[2]*fibre**powers[1] for powers,a in elim.dict().items())
        content = base.zero()
        for coefficient in ff: content = content.gcd(coefficient)
        assert content != 0
        assert all(coefficient % content == 0 for coefficient in ff)
        ff = PR([coefficient // content for coefficient in ff])
        assert ff.degree() == 4
        subres = fq.subresultants(fc,U)
        linear = next(f for f in reversed(subres) if f.degree(U) == 1)
        rational_polys = PolynomialRing(base.fraction_field(),"v")
        vv = rational_polys.gen()
        def poly_from_p3(f):
            assert f.degree(U) == 0
            return sum(a*target**powers[2]*vv**powers[1] for powers,a in f.dict().items())
        first = poly_from_p3(linear.coefficient({U:1}))
        zero = poly_from_p3(linear.subs({U:0}))
        quartic = rational_polys(ff)
        assert first.gcd(quartic) == 1
        quotient = rational_polys.quotient(quartic,"vbar")
        ux = -quotient(zero)*quotient(first).inverse_of_unit()
        vx = quotient.gen()
        assert fq(ux,vx,quotient(target)) == 0
        assert fc(ux,vx,quotient(target)) == 0
        print("PASS linear subresultant recovers the original affine curve from the quartic",flush=True)
        disc = ff.discriminant()
        assert disc != 0
        sqfree = disc.squarefree_decomposition()
        odd_degree = sum(f.degree() for f,e in sqfree if e % 2)
        assert odd_degree > 0
        print("QUARTIC_DISCRIMINANT squarefree degrees and multiplicities",
              [(f.degree(),e) for f,e in sqfree],flush=True)
        print("PASS quartic geometric monodromy S4: degree 4, inertia 3, nonsquare geometric discriminant",flush=True)
print("PASS_GLOBAL_OSCULATING_SECANT_SEPARATOR",flush=True)
print("Scope: a global geometric characterization certified on the exact models, not an equation-free existence theorem.")
