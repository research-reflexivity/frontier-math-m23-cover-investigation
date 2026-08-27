// Independent Magma certificate for the sextic Hurwitz Galois closure.

Q<x> := PolynomialRing(Rationals());
g := x^6 - 6*x^5 + 14*x^4 - 2*x^3 - 27*x^2 + 44*x - 44;
f := x^12 - 6*x^11 + 20*x^10 - 32*x^9 + 44*x^8 - 22*x^7
     + 6*x^6 - 22*x^5 + 44*x^4 - 32*x^3 + 20*x^2 - 6*x + 1;

assert IsIrreducible(g);
assert Discriminant(g) eq 2^22*11*23^4;

E<y> := NumberField(g);
OE := MaximalOrder(E);
assert Discriminant(OE) eq 2^4*11*23^4;

A<a> := NumberField(f);
OA := MaximalOrder(A);
assert Discriminant(OA) eq 2^8*11^2*23^8;
assert MinimalPolynomial(a + 1/a) eq g;

s := (-289*a^11 + 1647*a^10 - 5291*a^9 + 7715*a^8
      - 10615*a^7 + 3621*a^6 - 1155*a^5 + 6609*a^4
      - 10613*a^3 + 6621*a^2 - 4713*a + 887)/128;
assert s^2 eq -23;
assert Degree(MinimalPolynomial((a + 1/a) + s)) eq 12;

function FactorDegreesMod(poly, p)
    F := GF(p);
    R<z> := PolynomialRing(F);
    fac := Factorization(R!poly);
    degrees := &cat[[Degree(item[1]) : i in [1..item[2]]] : item in fac];
    Sort(~degrees);
    return degrees;
end function;

assert FactorDegreesMod(g, 3) eq [6];
assert FactorDegreesMod(g, 139) eq [1,5];
assert FactorDegreesMod(g, 2671) eq [1,1,1,1,2];

G, roots, galois_data := GaloisGroup(g);
assert #G eq 720;
number, degree := TransitiveGroupIdentification(G);
assert degree eq 6;
assert number eq 16;

print "PASS_MAGMA_HURWITZ_TRACE_SEXTIC_GALOIS_GROUP_S6";
print "PASS_MAGMA_HURWITZ_ABSOLUTE_FIELD_IS_TRACE_FIELD_TIMES_Q_SQRT_MINUS_23";
print "PASS_HURWITZ_GALOIS_CLOSURE_MAGMA_CERTIFICATE";
