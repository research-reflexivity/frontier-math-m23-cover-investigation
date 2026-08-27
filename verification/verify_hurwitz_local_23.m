// Independent Magma certificate for the local Hurwitz algebra at 23.

Q<x> := PolynomialRing(Rationals());
g := x^6 - 6*x^5 + 14*x^4 - 2*x^3 - 27*x^2 + 44*x - 44;
f := x^12 - 6*x^11 + 20*x^10 - 32*x^9 + 44*x^8 - 22*x^7
     + 6*x^6 - 22*x^5 + 44*x^4 - 32*x^3 + 20*x^2 - 6*x + 1;

F23 := GF(23);
R23<y> := PolynomialRing(F23);
assert R23!g eq (y - 5)^2*(y + 1)^4;

E<a> := NumberField(g);
OE := MaximalOrder(E);
decomposition_E := Factorization(23*OE);
data_E := [ <item[2], Valuation(Integers()!Norm(item[1]), 23)>
            : item in decomposition_E ];
Sort(~data_E);
assert data_E eq [ <2,1>, <4,1> ];

L<b> := NumberField(f);
OL := MaximalOrder(L);
decomposition_L := Factorization(23*OL);
data_L := [ <item[2], Valuation(Integers()!Norm(item[1]), 23)>
            : item in decomposition_L ];
Sort(~data_L);
assert data_L eq [ <2,2>, <4,2> ];

assert not IsSquare(F23!(-1));
assert not IsSquare(F23!5);
assert not IsSquare(F23!15);
assert F23!15/(F23!5) eq (F23!7)^2;

S6 := SymmetricGroup(6);
inertia := S6!(3,4)(5,6);
frobenius := S6!(1,2)(3,5)(4,6);
D := sub<S6 | inertia, frobenius>;
assert #D eq 4 and IsAbelian(D);
orbit_sizes := [ #orbit : orbit in Orbits(D) ];
Sort(~orbit_sizes);
assert orbit_sizes eq [2,4];
assert Order(inertia) eq 2 and #Support(inertia) eq 4;
assert Order(frobenius) eq 2 and #Support(frobenius) eq 6;

print "PASS_MAGMA_HURWITZ_LOCAL_PRIME_DATA_1_PLUS_2_PLUS_4";
print "PASS_MAGMA_HURWITZ_LOCAL_DECOMPOSITION_GROUP_V4";
print "PASS_HURWITZ_LOCAL_23_MAGMA_CERTIFICATE";
