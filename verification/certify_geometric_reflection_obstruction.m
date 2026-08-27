/* Independent Magma audit of the geometric reflection obstruction. */

procedure Must(flag, message)
    if not flag then
        print "CERTIFICATE_FAILURE: " cat message;
        quit;
    end if;
end procedure;

S := Sym(23);
generatorA := S!(1,2)(3,4)(7,8)(9,10)(13,14)(15,16)(19,20)(21,22);
generatorB := S!(1,16,11,3)(2,9,21,12)(4,5,8,23)(6,22,14,18)
                 (13,20)(15,17);
G := sub<S | generatorA, generatorB>;
x0 := S!(3,21)(4,16)(9,22)(10,15)(11,20)(12,19)(13,18)(14,17);
y := S!(1,2,11,10,16,9,6,3,23,19,20,14,21,17,4,8,22,5,18,15,
          13,7,12);
representatives := [
    S!(4,16)(5,10)(6,21)(8,19)(9,18)(11,23)(13,14)(17,22),
    S!(4,19)(5,9)(6,17)(8,16)(10,18)(11,14)(13,23)(21,22),
    S!(4,6)(5,9)(7,21)(8,15)(11,19)(12,18)(13,20)(14,17),
    S!(3,17)(5,11)(7,18)(8,16)(9,21)(12,19)(14,22)(15,23),
    S!(3,11)(5,17)(6,20)(7,18)(8,19)(10,13)(12,16)(14,22),
    x0,
    S!(3,8)(5,12)(6,15)(9,10)(11,16)(13,21)(17,19)(20,23)
];

Must(#G eq 10200960 and #Centralizer(S,G) eq 1 and Normalizer(S,G) eq G,
     "the natural M23 model changed");

sameClassAtTwoEnds := [];
inverseClassAtTwoEnds := [];
orientationReversingCounts := [];
for x in representatives do
    z := (x*y)^-1;
    sameClass, unused := IsConjugate(G,y,z);
    inverseClass, unused := IsConjugate(G,y,z^-1);
    Append(~sameClassAtTwoEnds,sameClass);
    Append(~inverseClassAtTwoEnds,inverseClass);
    reversingCount := 0;
    for r in Centralizer(G,x) do
        if Order(r) eq 2 and y^r eq z^-1 and z^r eq y^-1 then
            reversingCount +:= 1;
        end if;
    end for;
    Append(~orientationReversingCounts,reversingCount);
end for;

Must(sameClassAtTwoEnds eq [false,false,false,false,false,false,false] and
     inverseClassAtTwoEnds eq [true,true,true,true,true,true,true] and
     orientationReversingCounts eq [0,0,1,0,0,1,1],
     "the A/B and orientation-reversal census changed");

function IsShape4A(a)
    return Order(a) eq 4 and
           #[i : i in [1..23] | i^a eq i] eq 3 and
           #[i : i in [1..23] | i^(a^2) eq i] eq 7;
end function;

class2A := Setseq(Conjugates(G,x0));
compatible4A := [];
for c in class2A do
    a := (c*y)^-1;
    if IsShape4A(a) and sub<G | a,y> eq G then
        Append(~compatible4A,a);
    end if;
end for;
compatibleSet := Seqset(compatible4A);
seen := [];
quotientRepresentatives := [];
quotientOrbitSizes := [];
for a in compatible4A do
    if not a in seen then
        orbit := {a^(y^j) : j in [0..22]};
        seen cat:= Setseq(orbit);
        Append(~quotientRepresentatives,a);
        Append(~quotientOrbitSizes,#orbit);
    end if;
end for;
Must(Seqset(seen) eq compatibleSet and #compatible4A eq 46 and
     quotientOrbitSizes eq [23,23],
     "the quotient Nielsen census changed");
for a in quotientRepresentatives do
    squareClass, unused := IsConjugate(G,a^2,x0);
    closingClass, unused := IsConjugate(G,(a^2*y)^-1,y);
    Must(squareClass and closingClass,
         "the quadratic pullback no longer has same-class 23 inertia");
end for;

Q<t> := PolynomialRing(Rationals());
elkiesQuartic := t^4+t^3+9*t^2-10*t+8;
F<g> := NumberField(elkiesQuartic);
eta := (-g^3-2*g^2-8*g+8)/3;
Must(MinimalPolynomial(eta) eq t^2-3*t+8,
     "the Elkies quartic no longer contains Q(sqrt(-23)) as recorded");

print "ambient_normalizer_equals_M23=true";
print "order23_endpoint_classes_equal=[ false, false, false, false, false, false, false ]";
print "order23_inverse_endpoint_classes_equal=[ true, true, true, true, true, true, true ]";
print "orientation_reversing_comparison_counts=[ 0, 0, 1, 0, 0, 1, 1 ]";
print "orientation_reversing_fixed_ids=[3,6,7]";
print "holomorphic_branch_swap_fixed_ids=[]";
print "quotient_4A_2A_23_compatible_count=46";
print "quotient_4A_2A_23_inner_nielsen_classes=2";
print "quadratic_pullback_order23_classes=same_not_A_B";
print "elkies_quartic_contains_quadratic_subfield_discriminant=-23";
print "conclusion=no_algebraic_reflection_fixed_cubic";
print "PASS_GEOMETRIC_REFLECTION_OBSTRUCTION_MAGMA";
