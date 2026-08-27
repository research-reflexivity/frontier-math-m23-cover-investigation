# Distinguish the holomorphic base involution T -> -T from the
# orientation-reversing comparison induced by complex conjugation.

SizeScreen([100000,100000]);;

gA := (1,2)(3,4)(7,8)(9,10)(13,14)(15,16)(19,20)(21,22);;
gB := (1,16,11,3)(2,9,21,12)(4,5,8,23)(6,22,14,18)
      (13,20)(15,17);;
G := Group(gA,gB);;
ambient := SymmetricGroup(23);;
x0 := (3,21)(4,16)(9,22)(10,15)(11,20)(12,19)(13,18)(14,17);;
y := (1,2,11,10,16,9,6,3,23,19,20,14,21,17,4,8,22,5,18,
      15,13,7,12);;
representatives := [
    (4,16)(5,10)(6,21)(8,19)(9,18)(11,23)(13,14)(17,22),
    (4,19)(5,9)(6,17)(8,16)(10,18)(11,14)(13,23)(21,22),
    (4,6)(5,9)(7,21)(8,15)(11,19)(12,18)(13,20)(14,17),
    (3,17)(5,11)(7,18)(8,16)(9,21)(12,19)(14,22)(15,23),
    (3,11)(5,17)(6,20)(7,18)(8,19)(10,13)(12,16)(14,22),
    (3,21)(4,16)(9,22)(10,15)(11,20)(12,19)(13,18)(14,17),
    (3,8)(5,12)(6,15)(9,10)(11,16)(13,21)(17,19)(20,23)
];;

if Size(G)<>10200960 or Size(Centralizer(ambient,G))<>1
   or Normalizer(ambient,G)<>G then
    Error("the natural M23 model changed");
fi;

sameClassAtTwoEnds := [];;
inverseClassAtTwoEnds := [];;
orientationReversingCounts := [];;
ExactReversingComparisons := function(group,x,y0,z0)
    return Filtered(AsList(Centralizer(group,x)),r ->
        Order(r)=2 and y0^r=z0^-1 and z0^r=y0^-1);
end;;
for id in [1..7] do
    x := representatives[id];
    z := (x*y)^-1;
    Add(sameClassAtTwoEnds,IsConjugate(G,y,z));
    Add(inverseClassAtTwoEnds,IsConjugate(G,y,z^-1));
    reversing := ExactReversingComparisons(G,x,y,z);
    Add(orientationReversingCounts,Length(reversing));
od;

if sameClassAtTwoEnds<>[false,false,false,false,false,false,false]
   or inverseClassAtTwoEnds<>[true,true,true,true,true,true,true]
   or orientationReversingCounts<>[0,0,1,0,0,1,1] then
    Error("the A/B and orientation-reversal census changed");
fi;

# Audit the tempting genus-zero quotient passport.  There are two inner
# Nielsen classes with branch shapes (4A,2A,23A).  Ramified quadratic
# pullback squares the 4A element, but the two lifted order-23 branches
# remain in the same M23 class, so the result has type (2A,23A,23A), not
# the required (2A,23A,23B).
class2A := AsList(ConjugacyClass(G,x0));;
IsShape4A := function(a)
    return Order(a)=4
       and Number([1..23],i->i^a=i)=3
       and Number([1..23],i->i^(a^2)=i)=7;
end;;
compatible4A := [];;
for c in class2A do
    a := (c*y)^-1;
    if IsShape4A(a) and Group(a,y)=G then
        Add(compatible4A,a);
    fi;
od;
quotientOrbits := OrbitsDomain(Group(y),Set(compatible4A),OnPoints);;
quotientRepresentatives := List(quotientOrbits,o->o[1]);;
pullbackClosingCycles := List(quotientRepresentatives,a->(a^2*y)^-1);;

if Length(compatible4A)<>46 or Length(quotientOrbits)<>2
   or List(quotientOrbits,Length)<>[23,23]
   or not ForAll(quotientRepresentatives,a->IsConjugate(G,a^2,x0))
   or not ForAll(pullbackClosingCycles,z->IsConjugate(G,z,y)) then
    Error("the 4A quadratic-pullback obstruction changed");
fi;

Print("ambient_normalizer_equals_M23=true\n");
Print("order23_endpoint_classes_equal=",sameClassAtTwoEnds,"\n");
Print("order23_inverse_endpoint_classes_equal=",inverseClassAtTwoEnds,"\n");
Print("orientation_reversing_comparison_counts=",
      orientationReversingCounts,"\n");
Print("orientation_reversing_fixed_ids=[3,6,7]\n");
Print("holomorphic_branch_swap_fixed_ids=[]\n");
Print("quotient_4A_2A_23_compatible_count=",Length(compatible4A),"\n");
Print("quotient_4A_2A_23_inner_nielsen_classes=",Length(quotientOrbits),"\n");
Print("quadratic_pullback_order23_classes=same_not_A_B\n");
Print("conclusion=no_algebraic_reflection_fixed_cubic\n");
Print("PASS_GEOMETRIC_REFLECTION_OBSTRUCTION\n");
QUIT_GAP(0);
