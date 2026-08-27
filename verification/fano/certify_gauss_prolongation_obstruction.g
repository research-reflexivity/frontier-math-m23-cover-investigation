# Exact obstruction to selecting an unpointed order-23 inertia pair from the
# rational 2A inertia x and a chosen orientation-reversing comparison r alone.
#
# If I has order 23 in M23, its normalizer is D=23:11, so the conjugates of
# I form the 40320-point homogeneous space G/D.  On each comparison-
# fixed Nielsen class, A=C_G(x) cap C_G(r) is a 2-group
# of order 32.  Since D has odd order, A acts freely on G/D.  The central
# involution r therefore pairs the points, while A/<r> acts freely on the
# 20160 pairs.  There are 1260 residual gauge orbits, so (x,r) and the
# unpointed order-23 inertia do not canonically distinguish one pair.

SizeScreen([100000,100000]);;

generatorA :=
    (1,2)(3,4)(7,8)(9,10)(13,14)(15,16)(19,20)(21,22);;
generatorB :=
    (1,16,11,3)(2,9,21,12)(4,5,8,23)(6,22,14,18)
    (13,20)(15,17);;
G := Group(generatorA,generatorB);;
y :=
    (1,2,11,10,16,9,6,3,23,19,20,14,21,17,4,8,22,5,18,15,
     13,7,12);;
representatives := [
    (4,16)(5,10)(6,21)(8,19)(9,18)(11,23)(13,14)(17,22),
    (4,19)(5,9)(6,17)(8,16)(10,18)(11,14)(13,23)(21,22),
    (4,6)(5,9)(7,21)(8,15)(11,19)(12,18)(13,20)(14,17),
    (3,17)(5,11)(7,18)(8,16)(9,21)(12,19)(14,22)(15,23),
    (3,11)(5,17)(6,20)(7,18)(8,19)(10,13)(12,16)(14,22),
    (3,21)(4,16)(9,22)(10,15)(11,20)(12,19)(13,18)(14,17),
    (3,8)(5,12)(6,15)(9,10)(11,16)(13,21)(17,19)(20,23)
];;

if Size(G)<>10200960 or not IsSimpleGroup(G) then
    Error("the supplied generators do not give simple M23");
fi;

P := Group(y);;
D := Normalizer(G,P);;
gaussPointCount := Size(G)/Size(D);;
reflectionPairCount := gaussPointCount/2;;

if Size(P)<>23 or Size(D)<>253 or
        StructureDescription(D)<>"C23 : C11" or
        gaussPointCount<>40320 or reflectionPairCount<>20160 then
    Error("the characteristic-23 Gauss homogeneous space changed");
fi;

fixedIds := [];;
anchorOrders := [];;
pairGaugeOrbitCounts := [];;
xrPairs := [];;

ExactReflections := function(x,z)
    return Filtered(
        AsList(Centralizer(G,x)),
        reflection -> Order(reflection)=2 and
            y^reflection=z^-1 and z^reflection=y^-1
    );
end;;

Print("id\tdecomposition_order\tgauss_points\tanchor_gauge_order",
      "\treflection_pairs\tpair_gauge_orbits\n");
for id in [1..7] do
    x := representatives[id];;
    z := (x*y)^-1;;
    reflections := ExactReflections(x,z);;
    if Length(reflections)=1 then
        Add(fixedIds,id);;
        r := reflections[1];;
        Add(xrPairs,[x,r]);;
        anchorGauge := Intersection(Centralizer(G,x),Centralizer(G,r));;
        anchorOrder := Size(anchorGauge);;
        if not (r in anchorGauge) or anchorOrder<>32 or
                Set(FactorsInt(anchorOrder))<>[2] then
            Error("the (x,r)-anchor gauge is not the expected 2-group");
        fi;

        # Every point stabilizer is conjugate to odd-order D, so its
        # intersection with anchorGauge is trivial.  Thus anchorGauge acts
        # freely on G/D.  On the quotient by the central r-pairing the
        # stabilizer is exactly <r>, giving the following orbit count.
        pairGaugeOrbitCount :=
            reflectionPairCount/(anchorOrder/Size(Group(r)));;
        Add(anchorOrders,anchorOrder);;
        Add(pairGaugeOrbitCounts,pairGaugeOrbitCount);;
        Print(id,"\t",Size(D),"\t",gaussPointCount,"\t",anchorOrder,
              "\t",reflectionPairCount,"\t",pairGaugeOrbitCount,"\n");
    elif Length(reflections)<>0 then
        Error("branch reflection is not unique");
    fi;
od;

if fixedIds<>[3,6,7] or anchorOrders<>[32,32,32] or
        pairGaugeOrbitCounts<>[1260,1260,1260] then
    Error("the Gauss-pair obstruction profile changed");
fi;
if not ForAll(xrPairs,pair ->
        RepresentativeAction(G,xrPairs[1],pair,OnTuples)<>fail) then
    Error("the three (x,r) anchors are not simultaneously conjugate");
fi;

Print("# fixed_ids=3,6,7\n");
Print("# D=23:11_has_odd_order\n");
Print("# A=C_M23(x)_intersection_C_M23(r)_has_order_32\n");
Print("# A_acts_freely_on_M23/D_by_coprime_stabilizers\n");
Print("# A/<r>_acts_freely_on_the_20160_reflection_pairs\n");
Print("# 1260_pair_orbits_remain_for_every_fixed_class\n");
Print("# all_three_(x,r)_anchors_are_M23_conjugate\n");
Print("# additional_local_global_comparison_data_is_necessary\n");
Print("# PASS_GAUSS_PROLONGATION_OBSTRUCTION_CERTIFICATE\n");
QUIT_GAP(0);
