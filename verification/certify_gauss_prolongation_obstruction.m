// Exact Magma obstruction to selecting a Gauss-point pair from the
// rational 2A inertia x and a chosen orientation-reversing comparison r alone.
//
// At p=23 the completed Galois closure has decomposition group D=23:11,
// so its prolongations form the 40320-point homogeneous space G/D.  For
// every comparison-fixed Nielsen class this certificate checks that
// A=C_G(x) meet C_G(r) is a 2-group of order 32.  Coprimality with the
// odd-order point stabilizer D makes the A-action on G/D free.  Quotienting
// by the central involution r leaves 1260 A/<r>-orbits of unordered pairs,
// not a canonical pair.

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
y := S!(1,2,11,10,16,9,6,3,23,19,20,14,21,17,4,8,22,5,18,15,
          13,7,12);
representatives := [
    S!(4,16)(5,10)(6,21)(8,19)(9,18)(11,23)(13,14)(17,22),
    S!(4,19)(5,9)(6,17)(8,16)(10,18)(11,14)(13,23)(21,22),
    S!(4,6)(5,9)(7,21)(8,15)(11,19)(12,18)(13,20)(14,17),
    S!(3,17)(5,11)(7,18)(8,16)(9,21)(12,19)(14,22)(15,23),
    S!(3,11)(5,17)(6,20)(7,18)(8,19)(10,13)(12,16)(14,22),
    S!(3,21)(4,16)(9,22)(10,15)(11,20)(12,19)(13,18)(14,17),
    S!(3,8)(5,12)(6,15)(9,10)(11,16)(13,21)(17,19)(20,23)
];

Must(#G eq 10200960 and IsSimple(G),
     "the supplied generators do not give simple M23");
P := sub<G | y>;
D := Normalizer(G, P);
gaussPointCount := #G div #D;
reflectionPairCount := gaussPointCount div 2;
Must(#P eq 23 and #D eq 253 and gaussPointCount eq 40320 and
     reflectionPairCount eq 20160,
     "the characteristic-23 Gauss homogeneous space changed");

fixedIds := [];
anchorOrders := [];
pairGaugeOrbitCounts := [];
fixedXs := [ G | ];
fixedRs := [ G | ];

print "id\tdecomposition_order\tgauss_points\tanchor_gauge_order" cat
      "\treflection_pairs\tpair_gauge_orbits";
for id in [1..7] do
    x := representatives[id];
    z := (x*y)^-1;
    reflections := [ reflection : reflection in Centralizer(G, x) |
        Order(reflection) eq 2 and y^reflection eq z^-1 and
        z^reflection eq y^-1 ];
    if #reflections eq 1 then
        Append(~fixedIds, id);
        r := reflections[1];
        Append(~fixedXs, x);
        Append(~fixedRs, r);
        anchorGauge := Centralizer(G, x) meet Centralizer(G, r);
        anchorOrder := #anchorGauge;
        Must(r in anchorGauge and anchorOrder eq 32,
             "the (x,r)-anchor gauge is not the expected 2-group");

        // A point stabilizer in G/D is conjugate to odd-order D.  Its
        // intersection with the 2-group anchorGauge is therefore trivial,
        // so the action is free.  A pair stabilizer is exactly <r>.
        pairGaugeOrbitCount := reflectionPairCount div (anchorOrder div 2);
        Append(~anchorOrders, anchorOrder);
        Append(~pairGaugeOrbitCounts, pairGaugeOrbitCount);
        printf "%o\t%o\t%o\t%o\t%o\t%o\n", id, #D, gaussPointCount,
               anchorOrder, reflectionPairCount, pairGaugeOrbitCount;
    else
        Must(#reflections eq 0, "branch reflection is not unique");
    end if;
end for;

Must(fixedIds eq [3,6,7] and anchorOrders eq [32,32,32] and
     pairGaugeOrbitCounts eq [1260,1260,1260],
     "the Gauss-pair obstruction profile changed");

// The three ordered anchors (x,r) lie in one simultaneous-conjugacy orbit.
// Map x first, then adjust inside its centralizer to map r as well.
for index in [2..3] do
    xConjugate, transporter := IsConjugate(G, fixedXs[1], fixedXs[index]);
    Must(xConjugate, "the 2A anchors are not conjugate");
    rConjugate, adjustment := IsConjugate(
        Centralizer(G, fixedXs[index]), fixedRs[1]^transporter, fixedRs[index]
    );
    Must(rConjugate,
         "the three (x,r) anchors are not simultaneously conjugate");
end for;

print "fixed_ids=3,6,7";
print "D=23:11_has_odd_order";
print "A=C_M23(x)_intersection_C_M23(r)_has_order_32";
print "A_acts_freely_on_M23/D_by_coprime_stabilizers";
print "A/<r>_acts_freely_on_the_20160_reflection_pairs";
print "1260_pair_orbits_remain_for_every_fixed_class";
print "all_three_(x,r)_anchors_are_M23_conjugate";
print "additional_local_global_comparison_data_is_necessary";
print "PASS_GAUSS_PROLONGATION_OBSTRUCTION_MAGMA_CERTIFICATE";
