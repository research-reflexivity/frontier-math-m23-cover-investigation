# Certify that the Fano--affine relative-return fixed cycle is a literal
# intersection of two doubly marked endpoint traces.  The audit covers the
# full framed Nielsen locus (161 tuples = seven y-translation orbits).

SizeScreen([100000,100000]);;

gA := (1,2)(3,4)(7,8)(9,10)(13,14)(15,16)(19,20)(21,22);;
gB := (1,16,11,3)(2,9,21,12)(4,5,8,23)(6,22,14,18)
      (13,20)(15,17);;
G := Group(gA,gB);;
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
D := Normalizer(G,Group(y));;
class2A := AsList(ConjugacyClass(G,representatives[6]));;
compatible := Filtered(class2A,x ->
    Order(x*y)=23 and IsConjugate(G,x*y,y)
);;
translationOrbits := OrbitsDomain(Group(y),Set(compatible),OnPoints);;

# The returned pair (w,d) records the unique d in E with x^d=w.
EndpointReturns := function(E,x)
    local trace;
    trace := Set(Orbit(E,x,OnPoints));
    return List(trace,w -> [w,RepresentativeAction(E,x,w,OnPoints)]);
end;;

# For an x-orbit b on the sheets, lift the trace x^E by also remembering
# the transported orbit b^d.  These are actual pairs of sheet-orbit labels,
# not centralizer elements.
LiftedTrace := function(returns,b)
    return Set(List(returns,pair -> [pair[1],OnSets(b,pair[2])]));
end;;

Analyze := function(x)
    local z,Dz,returnsY,returnsZ,traceY,traceZ,overlap,returnY,returnZ,
          relativeReturns,H,A,B,liftedIntersectionCounts,liftedVector,
          fixedVector,bIndex,position,c,orbitRecordsY,orbitRecordsZ;
    z := (x*y)^-1;
    Dz := Normalizer(G,Group(z));
    returnsY := EndpointReturns(D,x);
    returnsZ := EndpointReturns(Dz,x);
    traceY := List(returnsY,pair -> pair[1]);
    traceZ := List(returnsZ,pair -> pair[1]);
    overlap := Intersection(traceY,traceZ);
    returnY := List(overlap,w -> returnsY[Position(traceY,w)][2]);
    returnZ := List(overlap,w -> returnsZ[Position(traceZ,w)][2]);
    relativeReturns := List([1..Length(overlap)],position ->
        returnY[position]*returnZ[position]^-1
    );

    H := Filtered([1..23],point -> point^x=point);
    A := List(Orbits(Group(x),Difference([1..23],H),OnPoints),Set);
    B := Concatenation(List(H,point -> [point]),A);

    liftedIntersectionCounts := [];
    liftedVector := [];
    fixedVector := [];
    for bIndex in [1..Length(B)] do
        orbitRecordsY := LiftedTrace(returnsY,B[bIndex]);
        orbitRecordsZ := LiftedTrace(returnsZ,B[bIndex]);
        Add(liftedIntersectionCounts,
            Length(Intersection(orbitRecordsY,orbitRecordsZ)));
        Add(liftedVector,liftedIntersectionCounts[bIndex] mod 2);
        Add(fixedVector,
            Number(relativeReturns,c -> OnSets(B[bIndex],c)=B[bIndex]) mod 2);
    od;

    if liftedVector<>fixedVector then
        Error("lifted-trace intersection differs from fixed-locus cycle");
    fi;
    if Sum(liftedVector) mod 2 <> Length(overlap) mod 2 then
        Error("odd fixed-point augmentation failed");
    fi;
    if ForAny(relativeReturns,c -> x^c<>x) then
        Error("a relative return does not centralize x");
    fi;

    return rec(
        overlapSize := Length(overlap),
        totalLiftedIntersection := Sum(liftedIntersectionCounts),
        fanoIntersection := Sum(liftedIntersectionCounts{[1..7]}),
        affineIntersection := Sum(liftedIntersectionCounts{[8..15]}),
        fanoHammingWeight := Sum(liftedVector{[1..7]}),
        affineHammingWeight := Sum(liftedVector{[8..15]}),
        vector := liftedVector
    );
end;;

allRecords := List(compatible,Analyze);;
orbitRecords := [];;
profiles := [];;
SameScalarProfile := function(records)
    local first;
    first := records[1];
    return ForAll(records,record ->
        record.overlapSize=first.overlapSize
        and record.totalLiftedIntersection=first.totalLiftedIntersection
        and record.fanoIntersection=first.fanoIntersection
        and record.affineIntersection=first.affineIntersection
        and record.fanoHammingWeight=first.fanoHammingWeight
        and record.affineHammingWeight=first.affineHammingWeight
    );
end;;
for id in [1..7] do
    orbitIndex := PositionProperty(
        translationOrbits,orbit -> representatives[id] in orbit
    );
    if orbitIndex=fail then
        Error("a Nielsen representative is missing from the framed locus");
    fi;
    positions := List(
        translationOrbits[orbitIndex],x -> Position(compatible,x)
    );
    records := List(positions,position -> allRecords[position]);
    first := records[1];
    if not SameScalarProfile(records) then
        Error("a scalar lifted-trace profile varies on a translation orbit");
    fi;
    Add(orbitRecords,records);
    Add(profiles,[
        first.overlapSize,
        first.totalLiftedIntersection,
        first.fanoIntersection,
        first.affineIntersection,
        first.fanoHammingWeight,
        first.affineHammingWeight
    ]);
od;

degreeParities := List(profiles,profile ->
    (profile[5]+profile[6]) mod 2
);;
overlapParities := List(profiles,profile -> profile[1] mod 2);;
if Size(G)<>10200960 or Size(D)<>253
        or Length(compatible)<>161
        or List(translationOrbits,Length)<>[23,23,23,23,23,23,23]
        or degreeParities<>overlapParities
        or overlapParities<>[1,1,1,1,1,0,1] then
    Error("lifted-trace Fano--affine incidence audit changed");
fi;

Print("framed_nielsen_tuple_count=",Length(compatible),"\n");
Print("translation_orbits=7x23\n");
Print("profile_fields=[overlap_size,total_lifted_intersection,",
      "fano_intersection,affine_intersection,fano_hamming_weight,",
      "affine_hamming_weight]\n");
for id in [1..7] do
    Print("id_",id,"=",profiles[id],"\n");
od;
Print("degree_parities=",degreeParities,"\n");
Print("overlap_parities=",overlapParities,"\n");
Print("PASS_LIFTED_TRACE_FANO_AFFINE_INCIDENCE\n");

QUIT_GAP(0);
