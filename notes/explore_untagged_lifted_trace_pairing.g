# Compare the tagged lifted equality incidence with the pairing obtained
# after the involution tag is forgotten.  This diagnoses whether the raw
# Hadamard node calculation can belong to the same generic carrier.

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
S := SymmetricGroup(23);;
N := Normalizer(S,D);;
n := Difference(AsList(N),AsList(D))[1];;

class2A := AsList(ConjugacyClass(G,representatives[6]));;
compatible := Filtered(class2A,x ->
    Order(x*y)=23 and IsConjugate(G,x*y,y)
);;

EndpointRecords := function(E,x,b)
    local records,d;
    records := [];
    for d in Elements(E) do
        Add(records,[x^d,OnSets(b,d)]);
    od;
    return records;
end;;

Analyze := function(x)
    local z,Dz,H,A,B,taggedTotal,untaggedTotal,b,recordsY,recordsZ,
          taggedSetY,taggedSetZ,targetsY,targetsZ,target,multiplicityY,
          multiplicityZ;
    z := (x*y)^-1;
    Dz := Normalizer(G,Group(z));
    H := Filtered([1..23],point -> point^x=point);
    A := List(Orbits(Group(x),Difference([1..23],H),OnPoints),Set);
    B := Concatenation(List(H,point -> [point]),A);
    taggedTotal := 0;
    untaggedTotal := 0;
    for b in B do
        recordsY := EndpointRecords(D,x,b);
        recordsZ := EndpointRecords(Dz,x,b);

        taggedSetY := Set(recordsY);
        taggedSetZ := Set(recordsZ);
        taggedTotal := taggedTotal+Length(Intersection(taggedSetY,taggedSetZ));

        targetsY := Set(List(recordsY,record -> record[2]));
        targetsZ := Set(List(recordsZ,record -> record[2]));
        for target in Intersection(targetsY,targetsZ) do
            multiplicityY := Number(recordsY,record -> record[2]=target);
            multiplicityZ := Number(recordsZ,record -> record[2]=target);
            untaggedTotal := untaggedTotal+multiplicityY*multiplicityZ;
        od;
    od;
    return [taggedTotal,untaggedTotal,taggedTotal mod 2,untaggedTotal mod 2];
end;;

records := List(representatives,Analyze);;
allRecords := List(compatible,Analyze);;

# The untagged value has a direct universal decomposition.  On a singleton
# source orbit, every 23:11 endpoint group sends the point to each of the
# 23 targets with multiplicity 11.  On a two-point source orbit it acts
# regularly on all C(23,2)=253 unordered pairs.
singletonContribution := 23*11^2;
pairContribution := 253;
universalUntaggedTotal := 7*singletonContribution+8*pairContribution;

# Diagnose the normalized and finite-node pieces for square and nonsquare
# returns.  The nonsquare element transports both the involution tag and
# the target orbit, while normalizing D.
x0 := representatives[6];;
B0 := Concatenation(
    List(Filtered([1..23],point -> point^x0=point),point -> [point]),
    List(Orbits(Group(x0),Filtered([1..23],point -> point^x0<>point),OnPoints),Set)
);;
baseRecords := function(b)
    return List(Elements(D),d -> [x0^d,OnSets(b,d)]);
end;;
returnedRecords := function(b)
    return List(Elements(D),d -> [x0^(d*n),OnSets(b,d*n)]);
end;;
PairTotals := function(recordFunction)
    local tagged,untagged,b,first,second,targetsFirst,targetsSecond,target;
    tagged := 0;
    untagged := 0;
    for b in B0 do
        first := baseRecords(b);
        second := recordFunction(b);
        tagged := tagged+Length(Intersection(Set(first),Set(second)));
        targetsFirst := Set(List(first,record -> record[2]));
        targetsSecond := Set(List(second,record -> record[2]));
        for target in Intersection(targetsFirst,targetsSecond) do
            untagged := untagged
                +Number(first,record -> record[2]=target)
                 *Number(second,record -> record[2]=target);
        od;
    od;
    return [tagged mod 2,untagged mod 2,(tagged+untagged) mod 2];
end;;
squareNormalized := PairTotals(baseRecords);;
nonsquareNormalized := PairTotals(returnedRecords);;

# At the finite node, forgetting the tag gives branch coefficients
# 77^2 and 8^2.  With the tag retained, a square return contributes 77
# and 8; a nonsquare return has no common involution tag.
squareNode := [
    (77-8) mod 2,
    (77^2-8^2) mod 2,
    ((77^2-8^2)-(77-8)) mod 2
];;
nonsquareNode := [0,(77^2-8^2) mod 2,(77^2-8^2) mod 2];;

# The actual finite conductor forgets both the source orbit b and the tag w,
# retaining only b^d.  Across all fifteen source orbits, its singleton
# buckets have size 77 and its two-set buckets have size 8.
ConductorBucketCounts := function(recordFunction)
    local targets,b,records,targetSet;
    targets := [];
    for b in B0 do
        records := recordFunction(b);
        Append(targets,List(records,record -> record[2]));
    od;
    targetSet := Set(targets);
    return List(targetSet,target -> [target,Number(targets,x -> x=target)]);
end;;
baseBuckets := ConductorBucketCounts(baseRecords);;
returnedBuckets := ConductorBucketCounts(returnedRecords);;
baseBucketSizes := Collected(List(baseBuckets,pair -> pair[2]));;
returnedBucketSizes := Collected(List(returnedBuckets,pair -> pair[2]));;
conductorPairingDegree := Sum(baseBuckets,pair ->
    pair[2]*returnedBuckets[Position(List(returnedBuckets,x -> x[1]),pair[1])][2]
);;

if Set(List(allRecords,r -> r[2]))<>[universalUntaggedTotal]
   or universalUntaggedTotal<>21505
   or squareNormalized<>[1,1,0]
   or nonsquareNormalized<>[0,1,1]
   or squareNode<>[1,1,0]
   or nonsquareNode<>[0,1,1]
   or baseBucketSizes<>[ [ 8, 253 ], [ 77, 23 ] ]
   or returnedBucketSizes<>baseBucketSizes
   or conductorPairingDegree<>23*77^2+253*8^2 then
    Error("untagged/tagged decomposition changed");
fi;

Print("fields=[tagged_integer,untagged_integer,tagged_parity,untagged_parity]\n");
for id in [1..7] do
    Print("id_",id,"=",records[id],"\n");
od;
Print("all_161_untagged_parities=",Set(List(allRecords,r -> r[4])),"\n");
Print("tagged_parities=",List(records,r -> r[3]),"\n");
Print("untagged_parities=",List(records,r -> r[4]),"\n");
Print("universal_untagged_decomposition=7*(23*11^2)+8*253=",
      universalUntaggedTotal,"\n");
Print("returned_fields=[tagged,untagged,off_diagonal]_mod_2\n");
Print("square_normalized=",squareNormalized,"\n");
Print("nonsquare_normalized=",nonsquareNormalized,"\n");
Print("square_node=",squareNode,"\n");
Print("nonsquare_node=",nonsquareNode,"\n");
Print("conductor_bucket_size_distribution=",baseBucketSizes,"\n");
Print("conductor_pairing_degree=23*77^2+253*8^2=",
      conductorPairingDegree,"\n");
Print("PASS_UNTAGGED_LIFTED_TRACE_PAIRING_DIAGNOSTIC\n");

QUIT_GAP(0);
