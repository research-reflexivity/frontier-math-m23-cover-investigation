# Legacy filename: certify the orientation-reversing comparisons supplied
# by complex conjugation, and distinguish them from nonexistent holomorphic
# lifts of the Q-defined base involution T -> -T.

SizeScreen([100000,100000]);;

gA := (1,2)(3,4)(7,8)(9,10)(13,14)(15,16)(19,20)(21,22);;
gB := (1,16,11,3)(2,9,21,12)(4,5,8,23)(6,22,14,18)
      (13,20)(15,17);;
G := Group(gA,gB);;
ambient := SymmetricGroup(23);;
if Size(G)<>10200960 or not IsSimpleGroup(G) then
    Error("generators do not give simple M23");
fi;

ambientCentralizer := Centralizer(ambient,G);;
ambientNormalizer := Normalizer(ambient,G);;
if Size(ambientCentralizer)<>1 or ambientNormalizer<>G then
    Error("the natural M23 action is not centerless and self-normalizing");
fi;

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

ExactReflections := function(x,z)
    return Filtered(
        AsList(Centralizer(G,x)),
        r -> y^r=z^-1 and z^r=y^-1
    );
end;;

reflectionCounts := [];;
fixedRecords := [];;
sameEndpointClasses := [];;
inverseEndpointClasses := [];;
for id in [1..7] do
    x := representatives[id];
    z := (x*y)^-1;
    reflections := ExactReflections(x,z);
    Add(sameEndpointClasses,IsConjugate(G,y,z));
    Add(inverseEndpointClasses,IsConjugate(G,y,z^-1));
    Add(reflectionCounts,Length(reflections));
    if Length(reflections)=1 then
        r := reflections[1];
        Add(fixedRecords,[
            id,
            Order(r),
            Size(Centralizer(ambient,Group(x,y,z))),
            Group(x,y,z)=G
        ]);
    fi;
od;

if reflectionCounts<>[0,0,1,0,0,1,1] or
   sameEndpointClasses<>[false,false,false,false,false,false,false] or
   inverseEndpointClasses<>[true,true,true,true,true,true,true] or
   List(fixedRecords,row -> row[1])<>[3,6,7] or
   not ForAll(fixedRecords,row -> row{[2..4]}=[2,1,true]) then
    Error("orientation-reversing comparison census changed");
fi;

PrintTo(
    "results/arithmetic_reflection_lift_summary.txt",
    "schema=m23.cover-investigation.geometric-reflection-obstruction.v2\n",
    "group_order=",Size(G),"\n",
    "ambient_degree=23\n",
    "ambient_centralizer_order=",Size(ambientCentralizer),"\n",
    "ambient_normalizer_order=",Size(ambientNormalizer),"\n",
    "ambient_normalizer_equals_M23=",ambientNormalizer=G,"\n",
    "reflection_counts=",reflectionCounts,"\n",
    "same_endpoint_classes=",sameEndpointClasses,"\n",
    "inverse_endpoint_classes=",inverseEndpointClasses,"\n",
    "fixed_record_columns=[id,lift_order,cover_automorphism_order,",
        "generates_M23]\n",
    "fixed_records=",fixedRecords,"\n",
    "fixed_ids=",List(fixedRecords,row -> row[1]),"\n",
    "orientation_reversing_comparison_is_unique_on_fixed_classes=true\n",
    "holomorphic_T_to_minus_T_lift_exists=false\n"
);;

Print("certified orientation-reversing comparisons and geometric obstruction\n");
Print("fixed_ids=",List(fixedRecords,row -> row[1]),"\n");
Print("PASS_GEOMETRIC_REFLECTION_OBSTRUCTION_LEGACY_CERTIFICATE\n");
QUIT_GAP(0);
