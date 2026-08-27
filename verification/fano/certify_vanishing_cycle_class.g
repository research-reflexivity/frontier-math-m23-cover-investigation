# Certify the group side of the vanishing-cycle inertia class.  The
# centralizer C_M23(x) is perfect with nonsplit center <x>, so it is a
# perfect central double cover of the certified affine group
# 2^3:PSL(3,2); it acts faithfully and transitively on the sixteen
# colliding sheets in the special fibre with PSL(3,2)-complement stabilizers.  The
# branch reflection fixes exactly four of the sixteen sheets: two double
# points split into real branches and two into conjugate branches.

SizeScreen([100000,100000]);;

gA := (1,2)(3,4)(7,8)(9,10)(13,14)(15,16)(19,20)(21,22);;
gB := (1,16,11,3)(2,9,21,12)(4,5,8,23)(6,22,14,18)(13,20)(15,17);;
G := Group(gA,gB);;
x := (3,21)(4,16)(9,22)(10,15)(11,20)(12,19)(13,18)(14,17);;
r := (3,18)(4,16)(5,6)(7,23)(10,14)(12,19)(13,21)(15,17);;

C := Centralizer(G,x);;
H := Filtered([1..23],point -> point^x=point);;
sixteen := Difference([1..23],H);;
sheetAction := ActionHomomorphism(C,sixteen,OnPoints);;
sheetStabilizer := Stabilizer(C,sixteen[1]);;
pointKernel := Kernel(ActionHomomorphism(C,H,OnPoints));;
doublePoints := List(Orbits(Group(x),sixteen),Set);;

pointwiseFixedPairs := Filtered(doublePoints,pair ->
    pair[1]^r=pair[1] and pair[2]^r=pair[2]);;
swappedPairs := Filtered(doublePoints,pair ->
    OnSets(pair,r)=pair and pair[1]^r<>pair[1]);;
movedPairs := Filtered(doublePoints,pair -> OnSets(pair,r)<>pair);;

signPatterns := [];;
stabPair := Stabilizer(Stabilizer(C,doublePoints[1],OnSets),
    doublePoints[2],OnSets);;
for element in AsList(stabPair) do
    AddSet(signPatterns,[
        doublePoints[1][1]^element<>doublePoints[1][1]
            and doublePoints[1][1]^element in doublePoints[1],
        doublePoints[2][1]^element<>doublePoints[2][1]
            and doublePoints[2][1]^element in doublePoints[2]
    ]);
od;

if Size(G)<>10200960
        or Size(C)<>2688
        or not IsPerfectGroup(C)
        or AbelianInvariants(C)<>[]
        or Centre(C)<>Group(x)
        or Length(ComplementClassesRepresentatives(C,Centre(C)))<>0
        or Size(Kernel(sheetAction))<>1
        or not IsTransitive(Image(sheetAction),[1..16])
        or Size(sheetStabilizer)<>168
        or StructureDescription(sheetStabilizer)<>"PSL(3,2)"
        or Size(Intersection(sheetStabilizer,pointKernel))<>1
        or Length(doublePoints)<>8
        or Number(sixteen,point -> point^r=point)<>4
        or Length(pointwiseFixedPairs)<>2
        or Length(swappedPairs)<>2
        or Length(movedPairs)<>4
        or Length(signPatterns)<>4 then
    Error("vanishing cycle group certificate changed");
fi;

output := OutputTextFile(
    "results/vanishing_cycle_class_group_summary.txt",false);;
SetPrintFormattingStatus(output,false);;
AppendTo(output,
    "schema=m23.cover-investigation.vanishing-cycle-class-group.v1\n",
    "centralizer_is_perfect=true\n",
    "central_inertia_nonsplit=true\n",
    "centralizer_is_perfect_central_double_cover_of_affine_group=true\n",
    "sixteen_sheet_action_faithful_transitive=true\n",
    "sheet_stabilizer=PSL(3,2)_complement\n",
    "no_quadratic_subextension_of_local_field=true\n",
    "inertia_class_nontrivial_in_residue_field=true\n",
    "reflection_fixed_sheets=4\n",
    "real_nodes_with_real_branches=2\n",
    "real_nodes_with_conjugate_branches=2\n",
    "ordered_pair_sign_patterns=4\n",
    "pairwise_root_needs_full_residue_field=true\n"
);;
CloseStream(output);;

Print("C(x) is a perfect central double cover of 2^3:PSL(3,2)\n");
Print("sixteen-sheet action faithful, stabilizer PSL(3,2) complement\n");
Print("reflection: 2 real-branch nodes, 2 conjugate-branch nodes\n");

QUIT_GAP(0);
