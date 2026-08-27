# Certify the group side of the field tower obtained from the special fibre.
# The quotient C_M23(x)/<x> acts on the eight ramified double points
# as the full affine group 2^3:PSL(3,2) of order 1344, with the line
# module K/<x> as regular normal subgroup; the 28 point pairs carry a
# canonical 4-per-line labelling.  The conjugacy classes give exactly ten
# admissible (octic,septic) Frobenius factorization pairs, exported for
# the PARI certificate.

SizeScreen([100000,100000]);;

gA := (1,2)(3,4)(7,8)(9,10)(13,14)(15,16)(19,20)(21,22);;
gB := (1,16,11,3)(2,9,21,12)(4,5,8,23)(6,22,14,18)(13,20)(15,17);;
G := Group(gA,gB);;
x := (3,21)(4,16)(9,22)(10,15)(11,20)(12,19)(13,18)(14,17);;

C := Centralizer(G,x);;
H := Filtered([1..23],point -> point^x=point);;
doublePoints := List(Orbits(Group(x),Difference([1..23],H)),Set);;
pairAction := ActionHomomorphism(C,doublePoints,OnSets);;
pairImage := Image(pairAction);;
pointAction := ActionHomomorphism(C,H,OnPoints);;
K := Kernel(pointAction);;
kernelImage := Image(pairAction,K);;

fanoLinesPos := Set(List(
    Filtered(AsList(Image(pointAction)),element -> Order(element)=2),
    element -> Set(Filtered([1..7],point -> point^element=point))
));;
KernelLine := function(k)
    local image;
    image := Image(pointAction,Centralizer(C,k));
    return Filtered([1..7],l -> ForAll(GeneratorsOfGroup(image),
        g -> OnSets(fanoLinesPos[l],g)=fanoLinesPos[l]));
end;;

if Size(G)<>10200960
        or Size(C)<>2688
        or Length(doublePoints)<>8
        or Size(pairImage)<>1344
        or Kernel(pairAction)<>Group(x)
        or Transitivity(pairImage,[1..8])<2
        or StructureDescription(pairImage)<>"(C2 x C2 x C2) : PSL(3,2)"
        or Size(kernelImage)<>8
        or not IsRegular(kernelImage,[1..8]) then
    Error("affine pair action changed");
fi;

pcgsK := Pcgs(K);;
kernelModule := GModuleByMats(List(GeneratorsOfGroup(C),
    c -> List(pcgsK,k -> ExponentsOfPcElement(pcgsK,k^c))*Z(2)^0),
    GF(2));;
lineModule := MTX.InducedActionFactorModule(kernelModule,
    MTX.BasesMinimalSubmodules(kernelModule)[1]);;
if not MTX.IsIrreducible(lineModule) then
    Error("line module is not irreducible");
fi;

# canonical 4-per-line labelling of the 28 double-point pairs
lineLabelCounts := List([1..7],l -> 0);;
for i in [1..8] do
    for j in [i+1..8] do
        transporters := Filtered(Difference(Elements(K),[()]),
            k -> OnSets(doublePoints[i],k)=doublePoints[j]);
        if Length(transporters)<>2
                or transporters[1]*x<>transporters[2] then
            Error("pair transporter is not an x-shift pair");
        fi;
        labels := Set(List(transporters,KernelLine));
        if Length(labels)<>1 or Length(labels[1])<>1 then
            Error("pair line label is not well-defined");
        fi;
        lineLabelCounts[labels[1][1]] :=
            lineLabelCounts[labels[1][1]]+1;
    od;
od;
if lineLabelCounts<>[4,4,4,4,4,4,4] then
    Error("pair line labelling changed");
fi;

# the Frobenius dictionary
dictionary := [];;
for class in ConjugacyClasses(C) do
    representative := Representative(class);
    AddSet(dictionary,[
        SortedList(List(OrbitsDomain(Group(representative),
            doublePoints,OnSets),Length)),
        SortedList(List(OrbitsDomain(Group(representative),
            H,OnPoints),Length))
    ]);
od;

expectedDictionary := [
    [[1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1]],
    [[1,1,1,1,2,2],[1,1,1,2,2]],
    [[1,1,2,4],[1,2,4]],
    [[1,1,3,3],[1,3,3]],
    [[1,7],[7]],
    [[2,2,2,2],[1,1,1,1,1,1,1]],
    [[2,2,2,2],[1,1,1,2,2]],
    [[2,6],[1,3,3]],
    [[4,4],[1,1,1,2,2]],
    [[4,4],[1,2,4]]
];;
if dictionary<>Set(expectedDictionary) then
    Error("Frobenius dictionary changed");
fi;

# export for the PARI certificate
gpFile := OutputTextFile("results/special_fibre_frobenius_dictionary.gp",false);;
SetPrintFormattingStatus(gpFile,false);;
AppendTo(gpFile,"allowedpairs = ",expectedDictionary,";\n");;
CloseStream(gpFile);;

output := OutputTextFile("results/special_fibre_field_tower_group_summary.txt",
    false);;
SetPrintFormattingStatus(output,false);;
AppendTo(output,
    "schema=m23.cover-investigation.special-fibre-field-tower-group.v1\n",
    "double_point_pair_count=8\n",
    "pair_action_order=1344\n",
    "pair_action_kernel=inertia\n",
    "pair_action_two_transitive=true\n",
    "pair_action_structure=2^3:PSL(3,2)\n",
    "line_module_regular_normal=true\n",
    "line_module_irreducible=true\n",
    "pair_line_labelling=4_per_line\n",
    "frobenius_dictionary_size=10\n",
    "dictionary_exported=results/special_fibre_frobenius_dictionary.gp\n"
);;
CloseStream(output);;

Print("pair action = 2^3:PSL(3,2), order 1344, kernel = inertia\n");
Print("line module regular on double points, irreducible\n");
Print("28 pairs labelled 4 per line; dictionary has 10 classes\n");

QUIT_GAP(0);
