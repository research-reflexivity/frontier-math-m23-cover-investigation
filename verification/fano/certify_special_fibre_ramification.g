# Certify the group side of the special-fibre ramification theorem.  The
# three arithmetic layers are separated by ramification behaviour:
#
#   p=23  the branch prime, where the two order-23 points collide.  In the
#         degree-1344 residue closure, inertia C7 is a Fano Singer cycle
#         and decomposition is its normalizer 7:3 (not C7 itself).  The
#         vanishing-cycle quadratic layer is ramified at 23, so in the
#         full degree-2688 field inertia is C14 and decomposition has
#         order 42 with quotient C3.  The Singer subgroup still permutes
#         the seven sheets and lines simply transitively and fixes the
#         affine origin.
#
#   p=3   an order-3 collineation, orbits 1,3,3 on sheets and 1,1,3,3 on
#         double points.
#
#   p=2   a unique conjugacy class of order-12 subgroups has orbits
#         1,3,3 on sheets and 4,4 on double points with its Sylow
#         2-subgroup inside the line module: an A4 whose Klein four-group
#         is wild and whose tame quotient is C3.

SizeScreen([100000,100000]);;

gA := (1,2)(3,4)(7,8)(9,10)(13,14)(15,16)(19,20)(21,22);;
gB := (1,16,11,3)(2,9,21,12)(4,5,8,23)(6,22,14,18)(13,20)(15,17);;
G := Group(gA,gB);;
x := (3,21)(4,16)(9,22)(10,15)(11,20)(12,19)(13,18)(14,17);;
r := (3,18)(4,16)(5,6)(7,23)(10,14)(12,19)(13,21)(15,17);;

C := Centralizer(G,x);;
H := Filtered([1..23],point -> point^x=point);;
doublePoints := List(Orbits(Group(x),Difference([1..23],H)),Set);;
pointHom := ActionHomomorphism(C,H,OnPoints);;
PA := Image(pointHom);;
K := Kernel(pointHom);;
lines := Set(List(
    Filtered(AsList(PA),element -> Order(element)=2),
    element -> Set(Filtered([1..7],point -> point^element=point))
));;
delta := Intersection(H,Filtered([1..23],point -> point^r=point));;

if Size(C)<>2688 or Size(PA)<>168 or Size(K)<>16
        or Length(lines)<>7 or Set(delta)<>[1,2,8] then
    Error("special-fibre ramification base data changed");
fi;

# ---- p = 23: the Fano Singer cycle ----

singerElements := Filtered(AsList(PA),element -> Order(element)=7);;
singer := singerElements[1];;
singerLift := PreImagesRepresentative(pointHom,singer);;
singerLift := singerLift^(Order(singerLift)/7);;
singerFixedDoublePoints := Filtered(doublePoints,
    pair -> OnSets(pair,singerLift)=pair);;
affineDecomposition23 := Normalizer(PA,Group(singer));;
fullInertia23 := Group(singerLift,x);;
fullDecomposition23 := Normalizer(C,fullInertia23);;

if Length(singerElements)<>48
        or not ForAll(singerElements,element ->
            CycleLengths(element,[1..7])=[7]
            and Length(Orbit(Group(element),lines[1],OnSets))=7)
        or Size(affineDecomposition23)<>21
        or Size(FactorGroup(affineDecomposition23,Group(singer)))<>3
        or Order(singerLift)<>7
        or SortedList(List(OrbitsDomain(Group(singerLift),doublePoints,
            OnSets),Length))<>[1,7]
        or Length(singerFixedDoublePoints)<>1
        or Size(fullInertia23)<>14
        or not IsCyclic(fullInertia23)
        or Size(fullDecomposition23)<>42
        or Size(FactorGroup(fullDecomposition23,fullInertia23))<>3 then
    Error("Singer cycle structure changed");
fi;

# Singer coordinates and the perfect difference set
cyclicIndex := [1];;
for step in [2..7] do
    Add(cyclicIndex,cyclicIndex[step-1]^singer);
od;
SingerCoordinate := point -> Position(cyclicIndex,point)-1;;
coordinateLines := Set(List(lines,line ->
    Set(List(line,point -> SingerCoordinate(point)))));;
baseSet := coordinateLines[1];;
translates := Set(List([0..6],shift ->
    Set(List(baseSet,entry -> (entry+shift) mod 7))));;
differenceMultiset := Collected(SortedList(Concatenation(List(
    baseSet,first -> List(Filtered(baseSet,second -> second<>first),
        second -> (first-second) mod 7)))));;
tripodCoordinates := Set(List(delta,sheet ->
    SingerCoordinate(Position(H,sheet))));;

if Length(coordinateLines)<>7
        or translates<>coordinateLines
        or Length(baseSet)<>3
        or differenceMultiset<>List([1..6],residue -> [residue,1])
        or not tripodCoordinates in coordinateLines then
    Error("Fano difference set structure changed");
fi;

# ---- p = 3: the order-three collineation ----

orderThree := First(AsList(PA),element -> Order(element)=3);;
orderThreeLift := PreImagesRepresentative(pointHom,orderThree);;
orderThreeLift := orderThreeLift^(Order(orderThreeLift)/3);;

if SortedList(CycleLengths(orderThree,[1..7]))<>[1,3,3]
        or SortedList(List(OrbitsDomain(Group(orderThreeLift),
            doublePoints,OnSets),Length))<>[1,1,3,3] then
    Error("order-three structure changed");
fi;

# ---- p = 2: the unique wild-in-the-line-module inertia ----

orderTwelveClasses := Filtered(ConjugacyClassesSubgroups(C),
    class -> Size(Representative(class))=12);;
inertiaCandidates := Filtered(orderTwelveClasses,class ->
    SortedList(List(OrbitsDomain(Representative(class),H,OnPoints),
        Length))=[1,3,3]
    and SortedList(List(OrbitsDomain(Representative(class),doublePoints,
        OnSets),Length))=[4,4]
    and IsSubgroup(K,SylowSubgroup(Representative(class),2)));;

if Length(orderTwelveClasses)<>6
        or Length(inertiaCandidates)<>1 then
    Error("order-twelve inertia classification changed");
fi;

wildInertiaGroup := Representative(inertiaCandidates[1]);;
wildPart := SylowSubgroup(wildInertiaGroup,2);;

if StructureDescription(wildInertiaGroup)<>"A4"
        or Size(wildPart)<>4
        or not IsElementaryAbelian(wildPart)
        or not IsSubgroup(K,wildPart)
        or Size(inertiaCandidates[1])<>28
        or not IsRegular(Image(ActionHomomorphism(K,doublePoints,OnSets)),
            [1..8]) then
    Error("wild inertia structure changed");
fi;

output := OutputTextFile("results/special_fibre_ramification_group_summary.txt",
    false);;
SetPrintFormattingStatus(output,false);;
AppendTo(output,
    "schema=m23.cover-investigation.special-fibre-ramification-group.v1\n",
    "order_seven_elements=48\n",
    "order_seven_is_singer_cycle=true\n",
    "singer_simply_transitive_on_sheets_and_lines=true\n",
    "residue_inertia_at_23=C7\n",
    "residue_decomposition_at_23=7:3\n",
    "residue_decomposition_quotient=C3\n",
    "full_inertia_at_23=C14\n",
    "full_decomposition_at_23=C2x(7:3)\n",
    "full_decomposition_order=42\n",
    "full_decomposition_quotient=C3\n",
    "singer_double_point_orbits=[1,7]\n",
    "singer_fixed_double_point_is_unique=true\n",
    "singer_fixed_double_point_is_affine_origin=true\n",
    "fano_lines_are_translates_of_a_difference_set=true\n",
    "difference_set_size=3\n",
    "difference_multiset_is_each_nonzero_residue_once=true\n",
    "distinguished_tripod_in_singer_coordinates=",tripodCoordinates,"\n",
    "distinguished_tripod_is_a_difference_set_translate=true\n",
    "order_three_sheet_orbits=[1,3,3]\n",
    "order_three_double_point_orbits=[1,1,3,3]\n",
    "order_twelve_subgroup_classes=6\n",
    "inertia_at_two_classes_matching_arithmetic=1\n",
    "inertia_at_two_structure=A4\n",
    "wild_part_order=4\n",
    "wild_part_elementary_abelian=true\n",
    "wild_part_inside_line_module=true\n",
    "tame_quotient_at_two=C3\n",
    "line_module_regular_on_double_points=true\n"
);;
CloseStream(output);;

Print("p=23: residue I=C7, D=7:3; full I=C14, |D|=42\n");
Print("Singer C7: simply transitive on 7 sheets and 7 lines, ",
    "fixes one double point\n");
Print("Fano lines = 7 translates of the difference set ",baseSet,"\n");
Print("class-6 tripod in Singer coordinates=",tripodCoordinates,"\n");
Print("unique order-12 inertia class at 2 is A4 with wild V4 in the ",
    "line module\n");

QUIT_GAP(0);
