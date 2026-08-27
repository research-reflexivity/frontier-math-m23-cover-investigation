# Identify the three-sheet Paley tripod as one Fano line inside the seven
# unramified sheets above the rational 2A branch point.  The full point and line
# schemes are the two degree-7 Gassmann actions of PSL(3,2).

SizeScreen([100000,100000]);;

gA := (1,2)(3,4)(7,8)(9,10)(13,14)(15,16)(19,20)(21,22);;
gB := (1,16,11,3)(2,9,21,12)(4,5,8,23)(6,22,14,18)
      (13,20)(15,17);;
G := Group(gA,gB);;
y := (1,2,11,10,16,9,6,3,23,19,20,14,21,17,4,8,22,5,18,15,
      13,7,12);;
x := (3,21)(4,16)(9,22)(10,15)(11,20)(12,19)(13,18)(14,17);;
r := (3,18)(4,16)(5,6)(7,23)(10,14)(12,19)(13,21)(15,17);;

heptad := Filtered([1..23],point -> point^x=point);;
delta := Intersection(
    heptad,
    Filtered([1..23],point -> point^r=point)
);;
centralizer := Centralizer(G,x);;
pointActionHom := ActionHomomorphism(centralizer,heptad,OnPoints);;
pointAction := Image(pointActionHom);;
pointKernel := Kernel(pointActionHom);;

actionInvolutions := Filtered(
    AsList(pointAction),element -> Order(element)=2
);;
fanoLines := Set(List(actionInvolutions,involution -> Set(Filtered(
    [1..7],point -> point^involution=point
))));;
sheetFanoLines := Set(List(fanoLines,line -> Set(List(
    line,index -> heptad[index]
))));;
deltaPosition := Position(sheetFanoLines,Set(delta));;
fanoPointDegrees := List(
    [1..7],point -> Number(fanoLines,line -> point in line)
);;
fanoPairs := Set(Concatenation(List(
    fanoLines,line -> Combinations(line,2)
)));;

lineActionHom := ActionHomomorphism(pointAction,fanoLines,OnSets);;
lineAction := Image(lineActionHom);;
pointStabilizer := Stabilizer(pointAction,1);;
lineStabilizer := Stabilizer(pointAction,fanoLines[1],OnSets);;

classes := ConjugacyClasses(pointAction);;
classRepresentatives := List(classes,Representative);;
pointFixedCounts := List(classRepresentatives,element -> Number(
    [1..7],point -> point^element=point
));;
lineFixedCounts := List(classRepresentatives,element -> Number(
    fanoLines,line -> OnSets(line,element)=line
));;
classOrders := List(classRepresentatives,Order);;
classSizes := List(classes,Size);;

tripodOrbit := Orbit(centralizer,Set(delta),OnSets);;

if Size(G)<>10200960
        or Length(heptad)<>7
        or Set(delta)<>[1,2,8]
        or Size(centralizer)<>2688
        or Size(pointKernel)<>16
        or not IsAbelian(pointKernel)
        or Exponent(pointKernel)<>2
        or Size(pointAction)<>168
        or StructureDescription(pointAction)<>"PSL(3,2)"
        or Length(actionInvolutions)<>21
        or Length(fanoLines)<>7
        or not ForAll(fanoLines,line -> Length(line)=3)
        or fanoPointDegrees<>[3,3,3,3,3,3,3]
        or Length(fanoPairs)<>21
        or deltaPosition=fail
        or Size(lineAction)<>168
        or Size(pointStabilizer)<>24
        or Size(lineStabilizer)<>24
        or IsConjugate(pointAction,pointStabilizer,lineStabilizer)
        or pointFixedCounts<>lineFixedCounts
        or Length(tripodOrbit)<>7
        or Set(tripodOrbit)<>sheetFanoLines then
    Error("special-fibre Fano geometry certificate changed");
fi;

output := OutputTextFile(
    "results/special_fibre_fano_geometry_summary.txt",false
);;
SetPrintFormattingStatus(output,false);;
AppendTo(output,
    "schema=m23.cover-investigation.special-fibre-fano-geometry.v1\n",
    "twoA_fixed_heptad=",heptad,"\n",
    "quartic_tripod_sheet_set=",Set(delta),"\n",
    "centralizer_order=",Size(centralizer),"\n",
    "centralizer_point_action_kernel_order=",Size(pointKernel),"\n",
    "centralizer_point_action=PSL(3,2)\n",
    "centralizer_point_action_order=",Size(pointAction),"\n",
    "action_involution_count=",Length(actionInvolutions),"\n",
    "fano_line_count=",Length(fanoLines),"\n",
    "fano_point_incidence_degrees=",fanoPointDegrees,"\n",
    "fano_pair_count=",Length(fanoPairs),"\n",
    "fano_plane_2_design=true\n",
    "fano_lines_on_special_fibre_sheets=",sheetFanoLines,"\n",
    "quartic_tripod_is_fano_line=true\n",
    "quartic_tripod_fano_line_position=",deltaPosition,"\n",
    "quartic_tripod_centralizer_orbit_size=",Length(tripodOrbit),"\n",
    "quartic_tripod_orbit_is_full_fano_line_scheme=true\n",
    "point_stabilizer_order=",Size(pointStabilizer),"\n",
    "line_stabilizer_order=",Size(lineStabilizer),"\n",
    "point_and_line_stabilizers_conjugate_in_PSL3_2=false\n",
    "conjugacy_class_orders=",classOrders,"\n",
    "conjugacy_class_sizes=",classSizes,"\n",
    "point_action_fixed_counts=",pointFixedCounts,"\n",
    "line_action_fixed_counts=",lineFixedCounts,"\n",
    "point_and_line_permutation_characters_equal=true\n",
    "degree7_point_and_line_fields_form_a_Gassmann_pair=true\n"
);;
CloseStream(output);;

Print("2A special-fibre heptad=",heptad,"\n");
Print("quartic tripod=",Set(delta)," is Fano line ",deltaPosition,"\n");
Print("PSL(3,2) point/line characters=",pointFixedCounts,"\n");
Print("point and line stabilizers are nonconjugate Gassmann twins\n");

QUIT_GAP(0);
