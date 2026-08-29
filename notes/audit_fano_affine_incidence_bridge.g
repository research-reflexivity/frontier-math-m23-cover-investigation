# Certify the intrinsic Fano-incidence bridge between the seven fixed sheets
# and the seven nonorigin blocks in the affine eight-point 2A fibre.

SizeScreen([100000,100000]);;

gA := (1,2)(3,4)(7,8)(9,10)(13,14)(15,16)(19,20)(21,22);;
gB := (1,16,11,3)(2,9,21,12)(4,5,8,23)(6,22,14,18)
      (13,20)(15,17);;
G := Group(gA,gB);;
x := (3,21)(4,16)(9,22)(10,15)(11,20)(12,19)(13,18)(14,17);;
C := Centralizer(G,x);;
K := FittingSubgroup(C);;
H := Filtered([1..23],point -> point^x=point);;
blocks := Set(Filtered(Orbits(Group(x),[1..23]),orbit -> Length(orbit)=2));;

phom := ActionHomomorphism(C,H,OnPoints);;
PA := Image(phom);;
fanoLinesPos := Set(List(
    Filtered(AsList(PA),element -> Order(element)=2),
    element -> Set(Filtered([1..7],point -> point^element=point))
));;
fanoLines := List(fanoLinesPos,line -> Set(List(line,i -> H[i])));;

DirectionLine := function(k)
    local image,fixedLines;
    image := Image(phom,Centralizer(C,k));
    fixedLines := Filtered([1..7],line -> ForAll(
        GeneratorsOfGroup(image),g -> OnSets(fanoLinesPos[line],g)=
                                      fanoLinesPos[line]
    ));
    if Length(fixedLines)<>1 then return fail; fi;
    return fanoLines[fixedLines[1]];
end;;

LeftDegrees := function(incidence,points)
    return List(points,point -> Number(
        incidence,pair -> pair[1]=point
    ));
end;;

RightDegrees := function(incidence,points)
    return List(points,point -> Number(
        incidence,pair -> pair[2]=point
    ));
end;;

records := [];;
for origin in blocks do
    nonorigin := Difference(blocks,[origin]);
    directions := [];
    incidence := [];
    for target in nonorigin do
        k := RepresentativeAction(K,origin,target,OnSets);
        line := DirectionLine(k);
        Add(directions,line);
        for point in line do
            Add(incidence,[point,target]);
        od;
    od;
    Add(records,rec(
        origin := origin,
        directionLines := Set(directions),
        incidence := Set(incidence),
        pointDegrees := LeftDegrees(incidence,H),
        targetDegrees := RightDegrees(incidence,nonorigin)
    ));
od;

if Size(G)<>10200960
   or Size(C)<>2688
   or Size(K)<>16
   or Length(H)<>7
   or Length(blocks)<>8
   or Size(Image(ActionHomomorphism(K,blocks,OnSets)))<>8
   or not IsTransitive(K,blocks,OnSets)
   or Length(fanoLines)<>7
   or not ForAll(fanoLines,line -> Length(line)=3)
   or not ForAll(records,record ->
       record.directionLines=Set(fanoLines)
       and Length(record.incidence)=21
       and Set(record.pointDegrees)=[3]
       and Set(record.targetDegrees)=[3]
   ) then
    Error("the Fano--affine incidence bridge audit changed");
fi;

Print("fixed_heptad_size=",Length(H),"\n");
Print("affine_block_count=",Length(blocks),"\n");
Print("translation_group_order=",Size(Image(
    ActionHomomorphism(K,blocks,OnSets))),"\n");
Print("origin_choices_certified=",Length(records),"\n");
Print("nonorigin_direction_count=7\n");
Print("Fano_incidence_count=21\n");
Print("incidence_degrees_point_to_direction=[3,3]\n");
Print("binary_cycle_identity=affine_eight_plus_incidence_pushforward_of_",
      "fixed_heptad_equals_origin\n");
Print("PASS_FANO_AFFINE_INCIDENCE_BRIDGE_AUDIT\n");

QUIT_GAP(0);
