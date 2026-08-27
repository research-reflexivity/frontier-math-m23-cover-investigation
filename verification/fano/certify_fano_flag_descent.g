# Certify the real-flag structure of the Fano plane on the special fibre and its
# reconstruction from the global y-frame Hurwitz classes.
#
# Four exact layers:
#   1. the full point--line--incidence envelope is internal to C_M23(x):
#      lines are the noncentral involutions of K=O_2(C) modulo the inertia
#      socle <x>, incidence flags are the 84 non-kernel involutions;
#   2. the branch reflection r marks the flag (8,{1,2,8}) as the mutual
#      transvection flag of the commuting 2A pair (x,r);
#   3. the point--line duality has no local realization: Out(C)=1,
#      N_S23(C)=C, and every twin gauge exits the Steiner system;
#   4. the seven y-frame Nielsen classes reconstruct the same flag: the
#      axis is the unique line avoided by all 161 compatible involutions,
#      and Nielsen class 6 is the unique class avoiding the center.

SizeScreen([100000,100000]);;

gA := (1,2)(3,4)(7,8)(9,10)(13,14)(15,16)(19,20)(21,22);;
gB := (1,16,11,3)(2,9,21,12)(4,5,8,23)(6,22,14,18)(13,20)(15,17);;
G := Group(gA,gB);;
y := (1,2,11,10,16,9,6,3,23,19,20,14,21,17,4,8,22,5,18,15,13,7,12);;
x := (3,21)(4,16)(9,22)(10,15)(11,20)(12,19)(13,18)(14,17);;
r := (3,18)(4,16)(5,6)(7,23)(10,14)(12,19)(13,21)(15,17);;

H := Filtered([1..23],point -> point^x=point);;
C := Centralizer(G,x);;
phom := ActionHomomorphism(C,H,OnPoints);;
PA := Image(phom);;
K := Kernel(phom);;

fanoLinesPos := Set(List(
    Filtered(AsList(PA),element -> Order(element)=2),
    element -> Set(Filtered([1..7],point -> point^element=point))
));;
fanoLines := List(fanoLinesPos,line -> Set(List(line,i -> H[i])));;
delta := Intersection(H,Filtered([1..23],point -> point^r=point));;
deltaIndex := Position(fanoLines,Set(delta));;
sourceSheet := 8;;
middleSheet := 1;;
targetSheet := 2;;

FlagOf := function(element)
    local axis,invariant;
    axis := Set(Filtered([1..7],point -> point^element=point));
    invariant := Filtered(fanoLinesPos,line -> OnSets(line,element)=line);
    return [Intersection(invariant)[1],Position(fanoLinesPos,axis)];
end;;

if Size(G)<>10200960
        or H<>[1,2,5,6,7,8,23]
        or Size(C)<>2688
        or Size(PA)<>168
        or StructureDescription(PA)<>"PSL(3,2)"
        or Size(K)<>16
        or Exponent(K)<>2
        or Length(fanoLines)<>7
        or Set(fanoLines)<>[[1,2,8],[1,5,23],[1,6,7],[2,5,7],
            [2,6,23],[5,6,8],[7,8,23]]
        or Set(delta)<>[1,2,8]
        or deltaIndex=fail
        or Set([sourceSheet,middleSheet,targetSheet])<>Set(delta) then
    Error("base envelope data changed");
fi;

# ---- layer 1: internal envelope ----

involutionsC := Filtered(AsList(C),element -> Order(element)=2);;
kernelInvolutions := Filtered(involutionsC,element -> element in K);;
flagInvolutions := Filtered(involutionsC,element -> not element in K);;
kernelVectors := Difference(Elements(K),[(),x]);;

kernelOrbits := OrbitsDomain(C,Set(kernelInvolutions),OnPoints);;

KernelLineData := function(k)
    local image;
    image := Image(phom,Centralizer(C,k));
    return [
        Size(image),
        Filtered([1..7],point -> ForAll(GeneratorsOfGroup(image),
            g -> point^g=point)),
        Filtered([1..7],l -> ForAll(GeneratorsOfGroup(image),
            g -> OnSets(fanoLinesPos[l],g)=fanoLinesPos[l]))
    ];
end;;
kernelLineTable := List(kernelVectors,KernelLineData);;
kernelLineOf := List(kernelLineTable,record -> record[3]);;

pcgsK := Pcgs(K);;
moduleMats := List(GeneratorsOfGroup(C),c -> List(pcgsK,k ->
    ExponentsOfPcElement(pcgsK,k^c))*Z(2)^0);;
kernelModule := GModuleByMats(moduleMats,GF(2));;
minimalSubmodules := MTX.BasesMinimalSubmodules(kernelModule);;
socleElement := PcElementByExponents(pcgsK,
    List(minimalSubmodules[1][1],IntFFE));;

flagTable := List(flagInvolutions,t -> FlagOf(Image(phom,t)));;
flagFibreSizes := Collected(flagTable);;
flagCentralizerData := List([1..Length(flagInvolutions)],position ->
    Filtered(Elements(K),k -> k^flagInvolutions[position]=k));;
flagHeptads := List(flagInvolutions,t ->
    Filtered([1..23],point -> point^t=point));;

flagOrbits := OrbitsDomain(C,Set(flagInvolutions),OnPoints);;

flagHeptadSet := Set(flagHeptads);;
heptadSharingOK := ForAll(flagHeptadSet,heptad ->
    Length(Filtered([1..84],i -> flagHeptads[i]=heptad))=3
    and Length(Set(List(Filtered([1..84],i -> flagHeptads[i]=heptad),
        i -> flagTable[i][2])))=1
    and Set(List(Filtered([1..84],i -> flagHeptads[i]=heptad),
        i -> flagTable[i][1]))
        =fanoLinesPos[flagTable[First([1..84],
            i -> flagHeptads[i]=heptad)][2]]
    and Size(Group(List(Filtered([1..84],i -> flagHeptads[i]=heptad),
        i -> flagInvolutions[i])))=4);;
heptadsPerLine := List([1..7],l -> Length(Set(List(
    Filtered([1..84],i -> flagTable[i][2]=l),
    i -> flagHeptads[i]))));;

if Length(involutionsC)<>99
        or Length(kernelInvolutions)<>15
        or Length(flagInvolutions)<>84
        or not ForAll(involutionsC,t -> IsConjugate(G,t,x))
        or not ForAll(kernelInvolutions,k ->
            Filtered([1..23],point -> point^k=point)=H)
        or SortedList(List(kernelOrbits,Length))<>[1,14]
        or not ForAll(kernelLineTable,record -> record[1]=12
            and record[2]=[] and Length(record[3])=1)
        or List([1..7],l -> Number(kernelLineOf,lines -> lines=[l]))
            <>[2,2,2,2,2,2,2]
        or not ForAll(kernelVectors,k -> KernelLineData(k)[3]
            =KernelLineData(k*x)[3])
        or SortedList(List(MTX.CompositionFactors(kernelModule),
            m -> m.dimension))<>[1,3]
        or Length(minimalSubmodules)<>1
        or Length(minimalSubmodules[1])<>1
        or socleElement<>x
        or not MTX.IsIndecomposable(kernelModule)
        or Length(ComplementClassesRepresentatives(C,K))<>1
        or Length(Set(flagTable))<>21
        or Set(List(flagFibreSizes,pair -> pair[2]))<>[4]
        or not ForAll(flagTable,flag -> flag[1] in fanoLinesPos[flag[2]])
        or not ForAll([1..84],position ->
            Length(flagCentralizerData[position])=4
            and Set(List(Difference(flagCentralizerData[position],[(),x]),
                k -> KernelLineData(k)[3]))
                =[[flagTable[position][2]]])
        or not ForAll([1..84],position ->
            Set(Intersection(H,flagHeptads[position]))
            =fanoLines[flagTable[position][2]])
        or Length(flagHeptadSet)<>28
        or not heptadSharingOK
        or heptadsPerLine<>[4,4,4,4,4,4,4]
        or List(flagOrbits,Length)<>[84] then
    Error("internal envelope certificate changed");
fi;

# ---- layer 2: the real flag ----

branchProduct := x*r;;
heptadR := Filtered([1..23],point -> point^r=point);;
heptadXR := Filtered([1..23],point -> point^branchProduct=point);;
kleinComplement := Difference([1..23],
    Union(H,heptadR,heptadXR));;

MutualFlag := function(involution,frame)
    local centralizer,hom,image,lines,axis,invariant;
    centralizer := Centralizer(G,frame);
    hom := ActionHomomorphism(centralizer,
        Filtered([1..23],point -> point^frame=point),OnPoints);
    image := Image(hom,involution);
    lines := Set(List(Filtered(AsList(Image(hom)),
        element -> Order(element)=2),
        element -> Set(Filtered([1..7],point -> point^element=point))));
    axis := Set(Filtered([1..7],point -> point^image=point));
    invariant := Filtered(lines,line -> OnSets(line,image)=line);
    return [axis,Intersection(invariant)];
end;;

frameFixed := [];;
frameHeptad := fail;;
frameFlag := fail;;
for frameData in [[r,x],[x,r],[branchProduct,x],[x,branchProduct]] do
    frameHeptad := Filtered([1..23],
        point -> point^frameData[2]=point);
    frameFlag := MutualFlag(frameData[1],frameData[2]);
    Add(frameFixed,[
        Set(List(frameFlag[1],i -> frameHeptad[i])),
        List(frameFlag[2],i -> frameHeptad[i])
    ]);
od;

rFlag := FlagOf(Image(phom,r));;
rKernelFix := Filtered(kernelVectors,k -> k^r=k);;

if not r in C
        or Order(branchProduct)<>2
        or not IsConjugate(G,r,x)
        or not IsConjugate(G,branchProduct,x)
        or Intersection(H,heptadR)<>delta
        or Intersection(H,heptadXR)<>delta
        or Intersection(heptadR,heptadXR)<>delta
        or Length(kleinComplement)<>8
        or rFlag[1]<>Position(H,sourceSheet)
        or rFlag[2]<>deltaIndex
        or Length(rKernelFix)<>2
        or Set(List(rKernelFix,k -> KernelLineData(k)[3]))
            <>[[deltaIndex]]
        or not ForAll(frameFixed,record -> record[1]=Set(delta)
            and record[2]=[sourceSheet]) then
    Error("real flag certificate changed");
fi;

# ---- layer 3: duality rigidity and the twin gauges ----

autC := AutomorphismGroup(C);;
pointStabPre := PreImage(phom,Stabilizer(PA,1));;
lineStabPre := PreImage(phom,
    Stabilizer(PA,fanoLinesPos[1],OnSets));;
normalizerC := Normalizer(SymmetricGroup(23),C);;

cyclePoints := [1];;
for index in [2..23] do
    Add(cyclePoints,cyclePoints[index-1]^y);
od;
CycleCoordinate := point -> Position(cyclePoints,point)-1;;
Gauge := j -> PermList(List([1..23],point ->
    cyclePoints[((j-CycleCoordinate(point)) mod 23)+1]));;

heptadOrbit := Set(Orbit(G,Set(H),OnSets));;
gaugeRecords := [];;
for gaugeIndex in [0..22] do
    gaugePerm := Gauge(gaugeIndex);
    gaugeImage := OnSets(Set(H),gaugePerm);
    gaugeCut := Intersection(Set(H),gaugeImage);
    Add(gaugeRecords,[
        gaugeIndex,
        cyclePoints[((12*gaugeIndex) mod 23)+1],
        gaugeCut,
        Length(gaugeCut),
        x^gaugePerm in G,
        gaugeImage in heptadOrbit,
        gaugeCut in List(fanoLines,Set),
        OnSets(gaugeCut,gaugePerm)=gaugeCut
    ]);
od;
cutSizeVector := List(gaugeRecords,record -> record[4]);;
tripleGaugeSheets := Set(List(Filtered(gaugeRecords,
    record -> record[4]=3),record -> record[2]));;
singletonRecords := Filtered(gaugeRecords,record -> record[4]=1);;
tripodGauge := Gauge(7);;
tripodCut := Intersection(Set(H),OnSets(Set(H),tripodGauge));;

if Size(autC)<>1344
        or Size(autC)<>Size(C)/Size(Centre(C))
        or Size(pointStabPre)<>384
        or Size(lineStabPre)<>384
        or IsConjugate(C,pointStabPre,lineStabPre)
        or normalizerC<>C
        or cutSizeVector<>[5,2,3,0,2,0,4,3,2,4,0,1,1,2,2,4,3,2,2,1,0,4,2]
        or ForAny(gaugeRecords,record -> record[5] or record[6]
            or record[7])
        or not ForAll(gaugeRecords,record -> record[8])
        or not ForAll(gaugeRecords,record ->
            (record[4] mod 2=1)=(record[2] in H))
        or tripleGaugeSheets<>[2,8,23]
        or not ForAll(singletonRecords,record ->
            record[3]=[record[2]])
        or Set(List(singletonRecords,record -> record[2]))<>[5,6,7]
        or gaugeRecords[1][2]<>1 or gaugeRecords[1][4]<>5
        or sourceSheet^tripodGauge<>sourceSheet
        or Number(gaugeRecords,record -> record[2]=sourceSheet)<>1
        or gaugeRecords[8][2]<>sourceSheet
        or tripodCut<>[2,6,8]
        or OnSets(tripodCut,tripodGauge)<>tripodCut
        or targetSheet^tripodGauge<>6 then
    Error("duality rigidity certificate changed");
fi;

# ---- layer 4: global cut reconstruction ----

classX := ConjugacyClass(G,x);;
compatible := Filtered(AsList(classX),element ->
    Order(element*y)=23 and IsConjugate(G,element*y,y));;
translationOrbits := OrbitsDomain(Group(y),Set(compatible),OnPoints);;
distinguishedIndex := First([1..Length(translationOrbits)],
    i -> x in translationOrbits[i]);;

orbitKeys := List(translationOrbits,orbit -> SortedList(
    List([1..22],n -> Order(orbit[1]*(orbit[1]^(y^n))))));;

lineMultTable := [];;
pointMultTable := [];;
nonlineTable := [];;
commutingCounts := [];;
for orbitIndex in [1..Length(translationOrbits)] do
    lineRow := List([1..7],l -> 0);
    pointRow := List([1..7],p -> 0);
    nonlineRow := [];
    for element in translationOrbits[orbitIndex] do
        elementCut := Intersection(H,
            Filtered([1..23],point -> point^element=point));
        if Length(elementCut)=3 then
            if Set(elementCut) in fanoLines then
                lineRow[Position(fanoLines,Set(elementCut))] :=
                    lineRow[Position(fanoLines,Set(elementCut))]+1;
            else
                Add(nonlineRow,Set(elementCut));
            fi;
        elif Length(elementCut)=1 then
            pointRow[Position(H,elementCut[1])] :=
                pointRow[Position(H,elementCut[1])]+1;
        elif Length(elementCut)=7 then
            if element<>x then
                Error("unexpected full heptad cut");
            fi;
        else
            Error("cut size outside the Steiner pattern");
        fi;
    od;
    Add(lineMultTable,lineRow);
    Add(pointMultTable,pointRow);
    Add(nonlineTable,Collected(nonlineRow));
    Add(commutingCounts,Number(translationOrbits[orbitIndex],
        element -> element*x=x*element));
od;

lineColumnSums := List([1..7],l -> Sum(List(lineMultTable,
    row -> row[l])));;
avoidedLines := Filtered([1..7],l -> lineColumnSums[l]=0);;
sourcePosition := Position(H,sourceSheet);;
sourceAvoidingClasses := Filtered([1..Length(translationOrbits)],
    i -> pointMultTable[i][sourcePosition]=0);;
distinguishedZeroPoints := Filtered([1..7],
    p -> pointMultTable[distinguishedIndex][p]=0);;
distinguishedLineSupport := Filtered([1..7],
    l -> lineMultTable[distinguishedIndex][l]>0);;
middlePencil := Filtered([1..7],l -> middleSheet in fanoLines[l]
    and l<>deltaIndex);;

commutingNontrivial := Filtered(compatible,element ->
    element<>x and element*x=x*element);;
commutingAxes := List(commutingNontrivial,element ->
    FlagOf(Image(phom,element))[2]);;

conjugationPairing := [];;
representative := fail;;
transporter := fail;;
for orbitIndex in [1..Length(translationOrbits)] do
    representative := translationOrbits[orbitIndex][1];
    transporter := RepresentativeAction(G,representative*y,y);
    Add(conjugationPairing,First([1..Length(translationOrbits)],
        i -> representative^transporter in translationOrbits[i]));
od;
fixedClasses := Filtered([1..Length(translationOrbits)],
    i -> conjugationPairing[i]=i);;

if Size(classX)<>3795
        or Length(compatible)<>161
        or Length(translationOrbits)<>7
        or Set(List(translationOrbits,Length))<>[23]
        or distinguishedIndex=fail
        or Length(Set(orbitKeys))<>7
        or Number(kernelVectors,k -> k in compatible)<>0
        or not r in C or r in compatible
        or avoidedLines<>[deltaIndex]
        or SortedList(lineColumnSums)<>[0,2,2,2,2,3,4]
        or Sum(lineColumnSums)<>15
        or sourceAvoidingClasses<>[distinguishedIndex]
        or distinguishedZeroPoints<>[sourcePosition]
        or pointMultTable[distinguishedIndex]<>[1,3,2,2,3,0,1]
        or distinguishedLineSupport<>middlePencil
        or List(middlePencil,l -> lineMultTable[distinguishedIndex][l])
            <>[1,1]
        or SortedList(List(lineMultTable,Sum))<>[1,1,2,2,2,2,5]
        or SortedList(List([1..7],i -> Sum(lineMultTable[i])
            +Sum(List(nonlineTable[i],pair -> pair[2]))))
            <>[10,13,13,13,13,13,13]
        or SortedList(commutingCounts)<>[0,0,0,0,1,1,3]
        or commutingCounts[distinguishedIndex]<>1
        or Length(commutingNontrivial)<>4
        or deltaIndex in commutingAxes
        or List(conjugationPairing,i -> conjugationPairing[i])
            <>[1..7]
        or Length(fixedClasses)<>3
        or not distinguishedIndex in fixedClasses then
    Error("global cut reconstruction certificate changed");
fi;

# ---- outputs ----

sortedOrder := SortingPerm(orbitKeys);;
canonicalOrder := List([1..7],i -> i/sortedOrder);;

output := OutputTextFile("results/fano_flag_descent_summary.txt",false);;
SetPrintFormattingStatus(output,false);;
AppendTo(output,
    "schema=m23.cover-investigation.fano-flag-descent.v1\n",
    "heptad=",H,"\n",
    "fano_lines=",fanoLines,"\n",
    "axis_line=",Set(delta),"\n",
    "axis_line_index=",deltaIndex,"\n",
    "flag_center_sheet=",sourceSheet,"\n",
    "involutions_in_centralizer=",Length(involutionsC),"\n",
    "kernel_involutions=",Length(kernelInvolutions),"\n",
    "flag_involutions=",Length(flagInvolutions),"\n",
    "all_centralizer_involutions_in_2A=true\n",
    "kernel_involutions_all_fix_heptad_pointwise=true\n",
    "kernel_orbit_sizes=",SortedList(List(kernelOrbits,Length)),"\n",
    "kernel_vector_stabilizer_image_order=12\n",
    "kernel_vector_fixes_no_point_one_line=true\n",
    "kernel_line_fibres_are_x_shift_pairs=true\n",
    "kernel_module_composition_dims=[1,3]\n",
    "kernel_module_unique_minimal_submodule=inertia_socle\n",
    "kernel_module_indecomposable=true\n",
    "line_scheme_is_K_mod_inertia=true\n",
    "flag_count=",Length(Set(flagTable)),"\n",
    "flag_lift_fibre_size=4\n",
    "flag_involution_heptads_distinct=",Length(flagHeptadSet),"\n",
    "flag_heptads_per_line=4\n",
    "flag_heptad_sharers_form_klein_four_with_axis_centers=true\n",
    "flag_heptad_cut_equals_axis=true\n",
    "flag_axis_equals_fixed_kernel_pair_label=true\n",
    "reflection_in_centralizer=true\n",
    "klein_four_heptads_pairwise_meet_in_axis=true\n",
    "klein_four_complement_size=",Length(kleinComplement),"\n",
    "klein_four_complement=",kleinComplement,"\n",
    "mutual_transvection_flags_all_equal=(8,[1,2,8])\n",
    "reflection_fixed_kernel_pair_is_axis_pair=true\n",
    "aut_C_order=",Size(autC),"\n",
    "out_C_trivial=true\n",
    "point_line_stabilizers_nonconjugate=true\n",
    "normalizer_in_S23_is_C=true\n",
    "gauge_cut_size_vector=",cutSizeVector,"\n",
    "gauge_images_never_G_heptads=true\n",
    "gauge_conjugates_of_x_never_in_G=true\n",
    "gauge_cuts_never_fano_lines=true\n",
    "gauge_cut_parity_law=fixed_sheet_in_heptad\n",
    "triple_cut_gauge_fixed_sheets=",tripleGaugeSheets,"\n",
    "singleton_cut_gauge_fixed_sheets=[5,6,7]\n",
    "tripod_gauge_fixes_center=true\n",
    "tripod_gauge_cut=",tripodCut,"\n",
    "tripod_gauge_cut_is_invariant_triangle=true\n",
    "compatible_involutions=",Length(compatible),"\n",
    "translation_orbits=7x23\n",
    "distinguished_orbit_scan_index=",distinguishedIndex,"\n",
    "line_column_sums=",lineColumnSums,"\n",
    "unique_universally_avoided_line=axis\n",
    "unique_center_avoiding_class=class_6\n",
    "distinguished_point_row=",pointMultTable[distinguishedIndex],"\n",
    "distinguished_line_support=punctured_middle_pencil\n",
    "line_cut_totals_sorted=",SortedList(List(lineMultTable,Sum)),"\n",
    "commuting_counts_sorted=",SortedList(commutingCounts),"\n",
    "no_compatible_involution_has_axis_delta=true\n",
    "reflection_not_y_compatible=true\n",
    "conjugation_fixed_classes=",Length(fixedClasses),"\n",
    "distinguished_class_conjugation_fixed=true\n",
    "remaining_step=galois_covariance_of_cut_correspondence_mod_frame_cocycle\n"
);;
CloseStream(output);;

tsv := OutputTextFile("results/fano_flag_descent_class_cuts.tsv",false);;
SetPrintFormattingStatus(tsv,false);;
AppendTo(tsv,"canonical_rank\tscan_index\tis_distinguished\t",
    "conjugation_partner_scan_index\tcommuting_count\t",
    "line_mults_by_fano_line\tpoint_mults_by_heptad_sheet\t",
    "nonline_triples\tproduct_order_key\n");;
for rankIndex in [1..7] do
    orbitIndex := canonicalOrder[rankIndex];
    AppendTo(tsv,
        rankIndex,"\t",
        orbitIndex,"\t",
        orbitIndex=distinguishedIndex,"\t",
        conjugationPairing[orbitIndex],"\t",
        commutingCounts[orbitIndex],"\t",
        lineMultTable[orbitIndex],"\t",
        pointMultTable[orbitIndex],"\t",
        nonlineTable[orbitIndex],"\t",
        orbitKeys[orbitIndex],"\n");
od;
CloseStream(tsv);;

Print("envelope: 99 involutions = 15 kernel + 84 flag, all 2A\n");
Print("lines = K/<x>: unique minimal submodule is the inertia socle\n");
Print("real flag: mutual transvection flags all (8,",Set(delta),")\n");
Print("rigidity: |Aut(C)|=",Size(autC),", N_S23(C)=C, gauges exit Steiner\n");
Print("cut table line column sums=",lineColumnSums,
    " avoided=",avoidedLines," (axis)\n");
Print("class-6 point row=",pointMultTable[distinguishedIndex],
    " avoids only the center\n");
Print("class-6 line support=",distinguishedLineSupport,
    " = punctured middle pencil ",middlePencil,"\n");
Print("conjugation pairing=",conjugationPairing,
    " fixed=",fixedClasses,"\n");

QUIT_GAP(0);
