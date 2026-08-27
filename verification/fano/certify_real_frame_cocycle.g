# Certify that the real structure acts on the x-framed Hurwitz classes as
# left multiplication by the inertia, that its frame cocycle at each
# real-fixed class is a unique internal involution of C_M23(x), and that
# at Nielsen class 6 this cocycle is exactly the branch reflection,
# whose transvection flag is the marked flag (8,{1,2,8}).
#
# The group-side output predicts the real-place arithmetic of the
# special-fibre fields: 3 real unramified sheets forming the axis, 3 real
# lines, and 4 real ramified double points.  The PARI companion
# certificate verifies these counts exactly.

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

FlagOf := function(element)
    local axis,invariant;
    axis := Set(Filtered([1..7],point -> point^element=point));
    invariant := Filtered(fanoLinesPos,line -> OnSets(line,element)=line);
    return [Intersection(invariant)[1],Position(fanoLinesPos,axis)];
end;;

# ---- the x-frame real structure ----

pairStabilizer := Intersection(Centralizer(G,x),Centralizer(G,y));;
inertiaCentralizer := Centralizer(G,y);;
conjugatePartner := x*y;;

if Size(G)<>10200960
        or Size(pairStabilizer)<>1
        or inertiaCentralizer<>Group(y)
        or Order(conjugatePartner)<>23
        or not IsConjugate(G,conjugatePartner,y)
        or x*conjugatePartner<>y then
    Error("x-frame real structure data changed");
fi;

# ---- the class-6 cocycle is the branch reflection ----

transporterBase := RepresentativeAction(G,y,conjugatePartner);;
transporterCoset := List(Elements(inertiaCentralizer),
    element -> element*transporterBase);;
cocycleSolutions := Filtered(transporterCoset,element -> element in C);;

if transporterBase=fail
        or Length(cocycleSolutions)<>1
        or cocycleSolutions[1]<>r
        or y^r<>conjugatePartner
        or x^r<>x
        or Order(r)<>2
        or r in K
        or FlagOf(Image(phom,r))<>[Position(H,8),deltaIndex]
        or Set(delta)<>[1,2,8] then
    Error("class-6 real frame cocycle changed");
fi;

# ---- group-side real-place predictions ----

ramifiedPairs := List(Orbits(Group(x),
    Filtered([1..23],point -> point^x<>point)),Set);;
fixedPairs := Filtered(ramifiedPairs,pair -> OnSets(pair,r)=pair);;
fixedSheets := Filtered(H,point -> point^r=point);;
fixedLines := Filtered(fanoLines,line ->
    OnSets(Set(line),r)=Set(line));;

if Length(ramifiedPairs)<>8
        or Length(fixedSheets)<>3
        or Set(fixedSheets)<>Set(delta)
        or Length(fixedLines)<>3
        or Length(fixedPairs)<>4 then
    Error("real-place prediction data changed");
fi;

# ---- the y-frame classes and the conjugation pairing ----

classX := ConjugacyClass(G,x);;
compatible := Filtered(AsList(classX),element ->
    Order(element*y)=23 and IsConjugate(G,element*y,y));;
translationOrbits := OrbitsDomain(Group(y),Set(compatible),OnPoints);;
distinguishedIndex := First([1..Length(translationOrbits)],
    i -> x in translationOrbits[i]);;

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

if Length(compatible)<>161
        or Length(translationOrbits)<>7
        or List(conjugationPairing,i -> conjugationPairing[i])<>[1..7]
        or Length(fixedClasses)<>3
        or not distinguishedIndex in fixedClasses then
    Error("conjugation pairing changed");
fi;

# ---- x-frame consistency and the cocycle trichotomy ----

muPairing := [];;
cocycleCounts := [];;
cocycleOrders := [];;
cocycleInKernel := [];;
cocycleFlagExamples := [];;
frameRep := fail;;
framePartner := fail;;
frameReturn := fail;;
for orbitIndex in [1..Length(translationOrbits)] do
    representative := translationOrbits[orbitIndex][1];
    frameRep := RepresentativeAction(G,representative,x);
    framePartner := y^frameRep;
    if Order(framePartner)<>23 or x*framePartner=framePartner then
        Error("x-frame transport failed");
    fi;
    frameReturn := RepresentativeAction(G,x*framePartner,y);
    Add(muPairing,First([1..Length(translationOrbits)],
        i -> x^frameReturn in translationOrbits[i]));
    transporter := RepresentativeAction(G,framePartner,x*framePartner);
    if transporter=fail then
        Error("partner transport failed");
    fi;
    cocycleSolutions := Filtered(List(
        Elements(Centralizer(G,framePartner)),
        element -> element*transporter),element -> element in C);
    Add(cocycleCounts,Length(cocycleSolutions));
    if Length(cocycleSolutions)=1 then
        Add(cocycleOrders,Order(cocycleSolutions[1]));
        Add(cocycleInKernel,cocycleSolutions[1] in K);
        Add(cocycleFlagExamples,FlagOf(Image(phom,
            cocycleSolutions[1])));
    else
        Add(cocycleOrders,0);
        Add(cocycleInKernel,fail);
        Add(cocycleFlagExamples,fail);
    fi;
od;

if muPairing<>conjugationPairing
        or List([1..7],i -> cocycleCounts[i]=1)
            <>List([1..7],i -> conjugationPairing[i]=i)
        or Filtered(cocycleOrders,order -> order<>0)<>[2,2,2]
        or Filtered(cocycleInKernel,value -> value<>fail)
            <>[false,false,false] then
    Error("cocycle trichotomy changed");
fi;

# ---- outputs ----

output := OutputTextFile("results/real_frame_cocycle_summary.txt",false);;
SetPrintFormattingStatus(output,false);;
AppendTo(output,
    "schema=m23.cover-investigation.real-frame-cocycle.v1\n",
    "pair_stabilizer_trivial=true\n",
    "x_frame_real_structure=left_multiplication_by_inertia\n",
    "conjugate_partner_valid=true\n",
    "distinguished_cocycle_solution_count=1\n",
    "distinguished_cocycle_is_branch_reflection=true\n",
    "distinguished_cocycle_order=2\n",
    "distinguished_cocycle_in_kernel=false\n",
    "distinguished_cocycle_flag=(8,[1,2,8])\n",
    "distinguished_cocycle_flag_is_marked_flag=true\n",
    "real_place_prediction_unramified_sheets=3\n",
    "real_place_prediction_real_sheets_form_axis=true\n",
    "real_place_prediction_lines=3\n",
    "real_place_prediction_ramified_pairs=4\n",
    "conjugation_pairing_fixed_classes=",Length(fixedClasses),"\n",
    "distinguished_class_conjugation_fixed=true\n",
    "mu_pairing_matches_conjugation_pairing=true\n",
    "cocycle_counts_by_class=",cocycleCounts,"\n",
    "cocycle_exists_iff_class_real_fixed=true\n",
    "all_cocycles_are_involutions=true\n",
    "no_cocycle_in_kernel=true\n",
    "cocycles_act_by_transvections=true\n",
    "frame_returns_confined_to_centralizer=true\n",
    "remaining_step=inertia_rigid_returns_for_all_galois_elements\n"
);;
CloseStream(output);;

Print("real structure in the x-frame: y' -> x*y', exact\n");
Print("class-6 frame cocycle = branch reflection r, flag (8,",
    Set(delta),")\n");
Print("cocycle counts by class=",cocycleCounts,
    " (1 iff real-fixed)\n");
Print("gauge-dependent cocycle flag examples=",
    cocycleFlagExamples,"\n");
Print("real-place predictions: 3 real sheets (the axis), 3 real lines, ",
    "4 real double points\n");

QUIT_GAP(0);
