# Exact relative-transporter invariants on all seven Nielsen IDs.

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

if Size(G)<>10200960 or not IsSimpleGroup(G) then
    Error("the supplied permutations do not generate M23");
fi;

CountPacketProductOrder := function(x,c,targetOrder)
    return Number([1..22],d -> Order(x*x^(c^d))=targetOrder);
end;;

Analyze := function(x)
    local z,nu,kappa2,kappa4,Ny,Nz,overlap,yReturns,zReturns,
          relativeTransporters,Cx,inertiaBlocks,BlockPermutation,blockImages,
          blockImage,multiplicities,mu,binaryWeight;

    z := (x*y)^-1;
    nu := Number(Cartesian([1..22],[1..22]),exponents ->
        Order(y^exponents[1]*z^exponents[2])=23
    );
    kappa2 := CountPacketProductOrder(x,y,2)
        + CountPacketProductOrder(x,z,2);
    kappa4 := CountPacketProductOrder(x,y,4)
        + CountPacketProductOrder(x,z,4);

    Ny := Normalizer(G,Group(y));
    Nz := Normalizer(G,Group(z));
    overlap := Intersection(
        Set(Orbit(Ny,x,OnPoints)),
        Set(Orbit(Nz,x,OnPoints))
    );
    yReturns := List(overlap,w -> RepresentativeAction(Ny,x,w,OnPoints));
    zReturns := List(overlap,w -> RepresentativeAction(Nz,x,w,OnPoints));
    if fail in yReturns or fail in zReturns then
        Error("an endpoint transporter is missing");
    fi;
    relativeTransporters := List([1..Length(overlap)],position ->
        yReturns[position]*zReturns[position]^-1
    );
    if Length(Set(relativeTransporters))<>Length(overlap)
            or ForAny(relativeTransporters,c -> x^c<>x) then
        Error("relative transporters changed");
    fi;

    Cx := Centralizer(G,x);
    inertiaBlocks := Orbits(Group(x),MovedPoints(x),OnPoints);
    if Length(inertiaBlocks)<>8
            or Set(List(inertiaBlocks,Length))<>[2] then
        Error("the inertia block system changed");
    fi;
    BlockPermutation := function(c)
        return PermList(List(inertiaBlocks,block -> Position(
            inertiaBlocks,Set(List(block,point -> point^c))
        )));
    end;
    blockImages := List(relativeTransporters,BlockPermutation);
    blockImage := Set(blockImages);
    multiplicities := List(blockImage,p -> Number(blockImages,q -> q=p));
    if ForAny(multiplicities,m -> not m in [1,2]) then
        Error("an unexpected affine-image multiplicity occurred");
    fi;
    mu := Length(blockImage);
    binaryWeight := Number(multiplicities,IsOddInt);

    return rec(
        nu := nu,
        kappa2 := kappa2,
        kappa4 := kappa4,
        overlapSize := Length(overlap),
        mu := mu,
        binaryWeight := binaryWeight,
        binaryAugmentation := Length(overlap) mod 2
    );
end;;

records := List(representatives,Analyze);;
nuValues := List(records,record -> record.nu);;
kappa2Values := List(records,record -> record.kappa2);;
kappa4Values := List(records,record -> record.kappa4);;
overlapSizes := List(records,record -> record.overlapSize);;
muValues := List(records,record -> record.mu);;
binaryWeights := List(records,record -> record.binaryWeight);;
binaryAugmentations := List(records,record -> record.binaryAugmentation);;

nuSelector := List(nuValues,value -> value=32);;
muSelector := List(muValues,value -> value=17);;
cyclicCountSelector := List([1..7],id ->
    kappa2Values[id]=0 and kappa4Values[id]=12
);;
binarySelector := List(binaryAugmentations,value -> value=0);;
distinguishedSelector := [false,false,false,false,false,true,false];;
sexticIds := [1,2,3,4,5,7];;

if nuValues<>[54,54,28,46,46,32,42]
        or kappa2Values<>[0,0,0,2,2,0,4]
        or kappa4Values<>[16,16,8,6,6,12,12]
        or overlapSizes<>[17,17,37,29,29,18,17]
        or muValues<>[16,16,31,28,28,17,14]
        or binaryWeights<>[15,15,25,27,27,16,11]
        or binaryAugmentations<>[1,1,1,1,1,0,1]
        or nuSelector<>distinguishedSelector
        or muSelector<>distinguishedSelector
        or cyclicCountSelector<>distinguishedSelector
        or binarySelector<>distinguishedSelector
        or Length(Set(List(sexticIds,id -> nuValues[id])))=1
        or Length(Set(List(sexticIds,id -> muValues[id])))=1
        or Length(Set(List(sexticIds,id ->
            [kappa2Values[id],kappa4Values[id]])))=1
        or Length(Set(List(sexticIds,id -> binaryWeights[id])))=1
        or Length(Set(List(sexticIds,id -> binaryAugmentations[id])))<>1 then
    Error("the full-seven relative-transporter data changed");
fi;

Print("nu_by_id=",nuValues,"\n");
Print("mu_by_id=",muValues,"\n");
Print("kappa2_by_id=",kappa2Values,"\n");
Print("kappa4_by_id=",kappa4Values,"\n");
Print("binary_theta_weight_by_id=",binaryWeights,"\n");
Print("binary_augmentation_by_id=",binaryAugmentations,"\n");
Print("all_four_boolean_selectors=degree_one_idempotent\n");
Print("raw_nu_mu_cyclic_counts_and_theta_do_not_descend_on_the_sextic_orbit\n");
Print("binary_augmentation_does_descend_on_the_1_plus_6_decomposition\n");
Print("PASS_HURWITZ_RELATIVE_TRANSPORTER_FINITE_DATA\n");

QUIT_GAP(0);
