# Explore the honest common-tag coefficient at the tame 2A branch.
# This deliberately intersects the endpoint trace label sets before
# pushing them to the degree-23 sheet-pair branches.

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

Analyze := function(x)
    local z,Dz,Ty,Tz,T,orbits,fixed,moved,diagCounts,swapCounts,
          nodeBits,fixedBits,totalNode,totalFixed,w,block,i;
    z := (x*y)^-1;
    Dz := Normalizer(G,Group(z));
    Ty := Set(Orbit(D,x,OnPoints));
    Tz := Set(Orbit(Dz,x,OnPoints));
    T := Intersection(Ty,Tz);
    orbits := Orbits(Group(x),[1..23]);
    fixed := Filtered(orbits,b -> Length(b)=1);
    moved := Filtered(orbits,b -> Length(b)=2);
    diagCounts := [];
    swapCounts := [];
    nodeBits := [];
    for block in moved do
        i := block[1];
        Add(diagCounts,Number(T,w -> i^w=i));
        Add(swapCounts,Number(T,w -> i^w=i^x));
        Add(nodeBits,(diagCounts[Length(diagCounts)]
                      -swapCounts[Length(swapCounts)]) mod 2);
    od;
    fixedBits := List(fixed,b -> Number(T,w -> b[1]^w=b[1]) mod 2);
    totalNode := Sum(nodeBits) mod 2;
    totalFixed := Sum(fixedBits) mod 2;
    return rec(
        overlap := Length(T),
        diagCounts := diagCounts,
        swapCounts := swapCounts,
        nodeBits := nodeBits,
        fixedBits := fixedBits,
        totalNode := totalNode,
        totalFixed := totalFixed
    );
end;;

records := List(representatives,Analyze);;
for i in [1..7] do
    Print("id_",i,
          " overlap=",records[i].overlap,
          " diag=",records[i].diagCounts,
          " swap=",records[i].swapCounts,
          " node_bits=",records[i].nodeBits,
          " fixed_bits=",records[i].fixedBits,
          " totals=[",records[i].totalNode,",",records[i].totalFixed,"]\n");
od;

QUIT_GAP(0);
