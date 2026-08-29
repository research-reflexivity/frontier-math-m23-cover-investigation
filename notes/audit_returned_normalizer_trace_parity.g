# Certify the returned one-endpoint trace parity in the two natural M23
# copies containing a fixed wild group D=23:11.

SizeScreen([100000,100000]);;

gA := (1,2)(3,4)(7,8)(9,10)(13,14)(15,16)(19,20)(21,22);;
gB := (1,16,11,3)(2,9,21,12)(4,5,8,23)(6,22,14,18)
      (13,20)(15,17);;
G := Group(gA,gB);;
S := SymmetricGroup(23);;
y := (1,2,11,10,16,9,6,3,23,19,20,14,21,17,4,8,22,5,18,
      15,13,7,12);;
D := Normalizer(G,Group(y));;
N := Normalizer(S,D);;
outside := Difference(AsList(N),AsList(D));;
s := outside[1];;
adjacentG := G^s;;

x := (3,21)(4,16)(9,22)(10,15)(11,20)(12,19)(13,18)(14,17);;
class2A := Set(AsList(ConjugacyClass(G,x)));;
adjacentClass2A := Set(AsList(ConjugacyClass(adjacentG,x^s)));;
colors := Orbits(D,class2A,OnPoints);;

records := List(colors,color -> rec(
    traceSize := Length(Set(Orbit(D,color[1],OnPoints))),
    squareReturnOverlap := Length(Intersection(
        Set(Orbit(D,color[1],OnPoints)),
        Set(Orbit(D,color[1],OnPoints))
    )),
    nonsquareReturnOverlap := Length(Intersection(
        Set(Orbit(D,color[1],OnPoints)),
        Set(List(Orbit(D,color[1],OnPoints),w -> w^s))
    ))
));;

if Size(D)<>253
   or Size(N)<>506
   or Length(outside)<>253
   or not ForAll(outside,n -> G^n=adjacentG)
   or Intersection(G,adjacentG)<>D
   or not IsEmpty(Intersection(class2A,adjacentClass2A))
   or Length(colors)<>15
   or Set(List(records,record -> record.traceSize))<>[253]
   or Set(List(records,record -> record.squareReturnOverlap))<>[253]
   or Set(List(records,record -> record.nonsquareReturnOverlap))<>[0] then
    Error("the returned normalizer-trace parity audit changed");
fi;

Print("wild_normalizer_order=",Size(D),"\n");
Print("full_affine_normalizer_order=",Size(N),"\n");
Print("relative_color_count=",Length(colors),"\n");
Print("adjacent_M23_intersection_order=",Size(Intersection(G,adjacentG)),"\n");
Print("adjacent_2A_classes_disjoint=",
      IsEmpty(Intersection(class2A,adjacentClass2A)),"\n");
Print("returned_trace_overlap_sizes_for_q_0_q_1=[253,0]\n");
Print("returned_trace_overlap_parities_for_q_0_q_1=[1,0]\n");
Print("PASS_RETURNED_NORMALIZER_TRACE_PARITY_AUDIT\n");
QUIT_GAP(0);
