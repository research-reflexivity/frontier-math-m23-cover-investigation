/* Independent Magma certificate for the Fano--affine odd fixed-point lemma. */

procedure Must(flag, message)
    if not flag then
        print "CERTIFICATE_FAILURE: " cat message;
        quit;
    end if;
end procedure;

S := Sym(23);
generatorA := S!(1,2)(3,4)(7,8)(9,10)(13,14)(15,16)(19,20)(21,22);
generatorB := S!(1,16,11,3)(2,9,21,12)(4,5,8,23)(6,22,14,18)
                 (13,20)(15,17);
G := sub<S | generatorA, generatorB>;
x := S!(3,21)(4,16)(9,22)(10,15)(11,20)(12,19)(13,18)(14,17);
C := Centralizer(G,x);

xOrbits := Orbits(sub<G | x>);
H := [Setseq(orbit)[1] : orbit in xOrbits | #orbit eq 1];
A := [orbit : orbit in xOrbits | #orbit eq 2];

Must(#G eq 10200960 and #C eq 2688 and #H eq 7 and #A eq 8,
     "the M23 tame-orbit data changed");

records := [];
for c in C do
    fixedH := #[point : point in H | point^c eq point];
    fixedA := #[block : block in A |
        {point^c : point in block} eq block];
    Must((fixedH+fixedA) mod 2 eq 1,
         "a centralizer element has even total fixed-point count");
    Append(~records,<fixedH,fixedA>);
end for;

profileCounts := [
    #[record : record in records | record eq <7,0>],
    #[record : record in records | record eq <7,8>],
    #[record : record in records | record eq <3,0>],
    #[record : record in records | record eq <3,4>],
    #[record : record in records | record eq <1,0>],
    #[record : record in records | record eq <1,2>],
    #[record : record in records | record eq <0,1>]
];

Must(profileCounts eq [14,2,252,84,784,784,768] and
     &+profileCounts eq #C,
     "the Fano--affine fixed-point profile changed");

print "centralizer_order=2688";
print "fixed_heptad_size=7";
print "affine_block_count=8";
print "profile_pairs=[<7,0>,<7,8>,<3,0>,<3,4>,<1,0>,<1,2>,<0,1>]";
print "profile_counts=[14,2,252,84,784,784,768]";
print "all_centralizer_elements_have_odd_total_fixed_points=true";
print "PASS_FANO_AFFINE_ODD_FIXED_POINT_LEMMA_MAGMA";
