# Exact M23 realization of the universal Fano--affine odd fixed-point
# lemma.  For every element of C_M23(x), the sum of the numbers of fixed
# points on the seven singleton x-orbits and the eight double x-orbits is
# odd.

SizeScreen([100000,100000]);;

gA := (1,2)(3,4)(7,8)(9,10)(13,14)(15,16)(19,20)(21,22);;
gB := (1,16,11,3)(2,9,21,12)(4,5,8,23)(6,22,14,18)
      (13,20)(15,17);;
G := Group(gA,gB);;
x := (3,21)(4,16)(9,22)(10,15)(11,20)(12,19)(13,18)(14,17);;
C := Centralizer(G,x);;
H := Set(Filtered([1..23],point -> point^x=point));;
A := List(Orbits(Group(x),Difference([1..23],H),OnPoints),Set);;
pointHom := ActionHomomorphism(C,H,OnPoints);;
pairHom := ActionHomomorphism(C,A,OnSets);;
F := Image(pointHom);;
Aff := Image(pairHom);;
K := Kernel(pointHom);;
translations := Image(pairHom,K);;

records := List(Elements(C),function(c)
    local fixedH,fixedA,linearOrder;
    fixedH := Number(H,point -> point^c=point);
    fixedA := Number(A,block -> OnSets(block,c)=block);
    linearOrder := Order(Image(pointHom,c));
    return [linearOrder,fixedH,fixedA,(fixedH+fixedA) mod 2];
end);;
profile := Collected(SortedList(records));;

if Size(G)<>10200960
        or Size(C)<>2688
        or Length(H)<>7
        or Length(A)<>8
        or Size(F)<>168
        or StructureDescription(F)<>"PSL(3,2)"
        or Size(Aff)<>1344
        or StructureDescription(Aff)<>"(C2 x C2 x C2) : PSL(3,2)"
        or Size(K)<>16
        or Size(translations)<>8
        or not IsRegular(translations,[1..8])
        or ForAny(records,record -> record[4]<>1) then
    Error("Fano--affine odd fixed-point lemma changed");
fi;

Print("centralizer_order=",Size(C),"\n");
Print("fixed_heptad_size=",Length(H),"\n");
Print("affine_block_count=",Length(A),"\n");
Print("linear_quotient=PSL3(2)\n");
Print("affine_quotient=2^3:PSL3(2)\n");
Print("record_fields=[linear_order,fixed_Fano_points,fixed_affine_points,",
      "sum_mod_2]\n");
Print("profile=",profile,"\n");
Print("all_2688_centralizer_elements_have_odd_total_fixed_points=true\n");
Print("PASS_FANO_AFFINE_ODD_FIXED_POINT_LEMMA\n");

QUIT_GAP(0);
