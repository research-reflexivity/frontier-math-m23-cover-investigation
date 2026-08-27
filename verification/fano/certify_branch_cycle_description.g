# Certify the group side of the branch cycle description: the class 2A
# has cycle shape 1^7 2^8, matching the branch fibre of the class-6
# curve; the triple (x,y,(xy)^-1) is a generating triple of the Nielsen
# class; the genus is four; and -- the point that explains the gate's
# square/nonsquare dichotomy -- the power y^lambda is conjugate to y
# exactly when lambda is a quadratic residue modulo 23.  Hence the two
# order-23 branch classes are exchanged precisely by the nonsquare part
# of the cyclotomic character, which is the Galois action on
# Q(sqrt(-23)), the quadratic subfield of Q(zeta_23) and the field of
# definition of the two order-23 branch points of the curve.

SizeScreen([100000,100000]);;

gA := (1,2)(3,4)(7,8)(9,10)(13,14)(15,16)(19,20)(21,22);;
gB := (1,16,11,3)(2,9,21,12)(4,5,8,23)(6,22,14,18)(13,20)(15,17);;
G := Group(gA,gB);;
y := (1,2,11,10,16,9,6,3,23,19,20,14,21,17,4,8,22,5,18,15,13,7,12);;
x := (3,21)(4,16)(9,22)(10,15)(11,20)(12,19)(13,18)(14,17);;

squares := Filtered([1..22],residue -> Jacobi(residue,23)=1);;
selfConjugatePowers := Filtered([1..22],residue ->
    IsConjugate(G,y^residue,y));;
partner := (x*y)^-1;;

ramificationSum := (23-15)+(23-1)+(23-1);;
genus := (-2*23+ramificationSum+2)/2;;

if Size(G)<>10200960
        or Collected(CycleLengths(x,[1..23]))<>[[1,7],[2,8]]
        or Order(x)<>2
        or Order(y)<>23
        or Order(x*y)<>23
        or x*y*partner<>()
        or Group(x,y)<>G
        or not IsConjugate(G,x*y,y)
        or IsConjugate(G,partner,y)
        or Jacobi(22,23)<>-1
        or selfConjugatePowers<>squares
        or Length(squares)<>11
        or ramificationSum<>52
        or genus<>4 then
    Error("branch cycle group certificate changed");
fi;

# the branch swap is exactly the nonsquare part of the cyclotomic action:
# y^lambda lands in the opposite order-23 class exactly for nonsquares
nonsquares := Difference([1..22],squares);;
if not ForAll(nonsquares,residue -> IsConjugate(G,y^residue,partner))
        or ForAny(squares,residue -> IsConjugate(G,y^residue,partner))
        or not ForAll(squares,residue -> IsConjugate(G,y^residue,y)) then
    Error("branch swap characterization changed");
fi;

output := OutputTextFile("results/branch_cycle_group_summary.txt",false);;
SetPrintFormattingStatus(output,false);;
AppendTo(output,
    "schema=m23.cover-investigation.branch-cycle-group.v1\n",
    "class_2A_cycle_shape=1^7_2^8\n",
    "special_fibre_matches_2A=true\n",
    "generating_triple_product_is_identity=true\n",
    "triple_generates_M23=true\n",
    "order23_classes_are_distinct=true\n",
    "self_conjugate_powers=quadratic_residues\n",
    "branch_swap_is_the_nonsquare_cyclotomic_action=true\n",
    "quadratic_subfield=Q(sqrt(-23))\n",
    "ramification_sum=52\n",
    "genus=4\n"
);;
CloseStream(output);;

Print("2A cycle shape 1^7 2^8 matches the special fibre\n");
Print("y^lambda ~ y exactly for the 11 quadratic residues mod 23\n");
Print("the two order-23 branch classes are swapped exactly by ",
    "nonsquares\n");
Print("Riemann-Hurwitz sum ",ramificationSum,", genus ",genus,"\n");

QUIT_GAP(0);
