# Certify the two-endpoint raw trace coefficient at the tame 2A node.
# The calculation is universal over the fifteen finite-inertia/wild-frame
# colors and also covers the adjacent-Mathieu nonsquare return.

SizeScreen([100000,100000]);;

gA := (1,2)(3,4)(7,8)(9,10)(13,14)(15,16)(19,20)(21,22);;
gB := (1,16,11,3)(2,9,21,12)(4,5,8,23)(6,22,14,18)
      (13,20)(15,17);;
G := Group(gA,gB);;
S := SymmetricGroup(23);;
y := (1,2,11,10,16,9,6,3,23,19,20,14,21,17,4,8,22,5,18,
      15,13,7,12);;
x := (3,21)(4,16)(9,22)(10,15)(11,20)(12,19)(13,18)(14,17);;
D := Normalizer(G,Group(y));;
C := Centralizer(G,x);;

IntegralSheetMatrix := function(packet)
    local matrix,w,point;
    matrix := NullMat(23,23,Integers);
    for w in packet do
        for point in [1..23] do
            matrix[point][point^w] := matrix[point][point^w]+1;
        od;
    od;
    return matrix;
end;;

ConstantOffDiagonalMatrix := function(diagonal,offDiagonal)
    local matrix,row,column;
    matrix := NullMat(23,23,Integers);
    for row in [1..23] do
        for column in [1..23] do
            if row=column then
                matrix[row][column] := diagonal;
            else
                matrix[row][column] := offDiagonal;
            fi;
        od;
    od;
    return matrix;
end;;

HadamardProduct := function(first,second)
    local product,row,column;
    product := NullMat(23,23,Integers);
    for row in [1..23] do
        for column in [1..23] do
            product[row][column] :=
                first[row][column]*second[row][column];
        od;
    od;
    return product;
end;;

expectedTrace := ConstantOffDiagonalMatrix(77,8);;
expectedProduct := ConstantOffDiagonalMatrix(5929,64);;
endpointDoubleCosets := DoubleCosets(G,D,C);;
traceMatrices := [];;
for color in endpointDoubleCosets do
    endpoint := D^Representative(color);
    trace := Set(Orbit(endpoint,x,OnPoints));
    Add(traceMatrices,IntegralSheetMatrix(trace));
od;

productMatrices := [];;
for first in traceMatrices do
    for second in traceMatrices do
        Add(productMatrices,HadamardProduct(first,second));
    od;
od;

# A nonsquare normalizer return carries the trace into the adjacent Mathieu
# copy.  Its ambient 23-sheet operator is unchanged.
N := Normalizer(S,D);;
n := Difference(AsList(N),AsList(D))[1];;
baseTrace := Set(Orbit(D,x,OnPoints));;
returnedTrace := Set(List(baseTrace,w -> w^n));;
returnedMatrix := IntegralSheetMatrix(returnedTrace);;
returnedProduct := HadamardProduct(expectedTrace,returnedMatrix);;

branchOrbits := Orbits(Group(x),[1..23]);;
fixedPoints := Filtered([1..23],point -> point^x=point);;
movedBlocks := Filtered(branchOrbits,orbit -> Length(orbit)=2);;

if Size(G)<>10200960 or Size(D)<>253 or Size(N)<>506
   or Length(endpointDoubleCosets)<>15
   or not ForAll(traceMatrices,matrix -> matrix=expectedTrace)
   or not ForAll(productMatrices,matrix -> matrix=expectedProduct)
   or returnedMatrix<>expectedTrace
   or returnedProduct<>expectedProduct
   or Length(fixedPoints)<>7 or Length(movedBlocks)<>8
   or expectedProduct[1][1]<>77^2
   or expectedProduct[1][2]<>8^2
   or (77^2-8^2) mod 2<>1
   or not ForAll(fixedPoints,point ->
       expectedProduct[point][point] mod 2=1)
   or not ForAll(movedBlocks,block ->
       expectedProduct[block[1]][block[1]] mod 2=1
       and expectedProduct[block[1]][block[2]] mod 2=0) then
    Error("raw lifted-trace node boundary audit changed");
fi;

Print("wild_endpoint_color_count=",Length(endpointDoubleCosets),"\n");
Print("endpoint_pair_count=",Length(productMatrices),"\n");
Print("trace_sheet_operator=8*J+69*I\n");
Print("two_endpoint_hadamard_operator=64*J+5865*I\n");
Print("diagonal_branch_coefficient=77^2=5929\n");
Print("transposition_branch_coefficient=8^2=64\n");
Print("node_branch_difference=5865=1_mod_2\n");
Print("ordinary_residual=fixed_heptad\n");
Print("tame_node_packet=full_affine_eight_set\n");
Print("adjacent_Mathieu_return_has_same_operator=true\n");
Print("PASS_RAW_LIFTED_TRACE_NODE_BOUNDARY\n");

QUIT_GAP(0);
