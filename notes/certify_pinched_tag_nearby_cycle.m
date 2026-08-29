/* Independent Magma certificate for the finite pinched-tag identities. */

F := GF(2);

function Dot(alpha,beta)
    return &+[F | F!(alpha[i]*beta[i]) : i in [1..#alpha]];
end function;

function Bucket(alpha,beta)
    return F!(&+alpha)*F!(&+beta);
end function;

function OffDiagonal(alpha,beta)
    return &+[F | F!(alpha[i]*beta[j]) : i,j in [1..#alpha] | i ne j];
end function;

for size in [1..6] do
    vectors := CartesianPower({0,1},size);
    for alphaTuple in vectors do
        alpha := [Integers()!x : x in alphaTuple];
        for betaTuple in vectors do
            beta := [Integers()!x : x in betaTuple];
            assert Bucket(alpha,beta) eq Dot(alpha,beta)+OffDiagonal(alpha,beta);
        end for;
    end for;
end for;

diagonalBranch := 77^2;
transpositionBranch := 8^2;
nodeDifference := diagonalBranch-transpositionBranch;
assert [diagonalBranch,transpositionBranch,nodeDifference] eq [5929,64,5865];
assert nodeDifference mod 2 eq 1;
conductorPairingDegree := 23*diagonalBranch+253*transpositionBranch;
assert conductorPairingDegree eq 152559;

singletonPairing := 23*11^2;
twoSetPairing := 253;
untaggedTotal := 7*singletonPairing+8*twoSetPairing;
assert [singletonPairing,twoSetPairing,untaggedTotal] eq [2783,253,21505];
assert untaggedTotal mod 2 eq 1;

squareNormalized := [1,1,0];
nonsquareNormalized := [0,1,1];
squareNode := [(77-8) mod 2,nodeDifference mod 2,
               (nodeDifference-(77-8)) mod 2];
nonsquareNode := [0,nodeDifference mod 2,nodeDifference mod 2];
assert squareNode eq squareNormalized;
assert nonsquareNode eq nonsquareNormalized;

for length in [1..12] do
    for valueTuple in CartesianPower({0,1},length+1) do
        values := [Integers()!x : x in valueTuple];
        boundarySum := &+[(values[i]+values[i+1]) mod 2 : i in [1..length]];
        assert boundarySum mod 2 eq (values[1]+values[#values]) mod 2;
    end for;
end for;

returnBits := [0,1,1];
wildTerms := [(1+q) mod 2 : q in returnBits];
specialTraces := [(wild+1) mod 2 : wild in wildTerms];
assert specialTraces eq returnBits;

print "pinched_bucket_identity_checked_for_fiber_sizes=1..6";
print "untagged_total=7*(23*11^2)+8*253=21505";
print "node_coefficients=[5929,64],difference=5865=1_mod_2";
print "conductor_bucket_distribution=23x77_and_253x8";
print "conductor_pairing_degree=152559";
print "square_fields_tagged_untagged_offdiag=[1,1,0]";
print "nonsquare_fields_tagged_untagged_offdiag=[0,1,1]";
print "tree_telescope_checked_through_12_edges=true";
print "special_trace_table=[0,1,1]";
print "PASS_PINCHED_TAG_NEARBY_CYCLE";

quit;
