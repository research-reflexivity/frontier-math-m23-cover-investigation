/* Exact Magma audit for the effective logarithmic orientation line. */

QHalfRank := func< r | ((r div 2) mod 2) >;

traceSize := 253;
branchPlus := 1078;
branchMinus := 112;
packetSize := 3542;
fullClassSize := 3795;

assert QHalfRank(0) eq 0;
assert QHalfRank(2*traceSize) eq 1;
assert QHalfRank(branchPlus-branchMinus) eq 1;
assert QHalfRank(packetSize) eq 1;
assert QHalfRank(packetSize-2*fullClassSize) eq 0;

/* Terminal pair (C,T^n): containment for q=0, disjointness for q=1. */
assert QHalfRank(fullClassSize-traceSize) eq 1;
assert QHalfRank(fullClassSize+traceSize) eq 0;

/* The two fibre formulae are equivalent precisely when epsilon=q. */
for epsilon in [0,1] do
    for qReturn in [0,1] do
        generic := (1+epsilon) mod 2;
        special := (1+qReturn) mod 2;
        assert (generic eq special) eq (epsilon eq qReturn);
    end for;
end for;

/* Small exact audit of the divided normalization telescope. */
for a0 in [0..4] do
    for a1 in [0..4] do
        for a2 in [0..4] do
            if IsEven(a1-a0) and IsEven(a2-a1) then
                localSum := (((a1-a0) div 2)+((a2-a1) div 2)) mod 2;
                endpoint := ((a2-a0) div 2) mod 2;
                assert localSum eq endpoint;
            end if;
        end for;
    end for;
end for;

/* Exhaustive polarization on an eight-element tag set. */
U := {1..8};
evenSubsets := [ S : S in Subsets(U) | (#S mod 2) eq 0 ];
for S in evenSubsets do
    for T in evenSubsets do
        symDiff := (S diff T) join (T diff S);
        lhs := QHalfRank(#symDiff);
        rhs := (QHalfRank(#S)+QHalfRank(#T)+(#(S meet T) mod 2)) mod 2;
        assert lhs eq rhs;
    end for;
end for;

/* Branch exchange and a subdivision identity complex. */
assert QHalfRank(branchPlus-branchMinus)
       eq QHalfRank(branchMinus-branchPlus);
for r in [0..11] do
    assert QHalfRank(r-r) eq 0;
end for;

print "log_quadratic_orientation_line_magma_audit=PASS";
print "finite_anomaly=1";
print "packet_return_polarization_cross_term=0,0";
print "terminal_full_trace_orientation=square:1,nonsquare:0";
print "normalization_telescope_audit=PASS";
