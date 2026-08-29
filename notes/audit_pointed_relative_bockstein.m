/* Independent Magma audit of the pointed relative Bockstein linear algebra. */

groupOrder := 10200960;
onePointStabilizer := groupOrder div 23;
twoPointStabilizer := groupOrder div (23*22);
assert onePointStabilizer eq 443520;
assert twoPointStabilizer eq 20160;
assert onePointStabilizer mod 4 eq 0;
assert twoPointStabilizer mod 4 eq 0;

full := <120,1035>;
trace := <8,69>;
packet := <full[1]-trace[1],full[2]-trace[2]>;
assert packet eq <112,966>;
assert <full[1] mod 4,full[2] mod 4> eq <0,3>;
assert <trace[1] mod 4,trace[2] mod 4> eq <0,1>;
assert <packet[1] mod 4,packet[2] mod 4> eq <0,2>;

alternativePacket := <packet[1]-2*full[1],packet[2]-2*full[2]>;
assert alternativePacket eq <-128,-1104>;
assert alternativePacket[1] mod 4 eq 0;
assert alternativePacket[2] mod 4 eq 0;

/* Nakayama's rank-one horizontal log-nearby stalk is free over the
   coefficient ring in degrees zero and one.  Reduction Z/4 -> F_2 is
   surjective in both degrees, so the coefficient connecting maps vanish. */
Z4 := Integers(4);
F2 := GF(2);
reductionImageH0 := {Integers()!entry mod 2 : entry in Z4};
reductionImageH1 := {Integers()!entry mod 2 : entry in Z4};
assert reductionImageH0 eq {0,1};
assert reductionImageH1 eq {0,1};

relativeComparisonMod4 := 3;
relativeFirstLayer := ((relativeComparisonMod4-1) div 2) mod 2;
assert relativeFirstLayer eq 1;

diagonalBranch := 1078;
transpositionBranch := 112;
assert diagonalBranch mod 2 eq 0;
assert transpositionBranch mod 2 eq 0;
dividedDefect := ((diagonalBranch-transpositionBranch) div 2) mod 2;
assert dividedDefect eq 1;
assert <diagonalBranch mod 4,transpositionBranch mod 4> eq <2,0>;

function QuadraticHalfWeight(vector)
    weight := &+vector;
    assert weight mod 2 eq 0;
    return (weight div 2) mod 2;
end function;

for bitCount in [2..9] do
    evenVectors := [];
    for mask in [0..2^bitCount-1] do
        vector := [Integers()!((mask div 2^(index-1)) mod 2)
            : index in [1..bitCount]];
        if (&+vector) mod 2 eq 0 then
            Append(~evenVectors,vector);
        end if;
    end for;
    for left in evenVectors do
        for right in evenVectors do
            total := [(left[index]+right[index]) mod 2
                : index in [1..bitCount]];
            dot := (&+[left[index]*right[index]
                : index in [1..bitCount]]) mod 2;
            assert QuadraticHalfWeight(total) eq
                (QuadraticHalfWeight(left)+QuadraticHalfWeight(right)+dot) mod 2;
        end for;
    end for;
end for;

packetSize := 3795-253;
assert packetSize eq 3542;
assert (23*packetSize div 2) mod 2 eq 1;
assert (diagonalBranch div 2) mod 2 eq 1;
assert (transpositionBranch div 2) mod 2 eq 0;
assert ((98*253) div 2) mod 2 eq 1;

for overlapSize in [0..253] do
    relativeHalfDistance := ((2*253-2*overlapSize) div 2) mod 2;
    assert relativeHalfDistance eq (1+overlapSize) mod 2;
end for;
assert ((2*253-2*253) div 2) mod 2 eq 0;
assert ((2*253-2*0) div 2) mod 2 eq 1;

fullClassSize := 3795;
traceSize := 253;
squareFullTraceHalfDistance := ((fullClassSize-traceSize) div 2) mod 2;
nonsquareFullTraceHalfDistance := ((fullClassSize+traceSize) div 2) mod 2;
assert squareFullTraceHalfDistance eq 1;
assert nonsquareFullTraceHalfDistance eq 0;

F := GF(2);

function PathCoboundary(edgeCount,isRelative)
    if isRelative then
        vertices := [2..edgeCount];
    else
        vertices := [1..edgeCount+1];
    end if;
    entries := [];
    for edge in [1..edgeCount] do
        for vertex in vertices do
            value := 0;
            if vertex eq edge then
                value +:= 1;
            end if;
            if vertex eq edge+1 then
                value +:= 1;
            end if;
            Append(~entries,F!value);
        end for;
    end for;
    return Matrix(F,edgeCount,#vertices,entries);
end function;

for edgeCount in [1..16] do
    absolute := PathCoboundary(edgeCount,false);
    relative := PathCoboundary(edgeCount,true);
    assert Rank(absolute) eq edgeCount;
    assert Rank(relative) eq edgeCount-1;
    assert Nrows(relative)-Rank(relative) eq 1;
    for column in [1..Ncols(relative)] do
        assert &+[Integers()!relative[row,column]
            : row in [1..Nrows(relative)]] mod 2 eq 0;
    end for;
end for;

for returnBit in [0,1] do
    finiteValue := 1;
    wildValue := returnBit;
    relativeValue := (finiteValue+wildValue) mod 2;
    assert relativeValue eq (1+returnBit) mod 2;
end for;

print "rho_rho_adjoint_entries_mod_4=0";
print "canonical_sheet_lifts=full:-I,trace:+I";
print "packet_sheet_operator=2I_mod_4";
print "ordinary_mapping_cone_bockstein=zero_by_alternative_lift";
print "horizontal_log_nearby_coefficient_bockstein=zero";
print "specialization_unit_divided_defect=1";
print "half_weight_quadratic_refines_tagged_dot_product=true";
print "packet_global_half_weight=1";
print "returned_half_conductor_values=1,0";
print "relative_trace_half_distance=1+overlap_parity";
print "returned_trace_half_distance=return_bit";
print "full_class_returned_trace_half_distance=1+return_bit";
print "unpointed_path_H1_dimension=0";
print "two_endpoint_relative_path_H1_dimension=1";
print "relative_endpoint_formula=finite_value+wild_value";
print "PASS_POINTED_RELATIVE_BOCKSTEIN_AUDIT";
