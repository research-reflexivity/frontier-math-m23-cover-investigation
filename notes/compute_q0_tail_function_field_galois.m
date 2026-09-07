/* Standalone group check for the primitive conductor-seven tail.
   No characteristic-zero M23 cover or Hurwitz polynomial is an input. */
F23 := GF(23);
K<t> := FunctionField(F23);
R<X> := PolynomialRing(K);
f := X^23 + 7*X^9 + 18*X^2 - t;
assert Degree(f) eq 23;
assert Derivative(f) ne 0;
assert IsIrreducible(f);
SetVerbose("GaloisGroup", 1);
time G, roots, data := GaloisGroup(f : ProofEffort := 10, ShortOK := false);
print "group_order=", #G;
print "group_degree=", Degree(G);
print "group_generators=", Generators(G);
print "splitting_prime=", data`Prime;
print "transitivity=", Transitivity(G);
print "simple=", IsSimple(G);
assert Degree(G) eq 23;
assert #G eq 10200960;
assert Transitivity(G) eq 4;
assert IsSimple(G);
print "PASS_Q0_TAIL_FUNCTION_FIELD_M23";
