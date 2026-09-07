/* Exact Magma computation of the rational characteristic-23 E8 tail group. */

F23 := GF(23);
K<z> := FunctionField(F23);
R<X> := PolynomialRing(K);
f := X^23 - X^13 + 5*X^8 - 2*X^3 + 6*z^3;

assert Degree(f) eq 23;
assert Derivative(f) ne 0;
assert IsIrreducible(f);

SetVerbose("GaloisGroup", 1);
time G, roots, data := GaloisGroup(f : ProofEffort := 10, ShortOK := false);

print "group_order=", #G;
print "group_degree=", Degree(G);
print "group_generators=", Generators(G);
print "splitting_prime=", data`Prime;
print "transitive=", IsTransitive(G);
print "primitive=", IsPrimitive(G);
print "transitivity=", Transitivity(G);
print "simple=", IsSimple(G);

assert Degree(G) eq 23;
assert #G eq 10200960;
assert Transitivity(G) eq 4;
assert IsPrimitive(G);
assert IsSimple(G);
print "PASS rational E8 tail function-field Galois group is M23";
