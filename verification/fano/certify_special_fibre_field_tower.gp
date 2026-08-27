\\ Certify the arithmetic side of the field tower obtained from the special fibre.  Every good
\\ prime below 20000 has (octic,septic) Frobenius factorization pair in
\\ the ten-class 2^3:PSL(3,2) dictionary computed by the GAP certificate,
\\ the two Gassmann septics factor identically everywhere, all ten
\\ classes occur, and the witness prime 599 splits both septics
\\ completely while the octic factors into four quadratics -- proving
\\ that the splitting field of the ramified octic strictly exceeds the
\\ septic closure.  With the certified irreducibility of the line module
\\ this forces the full affine group.

check(condition,message) = if(!condition,error(message));

read("results/special_fibre_frobenius_dictionary.gp");

shape(P,p) = vecsort(apply(poldegree, Vec(factormod(P,p)[,1]~)));

certify_special_fibre_field_tower() =
{
  my(H7, L7, R8, badDisc, counts, witness, identityPrime, total, out);

  H7 = V^7 - 312*V^6 + 43848*V^5 - 3914896*V^4
     + 241123008*V^3 - 9876755712*V^2 + 251799776256*V
     - 2694891036672;
  L7 = V^7 - 936*V^6 + 379728*V^5 - 86588640*V^4
     + 11922182400*V^3 - 982633068288*V^2 + 44572755885824*V
     - 859311381083136;
  R8 = V^8 + 156*V^7 - 1290*V^6 - 1217744*V^5 - 40611765*V^4
     + 2044635024*V^3 + 115083889782*V^2 + 968032514052*V
     - 6058197515007;

  check(polisirreducible(R8),"ramified octic is reducible");
  check(#allowedpairs == 10,"dictionary size changed");

  badDisc = poldisc(H7)*poldisc(L7)*poldisc(R8);
  counts = vector(10);
  witness = 0;
  identityPrime = 0;
  forprime(p = 5, 20000,
    if(badDisc % p == 0, next);
    my(s7 = shape(H7,p), s7L = shape(L7,p), s8 = shape(R8,p), idx = 0);
    check(s7 == s7L,Str("Gassmann violation at ",p));
    for(i = 1, 10, if(allowedpairs[i] == [s8,s7], idx = i));
    check(idx > 0,Str("dictionary violation at ",p));
    counts[idx]++;
    if(idx == 6 && witness == 0, witness = p);
    if(idx == 1 && identityPrime == 0, identityPrime = p);
  );
  total = vecsum(counts);

  \\ the witness prime: septics split completely, octic four quadratics
  check(witness == 599,"first witness prime changed");
  check(shape(H7,599) == [1,1,1,1,1,1,1],"witness septic shape changed");
  check(shape(L7,599) == [1,1,1,1,1,1,1],"witness line shape changed");
  check(shape(R8,599) == [2,2,2,2],"witness octic shape changed");
  check(vecmin(counts) >= 1,"some dictionary class never occurs");

  out = fileopen("results/special_fibre_field_tower_arithmetic_summary.txt","w");
  filewrite(out,"schema=m23.cover-investigation.special-fibre-field-tower-arithmetic.v1");
  filewrite(out,"prime_bound=20000");
  filewrite(out,Str("good_primes_scanned=",total));
  filewrite(out,"gassmann_violations=0");
  filewrite(out,"dictionary_violations=0");
  filewrite(out,Str("class_counts=",counts));
  filewrite(out,"all_ten_classes_realized=true");
  filewrite(out,Str("witness_prime=",witness));
  filewrite(out,"witness_pattern=septics_split_octic_four_quadratics");
  filewrite(out,Str("identity_prime=",identityPrime));
  filewrite(out,"octic_splitting_field_strictly_contains_septic_closure=true");
  filewrite(out,"with_irreducible_line_module_gives_full_affine_group=true");
  fileclose(out);

  print("scanned ",total," good primes: zero violations");
  print("witness prime 599: septics 1^7, octic 2^4");
  print("all ten Frobenius classes realized; identity prime ",
    identityPrime);
  print("the octic and septic realize the full affine Galois group")
};

certify_special_fibre_field_tower();
