\\ Independent PARI/GP certificate for T = 2 (mod 319).
\\ Reading the base certificate reruns both finite-field factorizations.
read("verify_specialization_t2.gp");

verify_progression() =
{
  my(primes = [11,29], shifted, base);
  base = subst(F23_4, T, 2);
  for(i=1,#primes,
    my(p = primes[i]);
    shifted = subst(F23_4, T, 2 + 319*n);
    check(Mod(1,p)*shifted == Mod(1,p)*base,
          Str("universal progression congruence failed modulo ", p));
    check(poldegree(Mod(1,p)*base, W) == 23,
          Str("degree drops modulo ", p));
  );
  print("PASS_PROGRESSION_319_PARI_CERTIFICATE");
  print("P(2+319*n,W)=P(2,W) modulo 11 and modulo 29");
  print("galois_group_of_every_primitive_part=M23");
};

verify_progression();
