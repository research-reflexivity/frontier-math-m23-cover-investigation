\\ Independent PARI/GP check of the specialization and modular factors.
read("../data/optimal_23_4_Z.gp");
read("../data/specialization_t2_Z.gp");

check(condition, message) = if(!condition, error(message));

verify_specialization() =
{
  my(f29, f11, q29, q2, q7, q14, maximal_orders);
  check(subst(P_t2, x, W) == subst(F23_4, T, 2),
        "stored polynomial is not P(2,x)");
  check(poldegree(P_t2, x) == 23, "specialization has wrong degree");
  check(content(P_t2) == 1, "specialization is not primitive");

  f29 = (Mod(1,29) / Mod(pollead(P_t2),29)) * P_t2;
  q29 = x^23 + 18*x^22 + 20*x^21 + 9*x^20 + 16*x^19 + 24*x^18
        + 27*x^17 + 5*x^16 + 28*x^15 + 20*x^14 + 14*x^13 + 3*x^12
        + 25*x^11 + 21*x^10 + 28*x^9 + x^8 + 9*x^7 + 17*x^6
        + 21*x^5 + 21*x^4 + 15*x^3 + 7*x^2 + 3*x + 5;
  check(f29 == Mod(1,29)*q29, "incorrect reduction modulo 29");
  check(polisirreducible(f29), "specialization is reducible modulo 29");

  f11 = (Mod(1,11) / Mod(pollead(P_t2),11)) * P_t2;
  q2 = x^2 + 8*x + 3;
  q7 = x^7 + 6*x^6 + 7*x^5 + 6*x^4 + 8*x^2 + 3*x + 1;
  q14 = x^14 + 8*x^13 + 6*x^12 + 8*x^11 + 6*x^10 + 6*x^9
        + 3*x^8 + 2*x^7 + 10*x^6 + 3*x^5 + x^4 + 8*x^3
        + 10*x^2 + 7*x + 6;
  check(f11 == Mod(1,11)*q2*q7*q14, "incorrect factorization modulo 11");
  check(polisirreducible(Mod(1,11)*q2), "quadratic factor is reducible");
  check(polisirreducible(Mod(1,11)*q7), "degree-7 factor is reducible");
  check(polisirreducible(Mod(1,11)*q14), "degree-14 factor is reducible");
  check(poldegree(gcd(f11, deriv(f11, x)), x) == 0,
        "reduction modulo 11 is not squarefree");

  maximal_orders = [443520,40320,40320,20160,7920,5760,253];
  check(sum(i=1,#maximal_orders,maximal_orders[i]%23 == 0) == 1,
        "unexpected maximal subgroup order divisible by 23");
  check(maximal_orders[7] == 23*11 && maximal_orders[7]%2 == 1,
        "23:11 order exclusion failed");

  print("PASS_SPECIALIZATION_T2_PARI_CERTIFICATE");
  print("degree=23 content=1 irreducible_mod_29=1");
  print("factor_degrees_mod_11=[2,7,14] frobenius_order=14");
};

verify_specialization();
