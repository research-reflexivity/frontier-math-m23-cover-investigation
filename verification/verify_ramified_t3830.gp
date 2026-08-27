\\ Exact branch-divisor and field-discriminant certificate at t=3830.
read("verify_progression_319.gp");
read("../data/Fint_coefficients_Z.gp");
y;

primitive_specialization(t0) =
{
  my(f = subst(F23_4,T,t0), c);
  c = abs(content(f));
  f/c;
};

monic_integral_model(f) =
{
  my(a = pollead(f,W));
  a^22 * subst(f,W,y/a);
};

verify_ramified_specialization() =
{
  my(t0 = 3830, f0, ft0, repeated, discriminant_factors, h62);
  my(f2, f3830, g2, g3830, disc2, disc3830);

  check(t0%319 == 2, "t=3830 is not in the M23 progression");
  check(t0%25 == 5 && valuation(t0,5) == 1,
        "t=3830 does not meet T=0 transversely at 5");

  \\ Exact discriminant factorization of the optimal plane equation.
  discriminant_factors = factor(poldisc(F23_4,W));
  check(matsize(discriminant_factors)[1] == 3,
        "unexpected number of discriminant factors");
  check(discriminant_factors[1,1] == T && discriminant_factors[1,2] == 8,
        "wrong T-factor in the discriminant");
  check(discriminant_factors[2,1] == T^2+23 && discriminant_factors[2,2] == 22,
        "wrong quadratic factor in the discriminant");
  h62 = discriminant_factors[3,1];
  check(poldegree(h62,T) == 62 && discriminant_factors[3,2] == 2,
        "wrong residual square factor in the discriminant");
  check(poldegree(gcd(h62,T*(T^2+23)),T) == 0,
        "residual discriminant factor meets the branch divisor");

  \\ Eight geometrically distinct transverse double roots at T=0 mod 5.
  f0 = Mod(1,5)*subst(Fhat,T,0);
  ft0 = Mod(1,5)*subst(deriv(Fhat,T),T,0);
  repeated = gcd(f0,deriv(f0,V));
  check(poldegree(repeated,V) == 8, "wrong repeated degree at T=0 mod 5");
  check(poldegree(gcd(repeated,deriv(repeated,V)),V) == 0,
        "repeated part is not squarefree mod 5");
  check(poldegree(gcd(repeated,ft0),V) == 0,
        "branch deformation is not transverse mod 5");

  f2 = primitive_specialization(2);
  f3830 = primitive_specialization(t0);
  check(polisirreducible(Mod(1,29)*f3830),
        "t=3830 specialization is reducible modulo 29");
  g2 = monic_integral_model(f2);
  g3830 = monic_integral_model(f3830);
  check(pollead(g2,y) == 1 && denominator(content(g2)) == 1,
        "failed to build the t=2 monic integral model");
  check(pollead(g3830,y) == 1 && denominator(content(g3830)) == 1,
        "failed to build the t=3830 monic integral model");

  disc2 = nfdisc(g2);
  disc3830 = nfdisc(g3830);
  check(valuation(disc2,5) == 0,
        "the t=2 root field should be unramified at 5");
  check(valuation(disc3830,5) == 8,
        "the t=3830 root field should have discriminant valuation 8 at 5");

  print("PASS_RAMIFIED_T3830_PARI_CERTIFICATE");
  print("disc_W(P)=unit*T^8*(T^2+23)^22*H62(T)^2");
  print("t=3830 congruences=2_mod_319,5_mod_25 galois_group=M23");
  print("root_field_discriminant_v5: t2=0 t3830=8");
};

verify_ramified_specialization();
