\\ Certify the ramification of the central vanishing-cycle layer.
\\
\\ Let E be the octic field cut out by R8 and c the square-class
\\ representative of the vanishing cycles.  The relative discriminant of
\\ E(sqrt(c))/E has norm 2^8*3^2*23^2.  Thus the central quadratic layer is
\\ unramified away from {2,3,23}, but it is ramified at every one of those
\\ rational primes.  In particular both octic primes above 23 ramify, so
\\ the full degree-2688 special-fibre field has central inertia at 23.

check(condition,message) = if(!condition,error(message));

default(parisizemax, 12*10^9);

certify_central_quadratic_ramification() =
{
  my(H7, R8, A, c, yv, nf, cLow, relativePolynomial,
     relativeDiscriminant, relativeFactorization, rationalSupport,
     primesAt23, normDiscriminant, out);

  read("data/Fint_coefficients_Z.gp");

  H7 = V^7 - 312*V^6 + 43848*V^5 - 3914896*V^4
     + 241123008*V^3 - 9876755712*V^2 + 251799776256*V
     - 2694891036672;
  R8 = V^8 + 156*V^7 - 1290*V^6 - 1217744*V^5 - 40611765*V^4
     + 2044635024*V^3 + 115083889782*V^2 + 968032514052*V
     - 6058197515007;

  check(subst(F,T,0) == 23^4*H7*R8^2,
    "special-fibre factorization changed");
  A = polcoef(F,1,T);
  c = lift(Mod(-A*H7,R8));

  yv = varlower("y",'V);
  nf = nfinit(subst(R8,'V,yv));
  cLow = subst(c,'V,yv);
  relativePolynomial = 'V^2 - Mod(cLow,nf.pol);
  relativeDiscriminant = rnfdisc(nf,relativePolynomial)[1];
  relativeFactorization = idealfactor(nf,relativeDiscriminant);
  normDiscriminant = idealnorm(nf,relativeDiscriminant);
  rationalSupport = vecsort(Set(vector(
    matsize(relativeFactorization)[1],
    i, relativeFactorization[i,1][1])));
  primesAt23 = select(
    i -> relativeFactorization[i,1][1] == 23,
    [1..matsize(relativeFactorization)[1]]);

  check(normDiscriminant == 2^8*3^2*23^2,
    "relative discriminant norm changed");
  check(rationalSupport == [2,3,23],
    "central layer ramifies outside the branch primes");
  check(#primesAt23 == 2
      && vector(#primesAt23,i,relativeFactorization[primesAt23[i],2])
          == [1,1],
    "central ramification above 23 changed");

  out = fileopen("results/central_quadratic_ramification_summary.txt","w");
  filewrite(out,
    "schema=m23.cover-investigation.central-quadratic-ramification.v1");
  filewrite(out,"relative_extension=E_sqrt_c_over_E_octic");
  filewrite(out,"relative_discriminant_norm=2^8*3^2*23^2");
  filewrite(out,"relative_discriminant_support=[2,3,23]");
  filewrite(out,"central_layer_unramified_outside_2_3_23=true");
  filewrite(out,"central_layer_ramified_at_2=true");
  filewrite(out,"central_layer_ramified_at_3=true");
  filewrite(out,"central_layer_ramified_at_23=true");
  filewrite(out,"octic_primes_above_23_ramified_in_c_layer=2");
  filewrite(out,"relative_discriminant_exponents_above_23=[1,1]");
  filewrite(out,"full_inertia_at_23_contains_central_involution=true");
  fileclose(out);

  print("central relative discriminant norm = 2^8*3^2*23^2");
  print("central layer ramifies at 2, 3 and 23, and nowhere else");
  print("both octic primes above 23 ramify in the vanishing-cycle layer")
};

certify_central_quadratic_ramification();
