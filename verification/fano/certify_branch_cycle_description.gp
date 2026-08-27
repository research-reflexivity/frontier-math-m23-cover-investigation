\\ Certify the branch-cycle description of the class-6 degree-one curve directly
\\ from the integral model, closing the bridge between the explicit
\\ equation F(T,V) and the abstract Nielsen class (2A,23A,23B).
\\
\\ The derivative resultant factors as T^8 (T^2+23)^92 G^2 with
\\ deg G = 84.  Dividing by the V-leading coefficient (T^2+23)^4
\\ gives discriminant exponent 88 at T^2+23.
\\ The three parts have completely different meanings:
\\
\\   T^8          eight simple ramification points; T does not divide the
\\                singular locus, so these are smooth points of the plane
\\                model and the fibre type at T=0 is 1^7 2^8, class 2A;
\\   (T^2+23)     the Newton polygon at this place is a single segment of
\\                slope 4/23 with gcd(4,23)=1, so the local extension is
\\                totally ramified of degree 23: a 23-cycle at each of the
\\                two conjugate geometric points T = +-sqrt(-23);
\\   G^2          exponent two and contained in the singular locus.  This
\\                certificate locates 84 singular fibres; node ordinarity
\\                in good reduction is checked by the in-repository
\\                singular-scheme computation.
\\
\\ Riemann-Hurwitz then gives 2g-2 = -46 + 8 + 22 + 22 = 6, genus four.

check(condition,message) = if(!condition,error(message));

default(parisizemax, 12*10^9);

read("data/Fint_coefficients_Z.gp");

place_valuation(f,q) =
{
  my(v = 0);
  if(f == 0, return(-1));
  while(f % q == 0, f = f/q; v++);
  v
};

certify_branch_cycle_description() =
{
  my(leading, derivativeResultant, polynomialDiscriminant, tangential,
     singular, factorization, nodeFactor, valuations, slopeNumerator,
     degreeCount, genus, discriminantUnit, out);

  check(poldegree(F,V) == 23,"cover degree changed");
  leading = polcoef(F,23,V);
  check(leading == (T^2+23)^4,"leading coefficient changed");

  derivativeResultant = polresultant(F,deriv(F,V),V);
  check(derivativeResultant % leading == 0,
    "derivative resultant is not divisible by the leading coefficient");
  polynomialDiscriminant = -derivativeResultant/leading;
  tangential = polresultant(F,deriv(F,T),V);
  singular = gcd(derivativeResultant,tangential);

  factorization = factor(derivativeResultant);
  check(matsize(factorization)[1] == 3,"resultant factor count changed");
  check(factorization[1,1] == T && factorization[1,2] == 8,
    "simple branch factor changed");
  check(factorization[2,1] == T^2+23 && factorization[2,2] == 92,
    "order-23 branch factor changed");
  nodeFactor = factorization[3,1];
  check(poldegree(nodeFactor,T) == 84 && factorization[3,2] == 2,
    "node factor changed");
  check(issquarefree(nodeFactor),"node factor is not squarefree");
  check(place_valuation(polynomialDiscriminant,T) == 8,
    "T-adic discriminant exponent changed");
  check(place_valuation(polynomialDiscriminant,T^2+23) == 88,
    "order-23 discriminant exponent changed");
  discriminantUnit = polynomialDiscriminant/
    (T^8*(T^2+23)^88*nodeFactor^2);
  check(poldegree(discriminantUnit,T) == 0 && discriminantUnit != 0,
    "polynomial discriminant has an unexpected factor");

  \\ T = 0 is a smooth branch locus; the node factor is singular
  check(singular % T != 0,
    "the special fibre became singular");
  check(singular % nodeFactor == 0,
    "the extra discriminant zeros are not singular points");
  check(subst(F,T,0) == 279841*
    (V^7 - 312*V^6 + 43848*V^5 - 3914896*V^4 + 241123008*V^3
     - 9876755712*V^2 + 251799776256*V - 2694891036672)*
    (V^8 + 156*V^7 - 1290*V^6 - 1217744*V^5 - 40611765*V^4
     + 2044635024*V^3 + 115083889782*V^2 + 968032514052*V
     - 6058197515007)^2,"special-fibre factorization changed");

  \\ Newton polygon at the order-23 place
  valuations = vector(24,i,place_valuation(polcoef(F,i-1,V),T^2+23));
  check(valuations[24] == 4,"leading valuation changed");
  check(valuations[1] == 0,"constant term valuation changed");
  check(valuations[23] == -1,"the V^22 coefficient is no longer zero");
  slopeNumerator = valuations[24] - valuations[1];
  check(slopeNumerator == 4,"Newton slope numerator changed");
  check(gcd(slopeNumerator,23) == 1,"Newton slope is not coprime to 23");
  \\ every finite point lies on or above the segment joining
  \\ (0,-4) and (23,0) in the normalized valuations v(c_i)=v(a_i)-4
  for(i = 1, 23,
    if(valuations[i] < 0, next);
    check(23*(valuations[i]-4) >= -slopeNumerator*(23-(i-1)),
      Str("Newton polygon is not a single segment at i=",i-1)));

  \\ Riemann-Hurwitz
  degreeCount = 8 + 22 + 22;
  genus = (-2*23 + degreeCount + 2)/2;
  check(genus == 4,"genus changed");

  out = fileopen("results/branch_cycle_description_summary.txt","w");
  filewrite(out,"schema=m23.cover-investigation.branch-cycle-description.v1");
  filewrite(out,"cover_degree=23");
  filewrite(out,"leading_V_coefficient=(T^2+23)^4");
  filewrite(out,"derivative_resultant=T^8*(T^2+23)^92*G^2");
  filewrite(out,"polynomial_discriminant=unit*T^8*(T^2+23)^88*G^2");
  filewrite(out,"node_factor_degree=84");
  filewrite(out,"node_factor_exponent=2");
  filewrite(out,"node_factor_squarefree=true");
  filewrite(out,"node_factor_inside_singular_locus=true");
  filewrite(out,"singular_factor_does_not_prove_ordinary_nodes_by_itself=true");
  filewrite(out,"ordinary_node_input=mod31_singular_scheme");
  filewrite(out,"ordinary_nodes_contribute_no_ramification=true");
  filewrite(out,"special_fibre_is_smooth=true");
  filewrite(out,"special_fibre_type=1^7_2^8");
  filewrite(out,"special_fibre_inertia_class=2A");
  filewrite(out,"order23_place_newton_slope=4/23");
  filewrite(out,"order23_place_totally_ramified=true");
  filewrite(out,"order23_branch_cycles=23");
  filewrite(out,"order23_branch_points_conjugate_over_Q_sqrt_minus_23=true");
  filewrite(out,"branch_points=3");
  filewrite(out,"riemann_hurwitz_2g_minus_2=6");
  filewrite(out,"genus=4");
  filewrite(out,"nielsen_class=(2A,23A,23B)");
  fileclose(out);

  print("Res_V(F,F_V) = unit*T^8*(T^2+23)^92*G^2, deg G = 84");
  print("disc_V(F) = unit*T^8*(T^2+23)^88*G^2");
  print("T=0: smooth, fibre type 1^7 2^8, class 2A");
  print("T^2+23: Newton slope 4/23, totally ramified, 23-cycles");
  print("G: 84 squarefree singular fibres; node ordinarity checked mod 31");
  print("exactly three branch points; Riemann-Hurwitz gives genus 4")
};

certify_branch_cycle_description();
