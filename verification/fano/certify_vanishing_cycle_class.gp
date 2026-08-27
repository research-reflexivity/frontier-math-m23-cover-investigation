\\ Certify the arithmetic side of the vanishing-cycle inertia class.  At
\\ each double root v of the ramified octic the colliding sheets are
\\ V = v +- sqrt(e(v) T) + O(T) with e(v) = -2 F_T(0,v)/F_VV(0,v); modulo
\\ squares e(v) is represented by c(V) = -F_T(0,V) H7(V) mod R8.  The
\\ certificate proves: c is not a square in the octic field even after
\\ any twist by squarefree divisors of 2*3*23; at every prime below
\\ 300000 splitting completely in the septic-octic closure the eight
\\ values c(v_i) have equal quadratic character, with both characters
\\ occurring -- the -1 primes prove the class is a nonsquare in the
\\ degree-1344 residue field, so the local decomposition group is the
\\ full perfect double cover; and exactly two of the four real nodes
\\ have positive class, matching the branch reflection.

check(condition,message) = if(!condition,error(message));

read("data/Fint_coefficients_Z.gp");

yv = varlower("y", 'V);

certify_vanishing_cycle_class() =
{
  my(H7, R8, A, c, cLow, nf, twists, badPrimes, identityPrimes,
     characters, realSigns, minusCount, plusCount, out);

  H7 = V^7 - 312*V^6 + 43848*V^5 - 3914896*V^4
     + 241123008*V^3 - 9876755712*V^2 + 251799776256*V
     - 2694891036672;
  R8 = V^8 + 156*V^7 - 1290*V^6 - 1217744*V^5 - 40611765*V^4
     + 2044635024*V^3 + 115083889782*V^2 + 968032514052*V
     - 6058197515007;

  check(subst(F,'T,0) == 279841*H7*R8^2,
    "special-fibre factorization changed");
  A = polcoef(F, 1, 'T);
  check(poldegree(A,'V) == 17,"T-linear coefficient changed");
  check(polresultant(H7,R8) != 0 && polresultant(R8,A) != 0
    && poldisc(R8) != 0,"vanishing cycle data degenerate");

  \\ the square-class representative: F_VV(0,v) = 2*23^4*H7(v)*R8'(v)^2
  \\ at octic roots, so e(v) = -A(v)*H7(v) modulo visible squares.
  c = lift(Mod(-A*H7, R8));
  check(issquare(content(c)),"class content is not a square");
  c = c/content(c);
  check(content(c) == 1,"class normalization failed");

  \\ real place: two positive and two negative nodes
  realSigns = vecsort(apply(v -> sign(subst(c,'V,v)),
    Vec(polrootsreal(R8))));
  check(realSigns == [-1,-1,1,1],"real node signs changed");

  \\ octic-level nontriviality, robust under rational square-class twists
  nf = nfinit(subst(R8,'V,yv));
  cLow = subst(c,'V,yv);
  twists = [1,-1,2,-2,3,-3,6,-6,23,-23,46,-46,69,-69,138,-138];
  for(i = 1, 16,
    check(#nffactor(nf, 'V^2 - twists[i]*Mod(cLow,nf.pol))~ == 1,
      Str("class became a square for twist ",twists[i])));

  \\ identity primes of the septic-octic closure: coherence of the
  \\ eight characters and both Frobenius classes of the double cover
  badPrimes = 2*3*23*poldisc(H7)*poldisc(R8)
    *polresultant(H7,R8)*polresultant(R8,A);
  identityPrimes = [];
  characters = [];
  forprime(p = 5, 300000,
    if(badPrimes % p == 0, next);
    if(#polrootsmod(H7,p) < 7, next);
    my(r8roots = polrootsmod(R8,p));
    if(#r8roots < 8, next);
    my(ch = Set(apply(a -> kronecker(lift(subst(Mod(c,p),'V,a)), p),
      Vec(r8roots))));
    check(#ch == 1,Str("incoherent characters at ",p));
    identityPrimes = concat(identityPrimes,[p]);
    characters = concat(characters,[ch[1]]);
  );
  check(#identityPrimes == 14,"identity prime count changed");
  check(identityPrimes[1] == 13049 && characters[1] == -1,
    "first identity prime changed");
  minusCount = #select(t -> t == -1, characters);
  plusCount = #select(t -> t == 1, characters);
  check(minusCount == 6 && plusCount == 8,
    "double cover Frobenius split changed");

  out = fileopen(
    "results/vanishing_cycle_class_arithmetic_summary.txt","w");
  filewrite(out,"schema=m23.cover-investigation.vanishing-cycle-class-arithmetic.v1");
  filewrite(out,"class_representative=-F_T(0,V)*H7(V)_mod_R8");
  filewrite(out,"class_content_square_root_removed=true");
  filewrite(out,"real_node_signs=[+,+,-,-]_sorted");
  filewrite(out,"real_branch_nodes=2_conjugate_branch_nodes=2");
  filewrite(out,"octic_twisted_nonsquare_tests=16");
  filewrite(out,"octic_class_nontrivial_mod_rational_squares=true");
  filewrite(out,"identity_prime_bound=300000");
  filewrite(out,Str("identity_primes=",identityPrimes));
  filewrite(out,Str("double_cover_characters=",characters));
  filewrite(out,"characters_coherent_at_every_identity_prime=true");
  filewrite(out,"minus_primes=6_plus_primes=8");
  filewrite(out,"class_nonsquare_in_residue_field=true");
  filewrite(out,"local_decomposition_group_is_full_double_cover=true");
  fileclose(out);

  print("class c = -F_T(0,V) H7(V) mod R8, content-free");
  print("real nodes: signs two + and two -, matching the reflection");
  print("sixteen twisted quadratics irreducible over the octic");
  print("14 identity primes, all coherent; characters 6 minus, 8 plus");
  print("the vanishing-cycle class realizes the central inertia")
};

certify_vanishing_cycle_class();
