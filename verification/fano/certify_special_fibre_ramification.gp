\\ Certify the arithmetic side of the special-fibre ramification theorem:
\\ all three constituent fields are unramified outside {2,3,23}, and the local
\\ behaviour at each of those primes matches exactly the group-theoretic
\\ inertia predicted by the companion GAP certificate.
\\
\\   p=23  totally tamely ramified with e=7 in both septics, and 1+7 in
\\         the octic.  These non-Galois factorization data identify the
\\         Singer inertia orbits, but not the decomposition group of the
\\         Galois closure; the companion GAP certificate gives D=7:3.
\\   p=3   wildly ramified with e=3, orbits 1,3,3 on sheets and 1,1,3,3
\\         on double points: an order-three collineation.
\\   p=2   tame of order 3 on the sheets but wild with e=4 on the double
\\         points: the wild inertia lies in the line module.

check(condition,message) = if(!condition,error(message));

ramification_data(P,p) =
{
  my(nf = nfinit(P));
  vecsort(apply(pr -> [pr.e,pr.f], idealprimedec(nf,p)))
};

certify_special_fibre_ramification() =
{
  my(H7, L7, R8, discs, out);

  H7 = V^7 - 312*V^6 + 43848*V^5 - 3914896*V^4
     + 241123008*V^3 - 9876755712*V^2 + 251799776256*V
     - 2694891036672;
  L7 = V^7 - 936*V^6 + 379728*V^5 - 86588640*V^4
     + 11922182400*V^3 - 982633068288*V^2 + 44572755885824*V
     - 859311381083136;
  R8 = V^8 + 156*V^7 - 1290*V^6 - 1217744*V^5 - 40611765*V^4
     + 2044635024*V^3 + 115083889782*V^2 + 968032514052*V
     - 6058197515007;

  \\ field discriminants: supported on the branch primes only
  discs = [nfdisc(H7), nfdisc(L7), nfdisc(R8)];
  check(discs[1] == 2^4*3^6*23^6,"point septic field discriminant changed");
  check(discs[2] == 2^4*3^6*23^6,"line septic field discriminant changed");
  check(discs[3] == 2^16*3^6*23^6,"octic field discriminant changed");
  for(i = 1, 3,
    check(vecsort(Vec(factor(abs(discs[i]))[,1]~)) == [2,3,23],
      "a special-fibre field ramifies outside the branch primes"));

  \\ p = 23: tame, totally ramified in the septics, 1+7 in the octic
  check(ramification_data(H7,23) == [[7,1]],
    "point septic ramification at 23 changed");
  check(ramification_data(L7,23) == [[7,1]],
    "line septic ramification at 23 changed");
  check(ramification_data(R8,23) == [[1,1],[7,1]],
    "octic ramification at 23 changed");
  check(23 % 7 != 0,"tameness at 23 lost");

  \\ p = 3: wild, order three
  check(ramification_data(H7,3) == [[1,1],[3,1],[3,1]],
    "point septic ramification at 3 changed");
  check(ramification_data(L7,3) == [[1,1],[3,1],[3,1]],
    "line septic ramification at 3 changed");
  check(ramification_data(R8,3) == [[1,1],[1,1],[3,1],[3,1]],
    "octic ramification at 3 changed");

  \\ p = 2: tame of order three on sheets, wild of order four below
  check(ramification_data(H7,2) == [[1,1],[3,1],[3,1]],
    "point septic ramification at 2 changed");
  check(ramification_data(L7,2) == [[1,1],[3,1],[3,1]],
    "line septic ramification at 2 changed");
  check(ramification_data(R8,2) == [[4,2]],
    "octic ramification at 2 changed");

  out = fileopen("results/special_fibre_ramification_arithmetic_summary.txt","w");
  filewrite(out,"schema=m23.cover-investigation.special-fibre-ramification-arithmetic.v1");
  filewrite(out,"point_septic_field_discriminant=2^4*3^6*23^6");
  filewrite(out,"line_septic_field_discriminant=2^4*3^6*23^6");
  filewrite(out,"octic_field_discriminant=2^16*3^6*23^6");
  filewrite(out,"unramified_outside_2_3_23=true");
  filewrite(out,"septic_ef_at_23=[[7,1]]");
  filewrite(out,"octic_ef_at_23=[[1,1],[7,1]]");
  filewrite(out,"ramification_at_23_is_tame=true");
  filewrite(out,"inertia_at_23_is_singer_cycle=true");
  filewrite(out,"non_galois_residue_degrees_do_not_determine_decomposition=true");
  filewrite(out,"residue_closure_decomposition_at_23=7:3");
  filewrite(out,"septic_ef_at_3=[[1,1],[3,1],[3,1]]");
  filewrite(out,"octic_ef_at_3=[[1,1],[1,1],[3,1],[3,1]]");
  filewrite(out,"ramification_at_3_is_wild=true");
  filewrite(out,"septic_ef_at_2=[[1,1],[3,1],[3,1]]");
  filewrite(out,"octic_ef_at_2=[[4,2]]");
  filewrite(out,"septics_tame_at_2_octic_wild_at_2=true");
  filewrite(out,"wild_inertia_at_2_lies_in_line_module=true");
  fileclose(out);

  print("all three constituent fields unramified outside {2,3,23}");
  print("p=23: Singer inertia e=7; residue closure decomposition is 7:3");
  print("p=3: wild, e=3 with orbits 1,3,3 and 1,1,3,3");
  print("p=2: septics tame e=3, octic wild e=4 -- wildness in the line module")
};

certify_special_fibre_ramification();
