\\ Certify the arithmetic realization of the local flag at every
\\ transvection place.  A prime p is a transvection place when H7 has
\\ factorization type 1^3 2^2 modulo p.  The group certificate predicts
\\ that the three split sheets are then the axis of a flag and that
\\ exactly one of them -- the center -- lies on both mixed lines.  Both
\\ predictions are verified here for every transvection prime below
\\ 10^6, after removing the finitely many primes dividing
\\ Res(L7,N28), where accidental collisions between line sums and
\\ nonline sums are possible.  The archimedean place is the same
\\ statement at infinity, already certified in REAL_FRAME_COCYCLE.md.

check(condition,message) = if(!condition,error(message));

wedge_triples() =
{
  my(triples = List());
  for(i = 1, 5, for(j = i+1, 6, for(k = j+1, 7,
    listput(triples,[i,j,k]))));
  Vec(triples)
};

third_additive_compound(companion,triples) =
{
  my(compound = matrix(35,35));
  for(col = 1, 35,
    my(base = triples[col]);
    for(pos = 1, 3,
      for(replacement = 1, 7,
        my(candidate = base, coefficient = companion[replacement,base[pos]]);
        if(coefficient == 0, next);
        candidate[pos] = replacement;
        if(#Set(Vec(candidate)) < 3, next);
        my(entries = Vec(candidate), inversions = 0, row = 0);
        for(a = 1, 3, for(b = a+1, 3,
          if(entries[a] > entries[b], inversions++)));
        my(sorted = vecsort(entries));
        for(t = 1, 35, if(triples[t] == sorted, row = t));
        compound[row,col] += (-1)^inversions*coefficient;
      )));
  compound
};

certify_local_flag_places() =
{
  my(H7, L7, companion, triples, resolvent, lineSeptic, nonline,
     collisionResultant, badPrimes, bound, placeCount, axisCount,
     centerCount, patterns, out);

  H7 = V^7 - 312*V^6 + 43848*V^5 - 3914896*V^4
     + 241123008*V^3 - 9876755712*V^2 + 251799776256*V
     - 2694891036672;
  L7 = Y^7 - 936*Y^6 + 379728*Y^5 - 86588640*Y^4
     + 11922182400*Y^3 - 982633068288*Y^2 + 44572755885824*Y
     - 859311381083136;

  \\ rebuild the 3-subset-sum resolvent and its nonline factor
  companion = matrix(7,7);
  for(i = 1, 6, companion[i+1,i] = 1);
  for(i = 1, 7, companion[i,7] = -polcoef(H7,i-1,V));
  triples = wedge_triples();
  resolvent = charpoly(third_additive_compound(companion,triples),Y);
  check(poldegree(resolvent,Y) == 35,"resolvent degree changed");
  check(resolvent % L7 == 0,"line septic does not divide the resolvent");
  nonline = resolvent/L7;
  check(poldegree(nonline,Y) == 28,"nonline factor degree changed");
  check(polisirreducible(nonline),"nonline factor is reducible");
  collisionResultant = polresultant(L7,nonline);
  check(collisionResultant != 0,"line and nonline factors share a root");

  bound = 1000000;
  badPrimes = 2*3*23*poldisc(H7)*poldisc(L7)*collisionResultant;
  placeCount = 0;
  axisCount = 0;
  centerCount = 0;
  patterns = List();
  forprime(p = 5, bound,
    if(badPrimes % p == 0, next);
    my(factorization = factormod(H7,p), splitSheets = [], pairTraces = []);
    if(vecsort(apply(poldegree,Vec(factorization[,1]~))) != [1,1,1,2,2],
      next);
    placeCount++;
    for(i = 1, matsize(factorization)[1],
      my(f = factorization[i,1]);
      if(poldegree(f) == 1,
        splitSheets = concat(splitSheets,[-polcoef(f,0)]),
        pairTraces = concat(pairTraces,[-polcoef(f,1)])));
    check(#splitSheets == 3 && #pairTraces == 2,
      Str("split structure failed at ",p));

    \\ the three split sheets form a line: the axis
    if(subst(L7*Mod(1,p),Y,splitSheets[1]+splitSheets[2]+splitSheets[3])
        == 0, axisCount++);

    \\ exactly one split sheet lies on both mixed lines: the center
    my(incidence = vector(3));
    for(i = 1, 3, for(j = 1, 2,
      if(subst(L7*Mod(1,p),Y,splitSheets[i]+pairTraces[j]) == 0,
        incidence[i]++)));
    if(vecsort(Vec(incidence)) == [0,0,2], centerCount++);
    listput(patterns,vecsort(Vec(incidence)));
  );

  check(placeCount > 9000,"transvection place count collapsed");
  check(axisCount == placeCount,
    "some transvection place has no axis line");
  check(centerCount == placeCount,
    "some transvection place has no unique center");
  check(Set(Vec(patterns)) == [[0,0,2]],
    "incidence pattern changed");

  out = fileopen("results/local_flag_places_arithmetic_summary.txt","w");
  filewrite(out,"schema=m23.cover-investigation.local-flag-places-arithmetic.v1");
  filewrite(out,Str("prime_bound=",bound));
  filewrite(out,Str("transvection_places=",placeCount));
  filewrite(out,"collision_primes_excluded_by_resultant=true");
  filewrite(out,"split_sheets_form_a_line_at_every_place=true");
  filewrite(out,Str("axis_confirmations=",axisCount));
  filewrite(out,"unique_center_at_every_place=true");
  filewrite(out,Str("center_confirmations=",centerCount));
  filewrite(out,"incidence_pattern=[0,0,2]_always");
  filewrite(out,"archimedean_place_is_the_same_statement=true");
  filewrite(out,"marked_flag_is_local_datum_of_a_place=true");
  fileclose(out);

  print(placeCount," transvection primes below ",bound);
  print("  split sheets form the axis line: ",axisCount,"/",placeCount);
  print("  unique center, pattern [0,0,2]: ",centerCount,"/",placeCount);
  print("every transvection place carries its own marked flag")
};

certify_local_flag_places();
