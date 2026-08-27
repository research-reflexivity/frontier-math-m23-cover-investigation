\\ Verify that the real place of the special fibre above the rational 2A point realizes the
\\ transvection predicted by the real frame cocycle: the unramified
\\ septic has exactly 3 real sheets forming one Fano line (the axis),
\\ the line septic has exactly 3 real lines, the ramified octic has
\\ exactly 4 real double points, and the flag center is the unique real
\\ sheet lying on both mixed real lines.

check(condition,message) = if(!condition,error(message));

certify_real_frame_cocycle() =
{
  my(H7, L7, R8, rootsH, realH, pairSums, realL, axisSum, axisGap,
     minSeparation, mixedCounts, centerIndex, tolerance, guard, out);

  default(realprecision, 200);

  H7 = V^7 - 312*V^6 + 43848*V^5 - 3914896*V^4
     + 241123008*V^3 - 9876755712*V^2 + 251799776256*V
     - 2694891036672;
  R8 = V^8 + 156*V^7 - 1290*V^6 - 1217744*V^5 - 40611765*V^4
     + 2044635024*V^3 + 115083889782*V^2 + 968032514052*V
     - 6058197515007;
  L7 = Y^7 - 936*Y^6 + 379728*Y^5 - 86588640*Y^4
     + 11922182400*Y^3 - 982633068288*Y^2 + 44572755885824*Y
     - 859311381083136;

  check(polisirreducible(H7),"unramified septic is reducible");
  check(polisirreducible(R8),"ramified octic is reducible");
  check(polisirreducible(L7),"line septic is reducible");

  \\ exact Sturm counts: the real-place cycle structure
  check(polsturm(H7) == 3,"unramified septic real root count changed");
  check(polsturm(L7) == 3,"line septic real root count changed");
  check(polsturm(R8) == 4,"ramified octic real root count changed");

  \\ numeric flag identification with explicit separation guards
  tolerance = 10.0^(-100);
  guard = 1.0;

  rootsH = polroots(H7);
  realH = [];
  pairSums = [];
  for(i = 1, 7,
    if(abs(imag(rootsH[i])) < tolerance,
      realH = concat(realH,[real(rootsH[i])]),
      if(imag(rootsH[i]) > 0,
        pairSums = concat(pairSums,[2*real(rootsH[i])])
      )
    )
  );
  realH = vecsort(realH);
  check(#realH == 3 && #pairSums == 2,"real sheet classification failed");

  realL = polrootsreal(L7);
  check(#realL == 3,"real line classification failed");
  minSeparation = vecmin(vector(#realL-1,k,abs(realL[k+1]-realL[k])));
  check(minSeparation > guard,"real line separation too small");

  \\ the three real sheets form one line: their sum is a real line root
  axisSum = realH[1]+realH[2]+realH[3];
  axisGap = vecmin(vector(#realL,k,abs(axisSum-realL[k])));
  check(axisGap < tolerance,"real sheet triple is not the axis line");

  \\ the center: the unique real sheet on both mixed real lines
  mixedCounts = vector(3);
  for(i = 1, 3,
    for(j = 1, 2,
      my(candidate = realH[i]+pairSums[j],
         gap = vecmin(vector(#realL,k,abs(candidate-realL[k]))));
      check(gap < tolerance || gap > guard,"ambiguous line match");
      if(gap < tolerance,mixedCounts[i]++)
    )
  );
  check(vecsort(mixedCounts) == [0,0,2],"center multiplicity changed");
  centerIndex = 0;
  for(i = 1, 3,if(mixedCounts[i] == 2,centerIndex = i));
  check(centerIndex == 3,"center is no longer the largest real sheet");

  out = fileopen("results/real_frame_cocycle_arithmetic_summary.txt","w");
  filewrite(out,"schema=m23.cover-investigation.real-frame-cocycle-arithmetic.v1");
  filewrite(out,"unramified_septic_real_roots=3");
  filewrite(out,"line_septic_real_roots=3");
  filewrite(out,"ramified_octic_real_roots=4");
  filewrite(out,"real_cycle_type_on_sheets=1^3_2^2");
  filewrite(out,"conjugation_is_a_transvection=true");
  filewrite(out,"real_sheet_triple_is_a_fano_line=true");
  filewrite(out,Str("axis_sum_gap<",1.0*tolerance));
  filewrite(out,Str("real_line_min_separation=",
    precision(minSeparation,20)));
  filewrite(out,"center_mixed_line_multiplicities=[0,0,2]");
  filewrite(out,"center_is_largest_real_sheet=true");
  filewrite(out,Str("center_real_sheet=",precision(realH[3],30)));
  filewrite(out,"real_place_matches_frame_cocycle_transvection=true");
  fileclose(out);

  print("real roots: septic 3, line septic 3, octic 4");
  print("the three real sheets form the axis line, gap < 1e-100");
  print("the center is the largest real sheet ",
    precision(realH[3],20));
  print("real place realizes the marked flag (P,Delta)")
};

certify_real_frame_cocycle();
