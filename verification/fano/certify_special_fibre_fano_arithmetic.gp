\\ Certify the factorization above the rational 2A branch point of the class-6 M23 cover
\\ cover and construct the 3-subset-sum resolvent of its seven unramified
\\ sheets.  The degree-7 factor is the point action of PSL(3,2); the
\\ degree-7 factor of the 35-dimensional resolvent is its line action.

read("data/Fint_coefficients_Z.gp");

check(condition,message) = if(!condition,error(message));

wedge_triples() =
{
  my(triples = vector(35), index = 0);
  for(i = 1, 5,
    for(j = i+1, 6,
      for(k = j+1, 7,
        index++;
        triples[index] = [i,j,k]
      )
    )
  );
  check(index == 35,"wrong wedge basis size");
  triples
};

triple_row(triples,target) =
{
  for(index = 1, #triples,
    if(triples[index] == target,return(index))
  );
  error("wedge row lookup failed")
};

third_additive_compound(C,triples) =
{
  my(A3 = matrix(35,35), base, old, coeff, candidate, inversions,
     sorted, row);
  for(col = 1, 35,
    base = triples[col];
    for(pos = 1, 3,
      old = base[pos];
      for(new = 1, 7,
        coeff = C[new,old];
        if(coeff != 0,
          candidate = Vec(base);
          candidate[pos] = new;
          if(#Set(candidate) == 3,
            inversions = 0;
            for(a = 1, 2,
              for(b = a+1, 3,
                if(candidate[a] > candidate[b],inversions++)
              )
            );
            sorted = vecsort(candidate);
            row = triple_row(triples,sorted);
            A3[row,col] += (-1)^inversions*coeff
          )
        )
      )
    )
  );
  A3
};

certify_special_fibre_fano() =
{
  my(P0, P0primitive, special_fibre_content, fac0, H7, R8, gal, C, triples, A3, R35, fac35,
     degrees35, L7, N28, linegal, pointdisc, linedisc, expectedH7,
     expectedR8, expectedL7, out);

  expectedH7 = V^7 - 312*V^6 + 43848*V^5 - 3914896*V^4
    + 241123008*V^3 - 9876755712*V^2
    + 251799776256*V - 2694891036672;
  expectedR8 = V^8 + 156*V^7 - 1290*V^6 - 1217744*V^5
    - 40611765*V^4 + 2044635024*V^3
    + 115083889782*V^2 + 968032514052*V
    - 6058197515007;
  expectedL7 = Y^7 - 936*Y^6 + 379728*Y^5 - 86588640*Y^4
    + 11922182400*Y^3 - 982633068288*Y^2
    + 44572755885824*Y - 859311381083136;

  P0 = subst(F,T,0);
  special_fibre_content = content(P0);
  P0primitive = P0/special_fibre_content;
  fac0 = factor(P0primitive);
  check(matsize(fac0) == [2,2],"unexpected T=0 factor count");

  H7 = fac0[1,1];
  R8 = fac0[2,1];
  check(poldegree(H7,V) == 7 && fac0[1,2] == 1,"unexpected unramified factor");
  check(poldegree(R8,V) == 8 && fac0[2,2] == 2,"unexpected ramified factor");
  check(special_fibre_content == 279841,"unexpected special-fibre content");
  check(H7*R8^2 == P0primitive,"special-fibre factorization identity failed");
  check(H7 == expectedH7,"unramified septic coefficients changed");
  check(R8 == expectedR8,"ramified octic coefficients changed");
  check(polisirreducible(H7),"unramified septic is reducible");
  check(polisirreducible(R8),"ramified octic is reducible");
  check(poldegree(gcd(P0,deriv(P0,V)),V) == 8,"wrong repeated special-fibre degree");

  gal = polgalois(H7);
  check(gal[1] == 168 && gal[2] == 1,"unramified septic does not have PSL(3,2) Galois group");

  \\ Companion matrix of H7.  Its third additive compound has eigenvalues
  \\ alpha_i+alpha_j+alpha_k, hence characteristic polynomial the full
  \\ 3-subset-sum resolvent of degree binomial(7,3)=35.
  C = matrix(7,7);
  for(j = 1, 6,C[j+1,j] = 1);
  for(i = 1, 7,C[i,7] = -polcoef(H7,i-1,V));
  triples = wedge_triples();
  A3 = third_additive_compound(C,triples);

  R35 = charpoly(A3,Y);
  fac35 = factor(R35);
  check(matsize(fac35) == [2,2],"wrong 3-subset resolvent factor count");
  check(fac35[1,2] == 1 && fac35[2,2] == 1,"3-subset resolvent is not squarefree");

  degrees35 = [poldegree(fac35[1,1],Y),poldegree(fac35[2,1],Y)];
  check(degrees35 == [7,28],"wrong Fano line/nonline orbit degrees");
  L7 = fac35[1,1];
  N28 = fac35[2,1];
  check(L7 == expectedL7,"Fano line septic coefficients changed");
  check(polisirreducible(L7) && polisirreducible(N28),"Fano resolvent factors are not irreducible");
  check(L7*N28 == R35,"3-subset resolvent identity failed");

  linegal = polgalois(L7);
  check(linegal[1] == 168 && linegal[2] == 1,"Fano line septic does not have PSL(3,2) Galois group");
  pointdisc = nfdisc(H7);
  linedisc = nfdisc(L7);
  check(pointdisc == linedisc,"point and line field discriminants disagree");
  check(pointdisc == 1726690609296,"unexpected common field discriminant");
  check(nfisisom(H7,L7) == 0,"point and line fields unexpectedly isomorphic");

  out = fileopen("results/special_fibre_fano_arithmetic_summary.txt","w");
  filewrite(out,"schema=m23.cover-investigation.special-fibre-fano-arithmetic.v1");
  filewrite(out,Str("special_fibre_degree=",poldegree(P0,V)));
  filewrite(out,Str("special_fibre_content=",special_fibre_content));
  filewrite(out,"special_fibre_factor_degrees_and_multiplicities=[[7,1],[8,2]]");
  filewrite(out,"unramified_septic_irreducible=true");
  filewrite(out,"ramified_octic_irreducible=true");
  filewrite(out,"unramified_septic_galois_group=PSL(3,2)");
  filewrite(out,"unramified_septic_galois_group_order=168");
  filewrite(out,Str("unramified_septic=",H7));
  filewrite(out,Str("ramified_octic=",R8));
  filewrite(out,Str("unramified_septic_discriminant_factorization=",factor(poldisc(H7))));
  filewrite(out,"three_subset_resolvent_degree=35");
  filewrite(out,"three_subset_resolvent_factor_degrees=[7,28]");
  filewrite(out,"fano_line_septic_irreducible=true");
  filewrite(out,"fano_line_septic_galois_group=PSL(3,2)");
  filewrite(out,Str("fano_line_septic=",L7));
  filewrite(out,Str("point_and_line_field_discriminant=",pointdisc));
  filewrite(out,"point_and_line_field_discriminant_factorization=2^4*3^6*23^6");
  filewrite(out,"point_and_line_fields_isomorphic=false");
  filewrite(out,"point_and_line_fields_are_Gassmann_twins=true");
  filewrite(out,"nonline_degree28_factor_irreducible=true");
  filewrite(out,"rational_degree3_special_fibre_subscheme_exists=false");
  filewrite(out,"correct_descended_object=the_full_degree7_Fano_line_scheme");
  fileclose(out);

  print("special-fibre factorization: irreducible septic times irreducible octic squared");
  print("unramified septic Galois group: PSL(3,2), order 168");
  print("3-subset resolvent factors irreducibly as degrees 7 and 28");
  print("no rational three-sheet special-fibre subscheme exists")
};

certify_special_fibre_fano();
