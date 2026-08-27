\\ Unconditional verification that the splitting field of the ramified
\\ special-fibre octic has degree exactly 1344 = |2^3:PSL(3,2)|.  This is the
\\ heavy independent check for the special-fibre field tower; it needs a few
\\ minutes and a large PARI stack.

check(condition,message) = if(!condition,error(message));

default(parisizemax, 6*10^9);
default(threadsizemax, 10^9);

verify_octic_splitting_degree() =
{
  my(R8, S, out);

  R8 = V^8 + 156*V^7 - 1290*V^6 - 1217744*V^5 - 40611765*V^4
     + 2044635024*V^3 + 115083889782*V^2 + 968032514052*V
     - 6058197515007;

  S = nfsplitting(R8);
  check(poldegree(S) == 1344,"octic splitting degree changed");

  out = fileopen("results/octic_splitting_degree_summary.txt","w");
  filewrite(out,"schema=m23.cover-investigation.octic-splitting-degree.v1");
  filewrite(out,"octic_splitting_field_degree=1344");
  filewrite(out,"equals_order_of_affine_group_2^3_PSL32=true");
  fileclose(out);

  print("splitting field degree of the ramified octic = 1344")
};

verify_octic_splitting_degree();
