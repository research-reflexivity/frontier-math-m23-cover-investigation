\\ Exact Newton-polygon and derivative-valuation certificate at D=T^2+23.
read("../data/Fint_coefficients_Z.gp");

check(condition, message) = if(!condition, error(message));

verify_boundary() =
{
  my(D = T^2 + 23, edge_points = 0, derivative_min = 10^9,
     derivative_attainers = 0, c, n, branch_order);
  for(i = 0, 23,
    c = polcoef(F, i, V);
    if(c != 0,
      n = valuation(c, D);
      check(23*n - 4*i >= 0, "coefficient lies below claimed Newton edge");
      if(23*n - 4*i == 0, edge_points++)
    )
  );
  check(edge_points == 2, "claimed Newton edge has unexpected lattice points");

  for(i = 1, 23,
    c = i*polcoef(F, i, V);
    if(c != 0,
      branch_order = 23*valuation(c,D) - 4*(i-1);
      if(branch_order < derivative_min,
        derivative_min = branch_order;
        derivative_attainers = 1,
        if(branch_order == derivative_min, derivative_attainers++)
      )
    )
  );
  check(derivative_min == 4, "unexpected valuation of Phi_V");
  check(derivative_attainers == 1, "possible cancellation in Phi_V leading term");

  print("PASS_BOUNDARY_VALUATIONS");
  print("ord_D_edge_slope=4/23 ord_V=-4 ord_dT=22 ord_PhiV=4 ord_dT_over_PhiV=18");
};

verify_boundary();
