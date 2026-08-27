\\ Independent PARI/GP check of the trace-sextic Galois group.

g = x^6 - 6*x^5 + 14*x^4 - 2*x^3 - 27*x^2 + 44*x - 44;

if (!polisirreducible(g), error("trace sextic is reducible"));
if (poldisc(g) != 2^22*11*23^4, error("wrong polynomial discriminant"));
if (nfdisc(g) != 2^4*11*23^4, error("wrong field discriminant"));

G = polgalois(g);
if (G[1] != 720, error("wrong Galois-group order"));
if (G[4] != "S6", error("wrong Galois group"));

print("PASS_PARIGP_HURWITZ_TRACE_SEXTIC_GALOIS_GROUP_S6");
quit
