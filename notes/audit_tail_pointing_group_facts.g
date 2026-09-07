# Group facts needed to transfer pointed quotient symmetries to G-covers.
SizeScreen([100000,100000]);;
G := MathieuGroup(23);;
P := SylowSubgroup(G,23);;
D := Normalizer(G,P);;
M := Stabilizer(G,1);;
if Size(G)<>10200960 or Size(D)<>23*11
   or Size(Centralizer(G,P))<>23 or Normalizer(G,D)<>D
   or Normalizer(G,M)<>M or Size(AutomorphismGroup(G))<>Size(G) then
    Error("tail-pointing group facts changed");
fi;
Print("D=23:11; C_G(P)=P; N_G(D)=D; N_G(M22)=M22; Out(G)=1\n");
Print("PASS_TAIL_POINTING_GROUP_FACTS\n");
QUIT_GAP(0);
