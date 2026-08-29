// Independent Magma certificate for normalized wild-parameter orientation.

p := 23;
F23 := GF(p);

function NormOneParity(value)
    assert value*value^p eq 1;
    sign := value^12;
    assert sign eq 1 or sign eq -1;
    return sign eq -1 select 1 else 0;
end function;

function NormalizedParameterTransport(raw, leading)
    // If delta^23=leading and u=s/delta, then
    // eta=raw*delta/delta^sigma.
    delta := leading^p;
    assert delta^p eq leading;
    return raw*delta/(delta^p);
end function;

// Unramified E8 tail.
PE<e> := PolynomialRing(F23);
kE<rE> := ext<F23 | e^2+18*e+1>;
baseUnits := [kE!j : j in [1..22]];
assert &and[IsSquare(b) : b in baseUnits];
frameOrientationChanges :=
    Seqset([b^11 eq -1 select 1 else 0 : b in baseUnits]);
assert frameOrientationChanges eq {0,1};
e8CurveA := 18*rE+5;
e8Leading := e8CurveA^8;
e8XTransport := 18*rE+1;
e8YTransport := rE+18;
e8Raw := e8YTransport/e8XTransport^2;
e8Target := 22*rE+5;
e8Eta := NormalizedParameterTransport(e8Raw,e8Leading);
assert e8Leading eq 15;
assert e8Leading^p*e8Raw^(-p) eq e8Target*e8Leading;
assert e8Eta eq e8Target;
assert e8Eta^(-p) eq e8Target;
assert NormOneParity(e8Eta) eq 1;

// Ramified A2 tail.
PA<a> := PolynomialRing(F23);
kA<rA> := ext<F23 | a^2+a+1>;
a2CurveA := 14*rA+15;
a2H10 := 18*rA+7;
a2Leading := a2H10*a2CurveA^(-11);
a2XTransport := 22*rA+22;
a2YTransport := 7*rA+11;
a2Raw := a2XTransport/a2YTransport;
a2Target := kA!1;
a2Eta := NormalizedParameterTransport(a2Raw,a2Leading);
assert a2Leading^p*a2Raw^(-p) eq a2Target*a2Leading;
assert a2Eta eq 1;
assert a2Eta^(-p) eq a2Target;
assert NormOneParity(a2Eta) eq 0;
assert NormOneParity(a2Raw) eq 1;
assert NormOneParity(a2Leading^22) eq 1;

// Inner A6 transport and the outer A*X^23+B*X^8 chart.
a6Target := 11*rA+7;
assert Order(a6Target) eq 24;
assert (-p) mod 24 eq 1;
assert NormOneParity(a6Target) eq 1;

PW<v> := PolynomialRing(kA);
kA4<w> := ext<kA | v^2-(6*rA+22)>;
r4 := kA4!rA;
outerA := 15*r4*w;
outerB := -4*r4-9;
outerASigma := outerA^p;
outerBSigma := outerB^p;
lambdaPower15 := outerA*outerBSigma/(outerB*outerASigma);
outerOverInnerPower15 :=
    (outerBSigma/outerB)^15*lambdaPower15^8/(kA4!a6Target)^15;
assert outerOverInnerPower15 eq r4;
assert r4^3 eq 1;

// Rational E8 is trivial.
assert NormOneParity(F23!1) eq 0;

print "local_degree23_orientation_identity=eta^(-23)=tau";
print "unramified_E8_normalized_parameter_transport=target_transport_order8";
print "A2_normalized_parameter_transport=target_transport_identity";
print "A6_outer_over_inner_transport_power15=r_of_order3";
print "A6_outer_and_inner_quadratic_orientations_agree=true";
print "Hilbert90_square_class_alone_does_not_orient_affine_F23_frame=true";
print "SCOPE=A6_target_transport_to_pointed_M23_node_realization_still_required";
print "PASS_WILD_PARAMETER_ORIENTATION";
quit;
