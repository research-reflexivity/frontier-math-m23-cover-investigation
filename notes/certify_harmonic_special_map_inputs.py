#!/usr/bin/env sage-python
"""Independent algebraic inputs to the harmonic special-map construction.

No characteristic-zero curve, map, Hurwitz polynomial or crosswalk is
read. The two existing Magma records are hash-checked, not rerun. The
geometric lifting and descent arguments are not certified by this script.
"""
import hashlib
import json
from math import gcd
from pathlib import Path

from sage.all import GF, QQ, PolynomialRing, divisors, lcm


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
k = GF(23)
R = PolynomialRing(k, "a")
a0 = R.gen()
K = R.fraction_field()
a = K(a0)

# Start with the logarithmic datum itself, not a Swan expression read
# off a characteristic-zero map. This normalization makes the residues
# above a^11=1 squares, and those above a^11=-1 nonsquares.
omega = -a**6/(a**22-1)
witness = K.one()
residue_classes = {1: set(), -1: set()}
for c in k:
    if not c:
        continue
    residue = -c**6/(k(22)*c**21)
    assert residue == c**7
    witness *= (a-c)**int(residue)
    orbit = 1 if c**11 == 1 else -1
    residue_classes[orbit].add(int(residue**11))
assert witness.derivative()/witness == omega
assert residue_classes == {1: {1}, -1: {22}}
lam = k.multiplicative_generator()**2
assert lam.multiplicative_order() == 11
assert omega(lam*a)*lam == lam**7*omega
assert pow(7, -1, 11) == 8
assert (lam**8)**7 == lam
assert omega.numerator().valuation(a0) == 6
assert omega.denominator().degree()-omega.numerator().degree()-2 == 14
signature = [(11, 7), (1, 0), (1, 0), (11, 15)]
sigmas = [QQ(h)/m for m, h in signature]
assert sum(s-1 for s in sigmas) == -2
assert sum(s < 1 for s in sigmas) == 3
assert all(s < 2 and s != 1 for s in sigmas)
print("PASS_HARMONIC_LOG_DATUM", "omega=-a^6 da/(a^22-1)",
      "signature", signature, "opposite residue classes", residue_classes)

# Target attachment positions on the deformation datum Z0. On the
# actual central source Y0 put a=b^23; its H-quotient has z=b^11.
# Keeping this Frobenius factor is essential: N -> X0 is inseparable.
target = 2*a**11/(a**11-1)
assert target(0) == 0
assert target.numerator().leading_coefficient()/target.denominator().leading_coefficient() == 2
assert (2*k(-1)/(k(-1)-1)) == 1
assert target(-a) == target/(target-1)
C = PolynomialRing(k, "b")
b = C.fraction_field().gen()
z = b**11
central_target = target(b**23)
assert central_target == 2*z**23/(z**23-1)
assert 1/(1-central_target) == (1-z**23)/(1+z**23)
print("PASS_HARMONIC_ATTACHMENTS", "target: primitive 0, wild 1/infinity, new 2")
print("PASS_CENTRAL_FROBENIUS_QUOTIENT", "a=b^23, z=b^11, beta0=(1-z^23)/(1+z^23)")

P = PolynomialRing(k, "X")
X = P.gen()
primitive = X**23+7*X**9+18*X**2
critical = X**8+17*X
simple = X**7+12
assert primitive == critical**2*simple
assert primitive.derivative() == 17*critical
assert critical.is_squarefree() and simple.is_squarefree()
assert critical.gcd(simple) == 1
assert critical[7] == 0  # Sum of the eight critical points is zero.
print("PASS_PRIMITIVE_TAIL", "finite fibre 2^8 1^7; degree 23; genus 0")

H = X**6+2*X
Q = X**5+17
plane = X**23-X**13+5*X**8-2*X**3
assert plane == H**3*Q
assert Q.is_squarefree() and H.gcd(Q) == 1
assert 3*H.derivative()*Q+H*Q.derivative() == 10
assert k(19) == -k(1)/6
assert (3-1)*(5-1)//2 == 4
assert 6*3+5 == 23
print("PASS_NEW_TAIL", "w^3=19*(X^5+17), T=(X^6+2*X)*w; etale over A1; genus 4")

# Scheme-theoretic identities for ALL pointed scaling symmetries.
A = PolynomialRing(k, "r")
r = A.gen()
A7 = A.quotient(r**7-1, "r7")
r7 = A7.gen()
T = PolynomialRing(A7, "x")
x = T.gen()
f = x**23+7*x**9+18*x**2
assert f(r7*x) == r7**23*f
assert gcd(23-9, 23-2) == 7
B = PolynomialRing(k, "r,s,X,Y")
r, s, Xb, Yb = B.gens()
relations = B.ideal(r**5-1, s**3-1)
curve = Yb**3-Xb**5-17
f = (Xb**6+2*Xb)*Yb
assert relations.reduce(curve(r,s,r*Xb,s*Yb)-curve) == 0
assert relations.reduce(f(r,s,r*Xb,s*Yb)-r*s*f) == 0
print("PASS_POINTED_SCALING_IDENTITIES", "orders 7 and 15; fullness is proved geometrically")

for genus, finite_different, h in [(0, 8, 7), (4, 0, 15)]:
    wild_different = 2*genus-2+2*23-finite_different
    assert wild_different % 22 != 0  # Excludes cyclic inertia 23:1.
    assert wild_different == 22+2*h
    target_length = QQ(23*11)/(22*h)
    source_length = target_length/23
    assert 2*source_length == QQ(1)/h
    print("PASS_TAIL_NUMERICS", "h", h, "different", wild_different,
          "source length in pi units", 2*source_length)

identified = 22*7*15
special_automorphisms = 11*7*15
assert identified == 22*lcm(7, 15) == 2310
assert identified//special_automorphisms == 2
assert [m for m in divisors(15) if (QQ(m)/15).denominator() == 1] == [15]
assert [e for e in range(15) if (QQ(1-e)/15).is_integer()
        and QQ(1-e)/15 >= 0] == [1]
assert (22+22+8-2*23+2)//2 == 4
print("PASS_LIFT_AND_DESCENT_NUMERICS", "2310 identified / 1155 automorphisms = 2",
      "minimal good degree 15; tangent exponent 1")

for name in ("q0_tail_function_field_galois_magma_summary.json",
             "rational_e8_tail_function_field_galois_magma_summary.json"):
    data = json.loads((HERE/name).read_text())
    digest = hashlib.sha256((ROOT/data["input"]).read_bytes()).hexdigest()
    assert digest == data["input_sha256"]
    assert data["engine"] == "Magma"
    assert data["group_degree"] == 23 and data["group_order"] == 10200960
    assert data["transitivity"] == 4 and data["simple"] is True
    print("PASS_EXISTING_MAGMA_RECORD_HASH", name, digest)
print("SCOPE: exact characteristic-23 identities and record hashes; no new Magma run; lifting/descent require the written proof")
