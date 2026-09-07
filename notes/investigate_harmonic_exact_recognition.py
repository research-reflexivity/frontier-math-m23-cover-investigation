#!/usr/bin/env sage-python
"""Independent local lifting and tentative rational recognition.

Input: harmonic jets only. No stored characteristic-zero model is read.
Reconstruction outputs are guesses until checked by exact identities.
The checkpoint contains computed residues, not exact algebraic values.
"""
import argparse
import json
import time
from pathlib import Path
from sage.all import GF, QQ, ZZ, matrix, vector
from harmonic_finite_equations import fixed_family, mixed
from harmonic_truncated_geometry import Arithmetic
from certify_harmonic_finite_complete_system import evaluate

JROWS = [
    [13,6,0,17,0,0,8,0,11,16,7],
    [0,0,8,0,0,1,0,1,0,14,14],
    [21,14,0,21,0,0,13,0,4,16,7],
    [0,0,10,0,3,2,0,18,0,19,19],
    [13,14,0,21,0,0,19,0,11,17,6],
    [0,0,5,0,15,5,0,0,0,20,20],
    [2,18,0,17,0,0,13,0,4,4,19],
    [0,0,20,0,6,8,0,12,0,3,3],
    [13,6,0,11,0,0,0,0,21,0,0],
    [0,0,21,0,15,22,0,17,0,7,7],
    [1,12,0,0,0,0,11,0,8,2,21],
]
ROOT = [0,0,-7,0,-1,-2,0,10,0,-7,-7]
DIGITS = [
    [6,16,0,14,0,0,0,0,7,19,4],
    [0,0,15,0,1,7,0,0,0,17,17],
    [8,20,0,2,0,0,19,0,19,21,2],
]


def selected(labels):
    wanted = ["H^7*V^3", "H^8*V^2", "H^9*V", "U*H^2*V^7"]
    return list(range(7))+[labels.index(("function",0,mon)) for mon in wanted]


def centered(a):
    a=int(a)%23
    return a if a<=11 else a-23


def encode(a):
    return [str(QQ(b)) for b in list(mixed.K(a))]


def decode(a):
    return mixed.K([QQ(b) for b in a])


def recognize(a,precision):
    coeffs=list(mixed.K(a))+[QQ(0),QQ(0)]
    result=[]
    for i,b in enumerate(coeffs[:2]):
        power=(precision+1-i)//2
        modulus=ZZ(23)**power
        residue=ZZ(b.numerator())*ZZ(b.denominator()).inverse_mod(modulus)%modulus
        try:
            r=residue.rational_reconstruction(modulus)
            height=max(abs(r.numerator()),r.denominator())
            # Leave substantial spare precision; a short reconstruction
            # alone is not evidence of an exact rational coefficient.
            result.append({"candidate":str(r),"height":str(height),
                           "spare_factor":str(modulus//max(1,2*height**2))})
        except (ValueError,ArithmeticError):
            result.append(None)
    return result


def main(args):
    checkpoint=Path(args.checkpoint)
    if checkpoint.exists():
        data=json.loads(checkpoint.read_text())
        assert data["scope"]=="finite residues only; exact recognition unproved"
        point=[decode(a) for a in data["point"]]
        precision=data["precision"]
    else:
        point=list(map(mixed.K,ROOT));precision=4
        for n,row in enumerate(DIGITS,1):
            point=[a+mixed.pi**n*centered(b) for a,b in zip(point,row)]
    inverse=matrix(GF(23),JROWS).inverse()
    while precision<args.digits:
        start=time.monotonic();n=precision
        labels,values,inc,g=evaluate(point,n+6)
        chosen=selected(labels);A=Arithmetic(mixed.K,n+1)
        correction=-inverse*vector(GF(23),[A.digit(values[j],n) for j in chosen])
        point=[a+mixed.pi**n*centered(b) for a,b in zip(point,correction)]
        precision+=1
        # The next evaluation verifies this correction. A checkpoint is
        # explicitly marked as pending its last correction until then.
        data={"scope":"finite residues only; exact recognition unproved",
              "precision":precision,"point":[encode(a) for a in point],
              "last_correction_verified":False}
        checkpoint.write_text(json.dumps(data,indent=2)+"\n")
        print("LIFT",precision,"correction",list(correction),
              "seconds",round(time.monotonic()-start,2),flush=True)
    labels,values,inc,g=evaluate(point,precision+5)
    assert all(a==0 for a in values)
    data={"scope":"finite residues only; exact recognition unproved",
          "precision":precision,"point":[encode(a) for a in point],
          "last_correction_verified":True,"incidence_div_pi3":[encode(a) for a in inc]}
    checkpoint.write_text(json.dumps(data,indent=2)+"\n")
    print("ALL_MAP_EQUATIONS_VANISH_MOD_PI",precision,flush=True)
    print("INCIDENCE_DIV_PI3",inc,flush=True)
    short=Arithmetic(mixed.K,precision+3)
    full=fixed_family(point,precision+5)
    print("TENTATIVE_COEFFICIENT_RECOGNITION",flush=True)
    for mon in (full.c-full.U*full.W**2+full.V**3).monomials():
        a=short.coefficient(full.c.monomial_coefficient(mon))
        print(str(mon),encode(a),recognize(a,precision+3),flush=True)
    for label,a in zip(["mark_A","mark_B"],full.marks):
        print(label,encode(short.coefficient(a)),recognize(a,precision+3),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--digits",type=int,default=14)
    p.add_argument("--checkpoint",default="/private/tmp/m23-harmonic-local-recognition.json")
    main(p.parse_args())
