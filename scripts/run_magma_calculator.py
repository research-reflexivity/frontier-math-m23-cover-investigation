#!/usr/bin/env python3
"""Run the public Magma inputs sequentially on Sydney's normal XML endpoint.

This transmits the selected source to the public calculator. It does not run
by default in verify-all. Exact request/response records are retained, including
failures; a printed PASS after a Magma error is never accepted as success.
"""
import argparse
import hashlib
import json
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ENDPOINT = "https://magma.maths.usyd.edu.au/xml/calculator.xml"
STANDARD = [
    "verification/verify_optimal_23_4.m",
    "verification/certify_gauss_prolongation_obstruction.m",
    "verification/certify_geometric_reflection_obstruction.m",
    "verification/verify_canonical_quadric.m",
    "verification/verify_hurwitz_degree23_branch.m",
    "verification/verify_hurwitz_degree23_geometry.m",
    "verification/verify_hurwitz_galois_closure.m",
    "verification/verify_hurwitz_local_23.m",
    "notes/certify_fano_affine_odd_fixed_point_lemma.m",
    "notes/certify_pinched_tag_finite_identities.m",
    "notes/certify_wild_parameter_orientation.m",
    "notes/audit_pointed_relative_bockstein.m",
    "notes/audit_log_quadratic_orientation_line.m",
]
EXTRA = [
    "verification/verify_harmonic_model_generated.m",
    "notes/compute_q0_tail_function_field_galois.m",
    "notes/compute_rational_e8_tail_function_field_galois.m",
]
ERROR_RE = re.compile(r"Runtime error|Syntax error|User error|Assertion failed|assertion failure|CERTIFICATE_FAILURE|Internal error|Segmentation fault|time limit|time exceeded|killed|interrupt", re.I)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def parse_response(raw, expected):
    tree = ET.fromstring(raw)
    if tree.tag != "calculator":
        raise ValueError("Unexpected response document")
    lines = [line.text or "" for line in tree.findall("./results/line")]
    headers = {child.tag: child.text for child in tree.findall("./headers/*")}
    warnings = [x.text for x in tree.findall("./headers/warning") + tree.findall("./headers/alert") + tree.findall("./offline") if x.text]
    output = "\n".join(lines)
    errors = [line for line in lines if ERROR_RE.search(line)]
    # Magma wraps long printed lines. Normalize whitespace, not mathematical
    # content, before matching the fixed completion messages.
    folded = " ".join(output.split())
    missing = [marker for marker in expected if " ".join(marker.split()) not in folded]
    passed = bool(expected) and not errors and not warnings and not missing and bool(headers.get("version"))
    return {"status": "PASS" if passed else "INCOMPLETE_OR_FAILED", "headers":headers,
            "lines":lines,"warnings":warnings,"errors":errors,"missing_markers":missing}


def expected_markers(source):
    # Literal, unconditional completion messages in the maintained scripts.
    # The full output is separately scanned for errors and server warnings.
    return re.findall(r'print\s+"([^"\n]*PASS[^"\n]*)"\s*;', source)


def run(relative, output_dir, seed=None):
    if relative not in STANDARD + EXTRA:
        raise ValueError("Only the explicitly listed public inputs may be submitted")
    raw = (ROOT/relative).read_bytes()
    if len(raw) > 50000:
        raise ValueError("Input exceeds the calculator's 50000-byte limit")
    source = raw.decode("utf8")
    markers = expected_markers(source)
    if not markers:
        raise ValueError("Input has no literal completion marker")
    started = datetime.now(timezone.utc)
    record = {"schema":"m23.magma-calculator-run.v1","started_at":started.isoformat(),
              "endpoint":ENDPOINT,"input_path":relative,"input_bytes":len(raw),
              "input_sha256":sha(raw),"input_source":source,"expected_markers":markers,
              "scope":"Execution of this exact arithmetic input; no certification of the written geometric or analytic arguments."}
    submitted = (f"SetSeed({seed});\n" if seed is not None else "") + source
    if len(submitted.encode('utf8')) > 50000:
        raise ValueError("Submitted wrapper exceeds the calculator input limit")
    record.update(seed_override=seed,submitted_source=submitted,
                  submitted_sha256=sha(submitted.encode('utf8')),
                  submitted_bytes=len(submitted.encode('utf8')))
    print(f"RUN {relative} ({len(raw)} bytes)",flush=True)
    before = time.monotonic()
    request = urllib.request.Request(ENDPOINT, data=urllib.parse.urlencode({"input":submitted}).encode("utf8"),
                                     headers={"Content-Type":"application/x-www-form-urlencoded"})
    try:
        with urllib.request.urlopen(request,timeout=85) as response:
            payload = response.read()
            record["http_status"] = response.status
        record["response_xml"] = payload.decode("utf8")
        record["response_sha256"] = sha(payload)
        record.update(parse_response(payload,markers))
    except Exception as error:
        record.update(status="TRANSPORT_OR_RESPONSE_ERROR",error=str(error))
    record["elapsed_wall_seconds"] = round(time.monotonic()-before,3)
    record["completed_at"] = datetime.now(timezone.utc).isoformat()
    stamp=started.strftime("%Y%m%dT%H%M%S%fZ")
    target=output_dir/f"{Path(relative).stem}--{stamp}.json"
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(record,indent=2)+"\n")
    print(f"{record['status']} {relative}: {record.get('headers',{})}",flush=True)
    if record["status"] != "PASS":
        print(json.dumps({k:record.get(k) for k in ['error','errors','warnings','missing_markers','lines']},indent=2),flush=True)
    print(f"RECORD {target.relative_to(ROOT)}",flush=True)
    return record


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--suite",choices=["all","standard","extra"],default="all")
    parser.add_argument("--only",choices=STANDARD+EXTRA)
    parser.add_argument("--seed",type=int,help="Explicit reproducible random seed; leaves the assertions unchanged")
    parser.add_argument("--output-dir",default="verification/magma_runs/2026-09-07")
    args=parser.parse_args()
    if args.seed is not None and (not args.only or not 0 <= args.seed < 2**32):
        parser.error("--seed requires --only and a 32-bit nonnegative seed")
    paths=[args.only] if args.only else (EXTRA+STANDARD if args.suite=="all" else STANDARD if args.suite=="standard" else EXTRA)
    output_dir=ROOT/args.output_dir
    records=[run(relative,output_dir,args.seed) for relative in paths]
    print(f"SUMMARY {sum(r['status']=='PASS' for r in records)}/{len(records)} complete successful executions",flush=True)
    raise SystemExit(0 if all(r["status"]=="PASS" for r in records) else 1)


if __name__ == '__main__':
    main()
