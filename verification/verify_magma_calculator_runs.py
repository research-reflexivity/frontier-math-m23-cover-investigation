#!/usr/bin/env python3
"""Audit retained calculator transcripts, current inputs and the compact index.

This validates provenance and parses execution results. It is not a new Magma
execution or a formal proof checker. It preserves unsuccessful attempts.
"""
import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT / "scripts"))
from run_magma_calculator import ENDPOINT, STANDARD, EXTRA, expected_markers, parse_response

RUN_DIR = ROOT / "verification/magma_runs/2026-09-07"
SUMMARY = ROOT / "verification/calculator_2026_09_07_magma_summary.json"


def require(value,message):
    if not value:
        raise ValueError(message)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def build_summary():
    records=[]
    for file in sorted(RUN_DIR.glob("*.json")):
        raw=file.read_bytes(); data=json.loads(raw)
        require(data["schema"]=="m23.magma-calculator-run.v1",f"Bad record schema: {file}")
        require(data["endpoint"]==ENDPOINT,f"Unexpected runner: {file}")
        relative=data["input_path"]
        require(relative in STANDARD+EXTRA,f"Unlisted input: {file}")
        source=data["input_source"].encode("utf8")
        require(sha(source)==data["input_sha256"],f"Input record digest mismatch: {file}")
        require(source==(ROOT/relative).read_bytes(),f"Record does not apply to current input: {relative}")
        require(len(source)==data["input_bytes"]<=50000,f"Input byte count mismatch: {file}")
        override=data.get("seed_override")
        require(override is None or (type(override) is int and 0<=override<2**32),f"Bad explicit seed: {file}")
        submitted=(f"SetSeed({override});\n" if override is not None else "")+data["input_source"]
        require(data.get("submitted_source",data["input_source"])==submitted,f"Unexpected source wrapper: {file}")
        if "submitted_sha256" in data:
            require(sha(submitted.encode())==data["submitted_sha256"],f"Submission digest mismatch: {file}")
            require(len(submitted.encode())==data["submitted_bytes"]<=50000,f"Submission size mismatch: {file}")
        expected=expected_markers(data["input_source"])
        require(expected==data["expected_markers"],f"Completion markers changed: {file}")
        parsed={"status":"TRANSPORT_OR_RESPONSE_ERROR"}
        if "response_xml" in data:
            response=data["response_xml"].encode('utf8')
            require(sha(response)==data["response_sha256"],f"Response digest mismatch: {file}")
            parsed=parse_response(response,expected)
            require(parsed["headers"]==data["headers"],f"Changed headers: {file}")
            require(parsed["lines"]==data["lines"],f"Changed transcript lines: {file}")
            if parsed["status"]=="PASS":
                require(data["http_status"]==200,f"Non-200 success: {file}")
                require(parsed["headers"]["version"]=="2.29-10",f"Unexpected version: {file}")
                require(float(parsed["headers"]["time"])>0,f"Missing execution time: {file}")
        # One first-pass geometry result was marked incomplete by the initial
        # reader because Magma wrapped a long line. Raw XML is immutable; the
        # current interpretation is separately derived rather than rewritten.
        if data["status"]!=parsed["status"]:
            require(data["status"]=="INCOMPLETE_OR_FAILED" and parsed["status"]=="PASS" and
                    not data["errors"] and not data["warnings"] and bool(data["missing_markers"]),
                    f"Unexplained status discrepancy: {file}")
        records.append({"record":file.relative_to(ROOT).as_posix(),"record_sha256":sha(raw),
                        "input":relative,"input_sha256":data["input_sha256"],
                        "input_bytes":len(source),"status":parsed["status"],
                        "initial_reader_status":data["status"],"started_at":data["started_at"],
                        "magma_version":parsed.get("headers",{}).get("version"),
                        "server_seed":parsed.get("headers",{}).get("seed"),"seed_override":override,
                        "seconds":parsed.get("headers",{}).get("time"),
                        "warnings":parsed.get("warnings",[])})
    require(records,"No calculator records found")
    chosen=[];missing=[]
    for relative in STANDARD+EXTRA:
        success=sorted([r for r in records if r["input"]==relative and r["status"]=="PASS"],key=lambda r:r["started_at"])
        if not success:missing.append(relative)
        else:chosen.append(success[-1])
    return {"schema":"m23.magma-calculator-suite.v1","date":"2026-09-07","magma_version":"2.29-10",
            "status":"PASS_ALL_PUBLIC_INPUTS" if not missing else "INCOMPLETE",
            "input_count":len(STANDARD+EXTRA),"successful_input_count":len(chosen),
            "attempt_count":len(records),"successful_runs":chosen,"all_attempts":records,
            "missing_successful_inputs":missing,
            "scope":"Exact executions of the 16 listed arithmetic scripts. No new Arb continuation, no local analytic Hensel audit, and no mechanical certification of the written geometric proofs or unresolved relative Fano-affine comparison."}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-summary",action="store_true")
    args=parser.parse_args()
    subprocess.run([sys.executable,str(ROOT/"scripts/export_harmonic_magma.py"),"--check"],check=True)
    summary=build_summary()
    if args.write_summary:
        SUMMARY.write_text(json.dumps(summary,indent=2)+"\n")
    else:
        require(json.loads(SUMMARY.read_text())==summary,"Calculator suite index is stale")
    require(not summary["missing_successful_inputs"],"No successful current run for: "+", ".join(summary["missing_successful_inputs"]))
    print(f"PASS {summary['successful_input_count']}/{summary['input_count']} current Magma inputs; {summary['attempt_count']} intact request/response records")
    print("SCOPE: recorded execution and provenance checks, not a fresh Magma run or a geometric proof checker")


if __name__ == '__main__':
    main()
