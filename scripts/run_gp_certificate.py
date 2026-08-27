#!/usr/bin/env python3
"""Run a PARI/GP certificate and reject errors that GP reports with status 0."""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ERROR_PATTERN = re.compile(
    r"(?:\*\*\*\s+(?:at top-level|in function|user error|syntax error)|"
    r"skipping file)",
    re.IGNORECASE,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("script", type=Path)
    parser.add_argument("--gp", default="gp")
    args = parser.parse_args()

    script = args.script
    if not script.is_absolute():
        script = ROOT / script
    if not script.is_file():
        raise SystemExit(f"GP certificate not found: {script}")

    completed = subprocess.run(
        [args.gp, "-q", "-f", str(script)],
        cwd=ROOT,
        text=True,
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    sys.stdout.write(completed.stdout)
    if completed.returncode != 0:
        raise SystemExit(completed.returncode)
    if ERROR_PATTERN.search(completed.stdout):
        raise SystemExit("PARI/GP reported a certificate error")


if __name__ == "__main__":
    main()
