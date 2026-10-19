#!/usr/bin/env python3
"""Correctness check for performance work: rebuild the Python package and run the unit tests.

  python3 check.py            # pip install -e . (rebuilds the C extension), then pytest
  python3 check.py --no-build

Skips tests that need network (Horizons) or matplotlib (plotting), and the four
module-level `test_method` helpers that pytest collects by mistake.
Prints one summary line: passed/failed/errors and runtime; exit code 0 if no failures.
"""
import argparse
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.abspath(os.path.join(HERE, "..", ".."))
DESELECT = "not horizons and not plotting and not test_method"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--no-build", action="store_true")
    p.add_argument("--python", default=sys.executable)
    a = p.parse_args()
    t0 = time.time()
    if not a.no_build:
        res = subprocess.run([a.python, "-m", "pip", "install", "-q", "-e", "."], cwd=SRC,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if res.returncode != 0:
            print(res.stdout[-3000:])
            print("CHECK: build FAILED")
            return 1
    t1 = time.time()
    res = subprocess.run([a.python, "-m", "pytest", "-q", "-p", "no:cacheprovider", "rebound/tests", "-k", DESELECT],
                         cwd=SRC, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    t2 = time.time()
    tail = res.stdout.strip().splitlines()[-1] if res.stdout.strip() else ""
    counts = {k: int(v) for v, k in re.findall(r"(\d+) (passed|failed|error|errors)", tail)}
    failed = counts.get("failed", 0) + counts.get("error", 0) + counts.get("errors", 0)
    if failed:
        print("\n".join(l for l in res.stdout.splitlines() if l.startswith(("FAILED", "ERROR"))))
    print("CHECK: %s | build %.0f s, tests %.0f s" % (tail, t1 - t0, t2 - t1))
    return 0 if (res.returncode == 0 and failed == 0) else 1


if __name__ == "__main__":
    sys.exit(main())
