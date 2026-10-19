#!/usr/bin/env python3
"""Append one iteration of performance work to perf_history.csv (next to the checkout).

  python3 history.py --iter 1 --item "EOS jerk loop bound" --decision kept \
      --tests "892 passed" --new RESULT.json --prev PREV.json --base BASE.json [--commit HASH]

Columns: iteration metadata, geomean speedup vs previous kept state and vs the original
baseline, then one speedup column per benchmark (vs the original baseline).
"""
import argparse
import csv
import datetime
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.abspath(os.path.join(HERE, "..", "..", "..", "perf_history.csv"))


def speedups(base, new):
    out = {}
    for name, nb in new["benchmarks"].items():
        bb = base["benchmarks"].get(name)
        if bb and "rate" in bb and "rate" in nb:
            out[name] = run.paired_speedup(bb, nb)[0]
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--iter", required=True)
    p.add_argument("--item", required=True)
    p.add_argument("--decision", required=True)
    p.add_argument("--tests", required=True)
    p.add_argument("--new", required=True)
    p.add_argument("--prev", required=True)
    p.add_argument("--base", required=True)
    p.add_argument("--commit", default="")
    p.add_argument("--csv", default=CSV)
    a = p.parse_args()
    new, prev, base = run.load(a.new), run.load(a.prev), run.load(a.base)
    s_prev, s_base = speedups(prev, new), speedups(base, new)
    names = list(base["benchmarks"].keys())
    row = {
        "iter": a.iter, "date": datetime.date.today().isoformat(), "item": a.item,
        "commit": a.commit or new["git"]["commit"], "decision": a.decision, "tests": a.tests,
        "geomean_vs_prev": "%.4f" % run.geomean(list(s_prev.values())),
        "geomean_vs_base": "%.4f" % run.geomean(list(s_base.values())),
        "result_file": os.path.basename(a.new),
    }
    for n in names:
        row[n] = "%.4f" % s_base[n] if n in s_base else ""
    fields = list(row.keys())
    exists = os.path.exists(a.csv)
    with open(a.csv, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        if not exists:
            w.writeheader()
        w.writerow(row)
    print("logged: iter %s, %s, vs prev %s, vs base %s" % (a.iter, a.decision, row["geomean_vs_prev"], row["geomean_vs_base"]))


if __name__ == "__main__":
    main()
