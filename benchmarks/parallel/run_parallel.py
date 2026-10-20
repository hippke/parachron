#!/usr/bin/env python3
"""Multicore benchmark runner for REBOUND (stdlib only).

Builds two libraries from one source tree:
  serial  - the stock build (`make`), the single-core reference users get by default
  openmp  - the same source with OPENMP=1 (built in a copy, <build>/omp)
and runs every case of pbench.c with the serial library and with the OpenMP library at
1, 2, ... threads. Reports the median wall time over --repeat runs, the speedup over the
serial build (the number that matters: does using more cores make the run finish sooner?),
and checks the final-state signatures:
  - "threads": the OpenMP result must not depend on the thread count (bitwise)
  - "vs serial": relative difference to the serial build (round-off level expected)

  python3 run_parallel.py                       # all cases, threads 1..nproc
  python3 run_parallel.py --only ias15_n512 --threads 1,2,4,8 --repeat 5
  python3 run_parallel.py --src ../../ --json out.json
  python3 run_parallel.py --compare-src ../rebound-base   # also interleaved A/B of the OpenMP builds of two trees
  python3 run_parallel.py --all-cores        # hybrid CPUs: always use the efficiency cores too
  python3 run_parallel.py --performance-cores   # hybrid CPUs: never use them (default: auto by N)
  python3 run_parallel.py --cflags="-DREB_OMP_MIN_WORK_GRAVITY=65536."   # try another gravity threshold

Run it on a quiet machine with as many physical cores as possible.
"""
import argparse
import json
import os
import shutil
import statistics
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_SRC = os.path.abspath(os.path.join(HERE, "..", ".."))
CORES = None   # None: auto (all cores for N_active >= 1024 on hybrid CPUs), "all", "performance"
CASES = ["ias15_hltau", "ias15_embryos", "ias15_n128", "ias15_n512", "ias15_n2048", "ias15_tp2000",
         "ias15_embryos+stages", "ias15_n128+stages", "ias15_n512+stages", "ias15_n2048+stages", "ias15_tp2000+stages",
         "whfast_ss", "whfast_tp10k", "whfast_n256"]


def sh(cmd, cwd=None, env=None):
    res = subprocess.run(cmd, cwd=cwd, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if res.returncode != 0:
        sys.stderr.write(res.stdout)
        raise SystemExit("command failed: %s" % " ".join(cmd))
    return res.stdout


def build(src, build_root, jobs, cflags=""):
    """Returns paths of the serial and OpenMP pbench binaries."""
    os.makedirs(build_root, exist_ok=True)
    pbench_c = os.path.join(HERE, "pbench.c")
    out = {}
    for kind in ("serial", "openmp"):
        tree = os.path.join(build_root, kind)
        # fresh copy of src/ so the user's tree is never touched and stale objects cannot leak in
        if os.path.exists(tree):
            shutil.rmtree(tree)
        shutil.copytree(os.path.join(src, "src"), os.path.join(tree, "src"),
                        ignore=shutil.ignore_patterns("*.o", "*.so", "*.gcda", "*.gcno"))
        srcdir = os.path.join(tree, "src")
        args = ["make", "-j%d" % jobs, "OPENGL=0", "SERVER=0", "librebound.so"]
        env = None
        if kind == "openmp":
            args.append("OPENMP=1")
            if cflags:
                env = dict(os.environ)
                env["OPT"] = (env.get("OPT", "") + " " + cflags).strip()
        sh(args, cwd=srcdir, env=env)
        exe = os.path.join(build_root, "pbench_" + kind)
        cmd = ["cc", "-O3", "-I" + srcdir, pbench_c, "-L" + srcdir, "-lrebound",
               "-Wl,-rpath," + srcdir, "-lm", "-o", exe]
        if kind == "openmp":
            cmd.insert(1, "-fopenmp")
        sh(cmd)
        out[kind] = exe
    return out


def run_once(exe, case, threads, scale):
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = str(threads)
    env.setdefault("OMP_PROC_BIND", "close")
    env.setdefault("OMP_PLACES", "cores")
    if CORES:
        env["REBOUND_OMP_CORES"] = CORES
    res = subprocess.run([exe, case, str(scale)], env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        sys.stderr.write(res.stderr)
        raise SystemExit("benchmark failed: %s %s" % (exe, case))
    return json.loads(res.stdout.strip().splitlines()[-1])


def measure(exes, cases, threads, repeat, scale):
    results = {}
    for case in cases:
        rows = {}
        # interleave configurations per repeat so slow drift of the machine affects all equally
        configs = [("serial", 1)] + [("openmp", t) for t in threads]
        if "openmp_b" in exes:
            configs += [("openmp_b", t) for t in threads]
        times = {c: [] for c in configs}
        sigs = {}
        for _ in range(repeat):
            for kind, t in configs:
                r = run_once(exes[kind], case, t, scale)
                times[(kind, t)].append(r["seconds"])
                sigs[(kind, t)] = r["signature"]
        for c in configs:
            rows["%s:%d" % c] = {"median_s": statistics.median(times[c]), "all_s": times[c], "signature": sigs[c]}
        results[case] = rows
    return results


def relerr(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


def report_ab(results, threads):
    print("\nA/B of the OpenMP builds (A = --src, B = --compare-src), median time ratio B/A (>1: A faster):")
    hdr = "%-15s" % "case" + "".join(" %8s" % ("T=%d" % t) for t in threads) + "  signatures"
    print(hdr)
    for case, rows in results.items():
        line = "%-15s" % case
        same = True
        for t in threads:
            a = rows["openmp:%d" % t]
            b = rows["openmp_b:%d" % t]
            line += " %7.3fx" % (b["median_s"] / a["median_s"])
            same = same and a["signature"] == b["signature"]
        print(line + ("  identical" if same else "  CHANGED"))


def report_stages(results, threads):
    """For CASE+stages also compare to the default IAS15 mode of the serial build."""
    lines = []
    for case, rows in results.items():
        base = case.replace("+stages", "")
        if case == base or base not in results:
            continue
        ser = results[base]["serial:1"]
        line = "%-20s" % case
        for t in threads:
            line += " %7.2fx" % (ser["median_s"] / rows["openmp:%d" % t]["median_s"])
        sigs = [rows["serial:1"]["signature"]] + [rows["openmp:%d" % t]["signature"] for t in threads]
        line += "  %-9s %.1e" % ("identical" if all(s == sigs[0] for s in sigs) else "DIFFER",
                                 relerr(sigs[0], ser["signature"]))
        lines.append(line)
    if lines:
        print("\nparallel_stages=1 vs the default IAS15 mode of the serial build (speedup at T threads;")
        print("'identical': serial build and all thread counts agree bitwise; last column: difference to default mode):")
        print("\n".join(lines))


def report(results, threads):
    hdr = "%-15s %9s" % ("case", "serial[s]")
    for t in threads:
        hdr += " %8s" % ("omp%d" % t)
    hdr += "  %-9s %s" % ("threads", "vs serial")
    print(hdr)
    print("-" * len(hdr))
    for case, rows in results.items():
        ser = rows["serial:1"]
        line = "%-15s %9.3f" % (case, ser["median_s"])
        for t in threads:
            line += " %7.2fx" % (ser["median_s"] / rows["openmp:%d" % t]["median_s"])
        s_omp = [rows["openmp:%d" % t]["signature"] for t in threads]
        same = all(s == s_omp[0] for s in s_omp)
        line += "  %-9s %.1e" % ("identical" if same else "DIFFER", relerr(s_omp[0], ser["signature"]))
        print(line)
    print("\nNumbers are speedups vs the serial (stock) build: >1 means the run finishes sooner.")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src", default=DEFAULT_SRC)
    ap.add_argument("--build", default=None, help="build directory (default: <src>/../build_parallel)")
    ap.add_argument("--only", default=None)
    ap.add_argument("--threads", default=None, help="comma list (default 1,2,4,... up to nproc)")
    ap.add_argument("--repeat", type=int, default=3)
    ap.add_argument("--scale", type=float, default=1.0)
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 2)
    ap.add_argument("--json", default=None)
    ap.add_argument("--cflags", default="", help="extra flags for the OpenMP build, e.g. "
                    "'-DREB_OMP_MIN_WORK_GRAVITY=65536.' (gravity threshold) or '-DREB_OMP_MIN_WORK=1024'")
    ap.add_argument("--all-cores", action="store_true", help="on hybrid CPUs, always use the efficiency cores too "
                    "(REBOUND_OMP_CORES=all)")
    ap.add_argument("--performance-cores", action="store_true", help="on hybrid CPUs, never use the efficiency "
                    "cores (REBOUND_OMP_CORES=performance). Default: auto, all cores only for N_active >= 1024")
    ap.add_argument("--compare-src", default=None, help="second source tree: interleaved A/B of the OpenMP builds")
    args = ap.parse_args()

    global CORES
    CORES = "all" if args.all_cores else ("performance" if args.performance_cores else None)
    hybrid = os.path.exists("/sys/devices/cpu_atom/cpus")
    if hybrid:
        print("Hybrid CPU (performance + efficiency cores): REBOUND uses %s." % {
            "all": "all cores (--all-cores)",
            "performance": "only the performance cores (--performance-cores)",
            None: "all cores for systems with N_active >= 1024, otherwise only the performance cores (auto)"}[CORES])
    src = os.path.abspath(args.src)
    build_root = args.build or os.path.join(os.path.dirname(src), "build_parallel", os.path.basename(src))
    if args.threads:
        threads = [int(x) for x in args.threads.split(",")]
    else:
        n = os.cpu_count() or 1
        threads, t = [], 1
        while t < n:
            threads.append(t)
            t *= 2
        threads.append(n)
    cases = args.only.split(",") if args.only else CASES

    exes = build(src, build_root, args.jobs, args.cflags)
    if args.compare_src:
        other = os.path.abspath(args.compare_src)
        exes_b = build(other, build_root + "_B", args.jobs)
        exes["openmp_b"] = exes_b["openmp"]
    results = measure(exes, cases, threads, args.repeat, args.scale)
    report(results, threads)
    report_stages(results, threads)
    if args.compare_src:
        report_ab(results, threads)
    if args.json:
        with open(args.json, "w") as f:
            json.dump({"src": src, "threads": threads, "repeat": args.repeat, "results": results}, f, indent=1)


if __name__ == "__main__":
    main()
