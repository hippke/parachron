#!/usr/bin/env python3
"""REBOUND use-case benchmark runner (stdlib only).

Builds librebound from a source tree, builds benchmarks/perf/bench.c against it,
runs every benchmark for a short fixed time, stores results as JSON and compares
them with a baseline: per-benchmark speedup plus one geometric-mean score.

  python3 run.py run                      # build, run, compare with baseline (if any)
  python3 run.py run --only ss_wh_gr,trace_tp --repeat 3 --scale 2
  python3 run.py run --src ../rebound-base --label base   # benchmark another worktree
  python3 run.py set-baseline results/<file>.json
  python3 run.py compare A.json B.json
  python3 run.py ab --a ../rebound-prev --repeat 3   # interleaved A/B of two trees (most robust)
  python3 run.py probe trace_tp --time 10  # how fast does the rate settle?

See README.md for the benchmark list and methodology.
"""
import argparse
import datetime
import json
import math
import os
import platform
import shutil
import statistics
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_SRC = os.path.abspath(os.path.join(HERE, "..", ".."))
PROJECT = os.path.abspath(os.path.join(DEFAULT_SRC, ".."))
DEFAULT_RESULTS = os.path.join(PROJECT, "bench_results")

# The amount of work per benchmark (warm-up and measured chunks) is defined in
# bench.c and calibrated to ~0.5 s per benchmark on the reference machine.


def sh(cmd, cwd=None, quiet=True):
    res = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if res.returncode != 0:
        sys.stderr.write(res.stdout)
        raise SystemExit("command failed: %s" % " ".join(cmd))
    return res.stdout


def git_info(src):
    def g(*a):
        try:
            return subprocess.run(["git", *a], cwd=src, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True).stdout.strip()
        except OSError:
            return ""
    commit = g("rev-parse", "--short=10", "HEAD")
    dirty = bool(g("status", "--porcelain", "--untracked-files=no", "--", "src"))
    subject = g("log", "-1", "--format=%s")
    return {"commit": commit, "dirty_src": dirty, "subject": subject, "branch": g("rev-parse", "--abbrev-ref", "HEAD")}


# Extra flags for "aligned" builds. Without them, an unrelated code change can move a hot loop
# across an alignment boundary and shift a benchmark by ~10% (seen on rings_sei_tree when a
# one-line change in gravity.c altered the layout of collision.c's loops). Aligning functions
# and loops removes most of that layout noise. Use --plain-build for the stock Makefile flags.
ALIGN_FLAGS = "-falign-functions=64 -falign-loops=32 -falign-jumps=32"


def build(src, build_dir, jobs, plain=False):
    """Build librebound.so in src/src (its own Makefile) and the bench binary in build_dir."""
    srcdir = os.path.join(src, "src")
    flags = "" if plain else ALIGN_FLAGS
    marker = os.path.join(srcdir, ".bench_build_flags")
    old = open(marker).read() if os.path.exists(marker) else None
    if old != flags:  # the Makefile does not track flags: rebuild everything when they change
        sh(["make", "clean"], cwd=srcdir)
    env = dict(os.environ)
    if flags:
        env["OPT"] = (env.get("OPT", "") + " " + flags).strip()  # Makefile.defs appends with OPT+=
    res = subprocess.run(["make", "-j%d" % jobs, "OPENGL=0", "SERVER=0", "librebound.so"], cwd=srcdir, env=env,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if res.returncode != 0:
        sys.stderr.write(res.stdout)
        raise SystemExit("building librebound failed")
    with open(marker, "w") as f:
        f.write(flags)
    os.makedirs(build_dir, exist_ok=True)
    lib_dir = os.path.join(src, "src")
    exe = os.path.join(build_dir, "bench")
    cc = os.environ.get("CC", "cc")
    sh([cc, "-O2", "-std=c99", "-D_GNU_SOURCE", "-Wall", "-I" + lib_dir, os.path.join(HERE, "bench.c"),
        "-L" + lib_dir, "-lrebound", "-Wl,-rpath," + lib_dir, "-lm", "-o", exe])
    cc_version = sh([cc, "--version"]).splitlines()[0]
    return exe, cc_version


def bench_cmd(exe, cpu, *args):
    cmd = [exe, *args]
    if cpu is not None and shutil.which("taskset"):
        cmd = ["taskset", "-c", str(cpu)] + cmd
    return cmd


def run_json(cmd, timeout=120):
    try:
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return {"error": "timed out after %d s" % timeout}
    if res.returncode != 0:
        return {"error": "exit code %d: %s" % (res.returncode, res.stderr.strip()[-500:])}
    try:
        return json.loads(res.stdout.strip().splitlines()[-1])
    except (ValueError, IndexError):
        return {"error": "bad output: %r" % res.stdout[-300:]}


def summarize_fixed(chunks, blocks=5):
    """Rate = total work / total time over the measured chunks.

    Spread: the measured chunks are split into `blocks` consecutive groups;
    spread = (max - min) / median of the group rates. It estimates timing noise
    within one run (the work itself is identical between runs).
    """
    work = sum(w for w, t in chunks)
    time = sum(t for w, t in chunks)
    rate = work / time if time > 0 else float("nan")
    n = len(chunks)
    k = min(blocks, n)
    rs = []
    for i in range(k):
        part = chunks[i * n // k:(i + 1) * n // k]
        tw, tt = sum(w for w, t in part), sum(t for w, t in part)
        if tt > 0:
            rs.append(tw / tt)
    spread = (max(rs) - min(rs)) / statistics.median(rs) if len(rs) > 1 else float("nan")
    return rate, spread, time


def geomean(xs):
    xs = [x for x in xs if x and x > 0 and math.isfinite(x)]
    return math.exp(sum(math.log(x) for x in xs) / len(xs)) if xs else float("nan")


def cmd_run(a):
    src = os.path.abspath(a.src)
    build_dir = os.path.join(a.results, ".build", os.path.basename(src))
    print("Building %s ..." % src, flush=True)
    exe, cc_version = build(src, build_dir, a.jobs, a.plain_build)
    listing = json.loads(sh([exe, "--list"]))
    names = [b["name"] for b in listing]
    if a.only:
        wanted = a.only.split(",")
        unknown = set(wanted) - set(names)
        if unknown:
            raise SystemExit("unknown benchmark(s): %s" % ", ".join(sorted(unknown)))
        names = [n for n in names if n in wanted]
    info = {b["name"]: b for b in listing}

    result = {
        "timestamp": datetime.datetime.now().isoformat(timespec="seconds"),
        "label": a.label,
        "src": src,
        "git": git_info(src),
        "machine": {"node": platform.node(), "cpu": cpu_model(), "python": platform.python_version(), "cc": cc_version},
        "settings": {"scale": a.scale, "repeat": a.repeat, "cpu": a.cpu,
                     "build_flags": "stock" if a.plain_build else ALIGN_FLAGS},
        "benchmarks": {},
    }
    total_t = 0.0
    for name in names:
        if not info[name]["available"]:
            result["benchmarks"][name] = {"unit": info[name]["unit"], "skipped": "not available on this CPU"}
            print("  %-16s skipped (not available on this CPU)" % name)
            continue
        meas = max(1, int(round(info[name]["meas_chunks"] * a.scale)))
        warm = info[name]["warm_chunks"]
        reps = []
        for _ in range(a.repeat):
            out = run_json(bench_cmd(exe, a.cpu, "--fixed", name, "--warm-chunks", str(warm), "--meas-chunks", str(meas)))
            if "error" in out:
                raise SystemExit("%s failed: %s" % (name, out["error"]))
            rate, spread, t = summarize_fixed(out["chunks"])
            total_t += t + out["warmup_s"] + out["setup_s"]
            reps.append({"rate": rate, "spread": spread, "measured_s": t, "warmup_s": out["warmup_s"],
                         "setup_s": out["setup_s"], "signature": out["signature"], "chunks": out["chunks"]})
        rate = statistics.median([r["rate"] for r in reps])
        spread = statistics.median([r["spread"] for r in reps])
        if a.repeat > 1:
            rr = [r["rate"] for r in reps]
            spread = max(spread, (max(rr) - min(rr)) / rate)
        sigs = {r["signature"] for r in reps}
        result["benchmarks"][name] = {
            "unit": info[name]["unit"], "usecase": info[name]["usecase"], "rate": rate, "spread": spread,
            "warm_chunks": warm, "meas_chunks": meas, "repeats": reps,
            "signature": reps[0]["signature"] if len(sigs) == 1 else None,
        }
        note = "" if len(sigs) == 1 else "  WARNING: non-deterministic result signature"
        print("  %-16s %12.4g %-9s %5.2f s%s" % (name, rate, info[name]["unit"], reps[0]["measured_s"], note), flush=True)
    print("  total benchmark time: %.1f s" % total_t)

    os.makedirs(a.results, exist_ok=True)
    tag = result["git"]["commit"] + ("-dirty" if result["git"]["dirty_src"] else "")
    fname = "%s_%s%s.json" % (datetime.datetime.now().strftime("%Y%m%d-%H%M%S"), tag, ("_" + a.label) if a.label else "")
    path = os.path.join(a.results, fname)
    with open(path, "w") as f:
        json.dump(result, f, indent=1)
    print("\nSaved %s" % path)

    baseline = a.baseline or os.path.join(a.results, "baseline.json")
    if os.path.exists(baseline) and os.path.abspath(baseline) != os.path.abspath(path):
        print()
        compare(load(baseline), result, baseline, path, a.threshold)
    elif not os.path.exists(baseline):
        print("No baseline yet. Make this run the baseline with:\n  python3 %s set-baseline %s" % (os.path.relpath(__file__), path))


def _new_result(src, label, a, cc_version):
    return {
        "timestamp": datetime.datetime.now().isoformat(timespec="seconds"),
        "label": label, "src": src, "git": git_info(src),
        "machine": {"node": platform.node(), "cpu": cpu_model(), "python": platform.python_version(), "cc": cc_version},
        "settings": {"scale": a.scale, "repeat": a.repeat, "cpu": a.cpu, "mode": "ab-interleaved",
                     "build_flags": "stock" if a.plain_build else ALIGN_FLAGS},
        "benchmarks": {},
    }


def _save(a, result):
    os.makedirs(a.results, exist_ok=True)
    tag = result["git"]["commit"] + ("-dirty" if result["git"]["dirty_src"] else "")
    fname = "%s_%s%s.json" % (datetime.datetime.now().strftime("%Y%m%d-%H%M%S"), tag, ("_" + result["label"]) if result["label"] else "")
    path = os.path.join(a.results, fname)
    with open(path, "w") as f:
        json.dump(result, f, indent=1)
    return path


def cmd_ab(a):
    """Interleaved A/B comparison: for every benchmark, alternate runs of build A and build B
    (A B A B ...), so slow drifts of the machine affect both equally."""
    srcs = [os.path.abspath(a.a), os.path.abspath(a.b)]
    labels = [a.label_a, a.label_b]
    exes, results = [], []
    for src, label in zip(srcs, labels):
        print("Building %s ..." % src, flush=True)
        exe, cc_version = build(src, os.path.join(a.results, ".build", os.path.basename(src)), a.jobs, a.plain_build)
        exes.append(exe)
        results.append(_new_result(src, label, a, cc_version))
    listing = json.loads(sh([exes[1], "--list"]))
    info = {b["name"]: b for b in listing}
    names = [b["name"] for b in listing]
    if a.only:
        wanted = a.only.split(",")
        names = [n for n in names if n in wanted]
    t0 = time.time()
    for name in names:
        if not info[name]["available"]:
            for res in results:
                res["benchmarks"][name] = {"unit": info[name]["unit"], "skipped": "not available on this CPU"}
            continue
        meas = max(1, int(round(info[name]["meas_chunks"] * a.scale)))
        warm = info[name]["warm_chunks"]
        reps = [[], []]
        for k in range(a.repeat):
            order = (0, 1) if k % 2 == 0 else (1, 0)  # ABBA... balances any within-pair trend
            for which in order:
                out = run_json(bench_cmd(exes[which], a.cpu, "--fixed", name, "--warm-chunks", str(warm), "--meas-chunks", str(meas)))
                if "error" in out:
                    raise SystemExit("%s (%s) failed: %s" % (name, labels[which], out["error"]))
                rate, spread, t = summarize_fixed(out["chunks"])
                reps[which].append({"rate": rate, "spread": spread, "measured_s": t, "warmup_s": out["warmup_s"],
                                    "setup_s": out["setup_s"], "signature": out["signature"], "chunks": out["chunks"]})
        for which in (0, 1):
            rr = reps[which]
            rate = statistics.median([r["rate"] for r in rr])
            sigs = {r["signature"] for r in rr}
            results[which]["benchmarks"][name] = {
                "unit": info[name]["unit"], "usecase": info[name]["usecase"], "rate": rate,
                "spread": statistics.median([r["spread"] for r in rr]),
                "warm_chunks": warm, "meas_chunks": meas, "repeats": rr,
                "signature": rr[0]["signature"] if len(sigs) == 1 else None,
            }
        print("  %-16s done" % name, flush=True)
    print("  wall time %.1f s" % (time.time() - t0))
    paths = [_save(a, res) for res in results]
    print("Saved %s\n      %s\n" % tuple(paths))
    compare(results[0], results[1], paths[0], paths[1], a.threshold)


def cpu_model():
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if line.startswith("model name"):
                    return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor()


def load(path):
    with open(path) as f:
        return json.load(f)


def best_chunks(bench):
    """Per-chunk [work, time] with the median time over repeats.

    The median (not the minimum) keeps runs with different --repeat counts comparable:
    a best-of-3 baseline would look systematically faster than a single run.
    """
    reps = [r["chunks"] for r in bench.get("repeats", []) if r.get("chunks")]
    if not reps or any(len(c) != len(reps[0]) for c in reps):
        return None
    return [[reps[0][k][0], statistics.median(c[k][1] for c in reps)] for k in range(len(reps[0]))]


def paired_speedup(bb, nb):
    """Median over chunks of (new chunk rate / baseline chunk rate), plus an uncertainty estimate.

    Both runs execute the same deterministic sequence of chunks, so chunk k does the same
    work in both; the median ignores chunks hit by a slow episode (other processes, VM steal).
    Falls back to the ratio of total rates when chunk sequences do not match.
    """
    cb, cn = best_chunks(bb), best_chunks(nb)
    if not cb or not cn or len(cb) != len(cn):
        return nb["rate"] / bb["rate"], float("nan"), "total"
    ratios = sorted((wn / tn) / (wb / tb) for (wb, tb), (wn, tn) in zip(cb, cn) if tb > 0 and tn > 0 and wb > 0 and wn > 0)
    if not ratios:
        return nb["rate"] / bb["rate"], float("nan"), "total"
    med = statistics.median(ratios)
    if len(ratios) >= 4:
        q = statistics.quantiles(ratios, n=4)
        err = 1.25 * (q[2] - q[0]) / 1.349 / math.sqrt(len(ratios))  # std. error of the median
    else:
        err = float("nan")
    return med, err, "paired"


def compare(base, new, base_name="baseline", new_name="new", threshold=0.03):
    """Print per-benchmark speedups and the geometric-mean score."""
    print("Baseline: %s (%s %s)" % (os.path.basename(base_name), base["git"]["commit"], base.get("label") or ""))
    print("New:      %s (%s%s %s)" % (os.path.basename(new_name), new["git"]["commit"], "-dirty" if new["git"].get("dirty_src") else "", new.get("label") or ""))
    if base.get("settings", {}).get("build_flags") != new.get("settings", {}).get("build_flags"):
        print("WARNING: different build flags (%s vs %s)" % (base.get("settings", {}).get("build_flags"), new.get("settings", {}).get("build_flags")))
    if base["machine"].get("cpu") != new["machine"].get("cpu"):
        print("WARNING: different CPUs (%s vs %s)" % (base["machine"].get("cpu"), new["machine"].get("cpu")))
    print()
    hdr = "%-16s %-9s %11s %11s %9s %8s %7s  %s" % ("benchmark", "unit", "baseline", "new", "speedup", "change", "+/-", "results")
    print(hdr)
    print("-" * len(hdr))
    ratios = []
    for name, nb in new["benchmarks"].items():
        bb = base["benchmarks"].get(name)
        if not bb or "rate" not in bb or "rate" not in nb:
            print("%-16s %-9s %11s %11s %9s" % (name, nb.get("unit", ""), "-", "%.4g" % nb["rate"] if "rate" in nb else "skipped", "-"))
            continue
        ratio, err, how = paired_speedup(bb, nb)
        ratios.append(ratio)
        tol = max(threshold, 3 * err) if math.isfinite(err) else max(threshold, 0.05)
        flag = "  "
        if ratio > 1 + tol:
            flag = " +"
        elif ratio < 1 - tol:
            flag = " -"
        if nb.get("signature") is None or bb.get("signature") is None:
            chk = "n/a"
        elif nb["signature"] == bb["signature"]:
            chk = "identical"
        else:
            den = abs(bb["signature"]) or 1.0
            chk = "CHANGED (rel %.1e)" % (abs(nb["signature"] - bb["signature"]) / den)
        print("%-16s %-9s %11.4g %11.4g %8.3fx %+7.1f%%%s %5.1f%%  %s%s" % (
            name, nb["unit"], bb["rate"], nb["rate"], ratio, 100 * (ratio - 1), flag,
            100 * err if math.isfinite(err) else float("nan"), chk, "" if how == "paired" else " (total-rate ratio)"))
    print("-" * len(hdr))
    gm = geomean(ratios)
    print("GEOMEAN speedup over %d benchmarks: %.3fx (%+.1f%%)" % (len(ratios), gm, 100 * (gm - 1)))
    print("speedup = median of paired per-chunk rate ratios; +/- = its standard error;")
    print("+/- flags: change beyond max(%.0f%%, 3 standard errors); results = fixed-work state signature vs baseline" % (100 * threshold))
    return gm


def cmd_compare(a):
    compare(load(a.base), load(a.new), a.base, a.new, a.threshold)


def cmd_set_baseline(a):
    os.makedirs(a.results, exist_ok=True)
    dst = os.path.join(a.results, "baseline.json")
    shutil.copyfile(a.file, dst)
    print("Baseline set: %s -> %s" % (a.file, dst))


def cmd_probe(a):
    """Run benchmarks for a long time and show how quickly the rate settles."""
    src = os.path.abspath(a.src)
    exe, _ = build(src, os.path.join(a.results, ".build", os.path.basename(src)), a.jobs, a.plain_build)
    listing = json.loads(sh([exe, "--list"]))
    names = a.names or [b["name"] for b in listing if b["available"]]
    print("%-16s %9s %10s  %s" % ("benchmark", "setup_s", "final", "rate relative to final median, per 0.5 s bin (first bin = 0-0.5 s)"))
    for name in names:
        out = run_json(bench_cmd(exe, a.cpu, "--run", name, "--time", str(a.time), "--window", "0.05"), timeout=10 * a.time + 30)
        if "error" in out:
            print("%-16s ERROR %s" % (name, out["error"]), flush=True)
            continue
        win = out["windows"]
        final = statistics.median([r for t, r in win if t > a.time / 2])
        bins = []
        for k in range(int(a.time / 0.5)):
            rs = [r for t, r in win if k * 0.5 < t <= (k + 1) * 0.5]
            bins.append("%.2f" % (statistics.median(rs) / final) if rs else " -- ")
        print("%-16s %9.3f %10.4g  %s" % (name, out["setup_s"], final, " ".join(bins)), flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results", default=DEFAULT_RESULTS, help="results directory (default: %(default)s)")
    p.add_argument("--threshold", type=float, default=0.03,
                   help="per-benchmark change flagged as real above this fraction (A/A runs stay within ~3%%; geomean within ~0.7%%)")
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="build and run benchmarks")
    r.add_argument("--src", default=DEFAULT_SRC, help="REBOUND source tree to benchmark (default: this checkout)")
    r.add_argument("--only", help="comma-separated benchmark names")
    r.add_argument("--repeat", type=int, default=1, help="independent runs per benchmark; median is reported")
    r.add_argument("--scale", type=float, default=1.0, help="multiply the measured work per benchmark (default ~0.5 s each)")
    r.add_argument("--cpu", type=int, default=None, help="pin to this CPU core with taskset")
    r.add_argument("--label", default="", help="free-form label stored with the results")
    r.add_argument("--plain-build", action="store_true", help="stock Makefile flags (no alignment flags; layout-noise prone)")
    r.add_argument("--baseline", help="baseline JSON (default: <results>/baseline.json)")
    r.add_argument("-j", "--jobs", type=int, default=os.cpu_count() or 4)
    r.set_defaults(func=cmd_run)

    c = sub.add_parser("compare", help="compare two result files")
    c.add_argument("base")
    c.add_argument("new")
    c.set_defaults(func=cmd_compare)

    ab = sub.add_parser("ab", help="interleaved A/B comparison of two source trees")
    ab.add_argument("--a", required=True, help="reference source tree (A)")
    ab.add_argument("--b", default=DEFAULT_SRC, help="candidate source tree (B, default: this checkout)")
    ab.add_argument("--label-a", default="A")
    ab.add_argument("--label-b", default="B")
    ab.add_argument("--only")
    ab.add_argument("--repeat", type=int, default=3, help="runs per build and benchmark (alternating)")
    ab.add_argument("--scale", type=float, default=1.0)
    ab.add_argument("--cpu", type=int, default=None)
    ab.add_argument("--plain-build", action="store_true")
    ab.add_argument("-j", "--jobs", type=int, default=os.cpu_count() or 4)
    ab.set_defaults(func=cmd_ab)

    s = sub.add_parser("set-baseline", help="make a result file the baseline")
    s.add_argument("file")
    s.set_defaults(func=cmd_set_baseline)

    pr = sub.add_parser("probe", help="long runs to check how fast rates settle")
    pr.add_argument("names", nargs="*")
    pr.add_argument("--src", default=DEFAULT_SRC)
    pr.add_argument("--time", type=float, default=8.0)
    pr.add_argument("--cpu", type=int, default=None)
    pr.add_argument("--plain-build", action="store_true")
    pr.add_argument("-j", "--jobs", type=int, default=os.cpu_count() or 4)
    pr.set_defaults(func=cmd_probe)

    a = p.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
