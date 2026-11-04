# Copyright (c) 2026 Michael Hippke. Part of Parachron (https://github.com/hippke/parachron), a modified version
# of REBOUND. Licensed under the GNU General Public License v3 or later; see LICENSE.
"""
Core-scaling benchmark of Parachron (WHFast512 and MEP) for the Solar System (Sun + 8 planets, dt = 6 d, GR,
17th-order corrector). Made for testers: one command, no options needed, one log file to send back.

    python3 bench_cores.py            # about 5-10 minutes
    python3 bench_cores.py --quick    # about 2 minutes, for a first check

Measures, for 1 core and then 2, 4, 8, ... up to the number of physical cores (plus the number of performance cores
on CPUs with performance and efficiency cores):
  1. single core: ns per WHFast512 step of one simulation (the speed of one core);
  2. throughput: n independent simulations at the same time, one per process (speed per core under load);
  3. MEP: one simulation on n cores, wall clock and speed-up against the serial run of the same length.
Writes everything to the screen and to parachron_bench_<computer>_<date>.txt in the current directory.
"""
import argparse
import math
import multiprocessing as mp
import os
import platform
import statistics
import subprocess
import sys
import time
import traceback


def sh(cmd):
    try:
        return subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=20).stdout.strip()
    except Exception:
        return ""


def sysctl(name):
    return sh("sysctl -n " + name)


def ps(cmd):
    return sh('powershell -NoProfile -Command "%s"' % cmd)


def loadavg():
    return " ".join("%.2f" % x for x in os.getloadavg()) if hasattr(os, "getloadavg") else "n/a"


class Log:
    def __init__(self, fn):
        self.f = open(fn, "w")

    def __call__(self, *a):
        s = " ".join(str(x) for x in a)
        print(s, flush=True)
        self.f.write(s + "\n")
        self.f.flush()


def physical_cores():
    """(physical cores, performance cores or None)"""
    if sys.platform == "darwin":
        n = int(sysctl("hw.physicalcpu") or os.cpu_count())
        p = sysctl("hw.perflevel0.physicalcpu")
        return n, (int(p) if p and int(sysctl("hw.nperflevels") or 1) > 1 else None)
    if sys.platform == "win32":
        n = ps("(Get-CimInstance Win32_Processor | Measure-Object -Property NumberOfCores -Sum).Sum")
        return (int(n) if n.isdigit() else os.cpu_count()), None
    cores = set()
    for line in sh("lscpu -p=core,socket").splitlines():
        if line and not line.startswith("#"):
            cores.add(line)
    return (len(cores) or os.cpu_count()), None


def core_counts(n, p):
    c = [1]
    while c[-1] * 2 <= n:
        c.append(c[-1] * 2)
    for extra in (p, n):
        if extra and extra not in c:
            c.append(extra)
    return sorted(c)


def machine(log):
    import parachron
    log("=== machine")
    log("date", time.strftime("%Y-%m-%d %H:%M %Z"))
    log("platform", platform.platform(), platform.machine())
    log("python", sys.version.split()[0], sys.executable)
    try:
        from importlib.metadata import version
        v = version("parachron")
    except Exception:
        v = "?"
    log("parachron", v, "built", getattr(parachron, "__build__", ""), parachron.__file__)
    if sys.platform == "darwin":
        log("cpu", sysctl("machdep.cpu.brand_string"))
        log("model", sysctl("hw.model"))
        log("macOS", sh("sw_vers -productVersion"), sh("sw_vers -buildVersion"))
        log("cores physical/logical", sysctl("hw.physicalcpu"), sysctl("hw.logicalcpu"))
        for i in range(int(sysctl("hw.nperflevels") or 0)):
            log("  perflevel%d" % i, sysctl("hw.perflevel%d.name" % i), sysctl("hw.perflevel%d.physicalcpu" % i),
                "cores, L2", sysctl("hw.perflevel%d.l2cachesize" % i))
        log("memory GB", round(int(sysctl("hw.memsize") or 0) / 2**30, 1))
        if sysctl("sysctl.proc_translated") == "1":
            log("WARNING: Python runs under Rosetta (x86 emulation). Use an arm64 Python; see the instructions.")
        log("power", sh("pmset -g batt | head -1"))
        log("low power mode", sh("pmset -g | grep -i lowpowermode") or "n/a")
    elif sys.platform == "win32":
        log("cpu", ps("(Get-CimInstance Win32_Processor).Name"))
        log("cores physical/logical", ps("(Get-CimInstance Win32_Processor | Measure-Object -Property NumberOfCores -Sum).Sum"),
            ps("(Get-CimInstance Win32_Processor | Measure-Object -Property NumberOfLogicalProcessors -Sum).Sum"))
        log("model", ps("(Get-CimInstance Win32_ComputerSystem).Model"))
        log("memory GB", round(int(ps("(Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory") or 0) / 2**30, 1))
        log("power plan", sh("powercfg /getactivescheme"))
        log("on battery", ps("(Get-CimInstance Win32_Battery).BatteryStatus") or "no battery")
    else:
        log("cpu", sh("grep -m1 'model name' /proc/cpuinfo | cut -d: -f2"))
        log(sh("lscpu | grep -E '^(Core|Thread|Socket|CPU\\(s\\)|CPU max MHz)'"))
    if sys.platform != "win32":
        log("compiler", sh("${CC:-cc} --version | head -1"))
    log("load average", loadavg())


def solar_system():
    import parachron
    sim = parachron.Simulation()
    sim.add("solar system")
    sim.integrator = "whfast512"
    sim.dt = 6. / 365.25 * 2. * math.pi
    sim.exact_finish_time = False
    sim.integrator.gr_potential = 1
    sim.integrator.corrector = 17
    return sim


def time_serial(steps):
    """Seconds for `steps` WHFast512 steps of one simulation."""
    sim = solar_system()
    sim.integrate(100 * sim.dt, exact_finish_time=0)              # initialisation outside the timing
    t0 = time.perf_counter()
    sim.integrate(sim.t + steps * sim.dt, exact_finish_time=0)
    return time.perf_counter() - t0


def worker(steps, barrier, q):
    sim = solar_system()
    sim.integrate(100 * sim.dt, exact_finish_time=0)
    barrier.wait()
    t0 = time.perf_counter()
    sim.integrate(sim.t + steps * sim.dt, exact_finish_time=0)
    q.put(time.perf_counter() - t0)


def throughput(n, steps):
    ctx = mp.get_context("spawn")
    barrier, q = ctx.Barrier(n), ctx.Queue()
    ps = [ctx.Process(target=worker, args=(steps, barrier, q)) for _ in range(n)]
    for p in ps:
        p.start()
    dts = [q.get() for _ in ps]
    for p in ps:
        p.join()
    return dts


def mep_run(cores, slices):
    import parachron
    sim = solar_system()
    m = parachron.MEP(sim, cores=cores)
    windows = max(1, math.ceil(slices / cores))
    t0 = time.perf_counter()
    m.integrate(windows * m.window_length)
    wall = time.perf_counter() - t0
    return wall, round(m.t / sim.dt), m.mean_iterations


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true", help="shorter runs (about 2 minutes)")
    ap.add_argument("--max-cores", type=int, default=None, help="test up to this many cores")
    args = ap.parse_args()
    host = (platform.node() or "computer").split(".")[0].replace(" ", "_")
    fn = "parachron_bench_%s_%s.txt" % (host, time.strftime("%Y%m%d_%H%M"))
    log = Log(fn)
    log("Parachron core-scaling benchmark", "(quick)" if args.quick else "")
    try:
        machine(log)
    except Exception:
        log(traceback.format_exc())
        log("Parachron is not installed correctly (see above). Stopping.")
        return
    n, p = physical_cores()
    if args.max_cores:
        n = min(n, args.max_cores)
    counts = core_counts(n, p if p and p < n else None)
    log("physical cores", n, "| performance cores", p if p else "n/a", "| core counts tested", counts)
    serial_steps = 3_000_000 if args.quick else 20_000_000
    tp_steps = 3_000_000 if args.quick else 30_000_000
    slices = 200 if args.quick else 1600

    log("\n=== 1. single core: one simulation, ns per step (5 runs of %d steps)" % serial_steps)
    sim = solar_system()
    sim.steps(1)
    log("kernel", sim.integrator.kernel_used)
    dts = [time_serial(serial_steps) for _ in range(5)]
    ns = [1e9 * d / serial_steps for d in dts]
    serial_ns = statistics.median(ns)
    log("ns/step runs", " ".join("%.1f" % x for x in ns))
    log("ns/step median %.2f  min %.2f" % (serial_ns, min(ns)))

    log("\n=== 2. throughput: n independent simulations at the same time (%d steps each)" % tp_steps)
    log("%6s %14s %14s %16s %12s" % ("cores", "ns/step mean", "ns/step worst", "steps/s total", "efficiency"))
    base = None
    for c in counts:
        try:
            d = throughput(c, tp_steps)
            mean = 1e9 * statistics.mean(d) / tp_steps
            worst = 1e9 * max(d) / tp_steps
            rate = c * tp_steps / max(d)
            base = base or rate
            log("%6d %14.2f %14.2f %16.4g %11.0f%%" % (c, mean, worst, rate, 100 * rate / (c * base)))
        except Exception:
            log("cores %d failed:\n%s" % (c, traceback.format_exc()))

    log("\n=== 3. MEP: one simulation on n cores, %d slices of 60000 steps; up and down (thermal drift check)" % slices)
    log("serial reference: %.2f ns/step (section 1)" % serial_ns)
    log("%6s %6s %12s %10s %10s %12s" % ("cores", "pass", "steps", "wall s", "speed-up", "iterations"))
    res = {}
    mep_counts = [c for c in counts if c >= 2]
    for pas, order in (("up", mep_counts), ("down", mep_counts[::-1])):
        for c in order:
            try:
                wall, steps, it = mep_run(c, slices)
                s = serial_ns * 1e-9 * steps / wall
                res.setdefault(c, []).append(s)
                log("%6d %6s %12d %10.2f %9.2fx %12.2f" % (c, pas, steps, wall, s, it))
            except Exception:
                log("cores %d failed:\n%s" % (c, traceback.format_exc()))

    log("\n=== summary")
    log("single core: %.2f ns/step" % serial_ns)
    for c in mep_counts:
        if c in res:
            log("MEP %3d cores: speed-up %.2fx (mean of %d passes)" % (c, statistics.mean(res[c]), len(res[c])))
    log("load average at the end", loadavg())
    log("\nDone. Please send the file %s" % os.path.abspath(fn))


if __name__ == "__main__":
    main()
