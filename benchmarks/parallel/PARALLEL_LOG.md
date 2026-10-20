# REBOUND multicore work log

Branch `parallel-multicore` of https://github.com/hippke/rebound, based on `main` (`62d6d93a`, the end of the
single-core `perf-bench` work). Goal: make IAS15 and/or WHFast(X) finish sooner when given several cores,
i.e. the wall time of the OpenMP build at T threads must be **below the wall time of the stock serial build**
(not merely below the OpenMP build at 1 thread).

## Summary (end of this round, commit `f0fe887`, measured on an Intel i5-13500)

| | stock OpenMP build (2 threads) | this branch, best thread count (P-cores only) |
|---|---|---|
| small systems (ias15_hltau N = 6, whfast_ss 8 planets, embryos N = 61) | 0.08–0.23x | 1.00x (never slower than serial) |
| IAS15 N = 128 / 512 / 2048, default mode | 0.45 / 0.47 / 0.70x | 1.01 / 3.15 / 3.85x |
| IAS15 `parallel_stages = 1` (vs default serial IAS15), N = 61 / 128 / 512 / 2048 | – | 1.89 / 2.89 / 4.51 / 5.69x |
| WHFast, 2 planets + 10,000 test particles | 0.89x | 2.28x |
| WHFast, 256 massive bodies | 0.42x | 1.93x (T = 4) |
| results | differ from serial (up to 5e-2 relative in whfast_ss) | bitwise identical to the serial build at every thread count |

Also fixed on the way: the OpenMP WHFast Jacobi μ bug, TRACE's OpenMP encounter-matrix indexing,
a data race (EOS jerk loop), thread-count-dependent reductions; the OpenMP test suite went from > 15
minutes (stopped) to 57 s, 956 tests passing. Not tested in parallel: WHFast512 (no OpenMP code paths,
~150 ns per step for 8 planets: expected 1.00x), MERCURIUS/TRACE performance (tests only), compensated
gravity performance.

## Process

Same loop as `PERF_LOG.md` of the single-core work, adapted to threads:

1. **Implement** one work item, commit.
2. **Correctness:**
   - `benchmarks/perf/check.py` (serial build, 953 tests) must pass: the default build users get must not change.
   - `REBOUND_OPENMP=1 OMP_NUM_THREADS=2 python3 benchmarks/perf/check.py` (OpenMP build, new switch in
     `setup.py`): every test that passes serially must pass with OpenMP (also with 3 threads).
   - `benchmarks/parallel/gravity_check.c`: bit-exact hash of all accelerations for 189 configurations
     (basic and compensated gravity, N = 3…1600, no / type-0 / type-1 test particles, all
     `gravity_ignore_terms`, ghost boxes with softening), serial library vs OpenMP library at T = 1, 2, 3, 5.
   - `run_parallel.py` signatures: results must be **independent of the thread count** (bitwise), and equal
     to the serial build where the change claims it.
3. **Performance:** `python3 benchmarks/parallel/run_parallel.py --repeat 5` (interleaved serial/omp1/omp2);
   for the serial build, `benchmarks/perf/run.py ab` against the base commit must show no regression.
4. **Decide:** keep if a targeted case gets >5% faster at T>1 relative to the serial build and nothing else
   regresses beyond noise; keep correctness fixes and small, safe changes; otherwise revert.
5. **Log** here.

Metric: speedup = t(serial build) / t(OpenMP build, T threads). Same fixed work in both.

Machine for this log: cloud VM, 2 vCPUs (Intel Xeon @ 2.8 GHz, AVX-512, 1 thread per core), gcc 13.3.
Two caveats about this machine, both measured:
- A parallel region (fork + join) costs ~1 µs with 2 threads (270 ns even with one thread). A throughput-
  bound kernel runs at ~1.8x on 2 threads (two copies of the serial gravity kernel side by side: 5.2 and
  6.0 ns/pair instead of 4.9).
- Noise: VM steal comes and goes; a 2-thread kernel with barriers suffers more from it than a serial one
  (the other thread waits at every barrier). Serial reference timings of the same binary varied by up to
  1.5x between runs an hour apart. Kernel comparisons were therefore made with interleaved repeats
  (min and median), and conclusions near 1.0x are not reliable here.
**Scaling beyond 2 cores could not be measured.** Run `python3 benchmarks/parallel/run_parallel.py` on a
workstation: it reports 1, 2, 4, ... nproc threads.

## Literature consulted (ideas tested in code are marked)

- Saha, Stadel & Tremaine (1997), parallel-in-time Picard iteration for planetary orbits: large speedups
  only with very many processors, at a much higher total cost. Here the cost has to stay near 1x, so only
  the cheap end of the idea is used: parallelism across the stages of one step (below), not across steps.
- **van der Houwen & Sommeijer, parallel iterated Runge-Kutta (PIRK) methods** (CWI,
  https://ir.cwi.nl/pub/1663/1663D.pdf): iterate an implicit Runge-Kutta method by Jacobi instead of
  Gauss-Seidel sweeps over the stages, so all s stage derivatives of one iteration are independent; the
  speedup approaches s when the right-hand side dominates. IAS15 is a Gauss-Radau collocation method with
  7 substeps iterated to convergence, i.e. a natural candidate. **Tested: iteration 6.**
- Stage-parallel implicit Runge-Kutta, arXiv:2505.17719: recent analysis of the same idea for IRK methods.
- Parallel direct N-body for planetary systems, arXiv:2208.13562: for small n, serial beats parallel
  because of synchronisation costs (agrees with iteration 0 and the thresholds of iteration 2).
- parallel-in-time.org (Parareal, PFASST): across-step parallelism; not tried here (needs many cores,
  and the coarse propagator of a chaotic system loses accuracy quickly).
- GENGA, GLISSE, Swiftest (GPU/OpenMP N-body), WHFast512 (Javaheri, Rein & Tamayo 2023: SIMD inside one
  step instead of threads), Jiménez-Pérez & Laskar (parallel symplectic integration): parallelism over
  particles or ensembles; the particle-level ideas are items P3-P6.

## Iteration 0 — baseline (commit `b25fbbd`)

New harness `benchmarks/parallel/` (`pbench.c`, `run_parallel.py`). Stock OpenMP build vs stock serial build:

| case | N | serial [s] | omp 1 thread | omp 2 threads | threads | vs serial |
|---|---|---|---|---|---|---|
| ias15_hltau | 6 | 0.471 | 0.33x | 0.08x | identical | 0 |
| ias15_embryos | 61 | 0.491 | 0.34x | 0.23x | identical | 0 |
| ias15_n128 | 128 | 0.511 | 0.44x | 0.45x | identical | 0 |
| ias15_n512 | 512 | 0.659 | 0.35x | 0.47x | identical | 0 |
| ias15_n2048 | 2048 | 1.478 | 0.43x | 0.70x | identical | 0 |
| whfast_ss | 9 | 0.262 | 0.48x | 0.13x | identical | 5e-2 |
| whfast_tp10k | 10,003 (3 active) | 0.420 | 0.62x | 0.89x | identical | 5e-5 |
| whfast_n256 | 256 | 0.674 | 0.33x | 0.42x | identical | 2e-4 |

**Using more cores made every case slower.** Causes, from reading the `#ifdef OPENMP` code paths:

1. The OpenMP build gave up the single-core optimisations: `REB_TARGET_CLONES` (AVX2 kernels) was disabled
   with OpenMP, the batched Kepler solver was replaced by scalar loops, and basic gravity switched from the
   half-pair (Newton's third law) loop to a full N² loop: twice the square roots and divisions.
2. Every small O(N) loop was a separate `#pragma omp parallel for` (fork/join each time), also for N = 6.
3. **Bug (also in upstream `hannorein/rebound`):** the OpenMP branch of the WHFast Jacobi Kepler step used
   μ = G(m₀ + Σ_{j<i} m_j), leaving out the planet's own mass; the serial branch includes it.
4. The stock OpenMP test-suite run took over 15 minutes (54 s serially) and had to be stopped.

## Iteration 1 — WHFast Jacobi μ (commit `18accf9`) — kept (correctness)

`integrator_whfast.c`: `j < MIN(i+1, N_active)`. OpenMP and serial WHFast now agree. Also `setup.py`:
`REBOUND_OPENMP=1` builds the Python extension with OpenMP so the test suite can run against it.

## Iteration 2 — serial code below a work threshold; deterministic reductions (commit `7c0788f`) — kept

- `reb_omp_parallel(work)` in `rebound_internal.h`: a loop runs in parallel only with > 1 thread and at least
  `REB_OMP_MIN_WORK` (2048) units of work; otherwise the serial code path runs, which is now compiled in
  OpenMP builds too (basic and compensated gravity, WHFast Kepler step and drifts). `if()` clauses on the
  remaining loops (transformations, boundary, collision trees, leapfrog, SEI, MERCURIUS, TRACE).
- `REB_TARGET_CLONES` also in OpenMP builds (the outlined parallel-region bodies use the default target).
- Thread-count independence: removed the OpenMP reductions in the WHFast jump step and the COM sums of
  `transformations.c` (order depended on the thread count; O(N_active)); `reb_simulation_two_largest_particles`
  is serial (its per-thread reduction could pick different indices on ties: all embryos have equal radii);
  removed the pragma of the EOS jerk loop (data race on `particles[j]`, code review item).
- Result (`p2.json`): omp1 0.91–1.12x everywhere, all results bitwise identical to the serial build.
  Still slower with 2 threads wherever the old one-sided gravity loop ran (0.31–0.70x).
- Serial build: bitwise identical, `benchmarks/perf` A/B geomean 1.000x (20 benchmarks; megno −5.7% in the
  first run, −1.3% on a 5-repeat recheck: noise).

## Iteration 3 — TRACE/MERCURIUS OpenMP branches (commit `20ed2de`) — kept (correctness)

With iteration 2 the OpenMP test suite ran in 60 s and showed **18 TRACE failures** (energy error 3e-4
instead of < 3e-10). TRACE's OpenMP loops read the encounter matrix at `[j*N+i]` for every pair, also for
j > i, i.e. the lower triangle, which is never set: pairs in a close encounter were not removed from the
kick. Fixed to `[min*N+max]` (also the test-particle back-reaction loops, Kepler mode). Remaining failures
were round-off: the OpenMP loops computed G·m/r³, the serial ones (G/r³)·m; with the serial expressions
the OpenMP branch is bitwise identical. TRACE and MERCURIUS now also take their serial branch below the
threshold (the last failing test was a wall-clock assertion, TRACE faster than IAS15).
**OpenMP build: 953 tests pass with 2 and with 3 threads.**

## Iteration 4 — parallel basic gravity, bitwise identical (commits `5a35d80`, `93d4a96`) — kept

The one-sided loop (each particle sums all partners: bitwise identical, embarrassingly parallel) does N²
pair evaluations instead of N²/2, so on 2 cores it cannot beat the serial loop (measured: 11.0 ns per pair
on one core vs 4.9 ns serial; 0.67–0.74x with 2 threads). Variants tried, all checked bit for bit:

| variant | idea | 1-core cost / serial | 2 threads vs serial | outcome |
|---|---|---|---|---|
| a | one-sided, vectorised in blocks | 2.2x | 0.67–0.74x | kept only as fallback (no shared data) |
| b | table of G/r³ (pass 1), then each particle's sum from the table in serial order (pass 2) | 2.0x | 0.70–0.92x | 56 vs 38 instructions/pair (callgrind) |
| c | b + structure-of-arrays copy, pass 2 vectorised over blocks of 8 particles | 1.4x | 0.98–1.10x | pass 2 latency-bound on misses (triangular row stride) |
| d | c + own-row terms fused into pass 1 | 1.4x | 0.73–0.94x | worse: pass 1 becomes one long dependency chain |
| e | column-panel table layout | 1.5–3.1x | 0.62–0.91x | pass 1 writes scatter over the table |
| f | **row blocks**: per block, off-diagonal rows (G/r³ once per pair, own terms, row-major table) + diagonal tile with the serial half-pair loop; pass 2 streams later rows vectorised over chunks of 64 particles | 1.4–1.6x | up to 1.29x (N = 2048) | **kept** |
| f' | f with vectorised square roots (see below) | 0.7–1.0x (N ≤ 1024), 1.3x (2048) | 1.18–1.22x (N ≥ 1024, min of 9) | **kept** |

Why it is exact: the serial loop adds to particle k its own-row terms (j < k, ascending), then one term per
later row (i > k, ascending), then the test-particle rows. In f, a particle in block b receives the
off-diagonal own-row terms (pass 1a), the diagonal tile's terms in the serial order (pass 1b: the serial
loop restricted to the block), then the later rows from the table in ascending order (pass 2), with the
same floating-point expressions. The block layout and thread count do not enter the arithmetic.
Negative control: the stock OpenMP library differs from the serial one in 61 of 189 configurations of
`gravity_check` (ghost boxes, small N); the serial library of this branch is identical to the base in all 189.

Found on the way (f'): the parallel kernels received `softening²` as a parameter; GCC cannot prove it
non-negative, so `sqrt` keeps its errno branch and the loop does not vectorise ("control flow in loop").
The serial loop squares `r->softening` itself. Squaring inside the kernel (same bits) vectorised all loops.
With that, the table kernel on **one** core costs about as much as the serial loop or less (structure of
arrays: 153 vs 172 µs at N = 256, 550 vs 654 µs at N = 512, 2213 vs 3177 µs at N = 1024), but 2 threads on
2 cores were slower than 2 threads on 1 core for N ≤ 512: pass 2 reads the table rows written by the other
core (1 MB per force evaluation at N = 512), and that transfer is slow on this VM.

Threshold: `reb_omp_parallel_gravity` uses the parallel kernel when N_active² × T ≥ 2 × 2^19 (N ≈ 724 on
2 cores, 362 on 8). Test-particle rows cost the same in parallel and serially, so they use the plain
threshold. Larger tables than 64 MB (N > ~4000) use the one-sided kernel.

## Iteration 5 — parallel batched WHFast Kepler drifts (in commit `5a35d80`) — kept

The OpenMP Kepler step used the scalar solver (1.2–1.5x slower than the batched one) and, for Jacobi
coordinates, recomputed each planet's mass sum from scratch (O(N²)). Now μ is computed once with the serial
running sum (same rounding), and threads drift contiguous chunks through the batched solver (each lane is
bitwise equal to the scalar solver). `reb_integrator_whfast_kepler_solver_range` (MERCURIUS/TRACE drifts)
likewise. Result: **whfast_tp10k 1.61x with 2 threads** (1.38–1.61x across runs), bitwise identical.
Gravity of the test particles: chunks of 64 per call (one AVX2-clone dispatch per chunk instead of per
particle).

## Iteration 6 — IAS15 stage-parallel mode, `ias15.parallel_stages = 1` — kept (opt-in)

Idea from the PIRK literature (above). In each predictor-corrector iteration, IAS15 normally evaluates the
force at substep n with positions predicted from the b coefficients already updated at substeps 1..n-1
(Gauss-Seidel). In stage mode, all seven substep positions are predicted from the b values at the start of
the iteration, and the seven force evaluations are independent (Jacobi). They run as one `omp parallel for`
over the stages (up to 7 threads), each calling the serial half-pair gravity loop on a structure-of-arrays
copy (`reb_gravity_basic_arrays`, the same per-pair expressions as the serial kernel). The b/g/e updates,
convergence test and step-size control are unchanged and serial.

- Exactness: the converged collocation solution is the same; only round-off differs from the default mode
  (1e-13 to 1e-16 relative in the final state, the same number of steps in all cases). With the flag on,
  results are **bitwise identical between the serial build and the OpenMP build at any thread count**
  (checked with pbench signatures at T = 1, 2, 3).
- Jacobi needs about the same number of iterations here: stage mode on one core costs 0.87–1.14x of
  the default (noise on this VM is of the same size).
- Saved in the simulation archive as an IAS15 setting, so restarts reproduce bit for bit. Used only with
  basic gravity, all particles active, no additional forces, no variational particles, no ghost boxes;
  otherwise the flag is silently ignored (test: identical to the default mode).
- New test `rebound/tests/test_ias15_parallel_stages.py`. `check.py`: 956 passed (serial build; one
  flaky `test_server` download failure that passes on rerun) and 956 passed (OpenMP build, 2 threads).
  `gravity_check`: 189/189 identical, serial = base = OpenMP at T = 1, 2, 3, 5.
- `pbench` cases `CASE+stages`; `run_parallel.py` reports them against the default mode of the serial build.

2-vCPU VM, `--repeat 5`, speedup vs the **default** IAS15 of the serial build:

| case | stages, 1 thread | stages, 2 threads |
|---|---|---|
| ias15_embryos (N = 61, collisions) | 0.87x | 1.02x |
| ias15_n128 | 1.09x | **1.71x** |
| ias15_n512 | 0.88x | **1.58x** |
| ias15_n2048 | 1.14x | **2.14x** |

This is the only change that helps IAS15 at N ~ 100: the work per parallel task is a whole force
evaluation instead of a slice of one, so one fork/join per iteration suffices and no data moves between
cores except the 7 small stage buffers. Limits: at most 7-way parallel (combine with the particle-parallel
kernel for more cores: nested parallelism, open); not useful for N ≲ 50, where a force evaluation costs
about as much as the fork/join.

### First workstation measurement (user, Intel Core i5-13500: 6 P-cores + 8 E-cores, 20 threads)

Commit before iteration 6 (no `+stages` cases), speedup vs the serial build:

| case | omp1 | omp2 | omp4 | omp8 | omp16 | omp20 |
|---|---|---|---|---|---|---|
| ias15_hltau | 0.99x | 0.98x | 0.99x | 0.99x | 1.00x | 1.00x |
| ias15_embryos | 1.00x | 1.00x | 1.00x | 1.00x | 1.00x | 1.00x |
| ias15_n128 | 1.01x | 1.01x | 1.01x | 1.01x | 1.01x | 1.01x |
| ias15_n512 | 1.01x | 1.02x | 2.70x | 2.87x | 2.40x | 2.62x |
| ias15_n2048 | 1.02x | 1.71x | 2.99x | 4.34x | **4.69x** | 3.29x |
| whfast_ss | 0.99x | 1.00x | 0.99x | 1.00x | 0.99x | 0.99x |
| whfast_tp10k | 0.99x | 1.49x | 2.06x | 1.99x | 2.25x | 2.29x |
| whfast_n256 | 1.01x | 1.01x | 1.02x | 0.98x | 1.66x | 1.86x |

All bitwise identical to the serial build at every thread count. Readings:
- The goal holds: nothing is slower than serial, and the large cases scale (n2048 4.7x, n512 2.9x).
- Threshold too conservative for a real CPU: n512 at T = 2 and whfast_n256 below T = 16 stay serial
  (N_active² × T < 2^20). The break-even of 2^19 was measured on the VM, whose cross-core transfer is slow.
  To be tuned with `run_parallel.py --cflags="-DREB_OMP_MIN_WORK_GRAVITY=65536."` (new option).
- T = 20 is slower than 16 for gravity: P- and E-cores mixed, and the static row-block split makes the
  P-cores wait at the barrier for the E-cores. More, smaller blocks with dynamic scheduling (item P11).

## Iteration 7 — gravity threshold from workstation data (commit `b301823`) — kept (side effect found in iteration 10: mixed team sizes, item P14)

Second i5-13500 run (stage cases included, and with `--cflags="-DREB_OMP_MIN_WORK_GRAVITY=65536."`):

| case | threshold | omp2 | omp4 | omp8 | omp16 | omp20 |
|---|---|---|---|---|---|---|
| ias15_n128 | 2^19·2/T (old) | 1.01x | 1.02x | 1.01x | 1.02x | 1.02x |
| ias15_n128 | 2^16·2/T | 1.01x | 1.01x | 1.02x | 0.93x | 0.80x |
| ias15_n512 | old | 1.02x | 2.73x | 2.82x | 2.46x | 2.57x |
| ias15_n512 | 2^16·2/T | **1.60x** | 2.60x | 2.87x | 2.50x | 1.70x |
| whfast_n256 | old | 1.03x | 1.02x | 0.99x | 1.59x | 1.79x |
| whfast_n256 | 2^16·2/T | 1.27x | **1.97x** | 1.61x | 1.60x | 1.42x |

The parallel kernel pays off from N ≈ 256 on, at any T ≥ 2; small N with many threads loses (n128 at
T = 16-20), and the best thread count grows with N: about N/64 (n256: 4, n512: 8, n2048: ≥ 16).
New rule: parallel when N_active² ≥ 2^16, independent of T, with at most max(2, N_active/64) threads
(all threads when test-particle rows dominate, e.g. whfast_tp10k). Bitwise results unchanged
(gravity_check 189/189 at T = 1, 2, 3, 5; OpenMP suite 956 passed at T = 3).

## Iteration 8 — IAS15 stage mode in row blocks (commit `921e2cc`) — kept (confirmed on the i5-13500 in iteration 10: 1.14–1.62x at T = 6)

Second i5-13500 run, stage mode vs the default serial IAS15:

| case | stages omp1 | omp2 | omp4 | omp8 | omp16 | omp20 |
|---|---|---|---|---|---|---|
| ias15_embryos+stages | 0.84x | 1.26x | **1.88x** | 1.62x | 1.56x | 1.45x |
| ias15_n128+stages | 0.95x | 1.52x | **2.54x** | 2.08x | 2.23x | 2.24x |
| ias15_n512+stages | 1.00x | 1.69x | **3.22x** | 2.44x | 2.72x | 2.34x |
| ias15_n2048+stages | 1.10x | 1.90x | **3.60x** | 2.62x | 3.05x | 3.58x |

Best at T = 4, worse at 8+. With 7 equal tasks, 4 threads need 2 rounds (≤ 3.5x, reached), 6 threads too,
and 8+ threads (`OMP_PROC_BIND=close`, `OMP_PLACES=cores`: P-cores first) put some stages on efficiency
cores, so every iteration waits for the slowest. Change: each stage is split into K row blocks of the
half-pair loop (K = N/64, at most 16 and 32 MB of scratch: a function of N only), each accumulated into
its own zeroed arrays, then summed per particle in ascending block order. 7K tasks with dynamic
scheduling. K = 1 (N < 128) keeps the previous fused loop and identical results. Results: identical in
the serial build and at T = 1, 2, 3, 5 (pbench signatures), round-off differences to iteration 6 for N ≥ 128.
2-vCPU VM A/B against iteration 7 (`--repeat 7`): 0.92–1.03x, i.e. neutral within noise. Two threads
gain little from balancing 7 tasks. The target is 6+ cores (6 P-cores: 7 tasks → 2 rounds, ≤ 3.5x;
7K tasks → ≤ 6x), to be measured with `--compare-src` on the workstation.

## Iteration 9 — use only performance cores on hybrid CPUs (user suggestion) — kept

User: "if parallelization is not optimal anyways (and scaling to many cores difficult), it seems
reasonable to restrict it to performance-cores". The data agrees: every case whose speedup drops at
T > 6 on the i5-13500 drops when efficiency cores join (barrier-synchronised loops run at the pace of the
slowest core). `reb_omp_threads()` = OpenMP thread count limited to the number of physical
performance cores; used by every REBOUND parallel region (`num_threads` clause on all 38) and by
`reb_omp_parallel`. Detection, once per process:
- Linux, Intel hybrid: `/sys/devices/cpu_core/cpus` (P-cores; physical cores via
  `topology/thread_siblings_list`) when `/sys/devices/cpu_atom/cpus` is non-empty.
- macOS, Apple silicon: `hw.nperflevels > 1` → `hw.perflevel0.physicalcpu`.
- Otherwise (and with `REBOUND_OMP_ALL_CORES=1`): no limit.
Tested on fake sysfs trees: i5-13500 layout → 6, non-hybrid → 0; this VM → 0 (no change). Only REBOUND's
loops are limited, not the user's own OpenMP code; results are independent of the thread count anyway.
For the threads to land on the P-cores, bind them (`OMP_PROC_BIND=close OMP_PLACES=cores`, as
`run_parallel.py` does; Linux usually prefers P-cores for busy threads without binding too).
`run_parallel.py --all-cores` measures the old behaviour for comparison.

## Iteration 10 — workstation measurement of iterations 7–9 (user, i5-13500, commit `f0fe887`) — measurement only

Raw output: `i5-13500_run3_commitf0fe887.txt` (work directory, `bench_results/parallel/`).

Default (P-cores only; T > 6 runs as 6), speedup vs the serial build:

| case | omp1 | omp2 | omp4 | omp8 | omp16 | omp20 |
|---|---|---|---|---|---|---|
| ias15_hltau | 0.97x | 0.96x | 0.99x | 0.99x | 0.99x | 0.98x |
| ias15_embryos | 1.01x | 1.01x | 1.00x | 1.01x | 1.01x | 0.99x |
| ias15_n128 | 1.01x | 1.01x | 1.01x | 1.01x | 1.01x | 1.01x |
| ias15_n512 | 0.99x | 1.62x | 2.65x | 3.15x | 3.11x | 3.13x |
| ias15_n2048 | 0.98x | 1.57x | 2.75x | 3.83x | 3.85x | 3.83x |
| ias15_embryos+stages (vs default serial) | 0.87x | 1.27x | 1.87x | 1.89x | 1.89x | 1.89x |
| ias15_n128+stages (vs default serial) | 0.94x | 1.62x | 2.45x | 2.88x | 2.89x | 2.86x |
| ias15_n512+stages (vs default serial) | 0.95x | 1.82x | 3.37x | 4.49x | 4.51x | 4.50x |
| ias15_n2048+stages (vs default serial) | 1.08x | 2.09x | 3.97x | 5.69x | 5.69x | 5.68x |
| whfast_ss | 1.00x | 0.99x | 1.00x | 1.00x | 1.00x | 0.99x |
| whfast_tp10k | 1.00x | 1.48x | 2.07x | 2.28x | 2.28x | 2.28x |
| whfast_n256 | 0.99x | 1.24x | 1.93x | 1.27x | 1.36x | 1.36x |

`--all-cores` (efficiency cores too), selected rows:

| case | omp4 | omp8 | omp16 | omp20 |
|---|---|---|---|---|
| ias15_n2048 | 2.85x | 3.95x | 4.42x | 4.30x |
| ias15_n512+stages (vs default serial) | 3.36x | 4.24x | 3.26x | 1.95x |
| ias15_n2048+stages (vs default serial) | 3.94x | 6.28x | 7.44x | **8.55x** |
| whfast_tp10k | 1.96x | 2.00x | 2.27x | 1.54x |
| whfast_n256 | 2.00x | 1.08x | **0.04x** | **0.03x** |

A/B at T = 6, iteration 9 tree vs iteration 6 tree (`--compare-src`, stage cases): n128 1.14x, n512
1.46x, n2048 1.62x faster. At T = 6 the P-core limit does not act and the gravity threshold does not
enter the stage mode, so this is the effect of the row blocks (iteration 8): **confirmed, kept.**

Readings:
- Iteration 9 (P-cores only) does what it should: no case gets slower from 6 to 20 threads; all results
  bitwise identical. Stage mode reaches 4.5x (N = 512) and 5.7x (N = 2048) vs serial IAS15, up from 3.2x
  and 3.6x at the best T before.
- Cost: large, well-balanced work loses the efficiency cores. ias15_n2048 default 3.85x (P-cores) vs
  4.42x (all cores), stage mode 5.69x vs 8.55x. Dynamic scheduling of 7×16 blocks absorbs slower cores;
  barrier-synchronised loops with few tasks do not.
- **Problem found: mixed team sizes.** whfast_n256 drops from 1.93x (T = 4) to 1.27x at T ≥ 8 (P-cores),
  and to 0.04x at T = 16–20 with all cores (25x slower than serial). Since iteration 7, gravity runs with
  N/64 = 4 threads while the other loops of the step (Kepler drift, jump) use all threads. Likely cause,
  not yet verified: libgomp keeps idle threads spinning between regions; with teams of different size the
  non-participating threads spin or sleep and wake in alternation, and with 16 threads on 14 places
  (`OMP_PLACES=cores`) spinning threads share cores with working ones. Run 1 (before iteration 7, one team
  size) had whfast_n256 1.66x at T = 16. To fix (P14): keep one team size per process and give the gravity
  work to the first N/64 threads inside it, or set `OMP_WAIT_POLICY=passive`-like behaviour for mixed
  regions; test also `OMP_WAIT_POLICY` on the workstation.

## Current state on the VM (commit `93d4a96` + threshold reorder)

`run_parallel.py --repeat 5` on the 2-vCPU VM (noisy; see above):

| case | serial [s] | omp 1 thread | omp 2 threads | threads | vs serial |
|---|---|---|---|---|---|
| ias15_hltau | 0.413 | 0.88–1.05x | 0.91–1.06x | identical | 0 (bitwise) |
| ias15_embryos | 0.509 | 0.90–1.02x | 0.98–1.00x | identical | 0 |
| ias15_n128 | 0.386 | 0.84–1.01x | 0.93–0.94x | identical | 0 |
| ias15_n512 | 0.641 | 1.03–1.05x | 1.00–1.01x | identical | 0 |
| ias15_n2048 | 1.140 | 0.95–0.97x | 0.77–0.92x (kernel alone: 1.18–1.22x) | identical | 0 |
| whfast_ss | 0.228 | 1.12x | 1.14x | identical | 0 |
| whfast_tp10k | 0.417 | 0.96–1.07x | **1.15–1.61x** | identical | 0 |
| whfast_n256 | 0.653 | 0.92–0.99x | 0.68–0.92x | identical | 0 |

(ranges: the last two full runs; the IAS15 N = 2048 run has only 3 adaptive steps, the kernel-only
measurement is `gravity_check`-style timing of one force evaluation, min of 9 interleaved runs.)

## Plan (work items)

Priority = expected speedup at T > 1 vs the serial build × breadth ÷ effort/risk.

| # | Item | Cases | Expected | Effort | Status |
|---|---|---|---|---|---|
| P1 | Fix the OpenMP Jacobi μ bug (correctness) | whfast_* | results = serial | trivial | **done (it. 1)** |
| P2 | Never slower than serial: serial code paths below a work threshold, AVX2 clones, no fork/join for small loops | all | omp ≥ 1.0x | medium | **done (it. 2)** |
| P3 | Basic gravity, parallel and bitwise identical to serial | ias15_n*, whfast_n256 | ~1.6–1.9x on 2 cores (estimate) | medium | **done (it. 4): 1.2x at N ≥ 1024 on this VM; single-core cost now ≤ serial for N ≤ 1024** |
| P4 | Half-pair loop with fixed block decomposition and per-block accumulators (deterministic, round-off-level differences to serial) | same | ~1.8x on 2 cores | medium | open: only if bitwise identity may be given up |
| P5 | WHFast test particles: parallel batched Kepler drift, parallel test-particle forces | whfast_tp10k | ~1.8x on 2 cores | low–medium | **done (it. 5): 1.4–1.6x** |
| P6 | WHFast massive bodies: Jacobi μ once, parallel batched drift | whfast_n256 | small (gravity dominates) | low | **done (with it. 5)** |
| P7 | One parallel region per IAS15 step, work-shared loops inside | ias15_embryos, n128 | reduces overhead at moderate N | medium–high | open: the O(N) IAS15 loops are < 5% of a step for N ≥ 128 |
| P8 | Direct collision search (serial O(N²)) | ias15_embryos | ? | low | open |
| P9 | TRACE/MERCURIUS OpenMP correctness | trace_*, mercurius_* | tests pass | low | **done (it. 3)** |
| P10 | Use the structure-of-arrays tile kernel in the serial build too (it is bitwise identical and faster on one core for 256 ≤ N ≤ 1024) | ias15_embryos(?), ias15_n*, whfast_n256 | 1.1–1.4x single-core | low | new idea (single-core, outside this round's goal) |
| L1 | IAS15 stage-parallel (Jacobi over the 7 substeps, PIRK) | ias15_* N ≳ 100 | up to 7x on the force part | medium | **done (it. 6): 1.6–2.1x on 2 cores, opt-in `parallel_stages`** |
| P11 | Hybrid CPUs: balance stage tasks; P-cores only | stage mode, gravity at T > 6 | no loss beyond 6 threads | low | **done (it. 8, 9); measured in it. 10** |
| P12 | Lower the gravity threshold after measurements on real CPUs | ias15_n512 at T = 2, whfast_n256 | 1.5x+ where now 1.0x | trivial | **done (it. 7)** |
| P13 | Stage mode × particle-parallel gravity (nested) for > 7 cores | ias15_n* | > 7x | medium | open |
| P14 | One team size per step (gravity on a subset of threads inside the full team), or a persistent parallel region per step (P7) | whfast_n256 T ≥ 8 (1.27x → ≥ 1.9x), all mixed-size cases | removes the 0.04x pathology | low–medium | open (high priority) |
| P15 | Use efficiency cores where tasks are many and dynamically scheduled (stage mode with K ≥ 8; test-particle loops) | ias15_n2048+stages | 5.7x → 8.5x | low | open |
| P16 | Stage mode for more configurations: test particles, compensated gravity, ghost boxes | IAS15 users with test particles | same speedups | medium | open |
| P17 | Speculative overlap of consecutive IAS15 steps (start step n+1 from the predicted b while step n converges) | IAS15 N ≲ 100 | 1.2–1.5x on top of stage mode (estimate) | high (research) | open |
| P18 | Ensembles: run independent simulations on separate cores (Python helper), WHFast512 several systems per vector | small systems (Solar System, WHFast512) | ~linear in cores | low | open (the way to use cores for N ≲ 10) |

Cases expected to stay at ~1.0x (too little work per step for thread-level parallelism: a fork/join costs
about as much as a whole step): ias15_hltau (N = 6, ~7 µs/step), whfast_ss (8 planets, ~1 µs/step),
WHFast512 (~150 ns/step). For these, the goal is not getting slower; ensembles remain the way to use cores.
