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

Later rounds (iterations 12–26) added:
- one OpenMP team size per step;
- automatic use of the efficiency cores for large systems;
- stage mode for test particles, compensated gravity and ghost boxes;
- deterministic parallel tree collisions (they depended on thread scheduling);
- bitwise MERCURIUS parallel loops;
- a parallel tree group walk;
- paired WHFast512 ensembles;
- the opt-in block-order basic gravity (`gravity_block_order = 1`).

User guide: `docs/openmp.md`.

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

## Merge with `main` (28 commits of the parallel single-core work, up to `870ca53`)

`main` (= `perf-bench`) had moved on from this branch's base `62d6d93a`: TRACE/MERCURIUS broadphase and
blocked loops, tree gravity group walk and cell reuse, WHFast512 NEON and cheaper restarts, an ensemble API
(`reb_simulations_integrate`/`steps`, Python `integrate_simulations`, a pthread pool), and `rebound.h`
mapping `_OPENMP` to `OPENMP` (so `-fopenmp` via CFLAGS selects the parallel code paths, which now are the
thread-threshold paths of this branch).
- One textual conflict, MERCURIUS encounter mode: `main` added an exact far-pair skip to the serial loop
  (`#ifndef OPENMP`), this branch had made serial/parallel a runtime choice. Resolved: the new serial
  loop is the serial branch of the runtime choice. TRACE merged automatically; checked that `main`'s new
  serial loops read the encounter flags at `[lower*N + higher]`, the convention the OpenMP branch uses.
- Interaction fixed in the merge: ensemble workers are threads of their own; in an OpenMP build each would
  have started full OpenMP teams inside REBOUND's loops (workers × threads on the same cores). Workers now
  run with one OpenMP thread (`omp_set_num_threads(1)` in the worker thread; the calling thread's setting
  is restored). Results do not depend on the thread count, so ensemble results are unchanged.
- Checks: test suites 966 passed (serial build) and 966 passed (OpenMP build, 2 threads) — `main`'s tests
  plus this branch's; `gravity_check` 189/189 bitwise: merged serial = `main` serial = merged OpenMP at
  T = 1, 2, 3, 5; pbench signatures (stage mode, whfast_n256, whfast_tp10k) identical in the serial build and
  at T = 2, 3; an ensemble of 4 stage-mode IAS15 simulations (N = 300) equals running them one by one.
  No performance runs (to be done on the workstation).

## Iteration 12 — one team size for all regions (item P14, commit `d50ee29`) — kept

Problem (iteration 10): since iteration 7, gravity ran with N/64 threads (4 at N = 256) and the other
loops of a step (Kepler drifts) with all threads. whfast_n256 fell from 1.93x (T = 4) to 1.27x
(T >= 8, 6 P-cores), and to 0.04x with all cores at T = 16-20. Likely mechanism: idle pool threads
are put to sleep and woken again whenever the team size changes, much more expensive when threads
outnumber cores (libgomp then throttles spinning).
Change: every parallel region runs with `reb_omp_threads()` threads. The gravity kernels let
threads 0..nw-1 (nw = the N/64 cap) work, through static ranges (`reb_omp_static_range`) and an
atomic block counter (`reb_omp_next`, `omp atomic capture`); the others wait at the barriers.
Results unchanged (gravity_check 189/189, pbench signatures). VM reproduction with 8 threads on
2 vCPUs (oversubscribed like T = 16-20 with all cores): whfast_n256 0.21 s (old) → 0.175 s (new);
at T = 2 equal within noise (min 0.187 vs 0.198 s). OpenMP suite 966 passed (T = 3).

## Iteration 13 — efficiency cores for large systems (item P15, commit `e8f5f50`) — kept

The core policy is chosen once per step and per force evaluation (`reb_omp_select_cores(r)` in
`reb_simulation_step` and `reb_simulation_update_acceleration`, thread-local flag), so the team size
stays constant within a step (iteration 12): all cores when N_active >= 1024 (`REB_OMP_ALL_CORES_N`),
otherwise only the performance cores. Basis (iteration 10, i5-13500): all cores won at N = 2048
(default 4.42x vs 3.85x, stage mode 8.04x vs 5.32x) and lost at N = 512 (default 2.75x vs 3.15x,
stage 4.53x vs 4.70x) and for whfast_tp10k (N_active = 3). Override: `REBOUND_OMP_CORES=all` or
`performance` (the old `REBOUND_OMP_ALL_CORES=1` still works); `run_parallel.py --all-cores`,
`--performance-cores`. Policy unit test: N = 500 → P-cores, N = 2000 → all, N_active = 3 with 10,000
test particles → P-cores, N_active = 1500 → all; overrides respected. Hybrid detection itself cannot
run on the VM (tested on fake sysfs trees in iteration 9). Tests 966 passed (serial and OpenMP).

## Iteration 14 — stage mode for test particles, compensated gravity, ghost boxes (item P16, commit `6c3dc56`) — kept

- Basic gravity: row blocks of the active particles (as before) plus row blocks of the test-particle
  rows (partners: the active particles; back-reaction only for `testparticle_type = 1`), about 8,192
  pair interactions per block; each task loops over the ghost boxes in the serial order; partial
  accelerations summed in a fixed task order. `gravity_ignore_terms` via the serial start indices.
- Compensated gravity: particle chunks, each particle summing its partners with the serial Kahan
  order (one-sided form of the serial loop), so the stage forces are bitwise equal to the serial
  compensated gravity at the same positions; no partial sums to combine.
- The plan depends on N, N_active and the gravity settings only. Previously supported configurations
  give unchanged results (pbench stage signatures).
- New `benchmarks/parallel/stage_check.c`: 14 configurations (all-active, test particles type 0/1
  incl. 3 active + 2,000 test particles, ignored terms, ghost boxes, compensated with and without test
  particles). Stage mode: same step counts as the default mode, final positions within 0–9e-15
  (relative), and bitwise identical between the serial build and the OpenMP build at T = 1, 2, 3.
  4 new Python tests. New pbench case `ias15_tp2000` (20 planets + 2,000 test particles).
- VM, 2 threads, ias15_tp2000 vs serial default IAS15: default 1.45x, stage mode **1.7x**.
- Still unsupported (flag ignored): additional forces (user callbacks read `r->particles`),
  variational particles, particle maps, MEGNO.

## Iteration 15 — stage mode beyond 7 cores: up to 64 blocks per stage (item P13, commit `61a3894`) — kept

Stage parallelism × row-block parallelism already existed since iteration 8 (7 × K tasks, K ≤ 16);
the i5-13500 reached 8.04x at T = 20 (N = 2048). Profile of stage mode (callgrind, one thread):
stage force kernel 86.6 % (N = 512) / 91.4 % (N = 2048), force at the step start 4.8 / 5.1 %
(parallel), prediction and reduction 3.5 / 1.3 % (parallel), IAS15's serial loops ~1 % / < 1 %: no
Amdahl limit below ~100 threads. The limit is granularity: 112 tasks are thin for many cores or mixed
core speeds. Change: K ≈ N/64 up to 64 (448 tasks per iteration), scratch ≤ 64 MB. Results
independent of T (serial and T = 1, 2, 3, 5 identical), round-off change for N ≥ 1024 in stage mode.
VM A/B (N = 2048, K 16 → 32): neutral at T = 1, 2 (min 0.542 vs 0.552 s, 0.288 vs 0.295 s). The gain
can only be measured on many cores. Tests 970 passed (serial and OpenMP).

## Iteration 16 — workstation measurement of iterations 12–15 (user, i5-13500, commit `a110755`) — measurement only

Raw output: `i5-13500_run4_commita110755.txt` (work directory, `bench_results/parallel/`). Best speedup
vs the serial build (stage rows: vs the default IAS15 of the serial build), per core policy:

| case | auto (default) | `--all-cores` | `--performance-cores` | before (it. 10, P-cores) |
|---|---|---|---|---|
| ias15_hltau / embryos / n128 | 0.92–1.01x | 0.96–1.00x | 0.95–1.01x | 0.96–1.01x |
| ias15_n512 | 3.18x | 2.87x | 3.17x | 3.15x |
| ias15_n2048 | **5.10x** (T = 20) | 5.00x | 3.94x | 3.85x |
| ias15_tp2000 (new) | 1.85x | 1.84x | 1.84x | – |
| embryos+stages | 1.91x | 1.80x | 1.89x | 1.89x |
| n128+stages | 2.89x | 3.11x (T = 8) | 2.88x | 2.89x |
| n512+stages | 4.60x | **5.79x** (T = 16) | 4.59x | 4.51x |
| n2048+stages | **8.87x** (T = 20) | 8.61x | 5.68x | 5.69x |
| tp2000+stages (new) | 2.45x | 2.66x (T = 16) | 2.45x | – |
| whfast_tp10k | 2.32x | 2.34x (T = 16), 1.27x (T = 20) | 2.29x | 2.28x |
| whfast_n256 | **1.99x, 1.94x at T ≥ 8** | 1.98x (T = 4), 1.37x (T = 20) | 2.02x, ~1.95x at T ≥ 8 | 1.93x (T = 4), 1.27x at T ≥ 8 |

All results bitwise identical at every thread count and to the serial build.

Readings:
- **Iteration 12 confirmed:** whfast_n256 holds 1.94–1.99x from 4 to 20 threads (was 1.27x at T >= 8),
  and the all-cores pathology is gone (1.37x at T = 20 instead of 0.04x).
- **Iteration 13 confirmed for the default mode:** auto picks the better policy in every default-mode
  case. ias15_n2048 5.10x (all cores, also better than 4.42x in it. 10 thanks to the single team size);
  N <= 512 and test-particle cases stay on the P-cores.
- **Stage mode wants the efficiency cores earlier:** n512+stages 5.79x with all cores at T = 16 vs
  4.60x auto. tp2000+stages 2.66x vs 2.45x. n128+stages 3.11x at T = 8, but worse at T = 20.
  Dynamically scheduled stage tasks absorb slower cores already at N ≈ 512. Follow-up P17b: in stage
  mode, use all cores from N_active ≈ 512 on.
- **Iteration 15 (64 blocks):** n2048+stages at T = 20 with all cores 8.61–8.87x vs 8.55x in
  iteration 10. Within noise on 20 threads, as expected. The finer split is meant for larger
  machines and stays (no cost at low T).
- **tp2000+stages at one thread costs 1.27x the default mode** (0.79x): the test-particle rows run
  the row kernel with only 20 partners per row (short inner loops), while the default mode's
  test-particle loop is cheaper. Still 2.45x vs serial at 8+ threads. Follow-up P18b: a test-particle
  kernel vectorised over test particles.
- ias15_hltau 0.92–0.97x in all three runs (N = 6, ~7 µs per step): a few percent of fixed per-step
  overhead in the OpenMP build. Earlier runs gave 0.95–1.00x, so this is close to noise. Item P19
  (check): the per-step `reb_omp_select_cores` call and `omp_get_max_threads` queries.

## Iteration 17 — stage mode uses all cores from N_active = 512 (item P17b, commit `efee47d`) — kept, to be measured

`reb_omp_select_cores`: when the simulation's integrator is IAS15 with `parallel_stages = 1`, the
all-cores threshold is `REB_OMP_ALL_CORES_N_STAGES` = 512 instead of 1024. Basis (iteration 16):
n512+stages 5.79x with all cores vs 4.60x on the P-cores; tp2000+stages 2.66x vs 2.45x; default mode
at N = 512 better on the P-cores (3.18x vs 2.87x), so it keeps 1024. Unit test: IAS15 N = 600 default
→ P-cores, stage mode → all, N = 300 stage → P-cores, WHFast N = 600 → P-cores. OpenMP suite 970 passed.

## Iteration 18 — test-particle rows ordered by partner; in-place writes (item P18b, commit `308a90a`) — kept

Profile (callgrind, one thread, ias15_tp2000): default mode 217 M instructions in gravity; stage mode
181 M in the row kernel plus 52 M in prediction, buffer zeroing and reduction (~1.7 MB of partial
buffers per iteration: 5 tasks × 3N per stage). The i5 showed 0.79x on one thread, the VM did not
(stage mode was already faster there on one thread), so memory traffic is the likely culprit on the i5.
Changes, both bitwise identical to the previous stage results:
- `reb_gravity_basic_tp_rows_arrays`: test-particle rows computed partner by partner over chunks of
  64 rows (each row has only N_active ≈ 20 partners); per particle the same terms in the same order.
- Tasks that are the only writers of their particles write the stage accelerations in place:
  test-particle blocks of type 0 always, the active block when there is only one and no back-reaction.
  With few massive bodies plus type-0 test particles, no partial buffers and no reduction remain.
Checks: stage_check 14/14 identical to iteration 14 (serial and T = 1, 2, 3), pbench signature
identical. VM, interleaved ×5: one thread median 0.297 s vs 0.314 s (min 0.293 vs 0.295); two threads
equal within noise. Tests 970 passed (serial and OpenMP). The effect on the i5's 0.79x is to be
measured.

## Iteration 19 — core policy once per step (item P19, commit `e3ab73b`) — kept

Profile of ias15_hltau (N = 6), OpenMP build at T = 1 vs serial build: +4.4 % instructions, of which
`reb_omp_select_cores` 2.1 % (called from `reb_simulation_update_acceleration`, i.e. ~20 times per
IAS15 step) and the gravity routine's thread checks ~1 %. The policy is now chosen only at the start
of each step (it only affects speed). Overhead 4.4 % → 1.5 % of instructions; results identical.
VM wall clock cannot resolve the difference (minimum 0.1120 s serial, 0.1127 s before, 0.1125 s
after). Tests 970 passed (serial and OpenMP).

## Iteration 20 — workstation measurement of iterations 17–19 (user, i5-13500, commit `0482a09`) — measurement only

Raw output: `i5-13500_run5_commit0482a09.txt`. Auto core policy, compared with iteration 16 (auto):

| case | it. 16 | it. 20 | reading |
|---|---|---|---|
| ias15_hltau (N = 6) | 0.92–0.97x | **0.98–1.00x** | it. 19 confirmed |
| ias15_n512+stages vs default serial, T = 20 | 4.60x (P-cores) | **6.02x** (all cores from N = 512) | it. 17 confirmed (+31 %) |
| ias15_tp2000+stages, one thread vs default serial | 0.79x | 0.81x (0.492 → 0.472 s) | it. 18: only ~4 % on the i5; the gap is not mainly the partial buffers |
| ias15_tp2000+stages, best | 2.45x | 2.51x | |
| ias15_n2048+stages vs default serial, T = 20 | 8.87x | 8.62x | within run-to-run noise |
| ias15_n2048 default | 5.10x | 5.21x | |
| whfast_n256 | 1.94–1.99x | 1.93–2.01x | stable from 4 to 20 threads |

All results bitwise identical at every thread count and to the serial build.

Open question from this run: on one thread, stage mode costs 1.06–1.23x of the default mode on the i5
(embryos 0.86x, n128 0.94x, tp2000 0.81x), but not on the VM. Instruction counts are within ~10 %, so
the likely difference is the number of predictor-corrector iterations per step: a Jacobi sweep
converges more slowly than IAS15's Gauss-Seidel sweep. Next diagnostic (cheap): count iterations per
step in both modes for these cases. If confirmed, options are a Gauss-Seidel first iteration or
stopping on the stage-mode convergence history.

## Iteration 21 — why stage mode is slower on one thread; sequential sweeps (diagnostic) — reverted

Instrumented copy of the library (iteration counter in `reb_integrator_ias15_step_try`; counts are
deterministic, so they hold for the i5 too). Predictor-corrector iterations per step attempt:

| case | default IAS15 (Gauss-Seidel over substeps) | stage mode (Jacobi) | 2 groups | 3 groups | 7 groups (sequential) |
|---|---|---|---|---|---|
| ias15_embryos, short / long run | 3.22 / 3.13 | 4.12 / 4.02 (+28 %) | 4.08 | 4.03 | 3.22 / 3.13 |
| ias15_n128 | 3.28 / 3.42 | 3.95 / 4.16 (+21 %) | 3.95 | 3.88 | 3.22 / 3.39 |
| ias15_tp2000 | 5.11 / 4.51 | 6.05 / 5.29 (+17 %) | 5.84 | 5.68 | 5.11 / 4.51 |
| ias15_n512 | 4.50 / 5.25 | 5.17 / 6.18 (+18 %) | 5.33 | 5.17 | 4.67 / 5.21 |
| ias15_n2048 | 4.67 / 3.43 | 5.33 / 4.11 (+20 %) | 5.00 | 5.00 | 4.67 |

This is the cause of the one-thread gap seen on the i5 (stage mode 0.81–0.94x). The force kernels
execute about the same number of instructions; the VM did not show the gap because its stage kernel
happens to run faster than the default kernel there. Jacobi iterations over the substeps converge
more slowly than IAS15's substep-by-substep sweep, as expected from the PIRK literature.

Tried (in the source, results independent of the thread count):
- substeps in 2 or 3 consecutive groups (each group from the b values updated by the earlier
  groups): saves almost nothing (0–7 %) and divides the parallel width;
- the first iteration of each step sequential (Gauss-Seidel), then Jacobi; used only with ≥ 32
  row blocks per substep, so a sequential iteration still has 32 parallel tasks. Iterations at
  N = 2048: 4.11 → 3.69 per step (long run). But wall time did not improve (VM, 107 steps, 3 runs
  each: one thread 23.2–23.9 s vs 23.3–24.5 s; two threads 12.8–13.0 s vs 12.7–13.6 s): a
  sequential iteration costs about 1.5x a Jacobi iteration there.
Reverted. Conclusion: the extra 17–28 % iterations are the price of evaluating the substeps in
parallel; stage mode pays off from about 2 threads (N ≥ 128) or 4 threads (N ≈ 60), and on 1–2
threads the default mode should be used (documented in the setting's description).

## Iteration 22 — new benchmark cases; tree collisions were not deterministic (commit `af5c99d`) — kept (correctness)

New pbench cases (setups from `benchmarks/perf/bench.c`): `tree_selfgrav` (leapfrog, tree gravity,
N = 10,000), `tree_selfgrav_group` (same with `tree_group_size = 32`), `tree_rings` (SEI shearing
sheet, ~5,000 particles, tree gravity + tree collisions), `mercurius_plts` (2 planets + 2,500
planetesimals, type 1), `mercurius_disc` / `trace_disc` (2 planets + 500 active planetesimals,
merging), `trace_tp` (2 planets + 2,000 test particles), and ensembles (`ensemble_ss512`,
`ensemble_3p512`, `ensemble_ss`; 32 clones through `reb_simulations_steps`, as many workers as
threads, also in the serial build).

First run on the VM found two correctness problems:
- **Tree collision search (TREE and LINETREE):** with 2+ threads the parallel search appends
  collisions (critical section) in thread order. Since the resolution order is a shuffle of that
  list, results changed from run to run (reproduced: two runs at T = 2 differ). Fix: after a
  parallel search, a stable counting sort by projectile restores exactly the serial order (each
  projectile is searched by one thread, in the serial order). O(N + N_collisions). Now identical to
  the serial build at T = 1, 2, 3 and on reruns, with and without ghost boxes and tree gravity.
  Regression test `test_openmp_determinism.py` (negative control: without the sort it fails with
  "1500 of 1500 particles differ between 1 and 2 threads").
  Side note: `sim.rand_seed` is random by default (time and PID), so runs with collisions are not
  reproducible unless it is set; that is upstream behaviour, not an OpenMP issue (documented).
- **MERCURIUS (iteration 23).**

## Iteration 23 — MERCURIUS and TRACE: bitwise parallel loops, used from 3 threads (commit `af5c99d`) — kept

Since the gravity threshold of iteration 7 (N ≈ 256), MERCURIUS's own OpenMP branches (WH-mode and
encounter-mode pair forces) engaged for the 500-body disc, and they used other expressions
(`-G*m*L/r^3`) and no far-pair test: results differed from the serial build at T = 2. Both are now
one-sided loops that reproduce the serial half-pair loop bitwise (per particle: star term,
partners in ascending order, then type-1 back-reaction; a partner-row term
`(prefact*m_i)*(x_i - x_k)` equals `(-prefact*m_i)*(x_k - x_i)` exactly; the same far-pair tests).
One-sided loops do twice the pair work; on 2 threads MERCURIUS and TRACE were slower than serial
(0.77x / 0.74x). New `reb_omp_parallel_onesided`: these loops run in parallel only from 3 threads
(the choice cannot change results). Checks: mercurius_disc, trace_disc, mercurius_plts, trace_tp
identical to the serial build at T = 1-4, also in a build with all thresholds forced to 4 (every
loop parallel).

## Iteration 24 — tree gravity: parallel group walk (commit `af5c99d`) — kept

Profile (callgrind, tree_selfgrav): tree walk 80 % (parallel), tree construction ~10 % (serial;
limits 2 threads to ~1.7x). The opt-in group walk (`tree_group_size`, added on main) was entirely
serial and used a shared in-group flag array. Now: group-id lookup (read-only), groups walked in
parallel with per-thread scratch, one group per call of an AVX2-cloned function (the first version
inside the OpenMP region ran at 0.72x on one thread: outlined region bodies get neither the clone
nor the same vectorisation). Bitwise identical to the serial group walk. VM: tree_selfgrav_group
T = 1 0.96-1.04x, T = 2 1.42-1.61x; tree_selfgrav T = 2 1.25-1.5x; tree_rings T = 2 1.2-1.6x.

## Iteration 25 — ensembles: WHFast512 simulations paired (commit `51083d4`) — kept

`reb_simulations_steps` / `rebound.steps_simulations` claim simulations in pairs. If both are
WHFast512 simulations that can share a paired kernel (C AVX-512; C AVX2 four-lane for <= 4
planets; C generic) and have no per-step extras (heartbeat is fine; no timestep modifications,
ODEs, variational particles, collisions, archive), each REBOUND step runs one paired kernel call
(`reb_integrator_whfast512_step_pair`, no extra synchronisation) and the rest of the step per
simulation exactly as `reb_simulation_steps`. `steps_pair` itself could not be used: it
synchronises after every call, and N REBOUND steps of C kernel steps are not bitwise equal to one
call of N*C steps (both checked). Test `test_ensemble_pairing.py`: bitwise equal to stepping each
simulation alone, kernels 1-4, encounter in one simulation of a pair. VM, 32 three-planet systems,
2 workers: 0.46 -> 0.36 s (AVX2 four-lane, 1.27x), 0.75 -> 0.57 s (C AVX-512, 1.31x).
Not paired: 8 planets on AVX2-only CPUs (no paired kernel, no gain), the assembly kernel; the
integrate mode (`integrate_simulations`, exact finish times).

## Iteration 26 — opt-in basic gravity in block order, `gravity_block_order = 1` (commit `384ee10`) — kept

The default parallel basic gravity reproduces the serial loop bitwise, at ~1.5x its work plus a
pair table. New saved setting `sim.gravity_block_order = 1`: row blocks (N_active/64, at most 64;
test-particle blocks of ~8,192 pairs; ghost boxes in the serial order), partial sums added in a
fixed order, in place where one block owns its particles; the same code in serial and OpenMP builds
(plain `omp for` with the full team). Results identical between builds and for any thread count
(`gravity_check blocks`: 126/126 at T = 1, 2, 3, 5), max difference to the default 7.7e-15
(relative to the largest acceleration). VM, speedup vs the default serial build:

| case | default T = 2 | block order T = 1 | block order T = 2 |
|---|---|---|---|
| whfast_n256 | 0.90x | 1.18x | 1.90x |
| ias15_n512 | 1.25x | 1.21x | 1.76x |
| ias15_n2048 | 1.04x | 1.24x | 2.33x |
| ias15_tp2000 | 1.36x | 1.31x | 1.53x |

Faster also on one core (structure-of-arrays kernel). New user documentation `docs/openmp.md`
(building, determinism, what runs in parallel, hybrid CPUs, stage mode, block order, ensembles);
`docs/gravity.md` corrected (it said OpenMP makes basic gravity O(N^2)).
Tests (whole round): 979 passed (serial build), 983 passed (OpenMP build, 3 threads).

## Iteration 27 — single-core regression in the group walk found and fixed (commit `639c71e`) — kept

Single-core A/B of the whole `benchmarks/perf` suite against `bbaacec` (the commit before this
round, `run.py ab --repeat 3`): geomean 0.983x, all 25 benchmarks bitwise identical, but
`selfgrav_group` 0.72x. Cause: the per-group function of iteration 24 received `softening^2` as a
parameter, so GCC kept the errno branch of `sqrt` and did not vectorise the cell and particle
loops (callgrind: 389 M vs 250 M instructions; the same pitfall as in iteration 4). The AVX2 clone
also cost more than it saved for groups of at most 32 particles. Fix: square the softening inside
the function, no clone. Targeted A/B (`--repeat 5`): selfgrav_group 0.983x, rings_poly_tree 1.041x,
selfgrav_tree 0.987x, ias15_embryos 1.006x, mercurius_disc 0.999x, trace_disc 0.977x; geomean
0.999x, all identical. The other benchmarks of the full run were within ±3 % (noise).

## Iteration 28 — ensembles did not scale with thread binding (commit `ca1564b`) — kept

The end-to-end `run_parallel.py` run showed ensembles at 1.0x on 2 threads (the direct run had
2.04x). With `OMP_PROC_BIND` set, the OpenMP runtime binds the calling thread to one core, and the
ensemble pool's pthreads inherited that binding. In OpenMP builds the pool is now an OpenMP
parallel region: the workers follow `OMP_PROC_BIND` / `OMP_PLACES`, and REBOUND's regions inside a
worker are nested (single-threaded). VM with binding: ensemble_ss512 1.97x, ensemble_3p512 1.88x,
ensemble_ss 1.65x at 2 threads, results identical.

VM run of all cases after iterations 22–28 (`run_parallel.py --repeat 1 --threads 1,2`, before the
ensemble fix; single repeat, so ±10–20 %): tree_selfgrav 1.86x, tree_selfgrav_group 1.48x,
tree_rings 1.48x, trace_tp 1.30x, mercurius_plts / mercurius_disc / trace_disc 0.96–1.03x (one-sided
loops engage only from 3 threads), block order 1.43–1.80x (vs default serial, at 2 threads:
ias15_n512 1.89x, ias15_n2048 2.21x, whfast_n256 1.82x), stage mode as before. All results
identical between builds and thread counts.
Tests at the end of the round: 979 passed (serial build), 983 passed (OpenMP build, 2 and 3
threads, also with binding), and 983 passed with every parallel path forced on (thresholds set to
4, 2 threads; 30 minutes because every tiny loop then runs in parallel).

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
| P13 | Stage mode × particle-parallel gravity (nested) for > 7 cores | ias15_n* | > 7x | medium | **done (it. 8, 15): 8.0x at T = 20 measured (it. 8 version); up to 448 tasks now** |
| P14 | One team size per step (gravity on a subset of threads inside the full team), or a persistent parallel region per step (P7) | whfast_n256 T ≥ 8 (1.27x → ≥ 1.9x), all mixed-size cases | removes the 0.04x pathology | low–medium | **done (it. 12)**; to be measured on the workstation |
| P15 | Use efficiency cores where tasks are many and dynamically scheduled (stage mode with K ≥ 8; test-particle loops) | ias15_n2048+stages | 5.7x → 8.5x | low | **done (it. 13)**; to be measured |
| P16 | Stage mode for more configurations: test particles, compensated gravity, ghost boxes | IAS15 users with test particles | same speedups | medium | **done (it. 14)** |
| P17 | Speculative overlap of consecutive IAS15 steps (start step n+1 from the predicted b while step n converges) | IAS15 N ≲ 100 | 1.2–1.5x on top of stage mode (estimate) | high (research) | open |
| P17b | Stage mode: use all cores from N_active ≈ 512 (it. 16: n512+stages 5.79x vs 4.60x) | n512+stages, tp2000+stages | 1.1–1.25x at T ≥ 16 | trivial | **done (it. 17)**; to be measured |
| P18b | Stage mode: test-particle rows vectorised over test particles | tp2000+stages (0.79x on one core) | 1.2x | low | **done (it. 18)**; to be measured on the i5 |
| P19 | Per-step OpenMP bookkeeping for tiny systems (hltau 0.92–0.97x) | ias15_hltau | back to 1.00x | trivial | **done (it. 19)**: 4.4 % → 1.5 % instructions |
| P18 | Ensembles: run independent simulations on separate cores (Python helper), WHFast512 several systems per vector | small systems (Solar System, WHFast512) | ~linear in cores | low | open (the way to use cores for N ≲ 10) |

Cases expected to stay at ~1.0x (too little work per step for thread-level parallelism: a fork/join costs
about as much as a whole step): ias15_hltau (N = 6, ~7 µs/step), whfast_ss (8 planets, ~1 µs/step),
WHFast512 (~150 ns/step). For these, the goal is not getting slower; ensembles remain the way to use cores.
