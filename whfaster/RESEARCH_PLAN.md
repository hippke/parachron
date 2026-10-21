# WHFaster: research plan for a faster-than-WHFast512 Solar System integrator

*Written Oct 20, 2026 for a future research agent. The work directory is `/home/michael/Downloads/rebound/whfaster/` on the user's workstation. The parent directory `/home/michael/Downloads/rebound/` holds the logs of two earlier rounds of work on the same fork.*

## 1. Goal

Integrate the Solar System (8 planets, general-relativistic correction) over Gyr timescales **at research-grade accuracy, faster than the current state of the art** described in:

> Dagli & Rein (2026), *A Very Fast Symplectic Integrator for Planetary Systems: WHFast512 in x86 Assembly*, arXiv:2609.25339 (`2609.25339v1.pdf` in this directory).

Start from **our fork's C implementation** (https://github.com/hippke/rebound, branch `main`, see §3).

- **Single-core and multi-core work are both allowed.** The user's example of an acceptable result: **2x faster wall-clock time using 5x as many CPU cores.**
- Always report both:
  - wall-clock time for one trajectory;
  - total core-hours, i.e. efficiency.

### 1.1 The SOTA benchmark to beat (from the paper)

| Item | Value |
|---|---|
| System | Sun + 8 planets, GR as an extra 1/r² potential on the central body (Nobili & Roxburgh 1986) |
| Splitting | Jacobi coordinates, symplectic correctors of order 17 |
| Timestep | dt = 6 days (the convergence test also uses 5 days) |
| Length | 5 Gyr ≈ 3.04 × 10¹¹ steps |
| Time | **9.0 h** on an AMD EPYC 9655 (Zen 5, 2.6 GHz), i.e. **≈ 106 ns per step** |
| Other results | 60 % faster than WHFast512 v1; > 8x faster than WHFast (-O3 -march=native); similar speedups on an Intel Xeon Gold 6148 |
| Accuracy | Median relative energy error < 1e-10 over 10 Gyr (80 runs); grows like 1e-16·√N_steps (Brouwer's law, unbiased); 5 orders of magnitude better than v1 at equal dt |
| Cost split | Kepler step 63 % (half of it in the iterations), interaction step 25 %, Jacobi⇄heliocentric transformations 11 % |
| Authors' claim | "almost optimal"; "speedups of ≳ 2× will be hard to achieve" on one core |

**Convert any measurement** to the 5 Gyr figure: `wall_5Gyr[h] = ns_per_step × 3.04e11 / 3.6e12`. Always measure the reference assembly kernel on the **same machine** in the same session, and report the ratio. Absolute numbers depend strongly on the CPU: our cloud Xeon runs the same assembly at ~160 ns/step, the paper's Zen 5 at ~106 ns/step.

### 1.2 What "research-grade" means here (the accuracy contract)

A faster method counts only if it keeps all of the following. Agree on this list before optimizing, and never relax it silently.

1. **Unbiased long-term energy error.** Over ≥ 10⁹ steps the ensemble-mean relative energy error is consistent with zero, and the median grows like √N_steps (paper Fig. 5 and Appendix C). A linear (biased) drift is a failure, even if it's small.
2. **Energy level.** Median |ΔE/E| comparable to the reference at the same physical setup: ≲ 1e-10 at 10 Gyr, and ~1e-13 over 10⁷ steps in our 16-run development ensembles.
3. **Correct secular dynamics.** The fundamental frequencies g₁…g₈ and s₁…s₈ (frequency analysis over 20–50 Myr) agree with a high-accuracy reference: WHFast512 at dt = 2–3 days with correctors, or IAS15 for short spans. Mercury's g₁ and its g₁–g₅ resonance drive the instability. Energy error alone is **not** a sufficient check (arXiv:2507.04987).
4. **Instability statistics.** Where affordable, Mercury's instability rate over 5 Gyr (~1 %, Laskar & Gastineau 2009; Abbot et al. 2023; Brown & Rein 2023) must not change beyond its statistical error. This is expensive and is the final validation, not a per-iteration test.
5. **Kepler solver quality.** The paper's grid tests (Fig. 2, 3, 7): relative a-error at machine precision with random sign over the relevant (a, e) region, convergence up to e ≈ 0.99, and a correct fallback for hyperbolic orbits.
6. **Self-contained steps.** No state carried between steps beyond positions and velocities, except where the method explicitly defines it. This keeps Simulationarchive restarts exact.

## 2. Leads and sources

### 2.1 Single-core leads (ranked by expected gain × probability)

| # | Lead | Why it might work | Expected | Risk / notes |
|---|---|---|---|---|
| S1 | **Larger timestep at equal accuracy:** higher-order splittings compatible with Jacobi coordinates, GR and correctors: SABA/SBAB, ABAH(8,6,4)/(10,6,4) (Laskar & Robutel 2001; Farrés et al. 2013), Blanes–Casas–Farrés families | Farrés et al. report "nearly an order of magnitude less CPU" for La2010-class accuracy. Step cost rises by roughly the number of extra stages; dt can rise more. | 1.3–3x | Must pass §1.2, especially secular frequencies. Wisdom (2015): ≥ 17 steps per pericenter timescale. arXiv:2507.04987: dt ≤ 6.5 d converged for instability statistics, larger dt changes rates. The most promising route to 2x on **one** core. |
| S2 | **Kepler solver without iterations, or with one fewer.** E.g. an explicit high-order series start for the universal variable X(dt), from r₀, η₀, ζ₀, β (series reversion of Kepler's equation in Stiefel form), then a single Newton step with the convergence check. Or direct f and g series (Lagrange/Danby): for Mercury n·dt ≈ 0.43 rad at dt = 6 d, so ~15–17 terms, branch-free Horner. | Kepler is 63 % of the time, half of that in iterations. The paper's ceiling for removing one iteration is ~20 %. A branch-free polynomial could also shorten the latency chain. | 1.1–1.3x | Must keep the random-sign round-off (Fig. 2, 3) and compensated summation of the final Stiefel evaluation (Appendix C). Our it. 48 (third-order start, skips one Halley step) gained +8 % on compact systems and ±0 on the Solar System: start there. |
| S3 | **"Momentum" Kepler iteration** from upstream PR #881 (WHFast512 asm ablation: momentum Kepler 1.077x; fused + momentum + fast rsqrt 1.08x) | Already measured upstream, never tried in the C kernels | ~1.05–1.08x | Changes results at round-off; check §1.2 items 1 and 5. |
| S4 | **Shorter dependency chain.** The kernel is latency-bound: our measurement shows cost almost independent of the number of occupied lanes (8 planets 161–168 ns, 4×2 planets 148 ns). Split the 8 lanes into two 4-lane chains that overlap (as the AVX2 kernel does), reorder so the interaction step's throughput-bound rsqrt work overlaps the tail of the Kepler step, and check the Zen 5 latency tables for `vsqrtpd`/`vdivpd` (256-bit vs 512-bit). | On AVX2, two 4-wide chains overlap and the C AVX2 kernel nearly matches AVX-512 | 1.05–1.2x | Hardware-specific. Measure on Zen 4/5 and Intel. |
| S5 | **Transformations (11 %)** | Already a triangular 8×8 multiply with 2 accumulators per coordinate | ≤ 1.1x | Low ceiling (Amdahl). |
| S6 | **Interaction step (25 %)** | rsqrt14 + Newton refinement is already near the hardware limit. Our C kernel uses a deterministic bit-trick seed + 3 Newton steps + 1 FMA. Two Newton steps failed the GR test (it. 41). | ≤ 1.05x | Don't spend much here. |

Rejected or not useful (measured or argued in earlier rounds; don't repeat without a new idea):
- single precision in any part (the interaction error is not small enough);
- `-ffast-math` / re-association;
- `-fno-math-errno` (measured slower);
- 7/9-term Halley series (accuracy worse, it. 43);
- two Newton steps in the planet-planet rsqrt (it. 41);
- warm starts from the previous step (break self-contained steps);
- multirate or RESPA between inner and outer planets (all planets share one SIMD vector, so the lanes cost the same);
- tables/gathers.

### 2.2 Multi-core leads

**The hard limit: one step takes ~100 ns.**
- An OpenMP fork/join costs ~1 µs (measured), and moving one cache line between cores takes ~50–100 ns.
- **Splitting a single step across cores cannot win.** The paper's footnote 6 makes the same point with the light-travel time.
- Every multi-core idea must therefore parallelize **across time** (or across trajectories, which doesn't shorten one trajectory).

| # | Lead | Idea | Expected | Risk |
|---|---|---|---|---|
| P1 | **Parareal with a symmetric/symplectic variant** (Lions, Maday & Turinici 2001; Dai, Le Bris, Legoll & Maday 2013; Gander & Hairer 2014) | Coarse propagator G (WHFast512 at large dt, lower order, or no correctors) runs serially over K time slices. Fine propagator F (the production WHFast512) runs on each slice in parallel from G's predicted start. Corrections iterate until the slice boundaries agree to a tolerance. Slices can be 10⁶–10⁹ steps long, so communication is negligible. | Speedup ≲ cores / iterations. **1.5–2x on 5–8 cores** is the realistic hope. | The Solar System is chaotic (Lyapunov time ~5 Myr): the coarse propagator diverges, and iterations grow with slice length × Lyapunov exponent. Results match the serial run only to a tolerance, never bitwise. Plain parareal destroys symplecticity (use the symmetric variants). **Not found in the planetary-dynamics literature** (§2.4): a genuine research contribution if it works. |
| P2 | **Parareal on short horizons, repeated** | Apply P1 only over windows shorter than a few Lyapunov times (e.g. 1–5 Myr), then restart from the converged end state. This keeps the iteration count small. | as P1 | The overall wall time still includes the serial coarse sweeps. |
| P3 | **Better coarse propagators** | A secular (averaged) model (Laskar's secular equations, or a first-order WH with dt = 30–60 d) as G. Errors in the fast phases (mean longitudes) must be corrected each iteration; consider parareal on slow variables only. | unknown | Research; the most novel direction. |
| P4 | **Saha, Stadel & Tremaine (1997) style waveform relaxation**, revisited | Their Picard iteration gave 50x speedup at 1000x cost. The Kepler/perturbation structure (WH) could make the iteration contract much faster (perturbation ~1e-3). | unknown | Earlier discussion in this project: the cost comes from many Picard sweeps over long windows. Shorter windows plus a WH-structured iteration might cut it. |
| P5 | **Ensembles and paired integration** (for throughput, not single-trajectory speed) | `rebound.integrate_simulations` (thread pool) + `steps_pair` (two systems per kernel call: 1.22–1.3x per system on AVX-512) + `N_systems` packing for ≤ 4 planets | ~linear in cores; 1.2–1.3x per core from pairing | Already implemented. Use it for the **validation ensembles** of §1.2 (80 runs × 10 Gyr is the expensive part of this project). |

What doesn't work, measured or argued in this fork:
- OpenMP for small N;
- multithreading inside one step;
- GPUs for a single small system;
- IAS15-style stage parallelism (WHFast is explicit, nothing to iterate);
- speculative execution of whole steps with bitwise verification. The state is 48 doubles, and a prediction is never bit-exact.

### 2.3 Speed budget

To reach **2x vs 106 ns on Zen 5** (≈ 53 ns per step-equivalent at the same accuracy):
- **Single core:** S1 is the only lead with a plausible 2x on its own. S2–S4 together might give 1.2–1.5x. S1 × (S2–S4) is the most likely single-core path.
- **Multi-core:** P1/P2 on 5–8 cores might add 1.5–2x on top, if the accuracy contract holds.

### 2.4 Sources

**The target and its background**
- Dagli & Rein 2026, WHFast512 v2 (this directory): arXiv:2609.25339. Draft: https://hanno-rein.de/DagliRein2026.pdf
- Javaheri, Rein & Tamayo 2023, WHFast512 v1: https://arxiv.org/abs/2307.05683
- Rein & Tamayo 2015, WHFast: MNRAS 452, 376
- Rein, Tamayo & Brown 2019, high-order symplectic integrators: https://arxiv.org/abs/1907.11335
- Rein & Tamayo 2019, splittings (DHC vs Jacobi): RNAAS 3, 16
- Rein, Dey & Tamayo 2026, DHC spurious precession: arXiv:2601.08019
- Wisdom, Holman & Touma 1996, symplectic correctors: Fields Inst. Commun. 10, 217
- Mikkola 1997, Kepler solver notation: CMDA 67, 145
- Wisdom 2015, ≥ 17 steps per pericenter timescale: AJ 150, 127
- History of 5 Gyr wall times (Fig. 1 of the paper): https://shadden.github.io/nbody_history/
- REBOUND: https://arxiv.org/abs/1110.4876 and https://rebound.hanno-rein.de

**Splittings for larger dt (lead S1)**
- Laskar & Robutel 2001, SABA/SBAB: CMDA 80, 39
- Farrés, Laskar, Blanes et al. 2013, high-precision splittings for the Solar System: CMDA 116, 141, doi:10.1007/s10569-013-9479-6
- Hernandez & Dehnen 2017, correctors and SIMD: MNRAS 468, 2614
- arXiv:2507.04987: WH step-size convergence for Solar System instability statistics (dt ≤ 6.5 d converged; up to 32 d "physical" but instability rate 1 % → 2 %; energy error not sufficient)

**Solar System long-term physics (validation, §1.2 items 3 and 4)**
- Laskar & Gastineau 2009, Nature 459, 817
- Abbot et al. 2023, ApJ 944, 190
- Brown & Rein 2023, MNRAS 521, 4349
- Hoang, Mogavero & Laskar 2022, MNRAS 514, 1342 (100 Gyr secular)
- Brouwer 1937, AJ 46, 149 (random-walk round-off)

**Parallel in time (leads P1–P4)**
- Saha, Stadel & Tremaine 1997: https://arxiv.org/abs/astro-ph/9605016
- Parareal overview: https://Www.wikipedia.org/wiki/Parareal (Lions, Maday & Turinici 2001)
- Dai, Le Bris, Legoll & Maday 2013, symmetric parareal for Hamiltonian systems: https://arxiv.org/pdf/1011.6222
- Gander & Hairer 2014, parareal for Hamiltonian ODEs: https://www.unige.ch/~hairer/preprints/parareal-gh.html
- Le Bris, Princeton talk on parareal and long-time Hamiltonian integration: https://www.math.princeton.edu/node/15074
- Maday, Fields Institute (Toronto) 2019 talk: https://www1.fields.utoronto.ca/talks/Parareal-Method-to-speed-time-dependent-simulations
- Kreienbuehl et al. 2015, parareal in numerical relativity: https://arxiv.org/pdf/1509.01572
- Karim 2014 (AUB thesis), parallel-in-time for a perturbed Kepler orbit: https://scholarworks.aub.edu.lb/items/9ac8506d-c3c8-4369-b9ab-bade4e500a45
- van der Houwen & Sommeijer, parallel iterated Runge-Kutta (stage parallelism, used for IAS15 in this fork): https://ir.cwi.nl/pub/1663/1663D.pdf
- Stage-parallel IRK: arXiv:2505.17719
- https://parallel-in-time.org
- Parallel direct N-body for planetary systems (serial wins at small N): https://arxiv.org/abs/2208.13562

**Hardware**
- Zen 5 AVX-512 teardown (latencies, 512-bit datapath): https://numberworld.org/blogs/2024_8_7_zen5_avx512_teardown/
- Upstream WHFast512 PRs: #881 (ablation, momentum Kepler), #913/#916/#921/#922/#932/#933 (asm refinements)

## 3. Our fork (the starting point)

**Repository:** https://github.com/hippke/rebound, branch **`main` at `2cfc0c4`** (Oct 20, 2026). It contains:
- upstream REBOUND `262a9cce` (Oct 16, 2026), including WHFast512 v2 assembly (`src/integrator_whfast512.s`);
- the single-core performance work: 2.43x geomean over upstream on 14 use-case benchmarks;
- the WHFast512 C kernels;
- the multicore (OpenMP) work and the ensemble API.

Full logs:
- `../PERF_LOG.md` (single-core, see the WHFast512 sections and iterations 37–48);
- `../PARALLEL_LOG.md` (multicore);
- `../reports/03_literature_review.md` (leads across fields).

### 3.1 WHFast512 in the fork

| File | Content |
|---|---|
| `src/integrator_whfast512.s` | Upstream v2 assembly (AVX-512 only). `integrator.kernel = 1`. Default on AVX-512 CPUs, so old runs reproduce. Fork fix: restores MXCSR and adds `vzeroupper` on exit (it. 37). |
| `src/integrator_whfast512_kernel.h` | **The C port of v2: the basis for this project.** One source written against 8-lane vector operations (add, fma, sqrt, masks, permutes), instantiated per instruction set. |
| `src/integrator_whfast512_c_avx512.c` | `kernel = 2` |
| `src/integrator_whfast512_c_avx2.c`, `_c_avx2x4.c` | `kernel = 3` (8 lanes as 2×256-bit; 4-lane variant for ≤ 4 planets) |
| `src/integrator_whfast512_c_generic.c`, `_c_genericx4.c` | `kernel = 4` (plain C + `fma()`; ARM) |
| `src/integrator_whfast512_neon_ops.h` | NEON instantiation (aarch64) |
| `src/integrator_whfast512.c/.h` | Setup, constants cache, transformation matrices, kernel dispatch, `steps_pair`, synchronize. Settings: `gr_potential`, `corrector` (order), `N_systems`, `concatenate_steps`, `kernel`, `kernel_used`. |

**Properties of the C kernels:**
- **Bitwise identical across all instruction sets and machines:** AVX-512, AVX2, generic, NEON; checked on x86 and under qemu for aarch64.
- Not bitwise equal to the assembly:
  - 1/sqrt: deterministic seed + Newton instead of `vrsqrt14pd`, whose bits are unspecified;
  - after Newton, a first-order Stiefel update (it. 40);
  - plain Horner inside Newton (it. 42);
  - third-order start with a per-lane safeguard (it. 48).
  
  Each change was validated with energy ensembles, eccentric/hyperbolic orbits and IAS15 comparisons.
- **Speed** on our cloud Xeon with AVX-512, Solar System, ns per system-step:
  - assembly 163;
  - **C AVX-512 149**;
  - C AVX2 164;
  - C AVX-512 with `steps_pair` **122** per system.
  
  On AVX-512 CPUs the C AVX-512 kernel is 1.04–1.17x faster than the assembly.
- **Hardware note:** the user's workstation (Intel i5-13500) has **no AVX-512**. The assembly can't run there, and the C AVX2 kernel is the reference. For headline numbers comparable to the paper, use a Zen 4/Zen 5 machine (the paper used an AMD EPYC 9655). Measure kernel 1 (assembly) and kernel 2 (C AVX-512) side by side.

**Rule from the earlier round** (keep it): changes may alter results at round-off level, but **the C kernels must stay bitwise identical to each other on every machine** (fixed operation sequences, no ISA-specific shortcuts), and accuracy must be shown unchanged (§1.2).

### 3.2 Build

```
git clone https://github.com/hippke/rebound && cd rebound   # branch main
cd src && make -j OPENGL=0 SERVER=0 librebound.so           # C library (assembly + all C kernels)
cd .. && pip install -e . --break-system-packages            # Python package (optional)
# OpenMP build (only for multi-core experiments): REBOUND_OPENMP=1 pip install -e .   or   make OPENMP=1
```

The examples `examples/whfast512_solar_system` and `examples/whfast512_2_planets` show the C API.

## 4. How to test

### 4.1 Correctness (every iteration)

1. **Full test suite:** `cd benchmarks/perf && python3 check.py`. It rebuilds and runs pytest; currently 966 tests pass, ~4 min including the build. WHFast512-specific: `python3 -m pytest rebound/tests/test_whfast512.py`. It covers 75+ tests: every kernel, encounter/escape stops, GR, correctors, the cross-kernel bitwise test, `steps_pair` == separate.
2. **Cross-kernel bitwise hashes:** kernels 2, 3, 4 (and 3 × four-lane, 1/2/4 systems) must give identical final states. The test suite checks this; for new code paths, add cases.
3. **Kepler grid tests** (paper Fig. 2, 3, 7): one test particle over a grid of a ∈ [0.1, 1] AU and e ∈ [0, 0.9], dt = 5 d. Record relative a-error and its **sign** after 1 step and after 1000 orbits, plus iteration counts; then the high-e grid (e up to 1 − 1e-8, dt/P from 10^-3.1 to 10^-0.1). Write `whfaster/tools/kepler_grid.c`. The paper's figures are the reference.
4. **Development energy ensemble** (as in the earlier round): 16 Solar System runs × 10⁷ steps (dt = 5–6 d, GR, correctors), initial conditions perturbed by meters. Compare the median and the spread of |ΔE/E| with the assembly kernel. Earlier values: median 3.2–5.9e-13.
5. **Bias test** (paper Appendix C, the most important accuracy test for any change to the Kepler solver or series): 32 runs × 3×10¹⁰ steps, perturbed at the 1e-15 level. The ensemble-mean ΔE/E must be consistent with zero (paper: (0.4 ± 1.2)·10⁻¹¹ with both corrections). Expensive: run as an ensemble with `rebound.integrate_simulations` on many cores, or on a cluster.
6. **Physics checks** for any change of dt or splitting (lead S1):
   - frequency analysis of g₁…g₈, s₁…s₈ over 20–50 Myr against a reference run (dt = 2 d, correctors);
   - positions against IAS15 over 10⁴–10⁵ yr;
   - finally, instability statistics.
7. **For parallel-in-time methods (P1–P4):** convergence of the slice boundaries per iteration; final state vs the serial fine run (expect agreement to the parareal tolerance, not bitwise); and the same energy/bias/frequency checks as for any new integrator. Report the iteration count and the total compute.

### 4.2 Speed

- **Micro-benchmark:**
  - set up the Solar System;
  - set `integrator.kernel` and `integrator.concatenate_steps = 1000` (or `reb_simulation_steps(r, N)` with N ≥ 10⁶);
  - time ≥ 10⁸ steps;
  - report ns/step and the projected 5 Gyr wall time.
  
  Always **interleave** candidate and reference (A, B, A, B, …; ≥ 5 repeats; report min and median). Fix the CPU frequency if possible (`cpupower frequency-set -g performance`), pin to one core (`taskset -c 2`), and use a quiet machine. Laptop/desktop turbo and hybrid cores distort single-run numbers by ±10 %.
- **Suite:** `cd benchmarks/perf && python3 run.py ab --a <reference tree> --repeat 5` (interleaved A/B of the whole suite). The WHFast512 benchmarks are `ss_whfast512` (assembly), `w512c_ss_avx512`, `w512c_ss_avx2`, `w512c_3p_avx2`, `w512c_4x2_avx2`, `w512c_pair_ss512`, `w512c_pair_3p_avx2` (`bench.c`). Add a dedicated `whfaster_ss_gr_corr_6d` benchmark that matches the paper's setup exactly (dt = 6 d, GR, order-17 correctors).
- **Multi-core:** report wall-clock time **and** core-hours for a fixed physical span. Use a span long enough that parallel-in-time slices are realistic (≥ 10 Myr), plus the projected 5 Gyr figure. The user's benchmark harness for threads is `benchmarks/parallel/run_parallel.py` (OpenMP builds); for parallel-in-time, write a new driver.

## 5. Process (carried over from the earlier rounds; it worked)

1. **First reproduce the SOTA** on the target hardware: assembly kernel 1, dt = 6 d, GR, correctors. Record ns/step for kernels 1–3 and log it as iteration 0.
2. **Write the accuracy contract** (§1.2) into the log with numbers before changing any code.
3. **One work item per iteration:** implement, run the correctness checks, interleaved A/B, then decide. Keep if > 5 % faster with unchanged accuracy, or if small and safe; otherwise revert.
4. **Log every iteration**, including reverted ones, in `whfaster/RESEARCH_LOG.md`, plus a CSV like `../perf_history_w512c.csv`. Commit kept changes on a branch `whfaster` of the fork; push.
5. **Separate result-changing from bitwise-safe changes.** For result-changing ones, run the development ensemble every time and the expensive bias test before calling the work done.
6. Long validations (ensembles, frequency analysis) are batch jobs: use the ensemble API across all cores, and checkpoint (Simulationarchive).

## 6. Advice

- **Start with S1 (larger dt at equal accuracy).** It's the only single-core lead with a plausible 2x, and it changes the algorithm, not just instruction scheduling. The paper's authors argue convincingly that instruction-level tuning of the current algorithm has < 2x left. Check that the splitting keeps:
  - Jacobi coordinates (needed for correctors with GR);
  - an exact Kepler drift;
  - symplectic correctors (or an equivalent processing step).
  
  Also make sure the speed comparison is at **equal accuracy** (§1.2 items 1–4), not equal dt.
- **The step is latency-bound, not throughput-bound.** Wider vectors don't help; shorter dependency chains and more independent work per cycle do. That's why two systems per call (`steps_pair`) gives 1.22–1.3x per system, and why the AVX2 kernel nearly matches AVX-512.
- **Round-off bias is the trap.** The paper spent most of its accuracy effort on making round-off unbiased: compensated summation, coefficient corrections, masked per-lane convergence. Any new Kepler formulation must reproduce the random sign of the error (Fig. 2/3) and pass the Appendix C ensemble test. A fast but biased solver is worse than useless for Gyr runs.
- **Keep lanes independent:** a planet's result must not depend on its neighbours in the vector (masked convergence). Our it. 40 first broke this; the test `independent_kepler_solvers` caught it.
- **Hardware matters:**
  - Zen 5 has a full 512-bit datapath, Intel client CPUs have no AVX-512, and older Intel server parts (e.g. Skylake-SP) downclock under heavy AVX-512 use.
  - Measure on the paper's class of machine (Zen 4/5 EPYC or Ryzen 9000) for the headline, and also on AVX2 hardware (the C kernels are what most users without AVX-512 get).
  - Use instruction latency tables (uops.info, the Zen 5 teardown) when shortening chains.
- **Parallel-in-time is a research gamble with a clear payoff if it works:** a symmetric parareal around WHFast512, validated for unbiased energy and correct secular frequencies over Gyr, appears to be unpublished in planetary dynamics. Budget it as an exploration: prototype in C with a simple driver, measure iterations-to-convergence vs slice length, and stop early if the iteration count approaches the number of cores.
- **Don't promise bitwise reproducibility across methods,** only within one method and kernel family. Users of REBOUND care about restarts and reproducibility: document what changes.
- **Upstream policy:** upstream REBOUND does not accept AI-generated code (PR #941 was closed for this reason), so this work stays in the fork. Write clear logs so a human can re-derive any result for a paper.
- **Cost of validation:**
  - Each 5 Gyr run at ~100 ns/step takes ~9 h on one core.
  - The 80-run, 10 Gyr ensemble of the paper is ~1,400 core-hours.
  - Plan validation compute early; use shorter spans (10⁹ steps) and 16–32 runs during development.

---

## 7. Status after round 2 (Oct 21, 2026) and next steps

Round 2 (`RESEARCH_LOG.md`, "Round 2", iterations 1–14; literature: `LITERATURE_PINT.md`) found a working
multi-core route: **parareal in mean elements** (lead P3, generalised). Fine F = WHFast512 (unchanged); coarse G =
inner planets by a numerical first-order mean-element theory (ln a short-period transform on 2-D FFT grids, rebuilt
from the input state, Kepler + trapezoidal first-order secular drift), giants by WHFast512 at 360 d; correction in
Jacobi elements; energy projection at the boundaries; windows of 10 × 986 yr. C code: `pint/` (coarse step
1.0 ms = 8 % of a fine slice; pipelined driver `pint_pipe`; how to run: `PINT_README.md` / `pint/README.md`).

Measured: K = 3 iterations to the fine solver's own noise level per window; energy unbiased; secular frequencies
over 10 Myr within 3e-6 (K = 3), 3e-7 (K = 4), 1e-7 (K = 5) ″/yr of the serial run. Modelled (from measured costs
and iteration counts) wall-clock speedups on 10 cores: 2.45× (K = 3), 1.94× (K = 4), 1.6× (K = 5); 20 cores 2.4–3.3×.

Next steps, ranked:
1. **Measure `pint_pipe` on a ≥ 10-core machine** (Zen 4/5 EPYC: also gives r on the paper's hardware).
2. **Better contraction** (each saved iteration is worth ~25 % in S): the oracle (iteration 8) shows that an exact
   ∂(ln a, λ)/∂(all) of the inner planets would gain about one iteration — second-order (in mass) mean-element
   theory for Mars/Venus/Earth, or the λ and z transforms; the 64² grid anomaly (iteration 13) is unexplained.
3. **Sliding/overlapping windows** to remove the N r start-up and the idle cores at the end of each window
   (bound ≈ N/(K(1+r)): 3.1× on 10 cores at K = 3).
4. Validation at scale: bias ensembles (≥ 1e9 steps), 50 Myr frequency analysis on several ICs, and finally
   instability statistics with parareal ensembles on a cluster.
5. A cheaper coarse step (analytic ∂E/∂v for the secular rates, vectorised pair loops): r 0.083 → ~0.05 (+10 %).

**Update (iteration 15):** parameters are on a flat optimum (slice length, grids, giants' step, pairs: ≤ ±10 %);
the gain is an adaptive iteration count per window (stop when the largest boundary update ≤ 3e-6 AU): mean K = 3.3
with frequency fidelity ≤ 1e-7 ″/yr over 10 Myr → modelled **2.3× on 10 cores** at research-grade accuracy
(implemented in `pint_pipe`). Fixed a silent step-rounding bug in the giants' WH run.
