# Multi-core speedup for one WHFast512 Solar System trajectory — summary (Oct 21, 2026)

This is the self-contained summary of the parallel-in-time work (whfaster round 2, iterations 1–16). Details,
all numbers and every rejected attempt: `RESEARCH_LOG.md` ("Round 2"); literature review: `LITERATURE_PINT.md`;
plans: `RESEARCH_PLAN.md` (rounds 1–2), `RESEARCH_PLAN_R3.md` (round 3: second-order theory); history table:
`pint_history.csv`; code: `pint/` (C), `parareal/` (Python prototypes and diagnostics).

## 1. Problem and constraint
* Target: one long trajectory (Sun + 8 planets, WHFast512, Jacobi coordinates, dt = 6 d, GR potential, order-17
  correctors), faster in wall-clock time by using several cores, at research-grade accuracy (accuracy contract,
  `RESEARCH_PLAN.md` §1.2: unbiased energy, correct secular frequencies, unchanged instability statistics).
* One step takes ~100–200 ns; splitting a step across cores cannot win (fork/join ≈ 1 µs; even with 20 ns
  same-CCD core-to-core latency a two-core split models at ≈ 0.7×). Parallelism must be across time.
* Round 1 had shown that parareal with a WH coarse solver fails (the coarse solver must resolve Mercury: cost
  ≥ 0.4 of the fine one; speedup < 1).

## 2. Why round-1 parareal failed (iteration 1)
The slice Jacobians (48 × 48, Jacobi elements) of fine F and coarse G predict parareal convergence exactly
(linearised recursion e^{k+1}_{n+1} = G' e^{k+1}_n + (F' − G') e^k_n). The killing mismatch is in the semi-major-axis
rows: osculating a depends on the planets' phases through short-period terms (∂ ln a/∂λ ~ 1e-6–1e-5) and Kepler
shear (∂λ/∂ ln a ≈ 4e4 per kyr for Mercury) amplifies any error there. A coarse model needs these rows to ≲ 1e-8
(absolute) and the rest to ≲ 1 %. Impossible in osculating elements with a cheap model; natural in **mean elements**.

## 3. The method that works
**Parareal in mean elements with a hybrid coarse propagator** (`pint/pint_coarse.c`, "config B"):
* Inner planets (Mercury–Mars): numerical first-order theory. For each planet pair (6 inner–inner pairs on 32²
  grids, 8 inner–Jupiter/Saturn pairs on 16² grids of mean longitudes) the first-order Jacobi perturbing
  accelerations give d ln a/dt = 2a/μ · v·a_pert; its 2-D FFT gives the short-period ln a displacement f_k/(i k·ν);
  GR as a 1-D term. Mean ln a = osculating − displacement. Tables are rebuilt from the actual input state at every
  call ("live": the slow-element dependence of G matters, k 6 → 3). Mean elements are drifted by Kepler mean motion
  plus first-order secular rates (grid means), trapezoidal (Heun) in time.
* Giants: WHFast512 at a 360-d step in a Sun(+inner masses) + 4-planet subsystem.
* Correction in osculating Jacobi elements (angles wrapped); energy projection of corrected boundary states
  (each slice's energy increment = one fine round-off realisation).
* Windows of N = 10 slices × 986 yr (360 000 d), chained; K iterations per window.
* Cost (measured, one core, 2.1 GHz Xeon): fine slice 11.5–12.5 ms, coarse step 1.0 ms → r = 0.083–0.093.

## 4. Accuracy (10 Myr validations against serial runs, two sets of initial conditions)
| iterations per window | worst inner-planet secular-frequency shift [″/yr] | energy |
|---|---|---|
| fixed K = 3 | 3.3e-6 (≈ 10 × the WH truncation error at 6 d) | unbiased |
| fixed K = 4 | 3.2e-7 (≈ WH error) | unbiased |
| fixed K = 5 | 1.2e-7 | unbiased |
| **adaptive: K ≥ 3, stop when max boundary update ≤ 3e-6 AU (mean K = 3.29)** | **8.6e-8** | unbiased |

References: asm vs C kernel (pure round-off) ≤ 1.5e-8 ″/yr; WH truncation error at 6 d (round 1) 3–9e-7 ″/yr.
Energy: parareal minus serial within the spread of two serial round-off realisations; boundary energy jumps at the
round-off level, no detectable bias.

## 5. Speedup (modelled from measured iteration counts and costs)
S = N / (N r + K̄ (1 + r)), pipelined schedule (validated on 2 threads: 0.77–0.86 measured vs 0.86 predicted).

| configuration | 10 cores | core-hours vs serial |
|---|---|---|
| adaptive rule (research-grade, recommended) | **≈ 2.3×** | ≈ 4.3× |
| fixed K = 4 | 1.94× | 5.2× |
| fixed K = 3 (frequencies only to ~3e-6 ″/yr) | 2.45× | 4.1× |
| 5 cores, adaptive | ≈ 1.4× | ≈ 3.6× |
| 20 cores (K + 1) | 2.4–3.3× | |

Not yet measured on a ≥ 10-core machine (this work ran on a 2-vCPU host). Program: `pint/pint_pipe` (README).

## 6. Bugs found and fixed
* Thread-unsafe static scratch arrays in the coarse step (races under the pipelined driver) → thread-local.
* Euler secular step: systematic bias of the inner planets' secular frequencies at K = 3 → Heun.
* Giants' WH step not dividing the slice was silently rounded (G integrated the wrong span) → step adjusted to
  divide the slice; drivers reject slices that are not multiples of 6 d.
* Prototype scripts: secant baseline slicing bug; garbled calibration update in `pint_long`.

## 7. What was tried and rejected (all measured)
Shorter/longer slices (0.25–1.5 kyr), previous-window calibration of G, secant-enhanced (quasi-Newton) parareal,
small-divisor cutoffs, larger grids (64²) everywhere or for close pairs only (better convergence, cost too high),
Uranus/Neptune pairs, giants' steps 180–720 d, an empirical Earth–Mars–Jupiter resonant term. Parameters are on a
flat optimum (±10 %).

## 8. What limits further gains (iteration 16)
Typical windows reach the fine solver's noise floor at k = 3; at k = 2 they are 1e2–1e4 × above it. Remaining errors
of the first-order theory: Venus–Earth high-harmonic aliasing (fixable, not worth its cost) and, for Mars, a 1836-yr,
2.2e-7 oscillation of the mean semi-major axis from the **Earth–Mars–Jupiter three-body near-resonance
4λ_E − 8λ_M + 3λ_J (+ϖ)**, second order in the masses. A substantial further gain requires K = 2 → a complete
second-order mean-element theory with this term as resonant dynamics: `RESEARCH_PLAN_R3.md`.

## 9. How to run
`pint/README.md`: build, `./pint_pipe 10 360000 3 <windows> 2 1 1 3e-6 6` (10 threads, adaptive), AVX2-only CPUs:
fifth argument 3. Long validation runs: `pint_long`, analysis `freq_ana.py`, `energy_ana.py`, `drift_ana.py`.
