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

## 9. Round 3 (iterations 17–26): a much better coarse step — but not cheaper
Goal: K = 2 (each saved iteration ≈ 25–30 % of S). Result (Python prototype, 10 test windows over 10 Myr):
* The round-2 "first-order" theory had first-order errors: small divisors from osculating mean motions (the
  Venus–Earth 8:13 term), a too coarse Venus–Earth secular grid, Uranus/Neptune missing from the secular rates,
  aliasing of the giants' 60-yr term in the secular quadrature, no transform of λ, k, h, q, p. Plus two genuine
  second-order effects: the Earth–Mars–Jupiter resonant term and the second-order secular rates (Mars–Jupiter).
* With all fixes: **k = 2 ≤ 2 × noise for 6/10 windows (median 1.45; round 2: 0/10, median 245) and k = 3 at the
  noise floor for every window (worst 1.6–2.7 × noise; round 2: 870)**. The three hardest windows need the complete
  second-order short-period theory (a partial one made things worse).
* Cost: one all-component transform pass costs 2.1 ms in C (the whole config-B step 1.0 ms); the improved step needs
  ≥ 2 passes → r ≈ 0.2–0.4 → no wall-clock gain at ≈ 10 cores (break-even r ≈ 0.12). Continuous pipelining
  instead of windows (simulated with the real numerics) is not better than the windowed pipeline at P = 10.
* Recommendation unchanged: config B with the adaptive rule (≈ 2.3× on 10 cores). Round-4 plan
  (`RESEARCH_PLAN_R4.md`): make the accurate coarse step cheap (reusable/Taylor-corrected tables, analytic
  derivatives, pruned evaluation).

## 10. How to run
`pint/README.md`: build, `./pint_pipe 16 360000 3 <windows> 1 1 1 3e-6 6 8 2` (16 cores, adaptive, 8 coarse helpers;
kernel 1 = asm is fastest on the 7950X), `PINT_OUT=prefix` writes the trajectory. AVX2-only CPUs: kernel argument 3.
Long validation runs: `pint_pipe` with `PINT_OUT` or `pint_long`; analysis `freq_ana.py`, `energy_ana.py`,
`ens_ana.py` (ensembles), `tau_ana.py` (vs time-step references), `gyr_ana.py` (Gyr statistics).

## 11. Round 4 — measured on 16 cores (AMD Ryzen 9 7950X, iterations 28–38, Oct 21–22, 2026)
**Driver engineering (bitwise identical results):** the coarse step runs the giants' WH integration on the calling
core while helpers build the tables in a staged task graph (it. 30; 2.7 → 2.85×); persistent slice threads hand
windows over without a join (it. 31; → 2.95–3.05×). Measured: fine slice 8.3 ms under full load (7.0 ms idle: clock),
coarse step 0.26 ms on the chain; window = 16 coarse hops + K̄ × (slowest fine slice + hop).

**Rejected (measured):** SMT (32 threads), 20/24/32 slices, other slice lengths, paired fine solves (`steps_pair`
gives no throughput on Zen 4), mean-n tables (no K̄ gain), 64² close-pair grids (K̄ −10 % but too slow on the chain),
giants' coarse step 180–1000 d (no gain / no convergence), sliding windows (round 3 simulation).

**What binds:** always the phase λ of Venus, Earth or Mars at the last boundaries (Kepler-shear image of 3e-10
errors in ln a); the stop rule needs the k = 2 iterate within τ at all 16 boundaries.

**Energy bias found and fixed (it. 36):** WHFast512's inertial ↔ Jacobi matrices are not exact inverses; every
integration call (synchronization) adds dE/E ≈ +7e-17 for the Solar System. The parareal driver makes 10 calls per
slice → linear drift +7.7e-13 per Myr (+1.3e-10 after 167 Myr, ≈ 4e-9 per 5 Gyr), invisible in 10 Myr tests. New
opt-in `jacobi_recursive` in WHFast512 (recursive transformations; default unchanged and bitwise identical to
upstream) removes it (per-call offset 0 within ±1e-17); the drivers use it by default (`PINT_JR`). This also affects
ordinary serial WHFast512 runs with frequent outputs.

**Accuracy yardstick (it. 35, 37):** frequency errors of serial runs vs dt = 3 d: dt 6 d up to 5.1e-6 ″/yr (Mars g4;
scales as dt², a systematic bias), dt 8 d 1.0e-5, dt 12 d 2.6e-5. Parareal is unbiased (7-realisation ensemble) and
adds, on top of the dt = 6 d error: τ 3e-6 → 1e-7, τ 1e-5 → 4e-7, τ 3e-5 → 7e-7, τ 1e-4 → 2e-6 (all with K_min 3).
Allowing a stop at k = 2 (K_min 2) gives K̄ < 3 but a systematic Saturn g6 error (dt 8–12 d level).
