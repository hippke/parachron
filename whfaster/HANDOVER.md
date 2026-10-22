# Hand-over: multi-core speedup of a WHFast512 Solar System integration (parallel-in-time / parareal)

*Written Oct 21, 2026 by the previous agent (cloud session, 2 vCPUs, no access to many-core hardware) for an agent
with server access. Repository: https://github.com/hippke/rebound, branch **`whfaster`**, directory `whfaster/`.
The user's local copy is `/home/michael/Downloads/rebound/whfaster`; the 16-core test machine is `host2`
(AMD Ryzen 9 7950X, 16 cores / 32 threads, AVX-512, repo at `/home/ubuntu/rebound`).*

## Update after session 2 (Oct 21–22, 2026; agent with direct access to host2) — read this before §0–7
* Machines: `MACHINES.md` (host2 = `ssh root@23.239.118.238`; it can pull but not push → edit and push locally, then
  `git pull` on host2; a second worktree `/home/ubuntu/rb2` was used to build without disturbing running jobs).
* Iterations 28–38 in `RESEARCH_LOG.md` ("Phase 9"), summary `PINT_SUMMARY.md` §11. Results:
  1. Driver: staged coarse step with the giants on the calling core (it. 30) and persistent slice threads (it. 31):
     **2.95–3.05×** at τ 3e-6 (was 2.71×), bitwise identical trajectory. `pint_pipe` now writes trajectories
     (`PINT_OUT`), see the header of `pint/pint_pipe.c` for all environment options.
  2. **Energy bug in WHFast512 found and fixed** (it. 36): +7e-17 per integration call from the inertial ↔ Jacobi
     matrices → linear drift in parareal (10 calls per slice) and in any serial run with frequent outputs. Opt-in
     `jacobi_recursive` (src/integrator_whfast512.c); drivers use it by default (`PINT_JR=1`). Worth reporting
     upstream (default kept bitwise identical to upstream on purpose).
  3. Accuracy yardstick (it. 35/37): the dt = 6 d WH error itself is up to 5e-6 ″/yr (Mars g4), a systematic bias;
     parareal adds 1e-7 (τ 3e-6) … 2e-6 (τ 1e-4) and is unbiased → larger τ is defensible; τ = 1e-5 validated by a
     7-realisation ensemble (3.24×).
  4. Ceiling of this scheme on 16 cores ≈ 16/K̄ × 0.84 − coarse chain; K̄ cannot go below 3 with K_min = 3. K_min = 2
     reaches K̄ 2.1–2.7 with systematic shifts at the dt 8–12 d level (giants' coarse step 180 d removes the Saturn
     part) — it. 37/38.
  5. **Recommended settings and speeds** (it. 39–40, `PINT_SUMMARY.md` §11): reference accuracy `K_min 3, τ 3e-5`:
     3.46× in the 100-window benchmark, **3.0× over a full 1 Gyr run** (36.2 min; 5 Gyr ≈ 3.1 h vs 9.1 h serial, asm
     kernel); dt 6–8 d class 3.5–3.8×; dt 12 d class 4.7× (benchmark figures).
  6. 1 Gyr runs: energy drift of parareal +7.5e-10/Gyr (old transformation) → −1.1e-10 (fix, 10 calls per slice) →
     −6.7e-11/Gyr (fix, 2 calls per slice, default now) — serial runs +2.9…+6.4e-11/Gyr. CCD load balancing
     (`PINT_BAL`) was tried and rejected (it. 40).
  7. Ensemble 16 × 200 Myr (it. 46): dynamics (Mercury e, secular frequencies) indistinguishable from serial; energy
     drift −1.25e-10 per Gyr relative to serial remains (3.5σ) — the first thing to investigate next (lead (0)).
* Open leads (ranked): (0) the residual parareal energy drift of −1.25e-10 per Gyr (it. 46; not τ-dependent, it. 47):
  next step an energy budget per slice in `pint_pipe` (projection target vs telescoping, per-call and conversion offsets); (a) K = 2 at reference accuracy needs the accurate coarse step ("config C",
  `RESEARCH_PLAN_R5.md`): it reaches K̄ ≈ 2.4 at N = 16 (it. 41), but its measured parallel cost under load makes it
  break-even at best on 16 cores (it. 45) — worth it only on machines with more cores;
  (b) an ensemble of Gyr runs (serial vs parareal) to test instability/secular statistics, which single runs cannot;
  (c) per-core fine-time imbalance (second CCD 3–5 % slower; the slowest slice sets every iteration);
  (d) the remaining per-call energy offset with corrector + GR (−0.9e-17; inverse corrector is algebraically exact,
  so it is a round-off asymmetry, probably Kepler solver or GR kick); (e) report the WHFast512 synchronization bias
  upstream (user decision).

## 0. Read this first
1. **Goal (user):** make one long Solar System trajectory (Sun + 8 planets, WHFast512, Jacobi coordinates,
   dt = 6 d, GR potential, order-17 correctors) faster in wall-clock time by using several cores, at research-grade
   accuracy (unbiased energy, correct secular frequencies, unchanged instability statistics).
2. **Working rules (user):** follow `/home/michael/Downloads/rebound/PERF_LOG.md` (performance-work guidelines);
   update the .md documents in `whfaster/` after every iteration (`RESEARCH_LOG.md` = chronological lab book,
   `pint_history.csv` = one row per experiment); commit and push every change to `hippke/rebound` branch
   `whfaster` (retry on connect failures); iterate until a speedup is found; if ideas run out, do a new literature
   search and write a new research plan.
3. **Status in one paragraph.** Parareal with a cheap *mean-element* coarse propagator ("config B") works and is
   validated over 10 Myr. On the 7950X with 16 threads it gives **2.71× wall-clock speedup** (100 windows, adaptive
   iterations, parallel coarse step with 8 helper threads), measured by the user. Round 3 built a much more accurate
   coarse step (2 instead of 3–4 iterations for most windows) but it is 2–5× more expensive and therefore does not
   pay yet → the main open lead is **"make the accurate coarse step cheap"** (§6, lead A). Second lead: remove the
   remaining overheads of the pipelined driver on 16 cores (measured/model = 0.83, lead B).

## 1. Documents (all in `whfaster/`)
| file | content |
|---|---|
| `PINT_SUMMARY.md` | self-contained summary of the method, accuracy, speedup, bugs, rejected ideas (rounds 2–3) |
| `RESEARCH_LOG.md` | every iteration with numbers: round 2 = it. 1–16, round 3 = it. 17–26, it. 27 + "Phase 8" = 7950X runs |
| `RESEARCH_PLAN.md`, `RESEARCH_PLAN_R3.md` | plans of rounds 1–2 and 3 (with status sections at the end) |
| `RESEARCH_PLAN_R4.md` | draft round-4 plan: make the accurate coarse step cheap |
| `RESEARCH_PLAN_R5.md` | round-5 plan (Oct 22): "config C" = minimal accurate coarse step for K = 2 at N = 16, cost gate |
| `LITERATURE_PINT.md` | literature review (parareal, MGRIT, waveform relaxation, mean-element theories, …) |
| `pint_history.csv` | experiment table |
| `pint/README.md` | build/run of the C code |
| `tools/bench_workflow.yml` | template for a GitHub Actions self-hosted-runner benchmark workflow (optional) |
| `HANDOVER.md` | this file |
| `WHFAST512_SYNC_BIAS.md` | draft upstream bug report: WHFast512 energy change per synchronization, fix `jacobi_recursive` |
| `MACHINES.md` | test machines: host2 (7950X, `ssh root@23.239.118.238`) and the local i5-13500 workstation (no AVX-512) |

## 2. Code map
**C (production path), `whfaster/pint/`**
* `pint_coarse.c/.h` — coarse propagator G ("config B") and element conversions. API: `pint_default_cfg`,
  `pint_coarse` (serial), `pint_coarse_par(cfg, E_in, E_out, par)` + `pint_set_threads(n)` (task-parallel version,
  **bitwise identical** to the serial one), `pint_lna_disp` (diagnostic: ln a transform + secular rates).
  Internals: `table_par()` = grids per planet → one task per inner–inner pair and per packed pair of giants → GR per
  planet → reduction in the serial order; `outer_wh()` = giants by WHFast512 (runs concurrently on a persistent
  worker when parallel); helper pool with optional spinning (`PINT_SPIN`, default 0 = sleep; 0 was best on 7950X).
* `pint_pipe.c` — **the benchmark/production driver**: pipelined parareal, one thread per slice, windows chained.
  `./pint_pipe N slice_days K windows [fine_kernel=2] [serial=1] [pin=1] [tol_AU=0] [Kmax=6] [helpers=0] [parmode=1]`.
  tol > 0 → adaptive: at least K iterations, stop at the first iteration whose largest boundary update ≤ tol AU (≤ Kmax).
  helpers > 0 → coarse tables on a helper pool; parmode 1 = only in the k = 0 sweep, **2 = every iteration (best)**.
  The coarse step is skipped for slices whose start is already exact (n ≤ k − 1).
* `pint_long.c` — long validation runs writing every boundary state (`<prefix>_pr.bin`, `<prefix>_ser<kernel>.bin`,
  records of (t, E_tot, E[8][6]) doubles); analysis: `freq_ana.py` (NAFF secular frequencies), `energy_ana.py`,
  `drift_ana.py`. `pint_run.c` (non-pipelined parareal with per-iteration errors), `pint_test.c`, `test_par.c`
  (bitwise check + timing of the parallel coarse step), `ser_pert.c` (IC-perturbation references).
* `pint_disp6.c` — round-3 cost study: C version of the **all-component** first-order transform (validated vs
  Python to 1e-10), optional mean motions for the divisors; 2.13 ms per pass (vs 0.45 ms for the ln a pass).
* Build: `cd ../../src && make -j OPENGL=0 SERVER=0 librebound.so; cd ../whfaster/pint && make REB=$(pwd)/../../src pint_pipe libpint.so`
  (the Makefile default `REB=/home/claude/rebound/src` is the old cloud path).

**Python prototypes and diagnostics, `whfaster/parareal/`** (need `rebound` of this fork importable and `libpint.so`)
* `ptlib.py` (masses, ICs, element conversions, `prop()` fine propagation), `ptc.py` (ctypes wrapper of the C coarse
  step: `CoarseC`, `lna_disp`, `lna_disp_R`, `disp6`), `meanel.py` (numerical first-order theory, all components,
  `n_override` for divisors, `secular_rates_fast` with per-pair grids), `mgc.py` (Python coarse step with all
  round-3 options: `meann`, `nsub`, `sec_ng`, `gavg`, `w2`, hooks), `so2.py` (second-order resonant term U2 for one
  harmonic), `g3.py` (E–M–J resonant term added to any base G; `freeze`, `ode`), `w2pair.py`, `w2triple.py`
  (second-order transform W2 and secular rates U2 on pair/triple tori; `W2Set`), `gfilt.py` (numerical filtering of
  mean elements), `seccal.py` (empirical secular-matrix calibration).
* Test harness: `ph1eval.py` — nonlinear parareal (N = 10 × 986 yr) from 10 epochs of the 10 Myr serial run
  (`EP = [0,1030,2000,3400,4500,5500,6000,7000,8160,9300]` slice indices), prints error/noise at k = 1..3;
  noise = asm-vs-C fine-kernel difference. `ph1mn.py` builds coarse variants from spec strings (token list at the
  top of the file), e.g. the best round-3 configuration:
  `python3 ph1mn.py eps=0,1030 all_mn1_g3_dt180_nsub2_un_ga120_w2k_frzs` (≈ 1.5 s per coarse call in Python;
  run two processes with different `eps=` lists in parallel; one full 10-epoch evaluation ≈ 30–40 min).
* Diagnostics: `phase0.py` (Jacobian oracle; data `phase0_jac.npz`), `jacblock.py`/`jacblock2.py`/`jac_new.py`
  (which Jacobian rows/columns limit convergence), `slice_quality.py`/`slice_err.py` (per-slice errors),
  `resid_sample.py`/`resid_at.py`/`resid_spec.py` (dense fine runs, residual spectra of mean elements),
  `term_search.py`, `so2_scan.py`, `mj3body.py`/`mj3w2*.py` (three-body validation of W2), `longwin.py` (long
  windows), `pint_sim.py`/`pint_sim2.py` (event-driven simulators of continuously pipelined parareal).
* Data: `pint/runs/m10_ser2.bin` (10 Myr serial run, C kernel, 10001 boundary records; committed),
  `parareal/phase0_jac.npz` (committed). Large or regenerable files are git-ignored (e.g. `resid10k.npz`, 258 MB,
  from `resid_sample.py`).

## 3. Method (config B) in brief
Parareal U^{k+1}_{n+1} = G(U^{k+1}_n) + F(U^k_n) − G(U^k_n) on windows of N slices × 986 yr (360 000 d), correction
in osculating Jacobi elements [ln a, λ, k = e cos ϖ, h = e sin ϖ, q = sin(i/2) cos Ω, p = sin(i/2) sin Ω], angles
wrapped, energy projection of corrected boundary states. G: inner planets by a numerical first-order mean-element
theory (pair accelerations on 32²/16² grids of mean longitudes → FFT → short-period ln a displacement; live tables
from the input state; Kepler + first-order secular drift, Heun), giants by WHFast512 at 360 d in a Sun(+inner
masses)+4-planet subsystem. Why mean elements: osculating a depends on the planets' phases and Kepler shear amplifies
any error in the ln a rows of the coarse Jacobian (round-1 failure, `RESEARCH_LOG.md` it. 1).

**Speedup model** (validated within 7–17 % on the 7950X): S = N / (N r + K̄ (1 + r)), r = coarse/fine cost per
slice, K̄ = mean iterations per window. The N r term is the serial k = 0 coarse sweep.

## 4. Results so far
**Accuracy (2.1 GHz Xeon runs, 10 Myr, two IC sets):** adaptive rule (K ≥ 3, stop when max boundary update
≤ 3e-6 AU) → inner-planet secular frequencies within 8.6e-8 ″/yr of serial (WH truncation error ~3–9e-7), energy
unbiased. Fixed K = 3: 3.3e-6 ″/yr (not research grade); K = 4: 3.2e-7.

**7950X (user measurements, `pint_pipe ... 3e-6 6` unless noted):**
| N | windows | variant | K̄ | speedup |
|---|---|---|---|---|
| 10 | 20 | serial coarse | 3.50 | 1.82 |
| 16 | 20 | serial coarse | 3.90 | 2.37 |
| 12 | 100 | serial coarse | 3.47 | 2.07 |
| 16 | 100 | serial coarse, tol 1e-5 | 3.45 | 2.51 |
| 16 | 100 | serial coarse (new build) | 3.80 | 2.36–2.38 |
| 16 | 100 | 8 helpers, parmode 1 | 3.80 | 2.37–2.42 |
| 16 | 100 | 8 helpers, parmode 2, spin 20000 | 3.80 | 2.67 |
| **16** | **100** | **8 helpers, parmode 2, PINT_SPIN=0** | 3.80 | **2.71** |
| 16 | 100 | 15 helpers, parmode 2 | 3.80 | 2.64 |
Fine slice 7.0 ms on the 7950X (12 ms on the Xeon); serial coarse step ≈ 0.84 ms → r ≈ 0.12 (the AVX-512 fine
kernel gains more from the fast core than the scalar coarse step). Parallel coarse (mode 2): coarse core time
5.17 → 2.45 core-s. Final position difference vs the serial run 2.8e-4 AU after 1.58 Myr (tol 3e-6) and 6.6e-4
(tol 1e-5): same order → round-off divergence, but **tol 1e-5 is not validated** (needs the 10 Myr frequency check).
Measured/model ≈ 0.83 with the parallel coarse step; ≈ 4.4 fine solves per slice incl. speculative ones → even at
r = 0 the bound is ≈ 16/4.4 ≈ 3.6 at N = 16 with this K̄.

**Not yet run (suggested to the user):** `./pint_pipe 14 360000 3 100 2 1 1 3e-6 6 4 2`,
`./pint_pipe 12 360000 3 100 2 1 1 3e-6 6 4 2`, `./pint_pipe 12 360000 3 100 2 1 1 3e-6 6 8 2`.

## 5. Knowledge that is not obvious (lessons, mostly from round 3)
* **Small divisors must use mean, not osculating, mean motions** (one fixed-point iteration: compute mean ln a,
  rebuild/re-evaluate the tables). The Venus–Earth 8:13 term (divisor 2π/240 yr) made the pointwise mean ln a wrong by
  up to 5.5e-7 with osculating n; with mean n 9e-9. Divisor-only re-evaluation differs from a full table rebuild at
  mean a by ~1e-8 in ln a (not tested in parareal).
* **Secular rates:** Venus–Earth grid needs ≥ 48² (32²: 0.5 % error), Earth–Mars ≥ 32², others 16–24²; **Uranus and
  Neptune must be included** in the inner planets' secular rates (they were missing; ~1e-3 of Jupiter's effect on q, p).
* **Giants' 60-yr term** (2λS − λJ) aliases into the inner secular drift if the giants are point-sampled; box-average
  their elements over ~120 yr around each quadrature node (`mgc gavg`). Giants' WH step 180 d instead of 360 d
  helped a lot once the inner planets were accurate (k2 median 10 → 3.5); 90 d = 180 d.
* **All six elements need a short-period transform** for K = 2 (ln a only: 1/10 windows at k = 2; all: 6/10).
* **Second-order effects that matter:** the Earth–Mars–Jupiter near-resonance −4λE + 8λM − 3λJ (P ≈ 1800–2000 yr,
  Mars ln a amplitude ≈ 1.5–2.2e-7; must be resonant dynamics, not transformed away; `so2.resonant_U2` reproduces the
  measured amplitude to 7 % but the phase only to ~30°) and **second-order secular rates U2** (dominated by the
  Mars–Jupiter same-pair term; U2·T equals the empirically fitted per-slice secular-matrix error δz exactly).
* **Mean ↔ osculating maps must be used in one consistent form.** Config B uses the "X-form" (W1 evaluated at the
  osculating angles both ways). The second-order W2 (`w2triple`) belongs to the "Y-form" (X = Y + W1(Y) + W2(Y)):
  switching forms with an incomplete W2 made Venus/Earth worse. U2 (secular part) alone is safe.
* **Exclude the giants' own displacements** from second-order products: G already uses the giants' osculating
  elements (their J–S great-inequality W1 produced huge spurious terms).
* **Finite-difference derivatives on top of finite-difference dE/dv are noisy**; G must be smooth at the 1e-10 level or
  the k = 3 floor rises (FD step 1e-4 in `so2.Pair`, not 1e-6).
* **Freezing:** second-order coefficients frozen per *window* are bad (U ∝ e·e^{−iϖ} rotates ~50° per window);
  frozen per *slice* (from the first coarse call of that slice, reused in later iterations) is as good as live.
* **The k = 0 prediction quality is essential**: using config B for the k = 0 sweep and the accurate G afterwards
  loses the K = 2 benefit (1/10).
* **Longer windows need more iterations** (binomial error growth): config B needs K = 4 for 40-slice windows.
* **Continuous pipelining** (no windows; three simulator variants with the real numerics) was not better than the
  windowed pipeline at P = 10 (S 1.4–2.2): every change of a slice's start after its iterations began costs
  iterations. Upper bound of any pipelined parareal: S ≤ min(P/(K(1 + r)), 1/r).
* Empirical calibration of G from previous windows, numerical filtering of mean elements, smaller/larger slices,
  secant parareal, 64² grids everywhere, small-divisor cutoffs: tried, no systematic gain (see the log).
* Practical: `pkill -f <pattern>` can kill its own shell; thread-local scratch is required in the coarse step
  (it is called concurrently by the slice threads); drivers reject slices that are not multiples of 6 d.

## 6. Open leads (ranked) with first steps
**A. Make the accurate coarse step cheap (round-4 plan, `RESEARCH_PLAN_R4.md`).** The accurate G (Python prototype
`all_mn1_g3_dt180_nsub2_un_ga120_w2k`) gives k2 ≤ 2 × noise for 6/10 windows and k3 at the noise floor everywhere
(config B: 0/10 and worst 870), i.e. K̄ ≈ 2.4 instead of 3.3–3.8 — worth ≈ 30–40 % in S *if* r stays ≲ 0.12 (at
10–16 cores; with the parallel coarse step the effective r is about half the serial one). Components and cost in C:
two all-component transform passes (start, end; 2.13 ms each serially), mean-n re-evaluation, giants at 180 d with
dense sampling and box averaging, RK4 secular drift (2 sub-steps) with the giants' samples and U/N terms, the E–M–J
term (`so2`), second-order secular rates U2 (`w2triple`, pairs + 3 triples) — the last two can be frozen per slice
(computed once at k = 0). Steps:
1. Instrument `pint_disp6.c` (grid / FFT / evaluation shares). Replace the finite-difference dE/dv by analytic Gauss
   equations; prune Fourier coefficients < 1e-12 before evaluation; pack two real grids per complex FFT for all six
   components; vectorise the pair loops. Target ≤ 0.5 ms per all-component pass.
2. Parallelise it with the existing task pool (`table_par` pattern; per-component tasks).
3. Test in Python (ph1eval) whether table reuse is allowed: tables per slice from the k = 0 call + first-order
   Taylor correction in the slow elements; divisor-only mean-n; then port the minimal set to C as a new cfg mode
   ("config C") behind a flag, keeping config B the default.
4. Port `so2.resonant_U2` (E–M–J) and the U2 pairs to C, computed once per slice (k = 0) on the helper pool.
5. Measure K̄ and r with `pint_pipe` on the 7950X; validate with `pint_long` + `freq_ana.py` over 10 Myr before
   recommending it.
**B. Remove the remaining driver overheads (measured/model = 0.83).** Profile `pint_pipe` on the 7950X (`perf`,
per-thread timelines). Candidates on the critical chain: energy projection and `energy_el` (simulation creation per
call), `poserr_el` for the stopping rule, thread creation per window (`pthread_create` for every window — use a
persistent pool), speculative fine chunks (10 chunks per slice; check that abandoned work does not delay the chain),
helper threads sharing physical cores with fine solves (try pinning helpers to SMT siblings explicitly, or N = 14 + 2
cores for helpers). Count fine solves per slice (≈ 4.4 now) and reduce speculation if it costs wall time.
**C. Window size and stopping rule on 16 cores.** Run N = 12, 14, 16 (and 20, 24, 32 with SMT) with helpers
(parmode 2); study the stopping rule (K ≥ 3 is forced; tol 1e-5 vs 3e-6; a per-window rule using the k = 2 update
may allow K = 2 for some windows). Every new setting must pass the 10 Myr frequency/energy contract
(`pint_long` + `freq_ana.py`, protocol of `RESEARCH_LOG.md` iterations 12/14).
**D. Cheaper config B without changing the numerics.** Analytic dE/dv; the giants' WH run is ≈ 1/3 of the step and
is now the critical part of the parallel step (test 720 d — round 2 found the optimum flat between 180 and 720 d,
but re-check the iteration counts); overlap the end-table inner–inner tasks with the giants' run (they only need
the inner mean elements).
**E. Complete second-order short-period theory** for the three hard test windows (Mars at 3.4/4.5 Myr,
Venus/Earth at 8.2 Myr) — only worth it if A succeeds.
**F. Many-core scaling study**: with the parallel coarse step, measure S(N) for N = 8…32 and core-hours ratio;
document in `pint_history.csv`.

## 7. Suggested first session on the server
1. `git pull` on `whfaster`, build (`make REB=... pint_pipe libpint.so test_par`), run `./test_par 8` (must print
   "bitwise identical: yes").
2. Reproduce the best number: `./pint_pipe 16 360000 3 100 2 1 1 3e-6 6 8 2` (expect ≈ 2.7×; serial ≈ 11.3 s).
3. Lead C quick scan (N = 12/14/16 with 4–8 helpers), then lead B profiling, then lead A step 1.
4. Log every run in `RESEARCH_LOG.md` (section "Phase 8") and `pint_history.csv`; commit and push.
