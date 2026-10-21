# WHFaster research log

Plan: `RESEARCH_PLAN.md`. Fork: hippke/rebound `main` @ `2cfc0c4`. Tools in `whfaster/tools/`.

## Machine used in this round

Cloud container: Intel Xeon @ 2.10 GHz (Sapphire-Rapids-class: AVX-512F/BW/DQ/VL/IFMA/VBMI), **2 vCPUs**, 8 GB RAM,
shared host (timing noise ±5–10 % between runs; always interleaved). gcc, `-O3 -march=native` for tools.
No Zen 4/5 machine was available in this session, so all speed numbers are ratios against kernel 1 (the paper's
assembly) measured in the same session on the same core. The user's i5-13500 has no AVX-512 (assembly can't run there).

## Iteration 0 — reproduce the SOTA (paper setup)

Setup: Sun + 8 planets (examples/whfast512_solar_system ICs), dt = 6 d, `gr_potential=1`, `corrector=17`,
`concatenate_steps = 10^7`, 5 interleaved repeats (`tools/bench_w512.c`).

| kernel | min ns/step | median ns/step | 5 Gyr projection (this CPU) | vs asm |
|---|---|---|---|---|
| 1 asm (paper v2) | 214.5 | 224.3 | 18.1 h | 1.00 |
| 2 C AVX-512 | 196.9 | 204.9 | 16.6 h | 0.92 (1.09x faster) |
| 3 C AVX2 | 213.5 | 216.8 | 18.0 h | 1.00 |

(The plan quotes 163/149/164 ns on "our cloud Xeon"; this instance is ~30 % slower in absolute terms, the ratios
agree. The paper's Zen 5 gives 106 ns for kernel 1. Conversion: `h_5Gyr = ns × 3.04e11 / 3.6e12` at dt = 6 d.)

## Accuracy contract (agreed before changing code; §1.2 of the plan, with the numbers used here)

1. **Unbiased energy:** ensemble mean of ΔE/E (incl. the GR potential −3 m0² m_i /(c² r_h²)) consistent with 0 over
   ≥ 10^9 steps; median ∝ √N.
2. **Energy level:** development ensemble (16 runs × 10^7 steps, perturbed ICs) median |ΔE/E| no worse than the
   reference kernel at dt = 6 d (earlier: 3.2–5.9e-13).
3. **Secular dynamics:** errors of the secular frequencies (Mercury g1, s1 and the others) relative to a converged
   reference must be **≤ those of the reference method (WH + corr17) at dt = 6 d**. Measured as the drift rate of the
   ϖ_i / Ω_i phase difference against a converged reference from identical ICs (arcsec/yr), plus mean-longitude
   phase errors. Energy alone is not accepted.
4. **Instability statistics:** final validation only (not affordable on 2 vCPUs in this session).
5. **Kepler solver quality:** grid tests (paper Fig. 2/3/7) for any Kepler change.
6. **Self-contained steps** (positions, velocities only).

Equal-accuracy speed comparisons use criterion 3 (the binding one for dt changes) and 1–2.

## Lead S1 — larger timestep with higher-order splittings (REJECTED, measured)

Tool: `tools/s1_acc.c` + `tools/s1_analyze.py`. Solar System (no GR; GR does not change integrator-error scaling),
100 kyr, identical ICs, samples every 360 000 d (exact step counts, so all dt hit the same times). Reference:
WHCKL (lazy-implementer kernel + corr17 + corrector2) at dt = 1 d (converged: WHCKM at 1–3 d agrees with it to ≤ 1e-10″/yr in g1).
Metrics: λ1 = Mercury mean-longitude error at 100 kyr; secular drift rates of ϖ1, ϖ2, Ω1, ϖ5 errors (″/yr, linear fit);
max |ΔE/E|. Cost = sub-steps (Kepler + interaction) per unit time relative to WHC @ 6 d.

| method @ dt | cost | λ1 err [rad] | dg1 [″/yr] | dg2 | ds1 | dg5 | max ΔE/E |
|---|---|---|---|---|---|---|---|
| WHC (= WHFast512) @ 5 d | 1.20 | 1.8e-4 | −5.9e-7 | −6.9e-7 | 2.6e-7 | 2.6e-7 | 2.1e-12 |
| **WHC @ 6 d (paper)** | **1.00** | **−2.9e-3** | **−8.8e-7** | −9.6e-7 | 3.4e-7 | 3.9e-7 | 2.4e-12 |
| WHC @ 7.2 d | 0.83 | −1.5e-3 | −1.2e-6 | −1.4e-6 | 5.1e-7 | 5.5e-7 | 7.3e-12 |
| WHC @ 8 d | 0.75 | 0.10 | −3.4e-7 | −3.3e-6 | 1.5e-6 | 3.9e-7 | 7.5e-11 |
| WHCKL @ 6 d | 1.36 | −2.9e-3 | −2.9e-8 | 3.0e-8 | −3.6e-8 | 5.7e-9 | 2.8e-12 |
| WHCKL @ 7.5 d | 1.09 | 1.5e-2 | 1.6e-7 | −1.9e-7 | 1.2e-7 | −4.0e-8 | 2.3e-11 |
| WHCKL @ 8 d | 1.02 | 0.10 | 1.2e-6 | −1.6e-6 | 8.6e-7 | −2.8e-7 | 7.5e-11 |
| SABA2 @ 12 d | 1.00 | −9.6e-3 | −5.7e-7 | −4.4e-7 | 1.2e-7 | 2.3e-7 | 1.2e-11 |
| SABA3 @ 18 d | 1.00 | 9.6e-3 | −5.2e-7 | −7.2e-7 | 3.7e-7 | 2.7e-7 | 1.1e-11 |
| SABA3 @ 12 d | 1.50 | 2.0e-4 | −2.8e-7 | −3.3e-7 | 1.2e-7 | 1.3e-7 | 8.3e-13 |
| SABA(10,6,4) @ 24 d | 2.00 | 9.3e-4 | 1.2e-8 | −1.7e-8 | −2.5e-10 | −4.4e-9 | 1.5e-12 |
| SABA(10,6,4) @ 36 d | 1.33 | 1.8e-2 | 2.1e-7 | 6.0e-8 | 1.6e-7 | −1.3e-8 | 2.1e-11 |
| SABA(10,6,4) @ 48 d | 1.00 | −0.25 | −2.5e-6 | 2.2e-6 | −2.5e-6 | 6.8e-7 | 2.1e-10 |
| SABAH(8,6,4), SABACL4 | (same pattern as SABA(10,6,4)) | | | | | | |

Findings:
1. **WHFast512 (C kernel) = WHC** in all metrics at 6 and 8 d (cross-check `w512_6`, `w512_8`).
2. **Beyond ≈ 6 d the error is dominated by step-size resonances with Mercury's eccentric orbit**, not by the
   ε²dt² term: the λ1 error is identical for WH and the lazy/modified-kick kernels, changes sign non-monotonically
   (5 d: 1.8e-4, 6 d: −2.9e-3, 6.4 d: 9.6e-3, 7.2 d: −1.5e-3, 8 d: 0.10 rad) and its envelope grows ~dt^10. Energy
   errors show the same jump (2.4e-12 → 7.5e-11 from 6 to 8 d). It also reaches g1 from 8 d on.
3. Better kernels (WHCKL/WHCKM) reduce the secular-frequency error 30× at fixed dt, but the resonance error limits
   dt to ~7.5 d at equal accuracy (1.25× fewer steps) for +36 % step cost (two interaction steps) → net loss.
4. Multi-stage SABA/ABA splittings sample the kicks at Gauss nodes and suppress the aliasing per *step*, but **per
   sub-step** (the unit of cost: one Kepler + one interaction) nothing beats WHC @ 6 d: at equal cost they trade
   λ1/energy error for g-error. SABA(10,6,4) @ 24 d is better than WHC @ 6 d on every metric but costs 2×.
5. Note for the plan: dt = 5 d has a 16× smaller λ1 error than 6 d over 100 kyr (Mercury aliasing), at 1.2× cost.

Conclusion: **no splitting gives ≥ 1.2× at equal accuracy; S1 cannot deliver 2×.** The binding limit is Mercury's
88-day orbit (resolution ≈ P1/15), consistent with arXiv:2507.04987.

## Leads P1/P2 — parareal (REJECTED, measured)

Prototype `parareal/pr.py`: F = WHFast512 C kernel, dt = 6 d, GR, corr17; G = WHFast512 at dt_G. Error = max position
error (AU) at slice boundaries vs the serial F run; tolerance = round-off noise of F itself (asm vs C kernel from the
same ICs: 1.3e-9 AU after 1 kyr, 7e-8 AU after 10 kyr, 7e-6 AU after 100 kyr).

* Plain parareal (Cartesian correction): **diverges** (errors up to 1e5 AU): G's errors are phase (mean-longitude)
  errors of O(0.1–1 rad), and additive Cartesian corrections across a phase error produce unphysical states.
* Corrections in Jacobi orbital elements (a, λ, e cos ϖ, e sin ϖ, sin(i/2) cos Ω, sin(i/2) sin Ω), angle-wrapped:
  bounded, and converges if G resolves Mercury:

| slices × length | dt_G | r = cost G/F | iterations to reach the F noise level |
|---|---|---|---|
| 10 × 986 yr | 24 d | 0.26 | no contraction (≈ 1 slice per iteration) |
| 10 × 986 yr | 16 d | 0.42 | 4 iterations: 9/11 slices |
| 10 × 986 yr | 12 d | 0.55 | 4 iterations: 10/11 slices |
| 10 × 986 yr | 8 d | 0.77 | 3 iterations: all |
| 10 × 9 856 yr | 12 d | 0.52 | 5 iterations: 7/11 (no contraction) |
| 10 × 9 856 yr | 16 d | 0.43 | none |

Speedup model S = N / (k + (k+1) N r) → best case here 0.3–0.5 (a slowdown). Parareal needs r ≲ 0.05 *and* k ≈ 2–3;
any WH-type coarse propagator must resolve Mercury's 88-d orbit (dt_G ≲ 12–16 d, r ≥ 0.4), and contraction degrades
with slice length. A secular coarse model (P3) would have r ~ 1e-3 but lacks the short-period terms that dominate the
Jacobian mismatch in λ; untested, research-level. **Multi-core (parallel-in-time) speedup of one trajectory: not
achievable with the methods tested.**

## Single-core work on the Kepler solver (leads S2/S3/S4)

### Method
`lab/`: variants of the C kernel (`kernel_*.h`, compiled with macros and their own suffix) are run on the same
`simd_data` (set up by the library: Solar System, dt = 6 d, GR, corrector 17), interleaved, 200–400 rounds of 3·10^5
steps, pinned to one core. Reported: median of the **paired ratio** (variant time / reference time in the same round),
which cancels the slow frequency drifts of this shared host (single-run noise ±10–30 %). Iteration counts via
counters (`-DW512_COUNT`).

### Where the time goes (stock C kernel, dt = 6 d)
* Kepler-only loop 0.61, interaction-only loop 0.47 of the full step (they overlap slightly).
* Every planet but Mercury is exact after the start + one Halley step; **Mercury needs a second Newton step in
  46 % of the steps** (post-Halley error up to 1e-8 > EPS = 1e-11). One Newton step costs ~15 % of a step.
* Forcing exactly one Newton step (timing only): 0.90.
* The interaction step is throughput-bound, not latency-bound: tree reductions of the force sums and of the 8×8
  transformations were *slower* (+0 % / +6 %, extra permutes on port 5); in-register broadcasts −0.9 % (noise).

### Iteration 1 (kept in modified form) — sixth-order series start
Series reversion of Kepler's equation in universal variables (`r0 X + η G2 + ζ G3 = dt`) to order τ^6, τ = dt/r0,
written as a bivariate polynomial in (η, v²) whose 12 coefficients depend only on r0 (available while the previous
interaction step runs), Estrin evaluation after η and v². Result: Mercury's post-Halley error ≤ 1e-11 at 5–6 d
(one Newton step always), 25 % of steps need two at 8 d (stock: 85 %).
Other things tried in this iteration: approximate 1/r for the start (unnecessary: the exact 1/r is ready early),
Newton pre-step instead of Halley (+1 %, less robust at 8 d: rejected), 9/13/15/17-term Halley series (no gain),
unconditional first-iteration bookkeeping (−1.2 %, bitwise identical), division-free Newton step with an approximate
1/f' from the Halley step (−5 %). Combined: **0.864 of the stock C kernel (1.157×), 1.24× faster than the assembly**.

### The round-off bias (Kepler grid test) — finding, and iteration 2
Tool `tools/kepler_grid.c` (and `lab/kgrid_lab.c` for lab variants): test particle, a ∈ [0.1, 1] AU, e ∈ [0, 0.9],
dt = 5 d, relative a-error after 1 step and after 1000 orbits, a evaluated from the double state in `__float128`.
Sign maps alone are ambiguous (the bias changes sign over the grid), so the decisive version runs **16 random phases
per (a, e)** (16×16 grid) and computes per point t = mean/sem of the 1000-orbit error; for an unbiased solver
⟨t²⟩ ≈ 15/13 = 1.15 (SE ≈ 0.22 over ~70 points). Restricted to dt/P < 0.07 (the Solar System regime; at larger
dt/P all kernels, the assembly included, have identical deterministic series-truncation errors, see below).

| kernel | ⟨t²⟩ e ∈ [0, 0.3) | ⟨t²⟩ e ∈ [0.3, 0.6) | median 1000-orbit error |
|---|---|---|---|
| assembly (paper v2) | 1.73 | 1.42 | 2.15e-14 |
| **C kernel as forked (items 40/42/46/48)** | 3.36 | **18.9** | 2.53e-14 |
| C, it. 1 (start + division-free Newton) | 4.29 | **26.8** | 2.87e-14 |
| C forked + compensated Newton + evaluation at final X (= assembly algorithm) | 0.77 | 1.09 | 2.10e-14 |
| it. 1 + compensated Newton + evaluation at final X, correction-form Newton | 0.99 | 1.86 | 2.35e-14 |
| it. 1 + compensated Newton + evaluation at final X, ratio-form Newton | 0.99 | 0.66 | 2.30e-14 |
| … same, plain Newton, compensated final evaluation | 0.89 | 3.16 | |
| … same, plain + constant correction in Newton, compensated final evaluation | 1.05 | 1.35 | 2.28e-14 |
| **it. 2: plain + constant correction everywhere, ratio form, evaluation at final X** | **1.21** | **1.17** | **2.22e-14** |
| (no evaluation at the final X: first-order update, any of the above) | 3.2–4.8 | 7.5–16 | |

Mechanism, established with a diagnostic that compares the updated G's against a compensated evaluation at the final
X over 2·10^6 steps: **G3 from plain Horner is systematically 5.5e-17 (relative) too small, in every orbit** — the
rounding error of the constant 1/3! (stored 0.5 ulp low); G1, G2 are unbiased at the 1e-18 level. In the Newton step
this shifts X (the effective step duration), in the f/g functions it biases g = dt − M G3. The fix costs ~nothing:
add the rounding errors of the leading constants (−E3 + zE5 − z²E7 for g3, zE4 − z²E6 for g2) before the final
rounding, by fused multiply-adds computed in parallel with the Horner chain. In addition the Stiefel functions have to
be evaluated at the final X (the first-order update of item 40 is biased even with exact G's at the old X, for reasons
not fully resolved) and the Newton update has to be the ratio form (the correction form X − f·(1/f') leaves a small
bias). **This bias is pre-existing in the fork's C kernels** (the default on AVX2 and ARM machines, e.g. the user's
i5-13500), not in the assembly. The Solar-System-level consequence is measured by the bias ensembles below.

Side finding (pre-existing, all kernels incl. the assembly): for dt/P ≳ 0.2 (z = βX² ≳ 2) the 19-term Stiefel series
is not converged; the Kepler step then violates energy conservation (e.g. a = 0.2, e = 0.5, dt = 6 d: 2.3 % of steps
off by > 1e-12; a = 0.1, e = 0.8: 62 %), identically for the assembly, the old and the new C kernels. This regime is
irrelevant for the Solar System (Mercury dt/P = 0.068, z ≈ 0.2) but matters for compact systems.

### Iteration 2 (KEPT) — final C kernel

`src/integrator_whfast512_kernel.h` (all C instantiations: AVX-512, AVX2, AVX2×4, generic, NEON; still bitwise
identical to each other — cross-kernel tests pass):
1. sixth-order bivariate series start (it. 1);
2. first Newton step in correction form X − f·RI with RI ≈ 1/f′ from the Halley step (exact division for lanes with
   |f″ΔX_h/f′| ≥ 1e-4 and for later iterations); NaN-safe convergence test (NaN in an occupied lane → bisection);
3. `w512_stiefel_Gs13_`: plain Horner + correction for the rounding of the constants 1/k! (−E3 + zE5 − z²E7,
   zE4 − z²E6, added by FMA before the final rounding);
4. Stiefel functions evaluated again at the final X (as in the assembly; replaces the first-order update of item 40,
   and with it the per-lane S-bookkeeping).

Kepler-step bias at the planets' own (a, e), dt = 6 d (64 random phases × 10^6 Kepler steps, t of the relative a-error;
`lab/kgrid_lab … ss`, figure `figures/kepler_bias_planets.png`):

| kernel | Mer | Ven | Ear | Mar | Jup | Sat | Ura | Nep |
|---|---|---|---|---|---|---|---|---|
| assembly | +1.0 | +0.5 | −0.3 | +0.4 | −1.0 | −0.1 | +0.2 | +0.2 |
| **C as forked** | **+59.1** | +0.9 | +1.1 | +0.4 | −0.6 | −1.9 | +0.2 | +0.2 |
| C final | −0.2 | +0.7 | +0.7 | +0.5 | −1.3 | −2.0 | +0.4 | −1.0 |
| C final, 3 more seeds | Σt² over 8 planets × 3 seeds = 28.6 (expected ≈ 25) | | | | | | | |

The forked C kernel drifts Mercury's semi-major axis by **+1.4e-18 per step**, i.e. Δa/a ≈ +4e-7 (+23 km) over 5 Gyr —
a change of g1 of order 1e-6 ″/yr, comparable to the WH truncation error at 6 d. The final kernel and the assembly
show no detectable bias. Kepler grid (16 phases, ⟨t²⟩ for e < 0.3 / 0.3–0.6): final 0.75–1.54 / 1.47–1.76 over 4 seeds; assembly 1.73 / 1.42;
forked C 3.4 / 18.9 (`figures/kepler_bias_tmap.png`). The final kernel is as unbiased as the assembly in this test. Standard and high-e grid maps:
`figures/kepler_grid_std.png`, `figures/kepler_grid_highe.png` (high-e: identical large-z failure region for all).

Other checks: full suite 966 passed (two round-off-sensitive WHFast512 tests adjusted with documented evidence:
one-step GR energy threshold 6e-15 → 7e-15 — over 300 perturbed ICs the new kernel exceeds 6e-15 in 0.3 % of cases,
the assembly in 1.0 %; independence test integrates 10^4 instead of 10^6 steps — an e = 0.8, dt ≈ P/2 test particle
can run off to infinity in the solver's invalid regime, identically fragile for old and new solvers, see above).
Development ensemble (16 × 10^7 steps, dt = 6 d, GR, corr17): median |ΔE/E| 1.19e-12 (assembly 1.41e-12, forked C
1.25e-12; the common offset ≈ −1.3e-12 is the deterministic truncation term, identical in all runs).
Secular phases at 6 d over 100 kyr unchanged (dg1 −8.82e-7 ″/yr, λ1 −2.90e-3 rad, as WHC).

Speed (paired ratios, this host; library path, paper setup dt = 6 d, GR, corr17, `tools/bench_w512.c`):

| kernel | time vs assembly | speedup | 5 Gyr projection on this host | projected on the paper's EPYC 9655 |
|---|---|---|---|---|
| 1 assembly | 1.000 | 1.00× | 13.6 h | 9.0 h (paper) |
| 2 C AVX-512 as forked | 0.926 | 1.08× | — | — |
| **2 C AVX-512 final** | **0.876** | **1.14×** | 12.3 h | ≈ 7.9 h (if the ratio transfers) |
| 3 C AVX2 final | 0.935 | 1.07× | 13.0 h | — |

Against the honest unbiased C baseline (the assembly's algorithm in C, `koub`): 1.24× (0.941 / 1.166, lab harness).
AVX2 instantiation: final = forked speed (1.025 vs 1.025 of forked AVX-512), now unbiased.
**Not measured on Zen 4/5**: the ratio vs the assembly may differ there (Zen 5 has 512-bit dividers with lower
latency, which shrinks the benefit of removing a division from the chain).

Solar System bias ensembles (paper Appendix C style, scaled to this machine; `tools/ens.c`, `tools/bias_summary.py`):
16 runs × 10^9 steps (dt = 6 d, GR, corr17, ICs perturbed at 1e-15), ΔE/E incl. the GR potential:

| kernel | N = 1e7 mean ± sem | N = 1e8 | **N = 1e9** | rms at 1e9 |
|---|---|---|---|---|
| assembly | −1.32e-12 ± 0.12e-12 | −1.01e-12 ± 0.41e-12 | **−1.18e-12 ± 1.2e-12** | 4.6e-12 |
| C final | −1.24e-12 ± 0.16e-12 | −1.17e-12 ± 0.40e-12 | **−0.09e-12 ± 1.3e-12** | 4.9e-12 |
| C as forked (3 runs only) | −1.33e-12 | −1.44e-12 | +1.3e-12 ± 3.9e-12 | 5.7e-12 |

(The common −1.3e-12 at 1e7 is the deterministic truncation term.) Both are consistent with zero drift at 10^9 steps
(|bias| ≲ 3e-21 per step at 2σ); the rms follows Brouwer's law (1e-16·√N = 3.2e-12). The forked kernel's Mercury bias
corresponds to ≈ +2.5e-12 in ΔE/E at 10^9 steps (Mercury holds ~1.8e-3 of the orbital energy) — at the edge of this
ensemble's sensitivity, which is why the per-planet Kepler test above is the decisive one. The paper's 32 × 3·10^10
test (≈ 5× more sensitive) was not affordable on 2 vCPUs and remains to be run on a cluster.

---

# Summary (Oct 20, 2026)

**Goal:** ≥ 2× faster Gyr Solar System integrations than WHFast512 v2 (paper: 106 ns/step on Zen 5) at research-grade
accuracy, single- or multi-core.

**Result:** 1.14× single-core (C AVX-512 kernel vs the paper's assembly, same machine, paper setup), with accuracy
equal to the assembly in every test of the contract that could be run here. **2× was not reached, and the measurements
indicate it is not reachable with the leads in the plan:**

* **S1 (larger dt):** above ≈ 6 d the error budget is dominated by step-size resonances with Mercury's eccentric orbit
  (λ1 error non-monotonic in dt, envelope ∝ dt^~10), which no splitting order removes; per unit cost nothing beats
  WH + correctors at 6 d (best alternative: equal cost, worse λ1/energy).
* **P1/P2 (parareal):** corrections must be applied in orbital elements to converge at all; then the coarse solver
  must still resolve Mercury (cost ≥ 0.4 × fine) and contraction degrades with slice length → speedup < 1.
* **S2 (Kepler solver):** the sixth-order series start and the division-free first Newton step give the measured
  1.14× (1.24× before the bias fix, which costs one more Stiefel evaluation per step).
* **S4/S5/S6:** the interaction step is throughput-bound; tree reductions and in-register broadcasts give nothing.

**Most important finding:** the fork's C kernels (default on AVX2/ARM machines) had a systematic round-off bias in the
Kepler step — Mercury's semi-major axis drifts by +1.4e-18 per step (+4e-7 over 5 Gyr), from plain Horner ignoring
the rounding of 1/3! and from the first-order update of the Stiefel functions. Fixed (constant correction + evaluation
at the final X), at no speed cost relative to the forked kernel.

**Open items / next steps for a human or the next round:**
1. Measure on Zen 4/Zen 5 (headline numbers; division latencies differ) and on the i5-13500 (AVX2 kernel 3).
2. Full bias test (32 × 3e10 steps) and instability statistics (§1.2 items 1, 4) on a cluster with `ens.c`.
3. Large-z Kepler regime (dt/P ≳ 0.2): all kernels violate energy conservation there (19-term series not converged);
   needs argument reduction (half-angle formulas for the Stiefel functions) if compact systems matter.
4. Robustness: a lane at infinity makes NaN in all lanes through 0·inf in the interaction step (and in the
   block-diagonal matrices for packed systems) — pre-existing.
5. P3 (secular coarse propagator for parareal) is the only untested parallel-in-time route; expect a research project.

Branch: pushed to github.com/hippke/rebound, branch `whfaster` (based on main @ 2cfc0c4; merges cleanly with the
later OpenMP/IAS15 commits on main, which touch none of the files changed here).

---

# Lessons, pitfalls and untested ideas (for the next round)

## Measurement methodology
* **Timing on shared cloud hosts:** single runs vary by ±10–30 % (host-level frequency/contention drifts, low steal
  time). Use many short interleaved rounds and the median of the *paired* ratio (variant/reference in the same round);
  differences of ~1 % are then resolvable. Minimum times are a usable second check; medians of unpaired runs are not.
* The two vCPUs of the cloud box are SMT siblings of one core: two simultaneous ensembles each run at ~200 ns/step
  instead of ~150. Never time anything while a background ensemble runs.
* "Calibration variants" that are deliberately wrong (force exactly 1/2/3 Newton iterations, skip the Kepler or the
  interaction step) are the quickest way to learn what a change can gain at most and where the critical path is.
* A latency model of the step was right for the Kepler step and wrong for the interaction step (throughput-bound):
  measure before restructuring.

## Accuracy methodology
* **Energy ensembles are a weak bias detector.** A per-step Kepler bias of 1.4e-18 on Mercury (the forked C kernel's
  defect) shows up as only ~2.5e-12 in ΔE/E after 10^9 steps — inside a 16-run ensemble's noise. The decisive tests:
  (1) Kepler-only, at the planets' own (a, e) and dt, 64 random phases × 10^6 steps, t-statistic of the relative
  a-error per planet (`lab/kgrid_lab <variant> ss`, ~4 s per kernel); (2) the (a, e) grid with 16 phases per point
  and ⟨t²⟩ ≈ 1.15 for an unbiased solver (SE ≈ 0.22 per bin). Use several seeds before believing a 2–3σ result:
  one variant's single-seed 1.86 turned out to be noise (1.47–1.76 over 3 more seeds).
* Sign maps after 1000 orbits (paper Fig. 2/3) are ambiguous: the bias changes sign across (a, e), so "fraction
  positive" and "mean/sem over the grid" disagree. Per-point t over several phases is unambiguous.
* For dt/P ≳ 0.07 (high e) or ≳ 0.2, every kernel (assembly included) has deterministic series-truncation errors;
  exclude that region when testing round-off bias, or the statistic is dominated by it.
* Measure the "reference noise level" of a chaotic integration by comparing two valid fine integrators (asm vs C)
  from identical ICs: 1.3e-9 AU after 1 kyr, 7e-8 AU after 10 kyr, 7e-6 AU after 100 kyr. Any parallel-in-time
  scheme only has to converge to that level.
* For dt comparisons, use exact step counts (`reb_simulation_steps`) with a sample interval divisible by all dt
  (360 000 d), so different runs are compared at identical times.

## Numerical facts learned
* Plain Horner with rounded constants 1/k! is systematically biased: the error of 1/3! (stored 0.5 ulp low) makes G3
  5.5e-17 too small *in every step*. Fixed cheaply by adding the constants' rounding errors (W512_IF_ERR) before the
  final rounding. Sign convention: W512_IF_ERR(k) = rounded − true (the compensated Horner *subtracts* it).
* Evaluating the Stiefel functions at the final X is necessary; the first-order update from the previous iterate is
  biased (residual Mercury t ≈ +4–7, Saturn ≈ −3) even with compensated or constant-corrected series and the ratio
  form. **Mechanism not understood** — worth a closer look, because it would save ~8 % (one series evaluation).
* After the constant correction, compensated Horner (paper Appendix C) is no longer needed for unbiasedness in these
  tests (plain + correction + evaluation at the final X: ⟨t²⟩ ≈ 1.2, planets clean).
* The Kepler solver's 19-term Stiefel series is not converged for z = βX² ≳ 2 (dt/P ≳ 0.2); the bisection fallback
  then returns wrong results too (the same for the assembly). Long-standing, undocumented.
* WHFast512's universal-variable start: a sixth-order series reversion (coefficients in the log, derived with sympy)
  makes one Halley step + one Newton step sufficient at dt ≤ 6 d; the start's latency matters only after η and v²
  (which need the new velocities), so position-only coefficients should be precomputed.
* The lazy-implementer kernel (WHCKL) is not exactly symplectic (kick ∇B(q + c∇B) is not a gradient at O(c²ε³));
  irrelevant here since S1 was rejected, but relevant if someone revisits it for Gyr runs.

## Practical pitfalls in this code base
* `reb_simulation_steps()` synchronizes at the end; WHFast512's default `concatenate_steps` is 10^6, so
  `sim.steps(1)` integrates 10^6 steps (tests rely on this).
* The Python package loads `librebound.cpython-*.so` built by `pip install -e .`, not `src/librebound.so`; rebuild
  with pip (≈ 70 s) or `benchmarks/perf/check.py` after any C change, or the tests run stale code.
* For ≤ 4 planets the C kernels dispatch to the four-lane instantiations (`*x4`); debug with those.
* Empty lanes are NaN inside the Kepler step; any new mask/convergence logic must restrict to `c->k1`, otherwise
  every step with < 8 planets falls into the bisection fallback.
* `integrator_saba.c`: synchronize with `keep_unsynchronized=1` before the first step crashed (fixed in this branch).
* `pkill -f <pattern>` kills the invoking shell when the pattern appears in its own command line.

## Untested ideas (ranked by my estimate)
1. Find why the first-order Stiefel update is biased; if fixable, the final evaluation (~8 %) can go again.
2. Shorter final evaluation: reuse the high-order tail of the Newton step's series and re-evaluate only the leading
   terms at the final X (needs a first-order correction of the tail; accuracy at the 1e-18 level must be shown).
3. A staggered/partitioned splitting that lets the next Kepler step start before the interaction step finishes
   (Kepler 0.5 and interaction 0.35 of the step would overlap); changes the method, needs a full S1-style analysis.
4. Argument reduction (half/quarter-angle identities) for the Stiefel functions for z ≳ 2 — correctness for compact
   systems, not speed.
5. P3: parareal with a secular coarse model in elements (only untested multi-core route; research project).
6. Measure on Zen 4/5: the division-free Newton step gains less where 512-bit division latency is short.

---

# Round 2 (Oct 21, 2026): parallel-in-time with a mean-element coarse propagator

Literature review of this round: `LITERATURE_PINT.md` (5 parallel searches; no measured single-trajectory multicore
speedup for a planetary integrator exists; the leading untested idea is parareal in mean elements). Process as in
`../PERF_LOG.md` (one work item per iteration, measure, keep/revert, log). Prototype code: `parareal/` (Python driving
the C kernel through the library; F = WHFast512 C kernel 2, dt = 6 d, GR, corrector 17, as in round 1).
Machine: same 2-vCPU cloud Xeon (one physical core; **multi-core wall-clock cannot be measured here** — speedups are
computed from measured iteration counts k and cost ratios r with the pipelined parareal model
S = N / (N r + k (1 + r)), N = slices = cores).
Convergence criterion: every boundary within the fine solver's own round-off noise (asm vs C kernel from identical ICs:
Mercury λ ≈ 1e-7 rad, ln a ≈ 1e-12 over 10 × 1 kyr; positions 1.3e-9 → 2e-8 AU).
History: `pint_history.csv`.

## Iteration 1 — slice-Jacobian test (`parareal/jac.py`, `oracle.py`)
Central-difference Jacobians (48 × 48, Jacobi elements ln a, λ, k, h, q, p) of one 986-yr slice of F and of coarse WH
(dt_G = 16, 24 d) at the 10 slice starts; linearised parareal recursion e^{k+1}_{n+1} = G'e^{k+1}_n + (F'−G')e^k_n
started from the true coarse errors.
* The linear model **reproduces the measured nonlinear parareal iteration by iteration** (dt_G = 16 d: Mercury λ
  4e-2 → 2e-2 → 9e-4 → 8e-5 → 1e-6 → 8e-8 (actual) vs 4e-2 → 2e-2 → 4e-4 → 6e-5 → 6e-8 (linear); 24 d: no contraction
  in both). So convergence is a property of F'−G' along the trajectory, and can be predicted cheaply.
* The mismatch is in the **ln a rows of Mercury–Mars**: |∂ln a_end/∂(·)| mismatch 1e-2…5e-2 (16 d), 0.8 (24 d). F' itself
  has ∂ln a_end/∂λ_j ≈ 1e-6…1e-5 (short-period terms) and ∂ln a_end/∂ln a_j up to 0.09 (short-period × shear, e.g.
  Earth ← Venus). Coarse WH gets these short-period derivatives wrong (aliasing at large dt).
* Oracle (coarse = F' + noise): k ≤ 3 to the noise level needs **ln a rows accurate to ≲ 1e-8 absolute and all other
  entries to ≲ 1e-2 relative**; with ln a rows at 1e-6 it takes 5–6 iterations; replacing only the ln a rows of the WH
  coarse Jacobian by the exact ones is not enough (λ rows are wrong too).
* Conclusion: no osculating-element coarse model can meet 1e-8 in the ln a rows (that is 1e-3 relative to the
  short-period terms themselves *and* requires correct end phases). In mean elements these rows are ≈ identity.

## Iteration 2 — numerical first-order mean-element theory (`parareal/meanel.py`, `test_meanel.py`)
For each planet pair (i<j), the first-order Jacobi perturbing accelerations (interior-on-exterior with the dipole
indirect term, exterior-on-interior heliocentric-like; checked against a direct N-body sum: 3–5e-5 relative for the
inner planets, 1–2e-3 for the giants = second-order terms) are evaluated on a 64 × 64 grid of (λ_i, λ_j) with frozen
slow elements; Gauss' equations (∂E/∂v by central differences) give the element rates; their 2-D FFT gives the
short-period displacements f_k/(i k·ν) plus the shear term −1.5 n·(ln a)_k/(i k·ν)² for λ. GR (the kernel's
−3 m0²/(c² r²) potential) as a 1-D single-planet term. Mean = osculating − Σ displacements. Set-up 0.5 s in numpy.
Residual short-period scatter of the mean elements along a 10-yr fine run (std after removing a quadratic trend):

| | Mer | Ven | Ear | Mar | Jup | Sat | Ura | Nep |
|---|---|---|---|---|---|---|---|---|
| osc ln a | 2.0e-6 | 6.3e-6 | 1.0e-5 | 3.3e-5 | 2.3e-5 | 5.5e-5 | 1.9e-5 | 2.2e-6 |
| mean ln a | **8.6e-10** | 7.2e-9 | 9.8e-9 | 3.0e-8 | 1.2e-7 | 3.2e-7 | 6.5e-8 | 1.3e-8 |
| osc λ | 2.8e-5 | 2.6e-5 | 2.4e-5 | 8.0e-5 | 3.9e-5 | 8.8e-5 | 2.7e-5 | 4.4e-6 |
| mean λ | 2.0e-8 | 6.4e-8 | 5.2e-8 | 9.1e-8 | 3.0e-6 | 7.5e-6 | 1.6e-7 | 1.7e-7 |

Inner planets: 300–2300× reduction, at the shear-scaled requirement (≈ 1e-8 × shear_Mer/shear_i). Giants: limited by
the near-resonant great inequality (2n_J − 5n_S, 880 yr) and Uranus–Neptune (≈ 2:1, 4200 yr), which a first-order
theory with small divisors cannot represent.

## Iteration 3 — mean-element parareal (`parareal/mpr.py`, `mgc.py`, `jac2.py`)
G = osculating → mean (start tables) → Kepler drift at n(ā) + first-order secular rates (k = 0 Fourier terms) →
osculating (end tables); tables cached per slice boundary. Correction in osculating elements as in round 1.
10 × 986 yr slices, max error over slices (ln a | λ) for Mercury … Neptune:

* **All planets in mean elements:** inner planets fine (k = 0: Mercury λ 2e-4 vs 4e-2 for WH 16 d), giants don't
  contract (λ ≈ 1e-2, Saturn ln a 2e-4): near-resonances, as predicted by iteration 2.
* **Hybrid** — giants by WHFast512 at dt = 180 d in a Sun(+inner mass) + 4-planet subsystem, inner planets in mean
  elements: giants converge in 1–2 iterations; inner planets: Mercury λ 2e-4, 2e-4, 2e-4, 4e-5, 4e-7, 2e-7, 8e-8;
  Mars slowest; **k = 6 to the noise level** (max position error 3e-8 AU). Linear prediction again exact.
* Dominant remaining mismatch (`jac2.py`, terms (F'−G')_rc e_c, ln a rows × shear): derivatives w.r.t. the
  **eccentricity vectors** (Mars λ ← Mars k, Mars ln a ← Mars k, Venus ln a ← Earth k/h): the frozen tables have no
  slow-element dependence. Next: tables rebuilt from the actual state (“live”), then a cheap equivalent.
