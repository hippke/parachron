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

## Iteration 4 — live tables (`parareal/prun.py --live 1`)
Tables (transform and secular rates) rebuilt from the actual input state at every G call, so G carries the
slow-element (e, i vector) derivatives. 10 × 986 yr: k = 1…4 max position errors **2.3e-4, 6.7e-6, 5.5e-8,
2.7e-8 AU** (noise 3.4e-8 AU): converged at **k = 3** (criterion from here on: max position error over all boundaries
≤ 3 × the asm-vs-C noise of the fine solver over the same span). Ablation: live secular rates only (frozen transform)
k = 4; frozen everything k = 6. Contraction ≈ 1e-2 per iteration.
10 × 9856 yr slices: Euler secular step unstable (Mars e grows 0.09 → 0.28 in 4 slices: 0.8 rad of g4 rotation per
step). With RK4 (2 substeps, averaged rates on 32² grids, giant elements sampled from their WH run): k = 5–6
(3.9e-7 AU over 100 kyr). One-slice diagnostic (`jac1.py`): Mars' **mean** a has a genuine slow phase dependence
(∂ ln ā/∂λ ≈ 3e-4 over 10 kyr; second-order/three-body near-resonances, not in a first-order theory); its effect grows
∝ T³ (mismatch × shear × e0), so slices of ~1 kyr are the natural size.

## Iteration 5 — the cheapest coarse model that keeps k = 3 (10 × 986 yr)
| variant | k |
|---|---|
| Ng = 48, all six element components transformed | 3 |
| Ng = 32 | 3 |
| Ng = 24 | 4 |
| Ng = 32, transform of ln a and λ only | 3 |
| **Ng = 32, transform of ln a only** (k, h, q, p, λ mean = osculating) | **3** |
| + pairs with Uranus/Neptune dropped | 3 |
| + pairs with Saturn dropped too | 4 |
| + inner–Jupiter/Saturn pairs on 16² grids (inner–inner 32²) | **3** (config A) |
| inner–inner 28² | 4 |
| α/e-based adaptive grids (8–32) | diverges at k = 2 (aliasing; not pursued) |

**Config A** (kept): ln a short-period transform only (d ln a/dt = 2a/μ · v·a_pert — one scalar per grid point),
6 inner–inner pairs on 32², 8 inner–{Jupiter, Saturn} pairs on 16² (8192 grid points per table), first-order secular
rates of λ, k, h, q, p from the same grid (mean over the grid, no FFT), Euler drift, giants by WH (180 or 360 d),
tables rebuilt per call.

## Iteration 6 — robustness of config A
| case | k |
|---|---|
| 10 × 986 yr from t = 0 | 3 |
| giants' WH step 360 d instead of 180 d | 3 |
| 10 × 986 yr starting at t = 1 Myr | 3 |
| 20 × 986 yr | 4 (k = 3: 5.9e-7 AU vs 3 × 1.4e-7) |

Model speedups (S = N/(N r + k(1 + r))): N = 10, k = 3: 2.7 at r = 0.05, 2.3 at r = 0.1; N = 20, k = 4: 3.8 at
r = 0.05; N = 5, k = 3: 1.5 at r = 0.05. The coarse cost r in C is the open number → iteration 7.

## Iteration 7 — coarse propagator and parareal driver in C (`pint/pint_coarse.c`, `pint_test.c`, `pint_run.c`)
`pint_coarse()` implements config A (giants' WH step 360 d, C kernel 2). Agreement with the Python prototype from the
same input: giants bitwise, inner ln a to 1e-15, λ to 1e-11 (1000 yr × n), k/h/q/p to 1e-10 (forward-difference
secular rates). Optimisation (one core, this host; per G call = 2 table builds + giants):

| version | G [ms] | tables | giants WH | fine slice [ms] | r |
|---|---|---|---|---|---|
| straight port (FD dE/dv for 6 components, naive FFT, complex division) | 2.39 | 2.12 | 0.28 | 13.2 | 0.18 |
| analytic ∂ln a/∂v, FD only for the secular components, converged Kepler, twiddle FFT, two real grids per complex FFT, division-free half evaluation | **0.85** | 0.61 | 0.24 | 11.1 | **0.071–0.076** |

Profile of the optimised table pair: grid r, v, ∂E/∂v 0.15 ms, pair loops 0.18, FFTs 0.15, series evaluation 0.12.

`pint_run 10 360000 2 5 360` (10 × 986 yr, 2 threads): k = 0…3 max position error 6.6e-3, 1.7e-4, 7.0e-6,
**2.7e-8 AU** (noise 1.35e-8 AU) → **k = 3**, same as Python. Measured: fine slice 12.7 ms, coarse step 0.93 ms,
r = 0.073. Pipelined model on 10 cores: **S = 2.53 (25 % efficiency)**; non-pipelined schedule (as implemented):
1.69. Measured wall on this 2-vCPU host: 0.60 × serial (expected for P < N; the host cannot show the 10-core number).
1.5-kyr slices (Python): k = 4 → keep ~1 kyr.

## Iteration 8 — what limits k (C coarse through ctypes: `parareal/ptc.py`; `oracle2.py`; `pint_run`)
* Jacobian diagnostic of config A (13 s with the C coarse step): linear prediction = measured iterations; dominant
  remaining mismatches ∂(ln a, λ)/∂λ and ∂(ln a, λ)/∂(k, h) of Mars and Venus.
* Oracle (coarse Jacobian blocks replaced by the exact ones), max a·|Δλ| at k = 1, 2, 3 [AU]:
  config A 2.3e-4, 6.5e-6, 3.0e-8; exact (ln a, λ) rows ← all columns 7.9e-5, 3.1e-7, 2.8e-9; exact inner-planet rows
  (all six elements) 1.3e-5, 3.4e-8, 4e-11. **k = 2 would need an essentially exact inner-planet Jacobian** (plus
  better e0); not reachable with cheap first-order improvements.
* Shorter slices do not help (the noise target shrinks with the window): 10 × 493 yr: k = 3, r = 0.154 (S10 = 2.0);
  10 × 246 yr: k = 3, r = 0.26 (S10 = 1.6). 1 kyr stays optimal.
* Scaling with the number of slices (= cores), 986-yr slices, C driver:

| N | k | r | modelled S (pipelined) | efficiency |
|---|---|---|---|---|
| 5 | 3 | 0.087 | 1.35 | 27 % |
| 10 | 3 | 0.073 | 2.53 | 25 % |
| 20 | 4 | 0.080 | 3.38 | 17 % |
| 40 | 4 | 0.078 | 5.38 | 13 % |

## Iteration 9 — long runs and the energy contract (`pint/pint_long.c`, `energy_ana.py`)
`pint_long`: consecutive windows of N slices, fixed k, optional serial runs (kernel 2 or 1) and boundary-jump
diagnostics. 1 Myr (100 windows × 10 × 986 yr, k = 3): all windows k = 3, last update 1e-6 – 3e-5 AU; parareal
energy history statistically equal to the serial one.
Boundary jumps U_n − F(U_{n−1}) (30 windows): **energy rms 5.8e-14, mean −2.3e-15 ± 3.3e-15**; the asm-vs-C
difference over one slice (pure round-off): rms 4.9e-14. Same at k = 4 and 5 (5.3e-14, 5.5e-14): the energy jumps
are round-off limited — the boundary state descends from F(U^{k−1}), i.e. another round-off realisation. Position
jumps: median 2.3e-9 AU (round-off 5e-10), growing along the window (max 3e-7 AU at slice 9).
Energy random walk (300 kyr, per-slice energy increments compared with the serial C run): serial asm 5.44e-14 rms,
parareal 6.25e-14, **parareal + energy projection 5.76e-14** — i.e. the parareal trajectory's own random-walk
amplitude is 1.28× (1.11× with projection) that of a serial run; no bias detectable (means ≤ 2e-15 ± 3e-15).
Energy projection (kept, default on): after each correction, scale barycentric velocities so that
E(U^k_{n+1}) = E(U^k_n) + [E(F(U^{k−1}_n)) − E(U^{k−1}_n)] (each slice's energy increment = one fine round-off
realisation, as in a serial run; leaves the parareal fixed point unchanged; one energy evaluation per boundary).

## Iteration 10 — pipelined driver (`pint/pint_pipe.c`)
One thread per slice; thread n waits for U^k_n (semaphore), computes G, forms and hands on U^k_{n+1} (with energy
projection), then runs F(U^k_n) for the next iteration; slices n < k are passed through. Schedule = the pipelined
model N T_G + K (T_F + T_G) per window; windows chained. Checked: 10 × 986 yr, K = 3 → final state within 1.7e-8 AU
of the serial run (3 × noise = 4e-8). Measures wall-clock and core-hours against a serial run of the same span —
the program to run on a ≥ 10-core machine (`pint_pipe 10 360000 3 <windows>`).

## Iteration 11 — previous-window calibration of G (rejected)
G_cal = G + mean defect (F − G) of the previous window (all elements, or λ only), N = 20: k = 4 in all variants
(errors at k = 3: 7.7e-7 uncalibrated, 6.9e-7 λ-calibrated, 9.2e-7 fully calibrated). The defects vary from slice to
slice (phase-dependent); a constant correction does not reduce the iteration error.

## Iteration 12 — 10 Myr validation: secular frequencies and drifts; Heun secular step (KEPT)
10 Myr = 1000 windows × 10 × 986 yr (energy projection on), compared with serial C (kernel 2), serial asm
(kernel 1, round-off reference) and serial C with all planets' x perturbed by 1e-12 and 1e-9 AU (chaos
reference). NAFF frequency analysis of k + ih, q + ip (`pint/freq_ana.py`), differences vs serial C [″/yr]:

| | asm | IC +1e-9 AU | **Euler K=3** | Euler K=4 | Euler K=5 | **Heun K=3** |
|---|---|---|---|---|---|---|
| g1 | −4e-10 | +1.1e-8 | **+1.7e-6** | −1.5e-7 | −1.6e-8 | **+3.3e-8** |
| g2 | −3e-11 | −1.5e-8 | +1.1e-6 | −2.9e-8 | −9.7e-9 | +1.1e-8 |
| g4 | −1e-10 | −8.3e-8 | **+1.05e-5** | +5.2e-8 | −1.0e-7 | −1.9e-7 |
| g6 | +6e-9 | −1.4e-7 | +4e-9 | +5e-9 | −6e-9 | −2.4e-8 |
| s1 | +6e-10 | −5.2e-8 | +9.5e-7 | +2.0e-7 | −7.1e-8 | −1.2e-7 |
| s2 = s3 | +3e-10 | +9.8e-8 | −5.0e-6 | −8.4e-8 | +6.4e-8 | +2.6e-7 |
| s4 | +1e-10 | +4.5e-8 | −3.6e-6 | −7.5e-8 | +4.5e-8 | +1.2e-7 |
| giants | ≤ 6e-9 | ≤ 1.4e-7 | ≤ 5e-9 | ≤ 5e-9 | ≤ 6e-9 | ≤ 2.4e-8 |

* **Euler K = 3 fails**: the inner planets' e and i vectors drift away from the serial run steadily (Mars |Δz| 4e-5,
  |Δζ| 1e-5 at 10 Myr, vs 3e-6 / 1e-7 for the 1e-9 AU IC perturbation), frequency shifts up to 1e-5 ″/yr, i.e.
  larger than the WH truncation error at 6 d (dg1 = −8.8e-7 ″/yr, round 1). Mechanism: G's secular step error is
  systematic (Euler over 0.03–0.09 rad of rotation per slice), so the residual after K iterations biases every
  boundary in the same direction. More iterations shrink it (K = 4: ≤ 2e-7; K = 5: ≤ 1e-7).
* **Heun (trapezoidal) secular step** (rates also at the predicted end, from the end table that is built anyway):
  K = 3 deviations ≤ 2.6e-7 ″/yr, at the level of the 1e-9 AU chaos reference and ≤ 4 % (g1) of the WH truncation
  error; e/i vector differences at 10 Myr chaos-like (Mars |Δz| 2e-6, |Δζ| 4e-7). Kept (config B).
* Energy (all runs): (E − E_serialC)/E at 10 Myr: asm −8.3e-12, IC 1e-9 +4.1e-12, Euler K3 −1.5e-12, Heun K3 +8.5e-12,
  K4 −8.5e-12, K5 +1.8e-12 — parareal within the spread of two serial round-off realisations; no drift beyond it.
* Convergence over 1000 windows (Heun K = 3): last update median 1.5e-6 AU, 89 windows > 1e-5 AU, max 2.5e-4 AU
  (clustered episodes, e.g. 3.3–3.6 Myr); with Euler K = 4: median 3e-8, max 1.1e-5.
* Cost (idle host, `pint_run`): Euler r = 0.077 (S10 = 2.50), **Heun r = 0.083 (S10 = 2.45)**.

**Config B = config A + Heun secular step + energy projection, K = 3 for 10 slices of 986 yr**: accuracy contract
items 1–3 met in these tests (energy unbiased, Brouwer-like; secular frequencies at chaos-noise level);
modelled wall-clock speedup 2.45 on 10 cores (25 % efficiency), 3.3 on 20 cores (K = 4).

## Iteration 13 — attempts to reach K = 2 (N = 10) or K = 3 (N = 20) (all rejected)
* Thread safety (bug found with `pint_pipe`): the coarse step used static scratch arrays; concurrent calls from the
  pipelined threads raced (REBOUND "e < 1" errors). Fixed: thread-local scratch, one-time FFT twiddles. The pipelined
  timing model checked on this host with N = 2, K = 2 (the only size 2 vCPUs allow): measured 0.77–0.86 × serial vs
  0.86 predicted from single-thread T_F, T_G (both vCPUs behave as independent cores for this load).
* **Secant-enhanced parareal** (`parareal/secant.py`; quasi-Newton/Krylov-parareal idea: correct (F' − G') along the
  stored iterate differences of the same slice, no extra F/G): rank 1 with damping 0.5: k = 2 error 1.2e-6 AU (plain
  6.9e-6), unweighted 2.4e-6; still K = 3 (needs 4e-8). With neighbour slices or full weight: worse. Not kept.
* Knob sweep for N = 20 (`parareal/sweepN.py`, error/noise at k = 3; plain config B = 1.7): 64² inner grids 15 (!),
  32² inner–giant grids 1.7, Uranus/Neptune pairs 1.5, giants at 180 d 1.6, Euler 3.8. Small-divisor cutoff
  (drop terms with period > P): P = 30 yr 320 (terms with 30–100 yr periods matter), 100 yr 2.1, 300 yr 1.7.
  The 64² failure is not a small-divisor effect (cutoff does not cure it); 32² is the sweet spot.
* The "≤ 3 × noise" criterion depends on the noise realisation: N = 20, k = 3 gives 1.7 × (noise from boundary-
  chained references) or 3.5 × (continuous serial references) → N = 20 is borderline between K = 3 and 4.

## Iteration 14 — second 10 Myr validation (different ICs: all planets' x + 1e-6 AU): K = 3 is not enough
Same protocol as iteration 12 (serial C reference; serial asm). Frequency differences vs serial C [″/yr]:

| | asm | Heun K=3 | Heun K=3 + defect calibration | Heun K=4 | Heun K=5 |
|---|---|---|---|---|---|
| g1 | −6e-10 | +2.5e-7 | +3.0e-7 | −1.6e-7 | +2.9e-8 |
| g4 | +3e-10 | **−2.2e-6** | −2.0e-6 | +3.1e-8 | −2.0e-9 |
| s1 | +6e-10 | +7.9e-7 | +4.9e-7 | +3.2e-7 | −1.2e-7 |
| s2 = s3 | −4e-10 | **+3.3e-6** | +3.0e-6 | −8.0e-9 | −1.5e-8 |
| s4 | −2e-10 | +1.7e-6 | +1.6e-6 | −4.3e-8 | +2.1e-8 |
| giants | ≤ 1.5e-8 | ≤ 1.9e-8 | ≤ 2.3e-8 | ≤ 1.8e-8 | ≤ 2.8e-8 |

* With these ICs Heun K = 3 shifts the inner planets' frequencies by up to 3.3e-6 ″/yr — 10× the WH truncation error
  at 6 d (round 1: dg1 −8.8e-7, ds1 +3.4e-7). Iteration 12's K = 3 result (≤ 2.6e-7) was a favourable realisation.
* Defect calibration (G + previous window's mean inner k, h, q, p defect) does not help → not a window-mean offset.
* Boundary-jump statistics in z = k + ih and ζ = q + ip (40 windows, components along/across the vector): no
  significant mean at K = 3 or 4 (e.g. Mercury ζ tangential +9e-14 ± 7e-14 per boundary at K = 4, which would move s1
  by only ~3e-10 ″/yr); rms z jumps 1e-11 (K = 3) → 2e-12 (K = 4) for Mercury, ζ jumps ~1e-11–1e-10 at both K
  (vs round-off 1e-15–6e-14).
* Chaos/IC-sensitivity reference (serial, all planets' x perturbed by δ): the 10 Myr frequencies depend linearly on
  δ: s2 shifts 1e-7, 9e-7, 9.6e-6, 9.4e-5 ″/yr for δ = 1e-9 … 1e-6 AU (all planets shift, giants included). The
  parareal pattern differs (inner planets only, giants unaffected).
* The deviations shrink with K and change sign between K = 4 and 5 (Mercury s1 +3.2e-7 → −1.2e-7): a
  convergence-level effect, not a fixed bias.

**Conclusion of round 2 on accuracy:** secular-frequency fidelity over 10 Myr (worst inner planet) ≈ 3e-6 ″/yr at
K = 3, ≈ 3e-7 at K = 4 (≈ WH truncation error at 6 d), ≈ 1e-7 at K = 5; energy unbiased and Brouwer-like for all K.

## Round 2 result (speedup at stated accuracy)
*Superseded by iteration 15 (adaptive stopping rule: 2.3× on 10 cores at K = 5-level fidelity); current summary in `PINT_SUMMARY.md`.*
Measured on one core of this host: fine slice T_F = 11.5–12.5 ms (986 yr, 6 d steps), coarse step T_G = 1.0 ms
(r = 0.083); iteration counts measured (10 × 986-yr windows); pipelined schedule validated at N = 2 (0.77–0.86
measured vs 0.86 predicted). Wall-clock speedup S = N/(N r + K(1 + r)), core-hours = N/S × serial:

| target accuracy | K | 10 cores | 20 cores (K+1) |
|---|---|---|---|
| frequencies ≲ 3e-6 ″/yr, energy unbiased | 3 | **2.45×** (4.1× core-h) | 3.3× |
| frequencies ≲ 3e-7 ″/yr (≈ WH 6 d error) | 4 | **1.94×** (5.2× core-h) | 2.8× |
| frequencies ≲ 1e-7 ″/yr | 5 | 1.61× (6.2× core-h) | 2.45× |

The user's example target ("2× with 5× the cores") is not reached at research-grade accuracy; "≈ 2× with 10× the
cores" is (K = 4). Not yet measured on a ≥ 10-core machine (this host has 2 vCPUs): `pint/README.md`.

## Iteration 15 — bug check, parameter search, adaptive stopping rule (KEPT)
**Bugs.** Fixed earlier and verified in the code: thread-unsafe static scratch in `pint_coarse` (thread-local now,
one-time FFT twiddles), Euler secular bias (Heun default), garbled calibration update in `pint_long`, slicing bug in
`secant.py`. Found and fixed now: a giants' step that does not divide the slice was silently rounded (`lround`), so G
integrated the giants over the wrong time span (the "540 d diverges" result of iteration 13 was this bug; correctly
done, 540 d converges but slowly: 5.6e-5 AU at k = 4). The step is now adjusted to divide the slice exactly; the
drivers reject slice lengths that are not multiples of 6 d.

**Parameter search** (`parareal/sweep2.py`): error/noise after k iterations for 10 × slice windows started at 10 epochs
of the 10 Myr serial run (0 … 9.3 Myr, including the slow windows at 1.03 and 3.4 Myr); measured r; modelled
10-core speedup at the mean (adaptive) and maximum (fixed) K needed to reach ≤ 2 × noise:

| configuration | r | K per epoch | S(mean K) | S(max K) |
|---|---|---|---|---|
| default (986 yr, giants 360 d, 32²/16², maxj 5, Heun) | 0.093 | 3,5,3,4,4,3,3,3,4,3 | 2.10 | 1.56 |
| 821-yr slices | 0.096 | mean 3.6 | 2.04 | 1.55 |
| 1183-yr slices | 0.080 | mean 4.2 | 1.87 | 1.20 |
| 1380-yr slices | 0.074 | mean 3.6 | 2.17 | 1.63 |
| 1183 yr, giants 300 d | 0.082 | mean 4.5 | 1.76 | 1.19 |
| 1183 yr, 64² inner grids | 0.201 | mean 3.6 | 1.58 | 1.25 |

One-at-a-time from the default (5 epochs, ≤ 1 × noise): inner–giant grids 8² worse (mean K 5.8), 32² no gain, inner
grids 16² much worse (6.8), 64² faster convergence (k = 3: 43 vs 870 × noise) but r = 0.23; Uranus/Neptune pairs,
giants' step 240/300/450 d, slice lengths 657–1380 yr: all within the epoch-to-epoch scatter (±5–10 % in S).
→ **The default is on a flat optimum; no parameter gives a significant gain.** The spread in K between windows (3 for
most, 5 for a few) is the dominant cost, i.e. the iteration count should be chosen per window.

**Stopping rule** (`parareal/stoprule.py`): the largest boundary update (known during the run) predicts the true
error: stop at the first k ≥ 3 with update ≤ τ: τ = 1e-5 AU → mean K 3.3 but one window at 14 × noise; τ = 3e-6 →
mean 3.7, worst 2.5 × noise; τ = 1e-6 → mean 4.0, worst 1.6 ×.
10 Myr validation on the hard ICs of iteration 14 (`pint_long … xtol`):

| run | mean K (histogram K=3/4/5) | max frequency shift [″/yr] | energy vs serial at 10 Myr |
|---|---|---|---|
| fixed K = 3 | 3 | 3.3e-6 | −3.0e-12 |
| fixed K = 4 | 4 | 3.2e-7 | +2.5e-13 |
| fixed K = 5 | 5 | 1.2e-7 | +8.0e-12 |
| **adaptive τ = 3e-6 AU** | **3.29** (710/285/5) | **8.6e-8** | +2.2e-12 |
| adaptive τ = 1e-6 AU | 3.64 (370/616/14) | 2.5e-7 | −1.3e-12 |

(asm vs C: frequencies ≤ 1.5e-8, energy +5.3e-12.) Shifts below ~2e-7 are at the scatter level of this measurement;
both rules give fixed-K = 4/5 fidelity at a lower mean cost. Kept: **τ = 3e-6 AU, K_min = 3, K_max = 6**.
Implemented in `pint_pipe` (args 8–9): the next iteration's fine slices run speculatively in 10 chunks and are
abandoned as soon as the last slice decides to stop; checked with N = 4 and 10, τ = 3e-6 and an unreachable τ (stops
at K_max, no deadlock).

Modelled wall-clock speedup at research-grade accuracy (frequency shifts ≲ 1e-7–3e-7 ″/yr, energy unbiased),
S = N/(N r + K̄(1 + r)), r = 0.083: **10 cores 2.3× (K̄ = 3.3)**, vs 1.94× (fixed K = 4) and 1.6× (fixed K = 5).

## Iteration 16 — diagnosis for a better mean-element theory (`parareal/resid_sample.py`, `resid_spec.py`, `res3b.py`)
Fine run of 10 kyr sampled every 6 d; first-order mean ln a of the inner planets (osculating − C transform, live
tables) analysed for residual frequencies, matched to combinations of mean motions.
* **Short-period residual** (rms over 300 yr: Mercury 1.7e-9, Venus 3.6e-8, Earth 2.6e-8, Mars 7e-8): Venus and
  Earth show a comb at (5+j) n_V − j n_E (P ≈ 0.09–0.12 yr), i.e. high harmonics of the Venus–Earth pair (α = 0.72)
  aliased by the 32² grid. New option `ng_hi`/`alpha_hi` (finer grid only for inner pairs with α ≥ 0.6, i.e. V–E and
  E–M): residual Venus 3.6e-8 → 5.2e-9, Earth 2.6e-8 → 9.2e-9 (Mars unchanged, see below). Parareal (10 epochs):
  worst error/noise at k = 3 870 → 42, k = 4 14 → 2.6, but r 0.093 → 0.154 and the per-window K at the 2 × noise
  target unchanged (floor-limited) → S(mean K) 2.11 → 1.76. Not kept (option left in, default off).
* **Long-period residual (the dominant one for Mars):** an oscillation of Mars' mean (and osculating) ln a with
  P = 1836 yr, amplitude 2.2e-7 (Earth 1.4e-8, giants none). The only low-order match is the **Earth–Mars–Jupiter
  three-body near-resonance 4λ_E − 8λ_M + 3λ_J (+ ϖ, first order in e)** — second order in the masses, frequency
  5.4e-4/yr = a small difference of large numbers. Over a 986-yr slice it changes Mars' mean ln a by up to 7e-7
  (→ λ error ~2e-3 rad, the bulk of Mars' coarse error) and its phase derivative (~8 × 2e-7) is the "genuine Mars
  mean-a phase dependence" of iteration 4. A first-order pairwise theory cannot contain it.
* Test with the term added to G (amplitude/phase fitted from the 30 preceding boundaries, frequency from the mean
  motions): mixed, no systematic gain (e.g. epoch 2000: k = 3 2.5 → 2.7; 4500: 10 → 5; 9300: 0.7 → 0.5). The
  frequency is extremely sensitive (osculating mean motions give a 100 % error; first-order mean ones ~5 %), and the
  fitted amplitudes scatter 7e-9 … 1.6e-7 between epochs — an empirical fit from 1-kyr samples is too crude.

**What the diagnosis implies.** Typical windows reach the noise floor at k = 3; at k = 2 they are 1e2–1e4 × noise.
A better theory therefore helps only (i) by ~10 %, by making the few hard windows converge like the others, or
(ii) substantially, if it makes k = 2 converge — which needs a 100–10⁴ × smaller k = 2 error, i.e. an inner-planet
coarse map exact to second order *including* the Earth–Mars–Jupiter near-resonant term as resonant dynamics (its
period ≈ 2 slices, so it cannot be averaged away) and with the frequency taken from second-order mean motions.
That is a full second-order (Lie-series) theory with three-planet terms — a multi-week project with uncertain
payoff; ranked as the main open research item in the plan.

---
# Round 3 — second-order mean-element theory (plan: `RESEARCH_PLAN_R3.md`)
Test bed for all round-3 iterations: `parareal/ph1eval.py` — nonlinear parareal (N = 10 × 986 yr) from 10 epochs of
the 10 Myr serial run (slices 0, 1030, …, 9300), error / noise (noise = asm vs C fine kernel) at k = 1, 2, 3.
Config B (round 2) at k = 2: median ≈ 245 × noise, worst ≈ 2e4; k2 ≤ 2 at 0/10 epochs.

## Iteration 17 — Phase 0: requirements oracle (`parareal/phase0.py`)
Slice Jacobians J_F, J_G and coarse errors e0 at the 10 epochs; linearised parareal with
J_G(ε) = J_F + ε (J_G − J_F) on selected rows and e0 → ε0 e0.
* K = 2 at all 10 epochs needs ε ≈ 0.03 and ε0 ≈ 0.03 on the inner-planet rows (all six elements); or ε 0.03 / ε0 0.1
  if all rows are improved. ln a and λ rows alone are not enough. **Gate 0 passed** (ε ≈ 0.03 is what a second-order
  theory can plausibly deliver).
* Row test (nonlinear, oracle rows replaced by F): Venus, Earth and Mars rows (all elements) are all needed
  (median k2 270 → 3.5); Mercury and giants irrelevant.

## Iteration 18 — Phase 1: completeness of the first-order theory (`ph1eval.py`, `ph1py.py`, `slice_quality.py`)
* Python transform of all six components (λ with shear, k, h, q, p) vs ln a only: median k2 245 → 225 (5 epochs) —
  little gain. Heun option added to the Python coarse step (`mgc.py`); the Python ln a-only step reproduces C exactly.
* Close-pair grid `ng_hi = 64`: worst epoch k3 870 → 40 but median k2 only 245 → 75–110; 128 = 64 (converged).
* Per-slice errors (`slice_quality.py`, 200 slices from boundary 1000): λ rms Venus 1.4e-3, Earth 8.6e-4,
  Mars 2.5e-4 rad (ng_hi 64: 1.0e-3 / 7.0e-4 / 2.5e-4); k, h errors smooth (≈ 2 % secular-rate error).

## Iteration 19 — Phase 3: second-order Earth–Mars–Jupiter term from first principles (`so2.py`, `g3.py`)
* `so2.resonant_U2`: second-order averaged field U2 = ⟨D V1·W1 − D W1·V̄1 + ½ D²V0[W1,W1]⟩_K for one harmonic K from
  the first-order pair Fourier tables (products of different pairs sharing one planet; slow derivatives by finite
  differences). For K = −4λE + 8λM − 3λJ it reproduces the measured Mars term (amplitude match with the true
  frequency, phase correlation 0.92).
* `g3.G3`: the term added to config B as resonant dynamics, integrated analytically over the slice (ln a, λ incl.
  shear, k, h, q, p). Epoch 0: k2 245 → 106; epoch 6000 unchanged. Mars per-slice λ error 2.5e-4 → 1.9e-4.
  Python cost 0.58 s per call.

## Iteration 20 — what limits Venus and Earth (`resid_at.py`, harmonic regressions)
* No single slow harmonic explains the Venus/Earth per-slice λ errors (best: 8λV − 13λE, R² = 0.05).
* The λ error is driven by the **pointwise error of the transformed mean ln a at the slice start**: regression of
  the λ error on the change of transformed mean ln a between boundaries has slope 0.75 n T (as expected), and
  replacing the start value by the slice average of a densely sampled fine run removes 2–5 × of the λ error.
* Dense fine run (6 d sampling, 3 kyr from boundary 1000, ng_hi 64): pointwise residual of the mean ln a 2.3e-8
  (Venus), 2.6e-8 (Earth), 4.0e-8 (Mars) — 4 × worse than at t = 0 — with a broadband line spectrum: high
  harmonics and three-body combinations around the Venus–Earth 3:5 near-resonance (P = 8 yr): 6λV − 10λE (4.1 yr),
  −3λV + 5λE − 2λJ (22 yr), 2λV − 5λM − 4λJ (3.9 yr), plus Mars–Jupiter/Saturn combinations. These are
  **second-order short-period terms** (squares and cross products of first-order terms with the 3:5 small divisor),
  i.e. W2 — not only the resonant U2. Long-period content of the mean ln a (1-yr binned) is only 1.3e-8 (235 yr) for
  Venus/Earth and 3.7e-8 (1500 yr, E–M–J) for Mars; its λ effect is small except for Mars.
* Grid resolution is converged (ng_ii up to 128, ng_io 64 give the same boundary jumps of the mean ln a).

## Iteration 21 — the real first-order errors: small divisors and missing secular terms
Diagnosis with oracles first (`parareal/oracle_lna.py`, `jacblock.py`):
* A perfect mean ln a (slice average of a densely sampled fine run, start and end) alone: median k2 245 → 85 — not
  sufficient. Linearised block test: setting the coarse error e0 of the V, E, M rows to zero gives median k2 7.8,
  but only if **all** elements are fixed (ln a + λ only: 123; k, h, q, p only: 131) → every element of the inner
  planets must improve together.
* Per-slice defects of k, h (config B) are white (no smoothness), ~4e-5: the untransformed first-order short-period
  terms of k, h (config B transforms ln a only).
* **Small divisors from osculating mean motions** were the main error of the first-order transform: tables rebuilt
  at the *mean* semi-major axes (one fixed-point iteration, `mgc.Coarse(meann=1)`) reduce the pointwise error of the
  mean ln a at the hard epoch 1030 from 5.5e-7 to 9e-9 (the Venus–Earth 8:13 term, divisor 2π/240 yr, is extremely
  sensitive to the 1e-5 short-period part of the osculating a). Per-slice λ error Venus 1.6e-3 → 2.6e-4.
* Python coarse step (`ph1mn.py`), parareal at the 10 epochs, error/noise:

| coarse step | k=2 median | k=2 worst | k=3 median | k=3 worst |
|---|---|---|---|---|
| config B (C, round 2) | 245 | 1.9e4 | 1.5 | 870 |
| ln a transform, mean-n tables | 47 | 400 | 1.0 | 11.8 |
| all-component transform, mean-n | 29 | 352 | 0.95 | 12.3 |
| + Earth–Mars–Jupiter 2nd-order term (`g3`, FD step 1e-4) | 10.2 | 48.7 | 1.04 | 2.14 |
| + giants' WH step 180 d | **3.5** | 49 | 0.76 | 1.59 |
| + giants' step 90 d | 3.9 | 49 | 0.87 | 1.82 |

* The E–M–J term first raised the k=3 floor (2.75): its slow-element derivatives by finite differences (h = 1e-6
  on top of the finite-difference dE/dv) were noisy (G not smooth at the 1e-10 level); h = 1e-4 fixed it.
* With all inner improvements, the giants (WH at 360 d) limited k=2 at several epochs; 180 d is converged.

## Iteration 22 — secular rates: grids and Uranus/Neptune
* Grid convergence of the secular rates per pair: Venus–Earth needs ≥ 48² (32²: 0.5 % error — the ~1e-5 per-slice
  k, h error), Earth–Mars ≥ 32², all others 16–24² suffice.
* **Uranus and Neptune were missing from the inner planets' secular rates** (maxj = 5 for the transform and the
  rates): their contribution (~1e-3 of Jupiter's) was the per-slice q, p error (Venus 9.7e-7 → 9.4e-8 with them).
* Secular drift by RK4 over sub-steps with the giants sampled from their WH run (`nsub`): Mars k error
  1.3e-5 → 5.7e-6 (Jupiter's short-period variation was aliased by the two-point Heun rule).
* Numerical filtering of the start mean elements (10-yr fine run, first-order transform point by point, Hann
  average; `gfilt.GF2`): Mars λ per-slice error 1.7e-4 → 8.9e-5, Venus/Earth only 10–15 %.
* **Parareal with all of it** (`all_mn1_g3_dt180_nsub10_un`: all-component transform at mean-n tables, E–M–J
  term, giants 180 d, RK4 secular drift with 10 sub-steps and Uranus/Neptune): error/noise at k = 1, 2, 3 —
  median 601 / **1.81** / 1.04, worst 2.0e3 / 13.2 / 1.59; **k2 ≤ 2 at 6/10 epochs** (config B: 0/10).
  Remaining hard epochs: 3400 (12.2), 4500 (13.2), 8160 (8.1), 1030 (4.6). Python cost ~7 s per call (not yet
  optimised; Phase 5).

## Iteration 23 — giants' short-period aliasing, filter, and the secular-matrix error
* Cheaper secular drift: RK4 with 2 sub-steps is as good as 10 (median k2 1.53, k2 ≤ 2 at 6/10).
* Mars (k, h) error at 1 Myr: a single line explained by **−λJ + 2λS** (the 60-yr Jupiter–Saturn term, R² = 0.73):
  the giants' 60-yr eccentricity oscillation was point-sampled (aliased) by the RK stages. Box-averaging the giants'
  elements over 120 yr around each stage (`gavg`, dense samples of the giants' WH run): Mars k error 1.2e-5 → 5.7e-6
  with 2 sub-steps.
* Filtered start/end mean elements (L = 10 yr) on top: helps some epochs (3400: 12 → 0.9; 8160: 8 → 2), hurts
  others (0: 0.6 → 3.5; 7000: 0.9 → 3.5); 5/10 at k2 ≤ 2 — not kept as default.
* **The remaining Mars/Venus/Earth (k, h) errors are a secular-matrix error**: at 4.5 Myr the per-slice errors are
  exactly linear in the eccentricity vectors (R² = 1.00): δz_M = (2e-5 − 2.1e-4 i) z_M per slice, i.e. Mars' g
  is off by 0.044″/yr (0.25 %), plus V–E couplings of the same size — the second-order (in the masses) secular terms
  missing from a first-order theory. λ and ln a errors are white (not predictable from slow variables).
* Empirical correction (`seccal.py`): a 4 × 4 complex matrix A (δz = A z_inner per slice) fitted from the coarse
  defects F − G of the preceding 50–150 slices predicts the next 20 slices' z errors 10–50 × better
  (1e-5 → 2–7e-7). Fits from other epochs transfer only 2–5 × (A depends on the slow state) and short fits
  (≤ 30 slices) are ill-conditioned. In production the defects of the previous converged windows are free
  (parareal computes F − G anyway); A is held fixed within a window so G stays one function per window.

## Iteration 24 — second-order theory proper (`w2pair.py`, `w2triple.py`)
* Jacobian attribution of the new G (`jac_new.py`, `jacblock2.py`): at 4500 the k = 2 limit (12) is the Mars ln a row
  (exact row → 2.5), at 8160 (9) the Venus and Earth rows (λ, e0 and J equally). The E–M–J term removes 60–100 % of
  the Mars ln a/ln a Jacobian mismatch; its amplitude agrees with a direct measurement (|A| 1.45e-7 vs 1.55e-7) but
  its phase differs by ~30°; integrating it numerically (pendulum feedback) changes nothing; no other three-planet
  term with P > 400 yr is significant (`so2_scan.py`). Empirical amplitude refits of the term: no gain.
* Pointwise residual of the mean ln a with **live** tables (tests with frozen tables over 200 yr were misleading):
  Mars 5.9e-8, Earth 8.9e-9, Venus 4.9e-9 at 4.5 Myr. Three-body test (Sun+Mars+Jupiter): the residual scales as
  m_J² — second order.
* `w2triple.TripleW2`: second-order transform W2 and secular rates U2 on the 2-D (pair) or 3-D (triple) torus by
  grid products f2 = Σ ∂V1/∂x · W1 + 1.125 n W1_lna² and FFT division; cross-only triples + pairs combine without
  double counting (verified: identical to the full triple); near-resonant harmonics (P > 500 yr) excluded (kept as
  resonant terms). **Displacements of the giants must be excluded**: G already uses the giants' osculating elements
  from their WH run, and the J–S great-inequality W1 gave wrong huge terms. Regression of the measured residual on
  each group's predicted W2: E–M–J cross 1.01, V–E–J cross 0.98–1.18 (correct), but they explain only part of the
  Mars residual (3.1e-8 of 5.9e-8).
* **The second-order secular rates U2 explain the secular-matrix error exactly**: U2·T per slice at 4.5 Myr
  = −(calibrated δz) for all four inner planets (Mars: −1.49e-5 + 5.8e-6 i vs 1.49e-5 − 4.0e-6 i), dominated by the
  Mars–Jupiter same-pair term. With W2 + U2 (14 pairs + V–E–J, E–M–J, V–E–M triples) in the coarse step the
  per-slice k errors drop 13–14 × (Venus 4.7e-6 → 3.7e-7, Earth 4.6e-6 → 3.5e-7, Mars 5.7e-6 → 4.1e-7), q 2–6 ×;
  λ errors only 5–20 %.
* Parareal with W2 + U2 in the coarse step: **worse** (k1 3–10 × larger; 1030: 4.9 → 40). Diagnosis: the mean ↔
  osculating maps must be used in one consistent form. The first-order step used the "X-form" (W1 evaluated at the
  osculating angles both ways); W2 is the second-order term of the "Y-form" (X = Y + W1(Y) + W2(Y)). Switching to the
  Y-form with an incomplete W2 (14 pairs + 3 triples) degraded the Venus/Earth λ (8:13, 3:5 small divisors).
* **U2 only** (second-order secular rates of k, h, q, p from the 14 pairs + 3 triples, transforms unchanged):
  k2 / noise per epoch 0.68, 1.49, 2.52, 11.4, 9.0, 1.36, 1.41, 0.87, 9.32, 1.22 — median 1.45, k2 ≤ 2 at 6/10
  (1030: 4.9 → 1.5; 2000 0.8 → 2.5), k3 worst 1.61. The three hard epochs (3400, 4500 Mars; 8160 Venus/Earth)
  remain at 9–11.
* Cost reduction of the second-order parts (E–M–J term U and U2): **frozen per window** (computed at the window
  start) is bad (k1 × 3–13, epoch 0 k2 0.7 → 27: U ∝ e e^{−iϖ} rotates by ~50° over a window); **frozen per slice**
  (computed at the first coarse evaluation of each slice, reused in later iterations) is as good as live
  (k2 per epoch 0.68, 1.7, 2.52, 11.5, 8.9, 1.25, 1.3, 0.89, 9.3, 0.9; k3 worst 2.7). So the expensive
  second-order coefficients cost only during the k = 0 coarse sweep.
* ln a-only transform with all other improvements: k2 ≤ 2 at 1/10 (worst 104) — the all-component transform is
  required for K = 2 (but k3 ≤ 2.1 × noise at all epochs even without it).

## Iteration 25 — Phase 7: continuous pipelining instead of windows (simulated with the real numerics)
* Long windows (config B, 40 slices, `longwin.py`): k = 3 converges only the first 10 boundaries, k = 4 all 40
  (≤ 1.5 × noise) — the binomial growth of the parareal error with distance from the last exact boundary. With the
  improved Python G (30 slices): k = 2 up to boundary 10–23, k = 3 all 30.
* Throughput bound of a continuously pipelined parareal: S ≤ min(P / (K(1 + r)), 1/r) — the coarse chain of the
  final iterate is serial (1/r ≈ 11 for config B).
* Event-driven simulators with the real F and G computed in simulated-time order (`pint_sim.py`: iterate-indexed
  sliding window, block entry, gated fine work; `pint_sim2.py`: asynchronous latest-iterate policy with a
  finalisation rule "predecessor final and the fine solve used was on an iterate within tol of its final value").
  Config B, P = 10, 80 slices: S = 1.4–2.2 at err ≤ 2–5 × noise (fine solves per slice 3.2–6, coarse updates per
  slice up to 37) — **not better than the windowed pipeline** (2.3 modelled). Every change of a slice's start state
  after its iterations began costs iterations; the windowed synchronous scheme is the efficient one at P ≈ 10.
* Consequence for the coarse step: in S = N/(N r + K(1 + r)) the improved G pays only if it is nearly as cheap as
  config B (r ≲ 0.12 for K̄ = 2.4). Cost-reduction options: run the k = 0 prediction sweep with config B and use the
  improved G only from k = 1 on (its evaluations on the k = 0 iterates run in parallel with the fine solves; the
  per-slice frozen second-order coefficients are computed there, off the critical path).
* Variant: k = 0 prediction sweep with config B (cheap), improved G from k = 1 on: k2 ≤ 2 at 1/10 (worst 443) — the
  k = 0 prediction quality is essential for K = 2; the improved G must be on the critical path.

## Iteration 26 — Phase 5: cost of the improved coarse step in C (`pint/pint_disp6.c`)
* C implementation of the all-component first-order transform (ln a, λ with shear, k, h, q, p; same grids as config
  B with ng_hi = 64; optional mean motions for the divisors), validated against the Python theory (relative
  agreement 1e-10 for ln a, λ, k, h; 1e-7 for q, p).
* Cost (one core, idle): config B G call 1.01 ms; one ln a transform pass 0.45 ms; **one all-component pass
  2.13 ms**. The improved G needs at least two such passes (start, end) plus the mean-n evaluation, the giants at
  180 d and the RK4 secular drift: ≥ 4.8 ms → **r ≥ 0.40** (config B: 0.084). Optimisations in reach (analytic
  dE/dv, packed real FFTs, pruned evaluation) might halve it → r ≈ 0.2–0.25.
* Speedup model S = N/(N r + K̄(1 + r)), N = 10: config B adaptive (K̄ 3.29, r 0.084) 2.3; improved G (K̄ ≈ 2.4)
  at r = 0.40: 1.4; at r = 0.20: 1.9; break-even needs r ≈ 0.12. **At ≈ 10 cores the round-3 coarse step does not
  pay in wall-clock time** — the accuracy gain is real, but the coarse cost enters N times through the serial
  k = 0 sweep.

## Round 3 result (iterations 17–26)
* The first-order mean-element theory of round 2 was far from first-order accurate. Fixed: small divisors from
  mean (not osculating) mean motions, all-component transform, Venus–Earth secular grid (≥ 48²), Uranus/Neptune in
  the secular rates, giants' 60-yr term aliasing (giants at 180 d, averaged forcing), Earth–Mars–Jupiter
  second-order resonant term, second-order secular rates (Mars–Jupiter dominated; explains the 0.25 % secular
  frequency error exactly). Python prototype: **k = 2 at ≤ 2 × noise for 6/10 test windows (median 1.45 × noise,
  config B: 0/10, median 245), k = 3 ≤ 1.6–2.7 × noise for all windows (config B worst 870)**.
* Not achieved: Gate 4 (8/10). The three hard windows (Mars at 3.4 and 4.5 Myr, Venus/Earth at 8.2 Myr) need the
  complete second-order short-period theory; the partial W2 (14 pairs + 3 triples, Y-form) made things worse
  because the mean ↔ osculating maps must be complete and consistent at second order.
* Not worth it at 10 cores: the improved coarse step costs r ≈ 0.2–0.4 in C (iteration 26); continuous
  pipelining (iteration 25) does not remove the startup term efficiently. **Config B with the adaptive rule remains
  the production configuration (≈ 2.3× on 10 cores).** Round-4 directions: `RESEARCH_PLAN_R4.md`.

## Phase 8 — first many-core measurement (user, AMD Ryzen 9 7950X, 16 cores, AVX-512; Oct 21, 2026)
`./pint_pipe N 360000 3 20 2 1 1 3e-6 6` (adaptive rule, 20 windows from t = 0):

| N (threads) | K̄ | wall | serial | speedup | fine busy | coarse busy | final Δpos |
|---|---|---|---|---|---|---|---|
| 10 | 3.50 (10× K=3, 10× K=4) | 0.772 s | 1.407 s | **1.82** | 5.06 core-s | 0.76 core-s | 1.2e-5 AU |
| 16 | 3.90 (2× K=3, 18× K=4) | 0.957 s | 2.273 s | **2.37** | 10.15 core-s | 1.32 core-s | 2.6e-5 AU |

Analysis:
* Fine slice 7.0 ms (serial 1.407 s / 200 slices) vs 12 ms on the 2.1 GHz Xeon; coarse ≈ 0.84 ms per call
  (0.76 core-s / ~20·10·4.5 calls) → **r ≈ 0.12** on Zen 4 (0.084 on the Xeon): the AVX-512 fine kernel gains more
  from the faster core than the scalar coarse step.
* Model S = N/(N r + K̄(1 + r)) with the measured K̄ and r = 0.12: N = 10 → 1.95 (measured 1.82, 93 %), N = 16 →
  2.54 (measured 2.37, 93 %). The pipelined model holds on real many-core hardware within ~7 %.
* K̄ is higher than in the 10 Myr runs (3.29): the first 20 windows from t = 0 are harder than average, and
  16-slice windows need more iterations (binomial growth, iteration 25).
* The coarse cost is now the main lever: N r is 1.2 (N = 10) and 1.9 (N = 16) of the 5–6 units of wall time.

### Phase 8, user run 2 (100 windows)
| N | tol | K̄ | wall | serial | speedup | model (r = 0.12) |
|---|---|---|---|---|---|---|
| 12 | 3e-6 | 3.47 (54/45/1 at K=3/4/5) | 4.09 s | 8.46 s | **2.07** | 2.26 |
| 16 | 1e-5 | 3.45 (57/41/2) | 4.49 s | 11.26 s | **2.51** | 2.78 |
Model/measured ≈ 0.90–0.92. The N r start-up term is 1.4–1.9 of ≈ 5.5 units of wall time. (tol 1e-5 final
difference 6.6e-4 AU after 1.58 Myr vs 7.7e-5 AU after 1.18 Myr at 3e-6 — needs the frequency check before use.)

## Iteration 27 — parallel coarse step (`pint_coarse_par`, `pint_set_threads`; `pint_pipe` args 10–11)
* The coarse step is split into tasks: grids per planet, one task per inner–inner pair and per packed pair of giants,
  GR per planet; the giants' WH run (≈ 1/3 of the step) runs concurrently with the start tables. Results are reduced
  in the serial order → **bitwise identical** to the serial step (`test_par.c`).
* Driver: `pint_pipe ... [helpers] [parmode]`: a pool of helper threads used in the k = 0 sweep (mode 1, when the
  other cores are mostly idle) or in all iterations (mode 2). Also: the coarse step is skipped for slices whose
  start is already exact (n ≤ k − 1) — it was computed and discarded on the critical path.
* Measured here only on 2 vCPUs (1 helper: 0.99 → 0.78–0.82 ms). Cost split (serial): two table passes 0.73 ms,
  giants' run 0.35 ms → with enough helpers the k = 0 step is bounded by the giants' run (≈ 0.35–0.45 ms, ≈ 2.3×
  faster). Expected on the 7950X (N = 16, K̄ 3.45): start-up term 1.9 → ≈ 0.8, model S ≈ 2.8 → ≈ 3.3.

### Phase 8, user run 3 (7950X, N = 16, 100 windows, tol 3e-6; K̄ 3.80: 26/68/6 at K = 3/4/5)
| variant | wall | speedup |
|---|---|---|
| baseline (new build, skips unneeded G) | 4.74–4.78 s | 2.36–2.38 |
| 8 helpers, k = 0 sweep only (mode 1) | 4.79 s | 2.37 |
| **8 helpers, all iterations (mode 2)** | **4.27 s** | **2.65** (coarse busy 5.17 → 3.19 core-s) |
| 15 helpers, mode 1 | 4.78 s | 2.37 |
Final position difference identical in all runs (2.84e-4 AU after 1.58 Myr; tol 1e-5 gave 6.6e-4 — same order, i.e.
round-off divergence, not the tolerance). Mode 1 gave nothing: sleeping helpers need ~0.1 ms to wake up, and the
step has 4 fork-join rounds plus a thread creation for the giants' run. Fix (iteration 27b): persistent giants
worker, lock-free task counter, workers spin `PINT_SPIN` pause instructions (default 20000) before sleeping.
Here (2 vCPUs, oversubscribed) spinning hurts (1.10 ms) and `PINT_SPIN=0` is best (0.68 ms vs 1.03 serial); to be
measured on the 7950X.

### Phase 8, user run 4 (7950X, N = 16, 100 windows, tol 3e-6, K̄ 3.80)
| variant | wall | speedup | coarse busy |
|---|---|---|---|
| 8 helpers, mode 2, spin 20000 | 4.22 s | 2.67 | 2.45 core-s |
| **8 helpers, mode 2, `PINT_SPIN=0`** | **4.17 s** | **2.71** | 2.45 |
| 8 helpers, mode 1 | 4.69 s | 2.42 | 4.61 |
| 15 helpers, mode 2 | 4.29 s | 2.64 | 2.57 |
Coarse step 2.1× faster (5.17 → 2.45 core-s) → speedup 2.37 → 2.71 (+14 %). Default `PINT_SPIN` set to 0.
Model with r_eff ≈ 0.057 predicts 3.25; measured/model = 0.83 — the remaining overheads (≈ 4.4 fine solves per slice
incl. speculative ones; energy projection and position checks on the chain; helpers sharing cores with fine solves)
now matter. Upper bound at r → 0 with 4.4 fine solves per slice: 16/4.4 ≈ 3.6.

## Phase 9 — new session with direct access to host2 (7950X), Oct 21, 2026
Machines are documented in `MACHINES.md`. From now on the agent runs the 7950X measurements itself (`ssh root@host2`,
repo `/home/ubuntu/rebound`, outputs in `whfaster/pint/runs/p9/`). numpy 1.26 / scipy 1.11 installed there (apt).
Build check: `./test_par 8` → bitwise identical: yes; coarse step serial 0.631 ms, parallel (8 helpers) 0.219 ms (idle machine).

### Iteration 28 — window size with helpers (the runs suggested at the end of Phase 8), tol 3e-6, Kmax 6, mode 2
| N | helpers | K̄ (K=3/4/5) | wall | serial | speedup | final Δx |
|---|---|---|---|---|---|---|
| 16 | 8 | 3.80 (26/68/6) | 4.14–4.21 s (4 runs) | 11.19–11.50 s | **2.66–2.74** | 2.84e-4 AU |
| 14 | 4 | 3.78 (29/64/7) | 3.95 s | 9.87 s | 2.50 | 1.55e-4 |
| 12 | 4 | 3.47 (54/45/1) | 3.55 s | 8.50 s | 2.39 | 7.7e-5 |
| 12 | 8 | 3.47 | 3.58 s | 8.48 s | 2.37 | 7.7e-5 |
N = 16 stays best: the shorter window saves iterations at N = 12 (K̄ 3.47), but not enough to pay for 25 % fewer
slices. Leaving physical cores free for the helpers does not help. Run-to-run spread of the speedup ≈ ±2 %.

### Iteration 29 — where the time goes on 16 cores (event trace), and the 10 Myr contract at N = 16
**Trace** (`PINT_TRACE=file ./pint_pipe 16 360000 3 100 2 1 1 3e-6 6 8 2`, new optional per-thread event log;
analysis: mean over 100 windows). Window = 41.1 ms:
| part | ms per window |
|---|---|
| k = 0 coarse sweep (16 × 0.347 ms) | 5.6 |
| iterations k = 1 … K: increments 9.79, 9.13, 8.88, 8.68 (74 %), 8.57 (6 %) | 34.4 |
| tail (decision → threads joined) + gap between windows | 0.7 + 0.2 |
* Fine slice **under full load 8.28 ms** vs 7.0 ms on an idle machine (16 concurrent serial runs pinned to cores 0–15:
  +11 % on CCD 0, +17 % on CCD 1). Clock effect, the serial reference runs at single-core boost. Governor: powersave
  (not changed). CPUs 0–15 are distinct physical cores, 16–31 their SMT siblings (pinning is correct).
* Coarse step on the chain 0.29–0.37 ms (0.22 ms idle): the helpers and the giants' worker run on SMT siblings of
  busy cores.
* Per-hop driver helpers (simulation creation in `project_energy`, `energy_el`, `poserr_el`) ≈ 30 µs → ≈ 1 %. Not worth it.
* Fine solves: 55.5 completed + 12.1 abandoned per window (3.5 + 0.8 per slice); 75 % core utilisation.
* SMT: 32 concurrent fine runs give +19 % total throughput but each is 1.68× slower → not useful for parareal latency.
→ K̄ × T_F is 77 % of the window; the coarse chain 14 %; overheads ≈ 5 %.

**10 Myr contract at N = 16** (hard ICs of iteration 14, `pint_long … 625 16 360000 16 3 1 τ 0 2 1 1 1e-6`; 7950X:
serial 71 s, parareal (non-pipelined driver) 43–50 s; `runs/p9/v16*`). Differences vs serial C [″/yr]:
| run | K̄ (K = 3/4/5/6) | max inner-planet shift | giants | energy vs serial C at 10 Myr |
|---|---|---|---|---|
| serial asm (round-off ref.) | — | 5.9e-10 | 1.5e-8 | +5.3e-12 |
| **τ = 3e-6** | 3.77 (171/437/16/1) | **9.3e-8** | 4.8e-8 | +6.1e-12 |
| τ = 1e-5 | 3.35 (412/206/7) | 3.8e-7 (s2) | 1.0e-7 | +6.5e-12 |
| τ = 3e-5 | 3.11 (557/66/2) | 7.3e-7 (s2) | 1.5e-7 | +2.5e-12 |
| K_min = 2, τ = 3e-6 (K_max 5) | 3.77 (172/435/18) | 3.2e-7 (s1) | 5.6e-8 | — |
* τ = 3e-6 is validated at N = 16 (same fidelity as at N = 10). Energy unbiased in all runs.
* **No window converges at k = 2** with config B (K_min = 2 changes nothing); but capping K at 5 instead of 6 in a
  few windows moved s1 by 4e-7 → the scatter of a single 10 Myr realisation is ≈ 3e-7 ″/yr, so τ = 1e-5 (3.8e-7)
  is not distinguishable from τ = 3e-6 with one realisation. Accepting τ = 1e-5 (K̄ 3.77 → 3.35, ≈ +11 % speed)
  needs an ensemble test. Kept: τ = 3e-6.

### Iteration 30 — staged parallel coarse step: giants on the calling core (KEPT, bitwise identical)
The giants' WH run (the longest single task, ≈ 0.17 ms) ran on a helper thread, i.e. on the SMT sibling of a core
busy with a fine slice (≈ 1.7× slower under load), and the step had four fork-join rounds plus the giants' worker.
New schedule in `pint_coarse_par` (`coarse_staged`): the **caller integrates the giants** on its own core while the
helpers run, as one staged job, start grids → start pairs/GR → reduction + predictor → end-table grids of the inner
planets → end-table inner–inner pairs and GR (these need only the inner planets' predicted elements, so they overlap
the giants' run); after the giants: their grids and the inner–giant pairs (second job), reduction and Heun corrector
in the caller. The pool is now a staged pool (per-stage task counters with a barrier between stages); the giants'
worker thread is gone. Serial path refactored into the same pieces (`table_setup/_reduce`, `inner_predict/_correct`).
* `test_par`: bitwise identical (4, 8, 15 helpers); idle coarse step 0.219 → 0.165 ms (4 helpers), 0.191 (8).
* `pint_pipe 16 360000 3 100 2 1 1 3e-6 6 h 2` (7950X): h = 8: **wall 3.98 s, speedup 2.84–2.86** (was 4.14–4.21 s,
  2.66–2.74); h = 4: 2.81–2.82; h = 2: 2.61. Coarse busy 2.43 → 2.06 core-s. Final Δx unchanged (2.84e-4 AU), i.e.
  the same trajectory bit for bit.

### Iteration 31 — persistent slice threads, windows handed over without a join (KEPT, bitwise identical)
`pint_pipe` created 16 threads per window and joined all of them, including threads still finishing an abandoned
speculative fine chunk (≤ 0.83 ms). Now the slice threads persist over all windows; the window state (U, semaphores,
stop flag, update maxima) is double-buffered by window parity. At the decision the last slice prepares the other
parity's state, posts `win_start` (slice 0 starts the next window immediately) and a `win_done` barrier (a thread may
finish k = K_max before the decision and must not touch the next window's semaphores before they are re-initialised).
Threads still in abandoned work join the next window's chain when they are done; slice n is needed only n coarse
steps after the window start.
* Checked: unreachable τ (stops at K_max every window), fixed K (τ = 0), N = 4, 6, 16: no deadlock; same K histogram
  and final Δx (2.84e-4 AU) as before → same trajectory.
* `pint_pipe 16 360000 3 100 2 1 1 3e-6 6 8 2`: **wall 3.75–3.77 s, speedup 2.97–3.05** (iteration 30: 3.98 s, 2.85).

### Iteration 32 — structure and coarse-model knobs at N = 16 (all rejected)
New trace after iteration 31: window 37.7 ms = k = 0 sweep 4.3 ms (16 × 0.26) + K̄ × ≈ 8.9 ms. The increment per
iteration equals the **slowest** slice's fine time + one hop (all slices are on the diagonal of the chain); cores
8–15 (second CCD) are 3–5 % slower than 0–7 under load (8.45–8.67 vs 8.04–8.22 ms).
| variant (`pint_pipe … 3e-6 6 8 2`) | K̄ | speedup |
|---|---|---|
| reference N = 16 × 986 yr | 3.80 | 2.92–3.05 |
| N = 32 × 493 yr on all 32 SMT threads (4 / 8 helpers) | 3.74 | 2.13 / 2.01 |
| N = 24 × 657 yr | 3.92 | 2.62 |
| N = 20 × 789 yr | 3.95 | 2.45 |
| N = 16 × 789 yr | 3.77 | 2.85 |
| N = 16 × 1183 yr | 4.27 | 2.80 |
| N = 16 × 1380 yr | 4.39 | 2.70 |
| **mean-n start tables** (`meann`, new option: start tables rebuilt at the inner planets' mean ln a) | 3.82 | 2.73 |
| 64² grids for V–E and E–M (`ng_hi` 64) | 3.67 | 2.36 |
| mean-n + 64² (≈ the Python "ln a, mean-n" row of iteration 21) | 3.43 | 2.04 |
| giants' coarse step 180 d | 3.82 | 2.75 |
| giants' coarse step 540 / 720 / 1000 d | 6.00 (never converges) | 1.96 |
* SMT does not pay (each fine solve 1.68× slower), more or shorter slices need as many iterations; 986 yr is optimal.
* The round-3 ingredients do not translate into fewer iterations at N = 16 with the update-based stop: the mean-n
  rebuild (free on the parallel critical path, 0.17 ms idle) changes K̄ by nothing; with 64² close-pair grids K̄
  drops by 10 % but the big pair tasks make the coarse step 2× slower. Stopping at K = 3 needs the k = 2 iterate
  within τ at all 16 boundaries — the binomial growth along the window dominates, not the per-slice coarse error.
* New `pint_cfg.meann` (default 0) and environment overrides for experiments in all drivers: `PINT_MEANN`,
  `PINT_NGHI`, `PINT_DTOUT` (`ptc.Cfg` updated). Bitwise serial = parallel checked for all four combinations.

### Iteration 33 — tolerance ensemble, paired fine solves, kernels, and what binds the stopping rule
**Ensemble** (`runs/p9/ens.sh`, `pint/ens_ana.py`): 7 IC realisations (all planets x + δ, δ = 5e-7 … 3e-6, −1e-6),
10 Myr, N = 16, `pint_long` with Heun + energy projection; secular-frequency shifts vs the serial run of the same ICs:
| τ | K̄ | inner rms / max [″/yr] | giants rms | mean shift per frequency | energy at 10 Myr vs serial |
|---|---|---|---|---|---|
| 3e-6 | 3.75–3.80 | 6.6e-8 / 2.4e-7 | 1.0e-8 | all within 1σ of 0 | −9e-12 … +6e-12 |
| 1e-5 | 3.35–3.42 | 1.4e-7 / 3.8e-7 | 2.7e-8 | ≤ 1.7σ (g1 +6.9e-8 ± 4e-8, s2 +1.2e-7 ± 8e-8) | −7e-12 … +8e-12 |
Both are far below the WH truncation error at 6 d (g1 ≈ 9e-7, s1 ≈ 3.4e-7 ″/yr) and show no significant bias.
→ **τ = 1e-5 is validated as a fast setting** (pipeline: K̄ 3.45, **speedup 3.24**); τ = 3e-6 stays the default.

**Paired fine solves** (`pint_pipe` arg 12 `M = 2`: a thread owns two consecutive slices and steps both in one
`reb_integrator_whfast512_steps_pair` call): on the 7950X the paired kernel gives no throughput (2 × 30000 steps
7.01 ms vs 60000 single 7.15 ms; 2 × 60000 → ×1.02), unlike the Xeon (×1.23). N = 32 × 493 yr, M = 2: 2.24–2.36.
Option kept (may pay on Intel AVX-512), not used.
**Kernels on the 7950X:** the assembly kernel (1) is 8 % faster than the C AVX-512 kernel (2): serial 10.49 vs
11.39 s; pipeline 3.53 vs 3.86 s (same speedup 2.97). Production on this machine: kernel 1.

**What binds the stopping rule** (`PINT_UPDDIAG=1/2`: per-planet and per-element maxima of the boundary updates,
N = 16, τ = 3e-6, 100 windows; medians):
| k | worst planet (count) | Venus λ·a | Venus |Δz|·a | Venus a·Δln a |
|---|---|---|---|---|
| 2 | Ven 45, Mar 35, Ear 11, Mer 8 | 2.5e-4 | 1.2e-6 | 1.4e-8 |
| 3 | Ven 42, Mar 34, Ear 18 | 3.7e-6 | 4.2e-8 | 2.9e-10 |
| 4 | Ven 37, Mar 22, Ear 13 | 1.3e-7 | 1.3e-9 | 1.1e-11 |
* Always the inner planets (giants ≤ 1e-7 at k = 3) and always **the phase λ**: it is the Kepler-shear image of the
  ln a error, δλ ≈ 1.5 n T δln a (Venus: 1.5 × 10.2/yr × 986 yr × 2.9e-10 ≈ 4.4e-6). The e/i vectors are 100× smaller,
  so a secular-matrix calibration (it. 23) cannot help. A relative (per-a) position criterion would not either.
* What helps is a better short-period ln a model of Venus–Earth: mean-n + 64² grid for **V–E only**
  (`PINT_ALPHAHI=0.7`) gives K̄ 3.80 → 3.43 (τ 3e-6) and 3.45 → 3.15 (τ 1e-5; 84 % of windows at the K_min = 3
  floor; 128² = 64²), but the single 64² task appears three times in sequence on the coarse critical path (parallel
  step 0.43 ms) → speedup 2.34. Would need the 64² pair task split into row/FFT/eval chunks; estimated net gain
  ≈ 4 % → deferred.
* Upper bound of this scheme on 16 cores: S ≈ 16/K̄ × (T_F idle/T_F loaded = 0.84) − coarse chain (≈ 12 %) ≈ 3.0
  (τ 3e-6) / 3.4 (τ 1e-5); measured 2.95–3.05 / 3.24 → the remaining gains need K̄ < 3, i.e. K = 2 windows.

### Iteration 34 — production output in the pipelined driver; ≈ 1 Gyr runs started
`pint_pipe` can now write the trajectory: `PINT_OUT=prefix` → accepted boundary states of every window in
`<prefix>_pr.bin` (and the serial run in `<prefix>_ser<kernel>.bin`), same records as `pint_long`, written by the
main thread from a ring filled by the deciding slice thread (off the critical path); `PINT_PERT` (IC shift),
`PINT_SERONLY` (serial reference only). New experiment override `PINT_ALPHAHI` (pair ratio for `ng_hi`).
* Production path validated on the hard ICs (10 Myr, kernel 2, vs `v16t3_ser2.bin`): K histogram 173/435/16/1
  (pint_long: 171/437/16/1); max inner frequency shift 1.1e-7 ″/yr (τ 3e-6) and 4.2e-7 (τ 1e-5); energy vs serial
  +1.0e-11 / +1.1e-11. **10 Myr in 24.0 s (τ 3e-6) / 22.2 s (τ 1e-5)** vs 71 s serial.
* Started on host2 (`runs/gyr/run.sh`): 62500 windows = 0.986 Gyr, asm kernel: serial (`s0`), serial with all planets
  shifted by 1e-12 AU (`s12`, chaos twin = null spread), parareal τ 3e-6 (`p3`) and τ 1e-5 (`p10`). The serial runs
  share cores with the parareal runs (SMT siblings 30, 31), so these are accuracy runs, not timings.
* Rejected without implementation: splitting the 64² V–E pair task (iteration 33) — even fully chunked, three
  sequential 64² passes (fill, FFT 64², evaluation) keep the helper chain at ≈ 0.35 ms, longer than the giants' run;
  net ≈ +1 %.

### Iteration 35 — the accuracy budget: parareal error vs the WH time-step error
Literature check (new): Rein, Brown & Kanda 2025 (arXiv 2507.04987, "On the statistical convergence of N-body
simulations of the Solar System"): fixed-step symplectic integrations give physical secular frequencies and
instability rates for dt up to 32 d in 5 Gyr ensembles. → The right yardstick for parareal is the integrator's own
time-step error, and the fair comparison is a serial run with a larger (cheaper) step. No 2024–26 parareal work on
planetary N-body found; GParareal/nnGParareal (Gattiglio, Grigoryeva, Tamborrino, arXiv 2405.12182) learn F − G from
past iterates (our secant attempt, it. 13, is a crude version).
**Serial references at dt = 3, 8, 12, 24 d** (hard ICs, 10 Myr, kernel 2, `pint_pipe_dt` = `pint_pipe` with
`PINT_SERONLY PINT_DT`; `runs/dt/`), frequencies vs dt = 3 d ["/yr]:
| run | g1 | g2 | g4 (Mars) | g6 | s1 | s2 | s4 | cost |
|---|---|---|---|---|---|---|---|---|
| serial dt 6 d | +1.9e-7 | −3.3e-7 | **−5.05e-6** | −1.55e-6 | +3.8e-7 | +1.10e-6 | +2.17e-6 | 1 |
| serial dt 8 d | +2.4e-7 | −7.1e-7 | −1.02e-5 | −3.23e-6 | +1.4e-6 | +2.15e-6 | +4.40e-6 | 0.75 |
| serial dt 12 d | +2.1e-6 | −1.5e-6 | −2.55e-5 | −7.93e-6 | −2.9e-6 | +5.73e-6 | +1.09e-5 | 0.5 |
| serial dt 24 d | +9.4e-3 | +2.6e-3 | +2.9e-4 | −3.1e-5 | −1.3e-2 | +5.4e-4 | −2.3e-4 | 0.25 |
| parareal τ 3e-6 | +1.8e-7 | −3.2e-7 | −5.01e-6 | −1.60e-6 | +2.8e-7 | +1.03e-6 | +2.16e-6 | 1/3.0 |
| parareal τ 1e-5 | +3.9e-7 | −3.0e-7 | −5.37e-6 | −1.71e-6 | +2.2e-7 | +1.52e-6 | +2.40e-6 | 1/3.24 |
| parareal τ 3e-5 | +1.8e-7 | −3.5e-7 | −5.50e-6 | −1.70e-6 | +4.9e-7 | +1.83e-6 | +2.63e-6 | — |
* The WH error at 6 d scales as dt² (8 d/6 d: 2.0 = (64−9)/(36−9)) and is larger than the round-1 estimate:
  Mars g4 ≈ −5e-6 (−6.7e-6 extrapolated to dt → 0), g6 ≈ −1.6e-6, s2 ≈ +1.1e-6. It is a systematic bias.
* Parareal adds, at τ = 3e-6 / 1e-5 / 3e-5, at most 6 % / 7 % / 9 % to the largest error (g4) and is unbiased (ensemble,
  iteration 33) — per frequency up to 0.6× of the dt = 6 d error for s2 at τ 3e-5. A serial run at 8 d (1.33× faster)
  doubles every error. → At "the accuracy of the dt = 6 d reference" there is room for larger τ.
* dt = 24 d is broken here (Mercury g1/s1 off by 1e-2: a step–orbit resonance), unlike the 32-d ensemble statistics of
  Rein et al. (single realisation, different ICs; not pursued).

### Iteration 36 — a linear energy drift in long parareal runs, traced to WHFast512's coordinate matrices (FIXED)
The first 111 Myr of the 1 Gyr runs (iteration 34): serial runs dE/E −2.7e-12 (s0) and +1.9e-11 (s12, IC + 1e-12 AU),
parareal τ 3e-6 **+8.6e-11, growing linearly** (+1.44e-11 at 10 Myr, +4.0e-11 at 61 Myr, +8.6e-11 at 111 Myr;
≈ 7.7e-13 per Myr → ≈ 4e-9 over 5 Gyr). The 10 Myr ensemble (it. 33) could not see it (±9e-12 scatter).
* Not the element ↔ Cartesian conversions of the driver: one sim → elements → sim round trip, mean over 20000 states
  along a trajectory: dE/E −2.8e-18 ± 3.1e-18 (rms 4.5e-16).
* **Every WHFast512 integration call adds a fixed energy offset**: serial runs with 100 calls per 986 yr drift
  linearly by +7.1e-17 (asm) / +7.7e-17 (C) per call. Same steps in 10 calls vs 1 call, averaged over 3000 start
  states (`cb.c`): +7.2e-17 (asm, corrector 17, GR), +1.1e-16 (C, corr. 17 or 3), +7.3e-17 (C, corr. 0), +7.4e-17
  (C, no GR), +6.5e-17 (C, corr. 0, no GR) ± 0.6–1.5e-17. The parareal driver makes 10 calls per fine slice
  (abandonment chunks) → 1.1e6 calls in 111 Myr × 7–8e-17 = 8e-11, as observed. The serial reference with one call
  per slice has the same bias 10× smaller; one long call (as in the paper) has none.
* **Cause:** at every call WHFast512 converts inertial → Jacobi and back (synchronization) with two precomputed 8 × 8
  matrices (`mat8_inertial_to_jacobi`, `mat8_jacobi_to_inertial`) that are not exact inverses in floating point. Their
  product (long-double evaluation) is diagonal to 1e-19 with residuals diag − 1 = [−1.37, −0.22, +0.07, −0.05, +0.45,
  −1.37, −0.17, −0.77] × 1e-16 (Mercury … Neptune), the same for every state. Scaling planet i's x and v by (1 + e_i)
  changes its energy by 2 e_i |U_i|; weighted by the planets' potential energies (Jupiter +0.45, Saturn −1.37
  dominate): predicted **+7.6e-17 per call**, measured +6.5 … 7.7e-17. This affects every WHFast512 run with frequent
  synchronizations (outputs): e.g. 5 Gyr with an output every kyr → ≈ 4e-10 of pure drift.
* **Fix** (`src/integrator_whfast512.c`, opt-in so the default stays bitwise identical to upstream): new state field
  `jacobi_recursive` = 1 → the conversions use REBOUND's recursive Jacobi transformations
  (`reb_transformations_inertial_to_jacobi_posvel` / `..._jacobi_to_inertial_posvel`, per system, COM frame), whose
  round-off depends on the state. Per-call offset (C AVX2 kernel, 15000 states): corr. 17 + GR −1.7e-17 ± 0.7e-17,
  corr. 0 no GR +0.3e-17 ± 0.2e-17 (matrices: +1.1e-16 / +6.5e-17). The small remaining part comes with the corrector
  (its apply/unapply pair is not an exact inverse either), ≈ 6× below the original.
* Drivers: `PINT_JR` (default **1**) sets `jacobi_recursive` for all fine solves and serial references in
  `pint_pipe` and `pint_long`; `PINT_CHUNKS` (default 10) = kernel calls per speculative fine slice.
* Checks on host2 (separate worktree `/home/ubuntu/rb2`, so the running jobs keep their library): with
  `jacobi_recursive = 0` the new library is **bitwise identical** to the old one (serial runs, kernels 1 and 2). Per-call
  offset with `jacobi_recursive = 1`, 6000 states: asm −0.3e-17 ± 1.1e-17, C AVX-512 −1.5e-17 ± 1.1e-17 (was +7.2e-17 /
  +11e-17). `test_jr.c`: 1 vs 0 agree to ≤ 8e-16 after one 10-step call for N_systems = 1, 2, 4 (C kernels 3, 4).
* Larger sample (local, C AVX2 kernel, 30000 states each, `call_bias`): with `jacobi_recursive` = 1 the per-call
  offset is +0.5e-17 ± 0.5e-17 (corrector 17, no GR), +0.1e-17 ± 0.2e-17 (no corrector, GR), −0.9e-17 ± 0.5e-17
  (corrector 17 + GR); with 0: +7.4e-17 ± 0.5e-17. → No significant residual (|offset| < 1.5e-17 at 2σ, ≥ 85 %
  removed); the "corrector residual" of the first sample was a fluctuation.
* Long runs re-planned: the τ = 1e-5 Gyr run with the biased transformation was cancelled; added a serial run with the
  fix (`s0jr`) and, after the τ = 3e-6 run finishes, the τ/K_min scan and a 1 Gyr parareal run with the fix (`p3jr`).
* 1 Gyr parareal τ 3e-6 (matrix transform, `p3`): 0.986 Gyr in 3152 s (52.5 min; three serial runs shared SMT
  siblings, so not a timing), K̄ 3.78 (K = 3/4/5/6: 19544/38282/3854/820).

### Iteration 37 — tolerance and K_min scan against the time-step error (`pint/tau_ana.py`, `runs/tau/`)
Hard ICs, 10 Myr, `PINT_JR=1`, kernel 2; frequency errors vs the dt = 3 d serial run; "extra" = largest difference
to the dt = 6 d serial run (whose own largest error vs dt 3 d is 5.1e-6 ″/yr, Mars g4). Single realisation.
| run | K̄ (K histogram) | g4 | g6 | s2 | max err | extra vs dt 6 |
|---|---|---|---|---|---|---|
| serial dt 6 / 8 / 12 d | — | −5.1e-6 / −1.0e-5 / −2.6e-5 | −1.6e-6 / −3.2e-6 / −7.9e-6 | +1.1e-6 / +2.1e-6 / +5.7e-6 | 5.1e-6 / 1.0e-5 / 2.6e-5 | 0 / 5.2e-6 / 2.0e-5 |
| K_min 3, τ 3e-6 (JR 0) | 3.75 | −5.1e-6 | −1.7e-6 | +1.1e-6 | 5.1e-6 | 1.2e-7 |
| K_min 3, τ 1e-5 (JR 0) | 3.35 | −5.4e-6 | −1.7e-6 | +1.5e-6 | 5.4e-6 | 4.2e-7 |
| K_min 3, τ 3e-5 | 3.11 (557/66/2) | −5.5e-6 | −1.7e-6 | +1.8e-6 | 5.5e-6 | 6.9e-7 |
| K_min 3, τ 1e-4 | 3.02 (612/12/1) | −6.4e-6 | −1.7e-6 | +3.1e-6 | 6.4e-6 | 2.0e-6 |
| K_min 3, τ 3e-4 | 3.00 | −8.4e-6 | −1.7e-6 | +6.3e-6 | 8.4e-6 | 5.2e-6 |
| K_min 3, τ 1e-3 | 3.00 | −3.6e-6 | −1.7e-6 | −9e-8 | 3.6e-6 | 1.5e-6 |
| K_min 2, τ 1e-4 | 3.01 (K=2: 9) | −7.2e-6 | −1.9e-6 | +4.2e-6 | 7.2e-6 | 3.1e-6 |
| K_min 2, τ 3e-4 | **2.62** (243/377/5) | −2.3e-6 | **−1.2e-5** | −1.3e-6 | 1.2e-5 | 1.1e-5 |
| K_min 2, τ 1e-3 | **2.11** (554/71) | −6.1e-6 | **−2.9e-5** | +2.4e-6 | 2.9e-5 | 2.8e-5 |
Energy vs the dt = 6 d serial run: within ±1e-11 in all parareal runs.
* With K_min = 3 the iteration count reaches its floor (3.0) at τ ≈ 1e-4; τ 3e-5 / 1e-4 add 14 % / 40 % to the
  dt = 6 d error (i.e. "between dt 6 and 8 d"); τ = 1e-3 at the floor shows a smaller error by chance (single run).
* K_min = 2 finally gives K̄ < 3 (2.62 at τ 3e-4, 2.11 at 1e-3), but a **Saturn g6 shift** appears (−1.2e-5,
  −2.9e-5) that grows with τ: the k = 2 iterate still carries part of the coarse step's own giants error (WH at 360 d,
  whose g6 error is large) — a systematic error, not noise. Accuracy ≈ dt 8–12 d, which a serial run gets for 0.75–0.5
  of the cost; parareal at K̄ 2.1–2.6 is still ≈ 2× faster than such a serial run. To be checked with an ensemble.

### Iteration 38 — bias ensemble of the candidate settings; giants' coarse step 180 d for K_min = 2
**1 Gyr with the fix** (`p3jr`, τ 3e-6, kernel 1, `PINT_JR=1`): 0.986 Gyr in 3046 s (shared machine), and
**K̄ 3.61** (K = 3/4/5/6: 25134/36360/994/12) vs 3.78 (19544/38282/3854/820) with the matrix transformation — the
per-call energy offsets also raised the convergence floor (the k ≥ 4 updates stalled at ≈ 1e-7 AU, it. 33), so the
fix is also a speed-up. Energy at 187 Myr: p3jr +1.2e-11 (serial runs +0.5e-11, −1.0e-11; old p3 +1.4e-10).
**Ensemble** (`runs/ens2/`, 6 IC realisations of it. 33, 10 Myr, JR 1, kernel 2; shifts vs the serial run of the same
ICs; `ens_ana.py … "runs/p9/e{p}_t3_ser2.bin"`):
| setting | K̄ (6 runs) | inner rms / max ["/yr] | inner means | Saturn g6 mean ± sem | giants rms |
|---|---|---|---|---|---|
| K_min 3, τ 3e-5 | 3.10–3.15 | 3.1e-7 / 7.3e-7 | all < 1σ | −7.0e-8 ± 0.8e-8 | 2.6e-8 |
| K_min 3, τ 1e-4 | 3.01–3.04 | 5.4e-7 / 1.5e-6 | g1 +6e-8 (2σ), g2 +1e-8 (2.7σ) | −8.7e-8 ± 1.0e-8 | 3.3e-8 |
| K_min 2, τ 3e-4 | 2.57–2.63 | 1.9e-6 / 5.1e-6 | < 1σ | **−9.0e-6 ± 2e-6** | 3.4e-6 |
| K_min 2, τ 3e-4, giants 180 d | **2.53–2.60** | 1.5e-6 / 3.3e-6 | s4 −8e-7 (1.7σ) | +4.7e-7 ± 0.7e-7 | 1.8e-7 |
Energy vs serial: within ±1e-11 everywhere.
* With K_min 3 the parareal error stays unbiased in the inner planets and adds a small systematic Saturn g6 shift
  (−7…−9e-8, i.e. 5 % of the dt = 6 d g6 error).
* The K_min 2 Saturn bias comes from the coarse giants' step: at 180 d it shrinks 20× (→ 30 % of the dt = 6 d g6
  error) and K̄ even drops a little (2.56). Inner rms 1.5e-6 ≈ 0.3 × the dt = 6 d error → **"between dt 6 and 8 d"
  accuracy with K̄ ≈ 2.56** (vs 3.78 at the reference setting). Its cost (giants' run twice as long on the coarse
  chain) is measured in the timing runs.
* Same ensemble, K_min 2 + giants 180 d at other τ:
  | τ | K̄ | inner rms / max | significant means (≥ 2σ) | Saturn g6 |
  |---|---|---|---|---|
  | 1e-4 | 2.90–2.95 | 8.3e-7 / 2.0e-6 | s2 = s3 +1.2e-6 ± 0.3e-6, s4 +7.6e-7 ± 2e-7, g4 −7.2e-7 ± 3e-7 | +9e-8 ± 2e-8 |
  | 3e-4 | 2.53–2.60 | 1.5e-6 / 3.3e-6 | none resolved (s4 −8e-7, 1.7σ) | +4.7e-7 ± 0.7e-7 |
  | 1e-3 | 2.11–2.16 | 1.1e-5 / 3.0e-5 | g4 −8.9e-6 ± 4e-6, g6 +1.2e-6 ± 0.1e-6 (s2 +1.1e-5 ± 7e-6) | +1.2e-6 |
  Windows that stop at k = 2 leave **systematic** inner shifts of ≈ 1e-6 ″/yr in s2/s4 (the size of going from
  dt 6 d to dt 8 d: s2 +1.1e-6 → +2.1e-6) — a legitimate accuracy class, but not the reference accuracy.
* **Accuracy classes** (errors relative to dt 3 d; dt 6 d reference ≈ 5e-6 max, systematic):
  (A) reference — K_min 3, τ 3e-6 … 3e-5: adds ≤ 7e-7 rms, unbiased inner planets (K̄ 3.6–3.1);
  (B) ≈ dt 6–8 d — K_min 3, τ 1e-4 (K̄ 3.02) or K_min 2 + giants 180 d, τ 1e-4 … 3e-4 (K̄ 2.9–2.56);
  (C) ≈ dt 12 d — K_min 2, τ 1e-3 (K̄ 2.1). Speeds: iteration 39.
* Tried: Uranus/Neptune paired with the inner planets in the coarse step (`PINT_MAXJ=7`, it. 22 found them missing
  from the inner secular rates; bitwise serial = parallel checked) for K_min 2 + giants 180 d: τ 1e-4 K̄ 2.92,
  s2 +9.7e-7 ± 4e-7, s4 +6.2e-7 ± 3e-7 (without: +1.2e-6 ± 3e-7, +7.6e-7 ± 2e-7); τ 3e-4 K̄ 2.55, inner rms
  3.1e-6. → The k = 2 shifts are a convergence-level effect, not this missing term. Not kept.

### Iteration 39 — clean timings of the accuracy classes; the 1 Gyr runs
**Timings** (idle host2, `runs/timing.sh`, `timing2.sh`; 100 windows = 1.58 Myr, `PINT_JR=1`, 8 helpers, mode 2;
speedup vs a serial run with the same kernel; two repetitions each):
| setting | class (it. 38) | K̄ | asm (kernel 1) | C (kernel 2) |
|---|---|---|---|---|
| K_min 3, τ 3e-6 | A reference | 3.81 | 2.89–2.91 (3.60 s) | 2.91–3.02 |
| K_min 3, τ 3e-5 | A reference | 3.18 | **3.44–3.48** (3.03 s) | 3.44 |
| K_min 3, τ 1e-4 | B dt 6–8 d | 3.05 | 3.53–3.60 | 3.51 |
| K_min 2, τ 3e-4, giants 180 d | B dt 6–8 d | 2.66 | **3.74–3.80** | 3.72 |
| K_min 2, τ 1e-4, giants 180 d | B | 3.00 | 3.42–3.43 | 3.35–3.39 |
| K_min 3, τ 3e-5, giants 180 d | A | 3.18 | — | 3.17–3.24 |
| K_min 2, τ 3e-4 (giants 360 d) | Saturn-biased | 2.68 | 3.96–4.03 | 4.05 |
| K_min 2, τ 1e-3 | C dt 12 d | 2.16 | **4.69–4.74** | 4.75 |
| K_min 2, τ 1e-3, giants 180 d | C | 2.14 | 4.43–4.45 | 4.38–4.47 |
| τ 3e-6, `PINT_JR=0` / `1` | | 3.80 / 3.81 | | 2.88–2.89 / 2.91–2.95 |
| τ 3e-6, `PINT_CHUNKS=2` / `1` | | 3.81 | | 3.03 / 2.78 |
* The serial asm reference runs at 6.6 ms per 986 yr: 5 Gyr ≈ 9.2 h (the paper's 9 h). Parareal (asm): class A
  τ 3e-5 ≈ **2.7 h**, τ 3e-6 ≈ 3.2 h; class B ≈ 2.4 h; class C ≈ 1.9 h. A serial run at the accuracy of class B
  (dt 8 d) would take 6.9 h, of class C (dt 12 d) 4.6 h → parareal gains 2.4–2.9× at equal accuracy.
* `jacobi_recursive` costs no speed. Two chunks per speculative fine slice are as fast as ten (one is slower: no
  early abandonment) → default `PINT_CHUNKS` 2.
**1 Gyr runs** (0.986 Gyr, 1e6 boundaries, asm kernel; `runs/gyr/`, `gyr_ana.py`): dE/E at the end — serial
`s0` +2.9e-11, `s12` (IC + 1e-12 AU) +6.8e-11 (nearly linear: one call per slice × 7e-17 with the matrices), serial
with the fix `s0jr` +3.5e-11; parareal τ 3e-6 with the matrices `p3` **+7.4e-10** (linear), with the fix `p3jr`
−8.8e-11 (wanders to −8.8e-11 after ≈ 400 Myr ≈ 1e7 calls × −0.9e-17, the small residual seen with corrector + GR in
it. 36 at 1.8σ) → the fix removes ≈ 90 % of the drift; with 2 chunks per slice (5× fewer calls) the remainder should
be at the level of the serial runs (run `p30c2`, τ 3e-5, started). Serial wall times on shared cores: 9707 s / 9066 s.
* Secular-frequency statistics over 49 windows of 20 Myr: the KS test rejects even between the two serial runs
  (p = 0.01 for g2, g3) — consecutive windows of one chaotic trajectory are strongly correlated (and NAFF switches
  between s2 and s3), so single Gyr realisations cannot test frequency or instability statistics; an ensemble is
  needed (as in Rein et al. 2025). Mercury's maximum eccentricity: 0.303–0.345 among the five runs.

### Iteration 40 — the recommended setting over 1 Gyr; CCD load balancing (rejected)
**`p30c2`** (K_min 3, τ 3e-5, asm, `PINT_JR=1`, 2 chunks; idle machine): **0.986 Gyr in 2174 s (36.2 min)**, K̄ 3.33
(K = 3/4/5/6: 45320/14419/2341/420) — higher than in the 100-window benchmark (3.18), so the long-run speedup is
6570 s (serial asm, 6.57 ms per slice) / 2174 s = **3.0×**; 5 Gyr ≈ 3.1 h vs 9.1 h serial. Energy (dE/E vs own start,
linear-fit slope per Gyr): serial s0 +4.0e-11, s12 +6.4e-11, s0jr +2.9e-11; parareal p3 (matrices, 10 calls/slice)
+7.5e-10, p3jr (fix, 10 calls) −1.1e-10, **p30c2 (fix, 2 calls) −6.7e-11** (end −5.5e-11, max 1.05e-10) → at the level
of the serial runs. Mercury e max 0.343 (serial runs 0.303–0.333), a within 0.387095–0.387104 AU in all runs.
**CCD load balancing** (`PINT_BAL=d`: slices on cores 0–7 get (1 + d) × the steps, cores 8–15 (1 − d) ×; 100 windows,
asm, τ 3e-6): d = 0: K̄ 3.81, 2.82–2.93×; d = 0.02 / 0.03 / 0.04: K̄ 4.02 / 4.06 / 4.07, 2.74–2.80× — unequal slices
cost more iterations than the balance gains. Rejected (option kept, off).
