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

