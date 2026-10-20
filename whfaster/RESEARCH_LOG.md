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

