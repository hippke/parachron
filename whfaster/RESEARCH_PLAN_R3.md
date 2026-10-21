# Round 3 research plan: a second-order mean-element coarse propagator (goal: K = 2)

*Written Oct 21, 2026, after round 2 (`PINT_SUMMARY.md`, `RESEARCH_LOG.md` iterations 1–16).*

## 1. Goal and why
Round 2 reached ≈ 2.3× wall-clock on 10 cores at research-grade accuracy (adaptive K, mean 3.3 iterations per
window, r ≈ 0.09). Typical windows reach the fine solver's noise floor at k = 3; at k = 2 they are 1e2–1e4 × above it.
**Each iteration saved is worth ~25–30 % in S**: K̄ = 2.2 would give ≈ 3.1× on 10 cores and ≈ 1.9× on 5 cores
(S = N/(N r + K̄(1 + r))). This requires a coarse map whose Jacobian and prediction error are 1–2 orders of
magnitude better for the inner planets — the first-order theory is accurate to ~1e-3 of the short-period terms; the
remaining error is second order in the masses plus one identified near-resonance (Earth–Mars–Jupiter,
4λ_E − 8λ_M + 3λ_J + ϖ, P = 1836 yr, 2.2e-7 in Mars' mean ln a).

Success criterion: k = 2 error ≤ 2 × noise for ≥ 90 % of windows (10 test epochs over 10 Myr, then 1000 windows),
with r ≤ 0.15, and the 10 Myr frequency/energy contract (iteration 12/14 protocol) met.

## 2. Phases, deliverables, go/no-go gates

### Phase 0 — requirements (oracle, linear model)
* Slice Jacobians J_F, J_G (current config B, C coarse) at the 10 test epochs; e0 = actual coarse errors.
* Model a better coarse map as J_G(ε) = J_F + ε (J_G − J_F) for the inner-planet rows (all six elements or
  subsets) and e0(ε0) = ε0 e0; find (ε, ε0) for which the linearised k = 2 error ≤ 2 × noise at every epoch.
* Deliverable: required improvement factors per block (ln a, λ, e/i rows; inner vs giants).
* **Gate 0:** if K = 2 needs ε < 1e-3 (i.e. better than a second-order theory can plausibly deliver), stop the theory
  track and go to Phase 7.

### Phase 1 — complete first-order theory
* Short-period transform for all six elements (λ with the shear term, k, h, q, p), live, consistent with the ln a
  transform (Python prototype in `meanel.py` already has them).
* Close-pair grid refinement (V–E, E–M at 64²) as an option.
* Second-order-consistent evaluation of the first-order displacement at mean (not osculating) elements
  (fixed-point inverse), which removes an O(ε²) inconsistency.
* Measure: residual spectrum (`resid_spec.py`), parareal error at k = 1, 2, 3 on the 10 epochs.
* **Gate 1:** quantify the gain; continue regardless (Phase 2 builds on it), but record ε achieved.

### Phase 2 — second-order averaging, two-planet part
* Use the averaging method for dz/dt = ε f1(z, θ), θ̇ = ω(z) (non-canonical, matches the Gauss-equation
  formulation of round 2): near-identity z = ζ + ε w1(ζ, θ) + ε² w2(ζ, θ); second-order averaged equations
  ζ̇ = ε⟨f1⟩ + ε²⟨Df1·w1 − Dw1·⟨f1⟩ (+ θ-derivative terms from the mean-motion dependence)⟩.
* Restrict to products of terms of the same pair (2-D tori): the second-order mean-motion corrections, the
  second-order secular rates and the w2 of each pair.
* Derivatives of f1 and w1 with respect to the slow elements and Λ by finite differences of the grid tables.
* Measure as in Phase 1; check that the constant part of Mars'/Earth's mean-motion error (the bulk of e0) is removed.

### Phase 3 — three-planet second-order terms
* Products of terms of different pairs that share a planet (i, j, l): functions on 3-D tori (λ_i, λ_j, λ_l).
  Select triplets by size: (E, M, J), (V, E, J), (M, J, S), (E, M, S), (V, E, M), (Mer, V, J) — confirm with the
  residual spectrum.
* Small divisors: combinations with |k·ν| < 2π / (2 slices) are **not** transformed away; they stay in the mean
  equations as resonant terms (amplitude function of the slow elements) and are integrated by the coarse step
  (their phase advances by up to ~π per slice). Frequencies from second-order mean motions (Phase 2).
* Deliverable: Mars' 1836-yr term reproduced in amplitude and phase (validation against the 10 kyr residual run).

### Phase 4 — the new coarse propagator (Python)
* G2 = second-order transform → second-order mean dynamics (secular + resonant, Heun or RK2) → inverse transform;
  giants unchanged (WH at 360 d).
* Measure k = 1…5 errors at the 10 epochs; adaptive K̄ with the round-2 stopping rule; compare with config B.
* **Gate 4:** k = 2 ≤ 2 × noise at ≥ 8 of 10 epochs → Phase 5. Otherwise analyse the remaining mismatch
  (`jac2.py`/`oracle2.py`), iterate once on Phases 2–3, then decide.

### Phase 5 — cost reduction and C port
* Second-order coefficients change slowly: compute them once per window (or per few windows), in parallel, and keep
  only the first-order part live; check that this keeps the convergence (live vs frozen test as in iteration 4).
* C implementation in `pint_coarse.c` (new cfg flags), bitwise-independent of the fine kernels; target r ≤ 0.15.

### Phase 6 — validation (accuracy contract)
* 10 Myr runs on both IC sets with the adaptive rule (and fixed K = 2 where it converges): secular frequencies vs serial,
  asm and IC-perturbation references (iteration 12/14 protocol); energy behaviour; boundary-jump statistics.
* Optional: 50 Myr frequency analysis; energy bias ensemble of 16 × 1e9 steps.

### Phase 7 — independent lever: overlapping/sliding windows (engineering)
* Start window w+1's coarse sweep from window w's provisional end state and refine as it converges; removes the
  N r start-up term and the idle cores of the last iterations. Bound S ≤ N/(K̄(1 + r)): ≈ 2.8× at K̄ = 3.3,
  ≈ 4.1× at K̄ = 2.2 on 10 cores. Can be done independently of Phases 1–6.

### Phase 8 — measurement on a many-core machine
`pint_pipe` on ≥ 10 physical cores (Zen 4/5 EPYC preferred) — needs hardware outside this environment.

## 3. Tools available (from round 2)
`parareal/ptc.py` (C coarse via ctypes), `jac2.py`/`oracle2.py` (Jacobians, oracles), `sweep2.py` (10-epoch
convergence and S), `stoprule.py`, `resid_sample.py`/`resid_spec.py` (residual spectra), `meanel.py` (Python
first-order theory incl. all components), `pint_long` + `freq_ana.py` (10 Myr validation).

## 4. Risks
* Second-order theory may still leave ε ≈ 1e-2 because of higher-order near-resonances (e.g. other three-body
  combinations, Venus–Earth 8:13) → Gate 0/4 catch this early.
* Cost: 3-D tori are expensive; mitigated by freezing second-order parts per window (Phase 5).
* Mercury (e = 0.2) may need more harmonics; GR × planet cross terms are negligible (≲ 1e-8 relative).
