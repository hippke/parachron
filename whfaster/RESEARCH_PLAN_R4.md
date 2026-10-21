# Round 4 research plan — make the accurate coarse step cheap (draft, Oct 21, 2026)

Round 3 (`RESEARCH_LOG.md` iterations 17–26) showed *what* the coarse step must contain for K = 2 (all-component
first-order transform at mean mean motions, E–M–J resonant term, second-order secular rates, giants at 180 d with
averaged forcing) and that a straightforward C implementation costs r ≈ 0.2–0.4, which cancels the gain at ≈ 10
cores (S = N/(N r + K̄(1 + r)); break-even r ≈ 0.12 for K̄ = 2.4). Round 4 attacks the cost.

## Ideas (ranked)
1. **Tables that do not have to be rebuilt per call.** The Fourier tables depend only on the slow elements, which
   change by < 1e-3 relative within a slice and by little between iterations. Options:
   (a) per-slice tables from the k = 0 evaluation, reused for k ≥ 1 together with a first-order Taylor correction
       in the slow elements (derivative tables computed once per slice, in parallel with the fine solves);
   (b) a per-window table at the window start plus Taylor corrections in the slow elements (secular drift over a
       window ≈ 1–5 %; check the Jacobian requirement found in round 2: "live" tables mattered);
   (c) literal (Laplace-coefficient × polynomial in k, h, q, p) expansions per pair, coefficients precomputed per
       window — evaluation cost only.
   Gate: per-call cost ≤ 0.5 ms with the round-3 accuracy (k2 ≤ 2 at ≥ 6/10 epochs).
2. **Cheaper evaluation**: analytic dE/dv (Gauss equations) instead of finite differences, packed real FFTs, keep only
   coefficients above 1e-12 (most of the 64² are negligible), SIMD over harmonics.
3. **Split the critical path**: only the k = 0 sweep is strictly serial (N r); later sweeps overlap with fine
   solves. Use the expensive G only where it is on the overlapped path and a cheaper but *equally predictive*
   G for k = 0 (round 3: config B for k = 0 destroys K = 2 — needs a better cheap predictor, e.g. the improved G
   with frozen tables from the previous iterate of the preceding window).
4. **Complete second-order short-period theory** (all pairs and triples, consistent Y-form) for the hard windows —
   only if 1–2 succeed (otherwise its cost is moot).
5. **Many-core measurement** (Phase 8 of round 3) — independent; needs a ≥ 10-core machine.

## First steps
* Instrument `pint_disp6.c` (grid, FFT, evaluation shares); implement analytic dE/dv and coefficient pruning;
  measure.
* Python test of idea 1(a)/(b) with the ph1eval harness (tables frozen per slice / per window + Taylor terms).
