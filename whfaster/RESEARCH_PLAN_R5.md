# Round 5 research plan — "config C": K = 2 at reference accuracy on 16 cores (Oct 22, 2026)

## Where we are (iterations 28–41)
* Config B on the 7950X (16 cores): reference accuracy (τ 3e-5, K_min 3) 3.46× in the 100-window benchmark, 3.0× over
  1 Gyr (K̄ 3.33). The scheme is at ≈ 90 % of its ceiling S ≈ 16/K̄ × 0.84 − coarse chain; K̄ cannot drop below 3
  with config B at reference accuracy (the first corrected iterate is 2e-4 AU off, the k = 2 update never ≤ τ).
* Iteration 41: the round-3 Python coarse step makes **6/10 windows of 16 slices stop at K = 2** at τ = 3e-5 (the rest
  at K = 3) → K̄ ≈ 2.4 with *better* accepted states than config B at K = 3–4. Ablation → minimal ingredient set
  ("config C"): mean-n tables, all-component first-order transform, Earth–Mars–Jupiter resonant term frozen per
  slice, giants' WH at 180 d sampled at slice start/middle/end, one RK4 step of the secular drift using the sampled
  giants, Uranus/Neptune in the rates. Not needed: second-order secular rates U2, giants' box averaging, 2nd RK4 step.

## Payoff and cost target
Window ≈ 16 T_G + K̄ (T_F + T_G), T_F ≈ 8.6 ms (loaded), config B T_G ≈ 0.26 ms → 31.5 ms (K̄ 3.18).
Config C with K̄ 2.4: T_G 0.3 ms → 26 ms (+20 %, ≈ 4.1×); 0.45 ms → 29 ms (+9 %); break-even ≈ 0.6 ms.
→ **Gate: T_G ≤ 0.35 ms on the critical path under load** (parallel staged schedule, 8 helpers), with K̄ ≤ 2.5.

## Steps
1. **C cost model first (1 day).** Port the pieces as separate C functions, each validated against the Python
   prototype to ≤ 1e-10 (relative) on the 10 epochs: (a) all-component transform (start from `pint_disp6.c`, which
   already exists, 2.13 ms per pass): analytic dE/dv (Gauss equations) instead of finite differences, packed real FFTs,
   pruning of coefficients < 1e-12; (b) secular rates only (grid means, no FFT) for the RK4 stages; (c) giants at
   180 d with the three samples (on the calling core, as now). Measure each serially and in the staged schedule.
   Stop if the sum cannot reach the gate.
   **Status (iteration 42):** (a) done in part — `pint_disp6k` (one kernel per pair, exact) 0.286 ms per pass on the
   7950X (was 0.640); estimate for the whole config C step ≈ 0.5 ms on the chain → still needs: dual-number/analytic
   dE/dv (0.100 ms of the pass), one packed FFT for both kernels, inner–inner rates frozen across RK4 stages.
   **Iteration 43:** dual-number dE/dv done (`pint_dual.h`, exact; grids 0.095 → 0.036 ms). Estimate for the whole
   step ≈ 2 ms serial / ≈ 0.55 ms on the chain as a straight port; ≈ 0.3 ms with inner–inner rates frozen across the
   RK4 stages, a divisor-only mean-n re-evaluation and one end evaluation (to be checked against the Python
   prototype's accuracy: the end map in Python is a 3-step fixed point at the osculating angles).
   **Iteration 45 (preliminary verdict):** task-parallel pass 0.09 ms idle (local) → config C chain ≈ 0.4–0.6 ms under
   load → about break-even on 16 cores. Pursue config C only together with pass fusion / task splitting, or for
   machines with more cores (≥ 24–32), where config B needs K ≥ 4.
2. **E–M–J term** (`so2.resonant_U2`, `g3.G3`): compute once per slice at k = 0 on a helper (off the chain; iteration
   24 showed per-slice freezing is as good as live), apply analytically in the drift.
3. **Config C in `pint_coarse.c`** behind a cfg flag (default stays config B), bitwise serial = parallel
   (`test_par`). Measure K̄ and speed with `pint_pipe` (`runs/timing.sh` pattern), then the accuracy contract:
   10 Myr ensemble (`ens_ana.py`, 6–7 ICs) against the dt-reference table (`tau_ana.py`), and a 1 Gyr run
   (`gyr_ana.py`, energy). Keep K_min 2 only if the k = 2 stops show no systematic shift (iteration 38 found
   ≈ 1e-6 ″/yr shifts in s2/s4 with config B — check whether config C removes them).
4. Fallbacks if the gate is missed: (a) table reuse between iterations (k ≥ 1 hops with the start tables of the
   previous iterate + first-order correction); (b) config C only for the k ≥ 1 hops would *not* work (iteration 25:
   the k = 0 prediction must be accurate too).

## Other open items
* Ensemble of Gyr runs (serial vs parareal) for secular and instability statistics — single runs cannot test them
  (iteration 39). E.g. 16 serial + 16 parareal runs of 200 Myr from perturbed ICs, compare distributions of g1, s1,
  Mercury's e_max.
* WHFast512 synchronization bias: report upstream (`WHFAST512_SYNC_BIAS.md`, user decision); the remaining
  −0.9e-17 per call with corrector + GR (round-off asymmetry) is now below the serial random walk.
