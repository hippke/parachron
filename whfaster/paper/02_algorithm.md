# Part 2 — Algorithm: the pipelined parareal scheme, its schedule and its cost

This part describes the parallel algorithm built around the mean-element coarse propagator of Part 1: windowing,
pipelining, the adaptive stopping rule, the parallel coarse step, and a quantitative account of where the wall-clock
time goes on 16 cores. Numbers are from `../RESEARCH_LOG.md` (iteration numbers “it.”), measured on an AMD Ryzen 9
7950X (16 cores / 32 hardware threads, AVX-512) unless stated otherwise.

---

## 1. Windowed parareal
* The trajectory is processed in consecutive **windows** of N slices of T = 986 yr (360 000 d = 60 000 fine steps);
  N = number of cores (16). Each window starts from the converged end state of the previous one.
* Within a window, iteration k: boundary U^k_{n+1} = G(U^k_n) + F(U^{k−1}_n) − G(U^{k−1}_n), formed in osculating Jacobi
  elements with wrapped longitudes, followed by the energy projection (Part 1 §4.2). Boundaries n ≤ k − 1 are exact
  copies of fine results; G is not evaluated for them.
* Why windows, and not one long parareal over the whole run: the parareal error grows binomially with the distance
  from the last exact boundary, so long windows need more iterations (Part 1 §4.3). A continuously sliding window was
  simulated with the real numerics (event-driven simulators with the actual F and G, round 3 it. 25): S = 1.4–2.2 at
  10 cores, below the windowed scheme, because every change of a slice's start state after its iterations have begun
  costs iterations.

## 2. Pipelining
One thread owns one slice. Thread n waits until U^k_n is available, computes G(U^k_n), forms and hands on U^k_{n+1},
and then immediately starts the fine solve F(U^k_n) that iteration k + 1 will need. The coarse sweeps of successive
iterations therefore follow each other through the window like a wave, and the fine solves of different iterations
overlap. With t(k, n) the time at which U^k_{n+1} is ready,

&nbsp;&nbsp;&nbsp;&nbsp;t(k, n+1) = max( t(k, n), t(k−1, n) + T_F(n) ) + T_hop,

so one window costs ≈ N T_hop (the k = 0 coarse sweep, inherently serial) + K × (slowest fine slice + one hop), and

&nbsp;&nbsp;&nbsp;&nbsp;S ≈ N T_F / (N T_G + K (T_F + T_G)),&nbsp;&nbsp;&nbsp;r = T_G / T_F.

The bound S ≤ min(N/K, 1/r) holds for any pipelined parareal (round 3 it. 25).

## 3. Adaptive stopping rule and speculative fine solves (round 2 it. 15; this work it. 29, 33, 37)
* A window stops at the first iteration k ≥ K_min (K_min = 3) whose **largest boundary update** max_n |x(U^k_n) −
  x(U^{k−1}_n)| (planet positions, AU) is ≤ τ, at the latest at K_max = 6. The update at iteration k is a good
  estimate of the error of iterate k − 1; the accepted iterate k is one contraction (≈ 10⁻¹–10⁻²) better.
* Default τ = 3 × 10⁻⁶ AU (validated in round 2 and at N = 16 here). Looser τ were validated by ensembles (Part 1 §6)
  and define faster accuracy classes (§6).
* Because the stop is decided only when the last slice has its update, every thread starts the fine solve for the
  next iteration speculatively; it runs as a few kernel calls (default 2) and is abandoned between calls when the
  window has stopped.

## 4. The parallel coarse step
The coarse chain (N hops per window) is serial, so the latency of one coarse call matters, not its CPU time.
* **Task decomposition (round 3 it. 27).** One coarse call = start tables (Fourier grids of each inner planet, 14 pair
  tasks: 6 inner–inner on 32², 8 inner–giant on 16², 4 GR terms) + predictor drift + the giants' integration (1000 steps
  of 360 d) + end tables + Heun corrector. Tasks of one stage are independent and are reduced in a fixed order, so
  the result is **bitwise identical** for any number of helper threads. Helpers sleep between calls; 8 helpers.
  Measured 2.37× → 2.71× on 16 cores.
* **Staged schedule (this work, it. 30).** The giants' integration is the longest single task (≈ 0.15 ms). It ran on
  a helper, i.e. on the SMT sibling of a core busy with a fine solve (≈ 1.7× slower). Now the **calling thread
  integrates the giants on its own core**, while the helpers run one staged job: start grids → start pairs/GR →
  reduction and predictor → end-table grids of the inner planets → end-table inner–inner pairs and GR (these need only
  the inner planets' predicted elements and therefore overlap the giants' run); after the giants: their grids and the
  inner–giant pairs. Coarse call 0.219 → 0.165 ms idle; speed-up **2.71 → 2.85×**, bitwise identical trajectory.
* **Exact derivatives (it. 43).** ∂E/∂v by forward-mode dual numbers instead of finite differences (which were only
  good to ≈ 10⁻⁶ in some components because of round-off): same iteration counts, serial coarse call −12 %, a few %
  overall at looser τ. Default.

## 5. Persistent threads and window hand-over (it. 31)
Creating and joining 16 threads per window also waited for threads finishing an abandoned speculative chunk
(≈ 0.7 ms tail + 0.2 ms gap per window). The slice threads now persist over all windows; the window state
(boundary arrays, semaphores, stop flag, update maxima) is double-buffered by window parity. At the decision the
last slice prepares the other parity, releases slice 0 for the next window, and posts a barrier that prevents a thread
that finished its last iteration early from touching re-initialised state. Speed-up **2.85 → 2.95–3.05×**, bitwise
identical trajectory.

## 6. Where the time goes, and the resulting speed (it. 29, 31, 39, 40, 44, 46)
**Event trace of one window** (16 × 986 yr, τ = 3e-6, after it. 31; mean of 100 windows), 37.7 ms:

| part | ms |
|---|---|
| k = 0 coarse sweep, 16 hops × 0.26 ms | 4.3 |
| K̄ = 3.8 iterations × (slowest fine slice ≈ 8.5–8.7 ms + one hop) | ≈ 33 |
| tail and hand-over | < 1 |

Facts that set these numbers: a fine slice takes 7.0 ms on an idle machine but **8.3 ms when all 16 cores run** (clock
reduction under all-core AVX-512 load; the serial reference runs at single-core boost); every iteration waits for the
slowest of the 16 fine slices, and the second CCD is 3–5 % slower; the coarse hop costs 0.26 ms on the chain (0.165 ms
idle) because its helpers share cores with fine solves.

**Speed-ups** (assembly kernel, 100-window benchmarks = 1.58 Myr; long runs ≈ 10 % lower because K̄ is somewhat
higher on average):

| setting | accuracy class | K̄ | speed-up | 5 Gyr wall-clock |
|---|---|---|---|---|
| serial, dt = 6 d | reference | — | 1 | 9.1 h |
| **K_min 3, τ 3e-6 (default)** | reference | 3.81 | **2.86–2.95×** | **≈ 3.2 h** |
| K_min 3, τ 3e-5 | reference (+≤ 14 % frequency error) | 3.18 (1 Gyr: 3.33) | 3.32–3.56× (1 Gyr: 3.0×) | 2.7–3.1 h |
| K_min 3, τ 1e-4 | between dt 6 and 8 d | 3.05 | 3.5–3.6× | 2.6 h |
| K_min 2, τ 3e-4, giants' coarse step 180 d | between dt 6 and 8 d | 2.66 | 3.7–3.8× | 2.4 h |
| K_min 2, τ 1e-3 | ≈ dt 12 d | 2.16 | 4.7× | 1.9 h |

At equal accuracy the fair competitor is a serial run with a larger step (dt = 8 d: 0.75 of the cost; 12 d: 0.5);
parareal is 2.4–2.9× faster than that. The ensemble of 16 × 200 Myr (Part 1 §6.3) ran at 405 s per member vs
1330 s serial (3.3×, τ = 3e-5).

**History of the speed-up on 16 cores (default accuracy):** coarse WH integrator < 1× → mean-element coarse step with
pipelined driver 2.37× → parallel coarse step 2.71× → staged coarse schedule 2.85× → persistent threads 2.95–3.0×.

## 7. Alternatives measured and rejected (it. 28–45)
| idea | result |
|---|---|
| fewer slices, cores left for coarse helpers (N = 12, 14) | 2.39–2.50× (fewer slices cost more than they save) |
| SMT: 32 half-length slices on 32 hardware threads | 2.0–2.1× (each fine solve 1.68× slower, total throughput only +19 %) |
| N = 20, 24 (shorter slices) | 2.45–2.62× |
| slice length 789 / 1183 / 1380 yr | 2.85 / 2.80 / 2.70× (986 yr optimal) |
| two slices per core in one paired kernel call | no gain on Zen 4 (paired throughput ×1.02; ×1.23 on Intel AVX-512) |
| mean-n tables, 64² Venus–Earth grid in config B | K̄ −10 % but the 64² task lengthens the coarse chain: 2.0–2.7× |
| giants' coarse step 180 / 400–1000 d | no gain / K̄ rises smoothly to non-convergence |
| per-CCD load balancing by unequal slice lengths | K̄ 3.81 → 4.02–4.07, slower |
| larger τ alone with K_min 3 | K̄ floors at 3.0 (τ ≥ 1e-4) |
| accurate coarse step only for the k = 0 prediction | helps the hardest windows only; no K = 2 |

## 8. Can K = 2 at reference accuracy pay? (it. 41–45, 48)
The accurate coarse step of Part 1 §5.2 gives K̄ ≈ 2.4 at N = 16, but its cost enters the serial coarse chain. With the
cost model of §2, one window costs ≈ 16 T_G + K̄ (T_F + T_G):
* break-even against config B requires **T_G ≲ 0.35 ms on the chain under load**;
* the all-component transform was reformulated so that one FFT per planet pair serves all six elements (the
  displacement Σ_k f_k e^{ik·θ}/ν_k equals a grid sum of f with the inverse-transformed kernel e^{ik·θ}/ν_k; exact,
  2.2× faster), with exact dual-number derivatives (grids 2.6× faster) and a task-parallel version: 0.068 ms per pass
  idle but **0.140 ms under the pipeline's load**;
* the minimal accurate step needs ≈ 5–9 dependent passes → 0.5–1.0 ms → −11 % … +3 % against config B on 16 cores;
* more cores do not change this by themselves: the coarse chain grows with N exactly like config B's iteration count
  (model at N = 40: config B ≈ 5.5×, accurate step ≈ 4.7×).

The accurate coarse step therefore pays only if a hop can be brought to ≈ 0.35 ms, e.g. by fusing the passes or
reusing tables across RK stages.

## 9. Scope and limitations
* One machine type (16-core desktop CPU); multi-node runs were not attempted.
* The coarse propagator is specific to the Solar System's architecture (four inner planets treated by perturbation
  theory, four well-separated giants integrated directly); compact or resonant systems would need a different split
  and re-validation.
* A small residual energy drift of the parareal trajectory (−1.25e-10 per Gyr relative to serial) is not yet explained
  (Part 1 §6.4).
