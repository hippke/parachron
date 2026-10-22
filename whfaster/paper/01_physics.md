# Part 1 — Physics: making parareal work for the long-term integration of the Solar System

This part explains *why* standard parallel-in-time integration fails for planetary systems and *which physical
insight* makes it work: performing the parareal correction with a coarse propagator built from a mean-element
(averaged) perturbation theory, so that the coarse and fine maps agree in the quantities that Kepler shear amplifies.
All numbers come from the lab book `../RESEARCH_LOG.md`; the iteration numbers (“it.”) refer to it.

Notation. Units G = 1, lengths in AU, time in yr/2π (REBOUND convention). Osculating Jacobi elements per planet:
**E = (ln a, λ, k, h, q, p)** with k + ih = e e^{iϖ} and q + ip = sin(i/2) e^{iΩ}; λ is the mean longitude.
Mean motion n = √(μ/a³). “Slice” = one parareal time sub-interval; “window” = N consecutive slices solved together.

---

## 1. The integration problem

* **System and fine propagator F.** Sun + 8 planets, symplectic Wisdom–Holman integrator in Jacobi coordinates
  (WHFast512, arXiv 2609.25339), time step dt = 6 d, 17th-order symplectic corrector, the GR potential term
  −3 m₀² m_i /(c² r_i²) of the kernel (heliocentric r_i, c = 10065.32 in these units). One step costs ≈ 100–120 ns on a
  modern AVX-512 core; 5 Gyr = 3.0 × 10¹¹ steps ≈ 9 h serially (paper: 9 h on Zen 5; this work: 9.1 h with the assembly
  kernel on an AMD Ryzen 9 7950X).
* **Why one trajectory cannot simply be split across cores.** Each step depends on the previous one. Splitting one
  step across cores cannot win: a fork/join costs ~1 µs, ten times a step (round 2 literature review,
  `LITERATURE_PINT.md`). The only available parallelism is *in time*.
* **Goal.** Reduce the wall-clock time of *one* trajectory on a many-core node at research-grade accuracy, defined
  before any code change (the “accuracy contract”, `RESEARCH_LOG.md` header): (1) unbiased energy, (2) energy error not
  worse than the serial integrator, (3) secular frequencies correct to within the serial integrator's own error at
  dt = 6 d, (4) unchanged instability statistics.

## 2. Parareal in one paragraph

The time span is cut into N slices of length T. A cheap coarse propagator G and the accurate fine propagator F (the
serial integrator applied to one slice) are combined in the iteration

&nbsp;&nbsp;&nbsp;&nbsp;U^{k+1}_{n+1} = G(U^{k+1}_n) + F(U^k_n) − G(U^k_n),&nbsp;&nbsp;&nbsp;U^k_0 = U_0,

where U^k_n is the state at the start of slice n in iteration k. All N fine solves F(U^k_n) of one iteration are
independent and run in parallel; only the cheap G is sequential. After k iterations the first k boundaries are exact
(they equal the serial fine chain), and if the iteration contracts fast, all boundaries reach the fine solver's own
accuracy after a few iterations. With pipelining (Part 2) the wall-clock speed-up of one window is

&nbsp;&nbsp;&nbsp;&nbsp;S ≈ N T_F / (N T_G + K (T_F + T_G)) ≤ N / K,

so what matters physically is the number of iterations K that G needs. Linearising around the fine trajectory gives
the error recursion e^{k+1}_{n+1} = G′ e^{k+1}_n + (F′ − G′) e^k_n: **convergence is governed by how well the coarse
Jacobian G′ reproduces the fine Jacobian F′ along the trajectory.** This was verified quantitatively: the linear
recursion with the measured slice Jacobians reproduces the nonlinear parareal iteration step by step (round 2 it. 1;
e.g. Mercury λ error 4e-2 → 2e-2 → 9e-4 → 8e-5 → 1e-6 measured vs 4e-2 → 2e-2 → 4e-4 → 6e-5 → 6e-8 predicted).

## 3. Why parareal fails for planetary orbits with a conventional coarse propagator

### 3.1 Measured failure (round 1)
The natural coarse propagator is the same symplectic integrator with a larger step. It must still resolve Mercury's
88-day orbit (dt_G ≲ 12–16 d), so its cost ratio r = T_G/T_F is ≥ 0.4, and even then it converges slowly:

| slices × length | dt_G | r | iterations to the fine solver's noise level |
|---|---|---|---|
| 10 × 986 yr | 24 d | 0.26 | no contraction |
| 10 × 986 yr | 16 d | 0.42 | 4 (9 of 11 boundaries) |
| 10 × 986 yr | 8 d | 0.77 | 3 |
| 10 × 9 856 yr | 16 d | 0.43 | no contraction |

Best modelled speed-up 0.3–0.5, i.e. a slowdown. Correcting in Cartesian coordinates is worse: it **diverges**
(errors up to 10⁵ AU), because the coarse errors are phase errors of 0.1–1 rad and additive Cartesian corrections
across a phase error produce unphysical states. Correcting in orbital elements with wrapped angles is bounded.

### 3.2 The mechanism: short-period terms × Kepler shear
The slice Jacobians (48 × 48, in the elements E) locate the mismatch (round 2 it. 1):

* **Kepler shear.** A relative error δ(ln a) in a planet's semi-major axis changes its mean motion by
  δn/n = −(3/2) δ(ln a), so its phase error grows linearly in time:

  &nbsp;&nbsp;&nbsp;&nbsp;δλ(T) = −(3/2) n T δ(ln a).

  Over one 986-yr slice the amplification 1.5 n T is ≈ 3.9 × 10⁴ for Mercury, 1.5 × 10⁴ for Venus,
  9.3 × 10³ for Earth, 4.9 × 10³ for Mars. A semi-major-axis error of 10⁻⁸ becomes a phase error of ~10⁻⁴ rad.
* **Phase dependence of the osculating semi-major axis.** The osculating a contains short-period terms with periods
  of the orbital and synodic periods; their amplitudes in ln a are 2e-6 (Mercury) to 3e-5 (Mars) for the inner
  planets (it. 2, table below). Hence the fine Jacobian has rows ∂(ln a_end)/∂λ_j ≈ 10⁻⁶…10⁻⁵ and, through shear,
  ∂(ln a_end)/∂(ln a_j) up to 0.09 (e.g. Earth ← Venus).
* **What the coarse model must get right.** An oracle experiment (coarse Jacobian = exact + controlled noise) shows
  that K ≤ 3 requires the **ln a rows accurate to ≲ 10⁻⁸ absolute** and the remaining entries to ≲ 1 %. A coarse WH
  integrator with a large step aliases the short-period terms and gets these rows wrong by 1–5 % (dt_G = 16 d) or
  ~80 % (24 d). No cheap model of the *osculating* elements can meet 10⁻⁸ in these rows: that is 10⁻³ relative to the
  short-period amplitudes and requires the correct phases at the slice end.

**Physical conclusion.** The quantity that Kepler shear amplifies, the osculating semi-major axis, is dominated by
fast, phase-dependent oscillations that a coarse model cannot reproduce cheaply. A coarse model in *mean* elements,
from which these oscillations are removed analytically, has ln a rows that are close to the identity, which is
exactly what the convergence of parareal needs.

## 4. The mean-element coarse propagator

### 4.1 Numerical first-order perturbation theory (round 2 it. 2)
For every planet pair (i < j) the first-order perturbing accelerations in Jacobi coordinates are evaluated on a
regular grid of the two mean longitudes (λ_i, λ_j), with the slow elements (a, e and i vectors) frozen at their
current values:
* interior planet perturbed by the exterior one: direct term plus the dipole indirect term;
* exterior planet perturbed by the interior one: the heliocentric-like Jacobi form (the planet sees the interior
  mass distribution);
* checked against a direct N-body sum: relative agreement 3–5 × 10⁻⁵ for the inner planets, 1–2 × 10⁻³ for the
  giants (the difference is second order in the masses).

Gauss' equations convert the accelerations into element rates, dE/dt = (∂E/∂v)·a_pert (for ln a analytically:
d ln a/dt = 2a v·a_pert/μ). The 2-D Fourier transform of each rate gives its harmonics f_k with wave vector
k = (k₁, k₂) and frequency ν_k = k₁ n_i + k₂ n_j. The **short-period displacement** of each element is

&nbsp;&nbsp;&nbsp;&nbsp;Δ E = Σ_{k ≠ 0} f_k e^{i k·θ} / (i ν_k)&nbsp;&nbsp;&nbsp;(θ = current mean longitudes),

plus, for λ, the shear term −(3/2) n Σ (ln a)_k e^{i k·θ}/(i ν_k)², and the **k = 0 term** of each rate is the
first-order **secular rate**. GR enters as a one-planet 1-D term. Mean elements = osculating − Σ_pairs Δ E. In effect
this is a first-order Lie-series (von Zeipel) transformation evaluated numerically on grids instead of by literal
expansions in eccentricities and inclinations, so it is exact in e and i at first order in the masses.

Quality (residual short-period scatter of the elements along a 10-yr fine run, std after removing a quadratic trend):

| | Mercury | Venus | Earth | Mars | Jupiter | Saturn | Uranus | Neptune |
|---|---|---|---|---|---|---|---|---|
| osculating ln a | 2.0e-6 | 6.3e-6 | 1.0e-5 | 3.3e-5 | 2.3e-5 | 5.5e-5 | 1.9e-5 | 2.2e-6 |
| **mean ln a** | **8.6e-10** | 7.2e-9 | 9.8e-9 | 3.0e-8 | 1.2e-7 | 3.2e-7 | 6.5e-8 | 1.3e-8 |
| osculating λ | 2.8e-5 | 2.6e-5 | 2.4e-5 | 8.0e-5 | 3.9e-5 | 8.8e-5 | 2.7e-5 | 4.4e-6 |
| mean λ | 2.0e-8 | 6.4e-8 | 5.2e-8 | 9.1e-8 | 3.0e-6 | 7.5e-6 | 1.6e-7 | 1.7e-7 |

For the inner planets the scatter drops 300–2300×, to the level the shear argument requires. For the giants it does
not: the **near-resonant great inequality 2n_J − 5n_S (period ≈ 880 yr)** and **Uranus–Neptune near 2:1 (≈ 4200 yr)**
have small divisors ν_k that a first-order theory cannot represent.

### 4.2 Hybrid coarse step (round 2 it. 3–6, 12)
* **Inner planets (Mercury–Mars):** osculating → mean (transform at the slice start), drift of the mean elements over
  the slice with the Kepler mean motion n(ā) plus the first-order secular rates, mean → osculating (transform at the
  slice end).
* **Giants (Jupiter–Neptune):** integrated by the same symplectic integrator with a **360-day step** in a reduced
  system (Sun with the inner planets' masses added, plus the four giants). The giants are far from the Sun, so a large
  step resolves them, and the near-resonant terms are followed dynamically instead of averaged. The inner planets'
  secular forcing uses the giants' elements from this run.
* **Live tables (it. 4).** The Fourier tables must be rebuilt from the actual input state at every coarse call.
  Frozen tables (computed once per slice boundary) lack the dependence of the short-period terms on the eccentricity
  and inclination vectors; the dominant remaining Jacobian mismatches were ∂(Mars λ)/∂(Mars k) and
  ∂(Venus ln a)/∂(Earth k, h). Iterations to the noise level (10 × 986 yr): frozen 6, live secular rates only 4,
  **fully live 3**; contraction ≈ 10⁻² per iteration.
* **Minimal sufficient model (“config B”).** Systematic ablation (it. 5, 6, 12) leaves: transform of ln a only (λ, k,
  h, q, p used as osculating), 6 inner–inner pairs on 32 × 32 grids and 8 inner–(Jupiter, Saturn) pairs on 16 × 16
  grids, secular rates from the same grids, **trapezoidal (Heun) secular step** using the rates at the start and at the
  predicted end, giants at 360 d. (The Euler secular step gave a *systematic* frequency bias: G's error is the same
  sign in every slice, so the parareal residual after a fixed number of iterations biases every boundary the same way;
  Heun removed it, it. 12.)
* **Correction variables.** The parareal correction F(U^k_n) − G(U^k_n) is formed in the osculating Jacobi elements
  with the mean longitudes wrapped to (−π, π]; this keeps corrections physically meaningful across phase errors.
* **Energy projection.** After each correction the barycentric velocities of the corrected boundary state are scaled
  so that its energy equals E(U^k_n) + [E(F(U^{k−1}_n)) − E(U^{k−1}_n)]: each slice then carries the energy increment of
  one fine solve, as in a serial run. This leaves the parareal fixed point unchanged. Boundary energy jumps are at the
  round-off level (rms 5.8e-14 vs 4.9e-14 for two serial kernels), and the per-slice energy random walk of the
  parareal trajectory is 1.11× that of a serial run (it. 9).
* **Slice length.** 986 yr (360 000 d) is optimal: the coarse error grows with T (the slow phase dependence of the
  mean semi-major axes acts ∝ T³: mismatch × shear × initial error), while shorter slices make the fine solves too cheap
  relative to the coarse step (it. 4, 8, 15, 32).

### 4.3 Convergence obtained
* 10 × 986 yr windows: **K = 3** to the fine solver's own round-off noise (max position error 2.7e-8 AU vs noise
  3.4e-8 AU after 10 kyr). The linear Jacobian prediction matches the measured iteration count.
* 16 × 986 yr windows (16 cores, this work): mean K̄ = 3.8 at the default tolerance (§6), 3.2 at a looser one.
* Longer windows need more iterations (the parareal error grows binomially with the distance from the last exact
  boundary): 20 slices K = 4, 40 slices K = 4 (round 2 it. 8; round 3 it. 25).

## 5. What limits the convergence — physics of the remaining coarse error

### 5.1 Which planet and which element (this work, it. 33)
Per-planet and per-element maxima of the boundary updates |U^k − U^{k−1}| over 100 windows of 16 slices:

| k | worst planet (number of windows) | Venus a·Δλ | Venus a·|Δz| | Venus a·Δln a |
|---|---|---|---|---|
| 2 | Venus 45, Mars 35, Earth 11, Mercury 8 | 2.5e-4 AU | 1.2e-6 | 1.4e-8 |
| 3 | Venus 42, Mars 34, Earth 18 | 3.7e-6 | 4.2e-8 | 2.9e-10 |
| 4 | Venus 37, Mars 22, Earth 13 | 1.3e-7 | 1.3e-9 | 1.1e-11 |

The stopping decision is always set by the **mean longitude of Venus, Earth or Mars**, and that phase error is the
Kepler-shear image of a semi-major-axis error: 1.5 n_Venus T × 2.9e-10 ≈ 4.4e-6 rad, matching the measured
3.7e-6 AU / 0.72 AU. Eccentricity and inclination vectors are two orders of magnitude smaller; the giants never bind.

### 5.2 Physical terms missing from the first-order theory (round 3, it. 16–24)
1. **Small divisors computed with osculating mean motions.** The Venus–Earth 8:13 term (divisor 2π/240 yr) is so
   sensitive to the 10⁻⁵ short-period part of the osculating a that the transformed mean ln a was wrong by up to
   5.5e-7. Using *mean* mean motions in the divisors (one fixed-point iteration) reduces this to 9e-9.
2. **Grid aliasing of Venus–Earth high harmonics** (near 3:5, α = 0.72): residual comb at (5 + j) n_V − j n_E; a 64 × 64
   grid for this pair cures it.
3. **Earth–Mars–Jupiter three-body near-resonance 4λ_E − 8λ_M + 3λ_J (+ ϖ)**, second order in the masses: an
   oscillation of Mars' mean semi-major axis with period ≈ 1836 yr and amplitude ≈ 2.2e-7. Its period is ≈ 2 slices,
   so it must be followed as resonant dynamics, not averaged. A second-order averaged term computed from the
   first-order pair tables reproduces the measured amplitude to 7 % (phase to ~30°).
4. **Second-order secular rates** (dominated by the Mars–Jupiter pair): the per-slice error of the eccentricity
   vectors is exactly linear in the eccentricity vectors (Mars' g off by 0.25 %); the second-order secular term U2
   reproduces it exactly.
5. **Uranus and Neptune** were missing from the inner planets' secular rates (≈ 10⁻³ of Jupiter's effect on q, p).
6. **The giants' 60-yr term (2λ_S − λ_J)** aliases into the inner secular drift if the giants are point-sampled.

Including items 1–6 in a Python prototype makes 6 of 10 test windows (N = 10) converge in K = 2, and all at the noise
floor at K = 3 (round 3). **At N = 16 (this work, it. 41)** the same accurate coarse step brings the error after two
iterations from 1.3e-6…3.7e-4 AU (config B) down to 3.5e-8…9.5e-7 AU, so 6 of 10 windows would stop at K = 2 at
τ = 3e-5 AU (K̄ ≈ 2.4). Ablation at N = 16 identifies the indispensable physics: the giants' coarse step must be
180 d (with 360 d no window reaches K = 2), the secular drift must use the giants' sampled positions within the
slice (one RK4 step), and then the Earth–Mars–Jupiter term, the all-component transform and the mean-n divisors;
second-order secular rates, the giants' 120-yr box averaging and a second RK4 sub-step are *not* needed.
(Cost is the obstacle, Part 2 §7.)

## 6. Accuracy: what “research grade” means and how parareal compares

### 6.1 The yardstick is the integrator's own time-step error (it. 35)
Secular frequencies (NAFF of k + ih and q + ip over 10 Myr; hard initial conditions with all planets shifted by
1e-6 AU) compared with a serial run at dt = 3 d:

| run | g4 (Mars) | g6 (Saturn) | s2 (Venus) | largest error | cost |
|---|---|---|---|---|---|
| serial dt = 6 d (reference setup) | −5.1e-6 | −1.6e-6 | +1.1e-6 | 5.1e-6 ″/yr | 1 |
| serial dt = 8 d | −1.0e-5 | −3.2e-6 | +2.1e-6 | 1.0e-5 | 0.75 |
| serial dt = 12 d | −2.6e-5 | −7.9e-6 | +5.7e-6 | 2.6e-5 | 0.5 |
| parareal τ = 3e-6 AU | −5.0e-6 | −1.6e-6 | +1.0e-6 | 5.1e-6 | 1/2.9 |
| parareal τ = 3e-5 AU | −5.5e-6 | −1.7e-6 | +1.8e-6 | 5.5e-6 | 1/3.4 |

The serial integrator's own error scales as dt² (8 d vs 6 d: ratio 2.0 = (64 − 9)/(36 − 9)) and is a **systematic
bias** of the secular frequencies, up to 5e-6 ″/yr for Mars at the standard 6-day step. Rein, Brown & Kanda (2025,
arXiv 2507.04987) show that ensembles of 5 Gyr integrations remain statistically physical up to dt ≈ 32 d.

### 6.2 Parareal errors are small and unbiased (it. 29, 33, 38)
* Ensemble of 7 initial-condition realisations, N = 16, 10 Myr: secular-frequency shifts relative to the serial run of
  the same ICs, inner planets: rms 6.6e-8 ″/yr (τ = 3e-6), 1.4e-7 (τ = 1e-5), 3.1e-7 (τ = 3e-5); all ensemble means
  consistent with zero (largest 1.7σ). Energy at 10 Myr within the scatter of two serial round-off realisations.
* The parareal error therefore adds a small *random* scatter on top of the serial integrator's *systematic* error:
  ≤ 6 % of the largest dt = 6 d error at τ = 3e-6, ≤ 14 % at τ = 3e-5.
* A small systematic Saturn g6 shift appears (−7…−9e-8 ″/yr, 5 % of the dt = 6 d g6 error) and becomes large
  (−9e-6, dt = 12 d level) if windows may stop after two iterations with the giants' coarse step at 360 d; a 180-d
  coarse step for the giants removes 95 % of it. Windows stopping at k = 2 also leave ≈ 1e-6 ″/yr shifts in the inner
  inclination frequencies s2, s4 (the size of going from dt = 6 d to 8 d).

### 6.3 Long-term statistics (it. 39, 40, 46)
* 1 Gyr single runs (asm kernel): Mercury's maximum eccentricity 0.303–0.345 in five runs (serial and parareal);
  single trajectories cannot test statistics (even two serial runs differ significantly in a KS test of consecutive
  20-Myr frequency windows, because the windows are correlated).
* **Ensemble of 16 IC realisations × 200 Myr, serial vs parareal (τ = 3e-5) from identical ICs:** Mercury's maximum
  eccentricity after decorrelation 0.3040 ± 0.0038 vs 0.3036 ± 0.0017 (KS p = 0.21), its mean 0.1876 vs 0.1884
  (p = 0.43); all eight inner secular frequencies indistinguishable (KS p = 0.06–0.94). **The parareal dynamics are
  statistically equivalent to serial integration** at this ensemble size.

### 6.4 Energy conservation and a numerical bias found on the way (it. 36, 46, 47)
* **A bias in the serial integrator itself.** Every WHFast512 integration call (i.e. every synchronisation with
  output) changes the total energy of the Solar System by a fixed dE/E ≈ +7.5 × 10⁻¹⁷, independent of the state
  (reproduced on unmodified upstream REBOUND: +7.54e-17 ± 0.49e-17, 15σ). Cause: the two precomputed 8 × 8 matrices
  for the inertial ↔ Jacobi transformation are not exact inverses in floating point. Their product is diagonal to
  10⁻¹⁹ with diagonal − 1 = [−1.37, −0.22, +0.07, −0.05, +0.45, −1.37, −0.17, −0.77] × 10⁻¹⁶ (Mercury … Neptune);
  scaling planet i's position and velocity by (1 + e_i) changes its energy by ≈ 2 e_i |U_i|, and the potential-energy
  weighted sum (Jupiter and Saturn dominate) predicts +7.6e-17, as measured. A 5 Gyr run with an output every kyr
  accumulates ≈ +4e-10 from this alone. Parareal, which restarts the integrator in every slice, showed it as a linear
  drift of +7.5e-10 per Gyr. Using REBOUND's recursive Jacobi transformation instead (state-dependent round-off)
  removes it (per-call offset consistent with zero, ≤ 1.5e-17 at 2σ).
* **Residual.** With the fix, the 16-member ensemble still shows a parareal drift of −1.25e-10 per Gyr relative to
  serial (≈ 3.5σ); it does not depend on the tolerance (τ = 3e-6: −1.15e-10 ± 0.62e-10; τ = 3e-5: −1.62e-10 ± 0.51e-10
  per Gyr, paired), which points to a per-slice systematic of the parareal driver rather than to incomplete
  convergence. Open.

## 7. Summary of the physical insight
1. Parareal for planetary systems is controlled by the coarse model's Jacobian in the **semi-major-axis rows**,
   because **Kepler shear** turns any error in a into a growing phase error; the osculating a is dominated by fast,
   phase-dependent short-period terms that a coarse integrator aliases.
2. A coarse propagator that works in **mean elements**, obtained by a numerical first-order perturbation theory with
   tables rebuilt from the current state, makes those Jacobian rows nearly exact; the near-resonant giants are better
   integrated directly with a large step than averaged. With this, parareal converges in ≈ 3–4 iterations for windows
   of 10–16 slices of ≈ 1 kyr.
3. The remaining coarse error is a phase error of the terrestrial planets caused by 10⁻¹⁰-level errors of their mean
   semi-major axes; its physical sources are small divisors near the Venus–Earth 8:13 / 3:5 commensurabilities, the
   Earth–Mars–Jupiter three-body near-resonance, second-order secular terms and the giants' coarse propagation. A
   second-order theory with these terms would allow two iterations.
4. The right accuracy yardstick is the serial integrator's own systematic time-step error; parareal adds an unbiased
   scatter well below it, and ensembles of 200 Myr show statistically identical dynamics.
