<p align="center"><img src="docs/img/parachron_logo.svg" alt="Parachron" width="560"></p>

<p align="center">
<a href="LICENSE"><img src="https://img.shields.io/badge/license-GPL--3.0-green.svg?style=flat" alt="GPL-3.0"></a>
<a href="https://github.com/hannorein/rebound"><img src="https://img.shields.io/badge/based%20on-REBOUND%205.2-blue.svg?style=flat" alt="based on REBOUND"></a>
</p>

# Parachron

**One planetary system, one trajectory, many cores.**
Parachron integrates a single long orbit of the Solar System in parallel in time. A 5 Gyr integration of the Sun and
the eight planets (WHFast512, general relativity, 6 d step) takes **3.0 h on 16 cores instead of 9.2 h** on one core, at
the accuracy of the serial integrator. With relaxed tolerances, which are still validated by ensembles, it takes
1.8 h (up to 5×).

Parachron is a fork of [**REBOUND**](https://github.com/hannorein/rebound) by Hanno Rein and contributors. Everything
REBOUND does, Parachron does too, with the same API and results. On top of it, it adds mean-element parareal and a
set of single- and multi-core speed-ups. If you do not need these, please use upstream REBOUND.

## The idea: mean-element parareal (MEP)

A step of a Wisdom–Holman integrator takes about 100 ns, too short to split across cores. Every step needs the
exact end state of the previous one, so extra cores have so far only been used to run more simulations, not to
finish one faster.

*Parareal* breaks this chain. It cuts the time axis into slices and integrates all slices at once, each from a guessed
start state. A cheap *coarse* propagator links the slices, and a few correction sweeps make the result converge to
the serial trajectory. For planetary orbits, standard parareal fails. An error δa in a semi-major axis grows into a
phase error δλ = −(3/2) n T δa/a through **Kepler shear**. A cheap coarse integrator gets the short-period terms of
*a* wrong by far more than the corrections can absorb.

MEP gives the coarse propagator what Kepler shear needs:

* **Mean elements.** The coarse step evaluates a numerical first-order perturbation theory on the fly: pair
  accelerations on grids of mean longitudes, a 2-D FFT and the small divisors. With it, the coarse step removes the
  short-period terms from the inner planets' semi-major axes and advances their secular evolution. The giants are
  integrated directly with a long step. The result is accurate where shear amplifies errors, and costs about 3 % of a
  fine slice.
* **Corrections in orbital elements.** The parareal correction is formed in Jacobi elements, with longitudes
  modulo 2π. An energy projection then gives each corrected state the energy of the symplectic fine solution.
* **Pipelined windows.** Slices are processed in windows of N slices, one core each. Each window iterates until no
  planet moves by more than τ (3×10⁻⁶ AU by default) between iterations, which takes typically 3 iterations.
  Iterations overlap like a pipeline, and the next window starts as soon as the current one has converged.

The results are statistically indistinguishable from serial integrations. The secular frequencies agree to within
10⁻⁷ ″/yr, far below the serial integrator's own time-step error. An ensemble of 16 runs of 200 Myr shows the same
dynamics as 16 serial runs, and the energy errors match. Details: Hippke et al. (in prep.), *Mean-Element Parareal:
Parallel-in-Time Symplectic Integration of the Solar System*.

| setting (16 cores, AMD Ryzen 9 7950X) | speed-up | 5 Gyr |
|---|---|---|
| serial WHFast512 (reference) | 1× | 9.2 h |
| MEP, τ = 3×10⁻⁶ AU (default, serial accuracy) | 3.1× | 3.0 h |
| MEP, relaxed tolerances (ensemble-validated) | 3.4–5.0× | 1.8–2.7 h |

## Example

```python
import math
import rebound          # Parachron keeps REBOUND's Python module name (see below)

sim = rebound.Simulation()
sim.add("solar system")                 # Sun + 8 planets
sim.integrator = "whfast512"
sim.dt = 6. / 365.25 * 2. * math.pi     # 6 days (units: AU, Msun, yr/2pi)
sim.exact_finish_time = False
sim.integrator.gr_potential = 1
sim.integrator.corrector = 17

mep = rebound.MEP(sim, checkpoint="run.mep")   # one thread per physical core
mep.integrate(1e8 * 2. * math.pi)              # 100 Myr; sim is updated to the end state
print(mep.mean_iterations, mep.walltime_last)

# later, possibly in another process: continues bitwise identically
mep = rebound.MEP.from_checkpoint("run.mep")
mep.integrate(2e8 * 2. * math.pi)
```

A C version is in [`examples/mep_solar_system`](examples/mep_solar_system/problem.c). MEP runs on Linux and needs a
central star with exactly 8 planets. It was designed and validated for the Solar System; other systems are on the
to-do list.

## What else is different from REBOUND

* **WHFast512 without AVX-512.** Portable C kernels (AVX2+FMA, AVX-512, generic, ARM NEON) next to the assembly. The
  kernels give bitwise identical results on every machine, and two simulations can be stepped together.
* **`jacobi_recursive`** (opt-in): removes a systematic energy offset of +7.5×10⁻¹⁷ per WHFast512 synchronization.
* **Faster on one core.** The use-case benchmarks run 2.4× faster than upstream (geometric mean). Results are
  bitwise identical, except for round-off-level changes in IAS15.
* **Deterministic OpenMP.** Results do not depend on the thread count and equal those of the serial build. There is
  also an opt-in stage-parallel IAS15 (`parallel_stages`) for N ≳ 100.

## Installation

```bash
git clone https://github.com/hippke/parachron && cd parachron && pip install -e .
```

The package is called `parachron`, but its Python module keeps the name `rebound`, so existing REBOUND scripts run
unchanged. It therefore replaces an installed REBOUND: use a separate virtual environment if you need both. The C
library is built as before (`make` in `src/`, or see the examples).

## Credit and citation

Parachron is built on REBOUND. Please cite REBOUND as described in its [documentation](https://rebound.hanno-rein.de/)
(`sim.cite()` lists the papers for the methods you use), and in addition the MEP paper if you use `rebound.MEP`.
REBOUND is developed by Hanno Rein and many contributors: see the
[upstream repository](https://github.com/hannorein/rebound) and its git history. All upstream copyright notices are
retained in the source files.

## License

Parachron is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
version, like REBOUND. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the [GNU General Public License](LICENSE)
for more details.

Copyright © 2011–2026 Hanno Rein and REBOUND contributors (REBOUND).
Modifications and additions © 2026 Michael Hippke. The changes relative to REBOUND are listed in
[`CHANGES_FROM_REBOUND.md`](CHANGES_FROM_REBOUND.md).
