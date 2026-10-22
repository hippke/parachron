# Changes from REBOUND

Parachron is a modified version of [REBOUND](https://github.com/hannorein/rebound) by Hanno Rein and contributors,
distributed under the GNU General Public License v3 (see `LICENSE`). This file lists the modifications, as required
by section 5(a) of the license.

* Upstream base: hannorein/rebound commit `262a9cce` (Oct 16, 2026; REBOUND 5.2.2).
* Modifications: 2026, by Michael Hippke. The full history of every change is in this repository's git log.
* All upstream copyright and license notices are retained. Modified upstream source files carry a "Modified in 2026"
  line below their copyright notice; new files carry their own copyright and license header.

Summary of the changes (details in `changelog.md`):
* Mean-element parareal (MEP): parallel-in-time integration of one WHFast512 planetary system (`rebound.MEP`, `reb_mep_*`).
* WHFast512 without AVX-512 (portable C kernels), unbiased C Kepler step, opt-in `jacobi_recursive`.
* Single-core speed-ups, deterministic OpenMP, ensembles, IAS15 `parallel_stages`, `gravity_block_order`, tree group walk.
* Bug fix in SABA. Continuous-integration workflows and their helper scripts removed (`.github/`).

## Modified files (46)
* `.github/ISSUE_TEMPLATE/advice.md`
* `.github/ISSUE_TEMPLATE/bug_report.md`
* `.github/ISSUE_TEMPLATE/feature_request.md`
* `.github/ISSUE_TEMPLATE/question.md`
* `.gitignore`
* `README.md` — Parachron README
* `changelog.md` — changes of Parachron
* `docs/gravity.md`
* `docs/index.md`
* `docs/make_c_examples.py`
* `mkdocs.yml`
* `pyproject.toml` — package name and metadata
* `rebound/__init__.py`
* `rebound/simulation.py` — Python bindings of the new options and APIs
* `rebound/tests/test_whfast512.py`
* `setup.py`
* `src/Makefile` — build of the new source files
* `src/binarydata.c`
* `src/boundary.c`
* `src/collision.c` — deterministic parallel tree collisions, root-box pruning
* `src/gravity.c` — blocked and parallel (deterministic) basic gravity, `gravity_block_order`, tree group walk
* `src/gravity.h`
* `src/integrator_eos.c`
* `src/integrator_ias15.c` — vectorised predictor/corrector, multiplications by reciprocals, opt-in `parallel_stages`
* `src/integrator_ias15.h`
* `src/integrator_leapfrog.c`
* `src/integrator_mercurius.c` — broadphase, blocked and deterministic parallel force loops
* `src/integrator_mercurius.h`
* `src/integrator_saba.c` — fix: no copy of the unsynchronized state when not needed
* `src/integrator_sei.c`
* `src/integrator_trace.c` — broadphase, blocked and deterministic parallel force loops
* `src/integrator_trace.h`
* `src/integrator_whfast.c` — batched and parallel Kepler drifts
* `src/integrator_whfast.h`
* `src/integrator_whfast512.c` — portable C kernels and kernel selection, paired steps, `jacobi_recursive`, constants cache
* `src/integrator_whfast512.h`
* `src/integrator_whfast512.s`
* `src/particle.c`
* `src/rebound.c` — OpenMP core policy (hybrid CPUs), ensemble thread pool
* `src/rebound.h` — new fields and API (MEP, ensembles, WHFast512 options)
* `src/rebound_internal.h`
* `src/simulation.c` — ensembles (`reb_simulations_steps`), OpenMP work thresholds
* `src/tools.c`
* `src/transformations.c`
* `src/tree.c`
* `src/tree.h`

## New files (27)
* `CHANGES_FROM_REBOUND.md`
* `docs/img/parachron_icon.svg`
* `docs/img/parachron_logo.svg`
* `docs/openmp.md`
* `examples/mep_solar_system/Makefile`
* `examples/mep_solar_system/problem.c`
* `ipython_examples/MEP.ipynb`
* `rebound/mep.py` — Python interface of MEP
* `rebound/tests/test_ensemble.py`
* `rebound/tests/test_ensemble_pairing.py`
* `rebound/tests/test_gravity_block_order.py`
* `rebound/tests/test_ias15_parallel_stages.py`
* `rebound/tests/test_mep.py`
* `rebound/tests/test_openmp_determinism.py`
* `rebound/tests/test_tree_group.py`
* `src/integrator_whfast512_c_avx2.c`
* `src/integrator_whfast512_c_avx2x4.c`
* `src/integrator_whfast512_c_avx512.c`
* `src/integrator_whfast512_c_generic.c`
* `src/integrator_whfast512_c_genericx4.c`
* `src/integrator_whfast512_data.h`
* `src/integrator_whfast512_kernel.h` — portable C kernel of WHFast512 (unbiased Kepler step)
* `src/integrator_whfast512_neon_ops.h`
* `src/mep.c` — MEP driver
* `src/mep.h`
* `src/mep_coarse.c` — MEP coarse propagator
* `src/mep_dual.h`

## Removed files (13)
* `.github/codecov.yml`
* `.github/function_prefix_checker.bash`
* `.github/lint.bash`
* `.github/workflows/build_wheels.yml`
* `.github/workflows/c.yml`
* `.github/workflows/cmpi.yml`
* `.github/workflows/coverage.yml`
* `.github/workflows/emscripten.yml`
* `.github/workflows/function_prefix_checker.yml`
* `.github/workflows/ipynb.yml`
* `.github/workflows/lint.yml`
* `.github/workflows/python.yml`
* `.github/workflows/reboundx.yml`
