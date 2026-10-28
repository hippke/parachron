# Changes from REBOUND

Parachron is a modified version of [REBOUND](https://github.com/hannorein/rebound) by Hanno Rein and contributors,
distributed under the GNU General Public License v3 (see `LICENSE`). This file lists the modifications, as required
by section 5(a) of the license.

* Upstream base: hannorein/rebound commit `262a9cce` (Oct 16, 2026; REBOUND 5.2.2).
* Modifications: 2026, by Michael Hippke. The full history of every change is in this repository's git log.
* All upstream copyright and license notices are retained. Modified upstream source files carry a "Modified in 2026"
  line below their copyright notice; new files carry their own copyright and license header.

Summary of the changes (details in `changelog.md`):
* Mean-element parareal (MEP): parallel-in-time integration of one WHFast512 planetary system (`parachron.MEP`, `reb_mep_*`).
* WHFast512 without AVX-512 (portable C kernels), unbiased C Kepler step, opt-in `jacobi_recursive`.
* Single-core speed-ups, deterministic OpenMP, ensembles, IAS15 `parallel_stages`, `gravity_block_order`, tree group walk.
* Python module renamed from `rebound` to `parachron` (compiled library `libparachron`); the C API is unchanged.
* Bug fix in SABA. Continuous-integration workflows and their helper scripts removed (`.github/`).

## Modified files (118)
* `.github/ISSUE_TEMPLATE/advice.md`
* `.github/ISSUE_TEMPLATE/bug_report.md`
* `.github/ISSUE_TEMPLATE/feature_request.md`
* `.github/ISSUE_TEMPLATE/question.md`
* `.gitignore`
* `MANIFEST.in` — include this file in source distributions
* `README.md` — Parachron README
* `changelog.md` — changes of Parachron
* `docs/binaryformat.md`
* `docs/boundaryconditions.md`
* `docs/chaos.md`
* `docs/collisions.md`
* `docs/gravity.md`
* `docs/index.md`
* `docs/make_c_examples.py`
* `docs/make_integrators.py`
* `docs/miscellaneous.md`
* `docs/namingconvention.md`
* `docs/particles/add.md`
* `docs/particles/index.md`
* `docs/particles/operators.md`
* `docs/particles/orbitalelements.md`
* `docs/particles/remove.md`
* `docs/quickstart/firstexample.md`
* `docs/quickstart/installation.md`
* `docs/simulation/binaryfiles.md`
* `docs/simulation/diagnostics.md`
* `docs/simulation/index.md`
* `docs/simulation/operators.md`
* `docs/simulation/referenceframes.md`
* `docs/simulation/timestepping.md`
* `docs/simulation/variables.md`
* `docs/simulationarchive.md`
* `docs/units.md`
* `docs/visualization.md`
* `ipython_examples/AdvWHFast.ipynb`
* `ipython_examples/ChaoticHyperion.ipynb`
* `ipython_examples/Cheartbeat.ipynb`
* `ipython_examples/Checkpoints.ipynb`
* `ipython_examples/Churyumov-Gerasimenko.ipynb`
* `ipython_examples/CloseEncounters.ipynb`
* `ipython_examples/EccentricComets.ipynb`
* `ipython_examples/EmbeddedOperatorSplittingMethods.ipynb`
* `ipython_examples/EscapingParticles.ipynb`
* `ipython_examples/Forces.ipynb`
* `ipython_examples/FourierSpectrum.ipynb`
* `ipython_examples/FrequencyAnalysis.ipynb`
* `ipython_examples/HighOrderSymplectic.ipynb`
* `ipython_examples/Holmberg.ipynb`
* `ipython_examples/Horizons.ipynb`
* `ipython_examples/HybridIntegrationsWithTRACE.ipynb`
* `ipython_examples/HyperbolicOrbits.ipynb`
* `ipython_examples/IntegratingArbitraryODEs.ipynb`
* `ipython_examples/Megno.ipynb`
* `ipython_examples/OrbitPlot.ipynb`
* `ipython_examples/OrbitalElements.ipynb`
* `ipython_examples/PoincareMap.ipynb`
* `ipython_examples/PoincareSurfaceOfSection.ipynb`
* `ipython_examples/PrimordialEarth.ipynb`
* `ipython_examples/RadialVelocity.ipynb`
* `ipython_examples/RealtimeVisualizations.ipynb`
* `ipython_examples/Resonances_of_Jupiters_moons.ipynb`
* `ipython_examples/Rotations.ipynb`
* `ipython_examples/SaturnsRings.ipynb`
* `ipython_examples/Simulationarchive.ipynb`
* `ipython_examples/SimulationarchiveRestart.ipynb`
* `ipython_examples/Starman.ipynb`
* `ipython_examples/Testparticles.ipynb`
* `ipython_examples/TransitTimingVariations.ipynb`
* `ipython_examples/UniquelyIdentifyingParticlesWithNames.ipynb`
* `ipython_examples/Units.ipynb`
* `ipython_examples/User_Defined_Collision_Resolve.ipynb`
* `ipython_examples/VariationalEquations.ipynb`
* `ipython_examples/VariationalEquationsWithChainRule.ipynb`
* `ipython_examples/WHFast.ipynb`
* `ipython_examples/WHFast512.ipynb`
* `mkdocs.yml` — Parachron name, logo and links
* `pyproject.toml` — package name and metadata
* `python_examples/dragforce/problem.py`
* `python_examples/horizons/problem.py`
* `python_examples/longtermtest/problem.py`
* `python_examples/megno/problem.py`
* `python_examples/megno_simple/problem.py`
* `python_examples/orbital_elements/problem.py`
* `python_examples/outersolarsystem/problem.py`
* `python_examples/simple_orbit/problem.py`
* `python_examples/simulationarchive/problem.py`
* `setup.py` — compiled library named libparachron
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

## Renamed files (65)
* `rebound/__init__.py` → `parachron/__init__.py` (renamed and modified) — module name parachron, compiled library libparachron
* `rebound/binarydata_field_descriptor.py` → `parachron/binarydata_field_descriptor.py` (renamed only)
* `rebound/citations.py` → `parachron/citations.py` (renamed only)
* `rebound/frequency_analysis.py` → `parachron/frequency_analysis.py` (renamed only)
* `rebound/horizons.py` → `parachron/horizons.py` (renamed and modified)
* `rebound/integrator.py` → `parachron/integrator.py` (renamed only)
* `rebound/ode.py` → `parachron/ode.py` (renamed only)
* `rebound/orbit.py` → `parachron/orbit.py` (renamed and modified)
* `rebound/particle.py` → `parachron/particle.py` (renamed and modified)
* `rebound/particles.py` → `parachron/particles.py` (renamed only)
* `rebound/plotting.py` → `parachron/plotting.py` (renamed and modified)
* `rebound/rotation.py` → `parachron/rotation.py` (renamed and modified)
* `rebound/simulation.py` → `parachron/simulation.py` (renamed and modified) — Python bindings of the new options and APIs
* `rebound/simulationarchive.py` → `parachron/simulationarchive.py` (renamed and modified)
* `rebound/tests/__init__.py` → `parachron/tests/__init__.py` (renamed only)
* `rebound/tests/test_additional_forces.py` → `parachron/tests/test_additional_forces.py` (renamed and modified)
* `rebound/tests/test_binary_fields.py` → `parachron/tests/test_binary_fields.py` (renamed and modified)
* `rebound/tests/test_boundary.py` → `parachron/tests/test_boundary.py` (renamed and modified)
* `rebound/tests/test_bs.py` → `parachron/tests/test_bs.py` (renamed and modified)
* `rebound/tests/test_collisions.py` → `parachron/tests/test_collisions.py` (renamed and modified)
* `rebound/tests/test_copy.py` → `parachron/tests/test_copy.py` (renamed and modified)
* `rebound/tests/test_data.py` → `parachron/tests/test_data.py` (renamed and modified)
* `rebound/tests/test_eos.py` → `parachron/tests/test_eos.py` (renamed and modified)
* `rebound/tests/test_fpcontract.py` → `parachron/tests/test_fpcontract.py` (renamed and modified)
* `rebound/tests/test_frequency_analysis.py` → `parachron/tests/test_frequency_analysis.py` (renamed and modified)
* `rebound/tests/test_gravity.py` → `parachron/tests/test_gravity.py` (renamed and modified)
* `rebound/tests/test_horizons.py` → `parachron/tests/test_horizons.py` (renamed and modified)
* `rebound/tests/test_import.py` → `parachron/tests/test_import.py` (renamed and modified)
* `rebound/tests/test_integrator.py` → `parachron/tests/test_integrator.py` (renamed and modified)
* `rebound/tests/test_janus.py` → `parachron/tests/test_janus.py` (renamed and modified)
* `rebound/tests/test_leapfrog.py` → `parachron/tests/test_leapfrog.py` (renamed and modified)
* `rebound/tests/test_megno.py` → `parachron/tests/test_megno.py` (renamed and modified)
* `rebound/tests/test_mercurius.py` → `parachron/tests/test_mercurius.py` (renamed and modified)
* `rebound/tests/test_modify_orbital_parameters.py` → `parachron/tests/test_modify_orbital_parameters.py` (renamed and modified)
* `rebound/tests/test_name.py` → `parachron/tests/test_name.py` (renamed and modified)
* `rebound/tests/test_orbital_elements.py` → `parachron/tests/test_orbital_elements.py` (renamed and modified)
* `rebound/tests/test_particle.py` → `parachron/tests/test_particle.py` (renamed and modified)
* `rebound/tests/test_pickle.py` → `parachron/tests/test_pickle.py` (renamed and modified)
* `rebound/tests/test_plotting.py` → `parachron/tests/test_plotting.py` (renamed and modified)
* `rebound/tests/test_post_timestep_modifications.py` → `parachron/tests/test_post_timestep_modifications.py` (renamed and modified)
* `rebound/tests/test_restart_mercurius.py` → `parachron/tests/test_restart_mercurius.py` (renamed and modified)
* `rebound/tests/test_rotations.py` → `parachron/tests/test_rotations.py` (renamed and modified)
* `rebound/tests/test_saba.py` → `parachron/tests/test_saba.py` (renamed and modified)
* `rebound/tests/test_serialize.py` → `parachron/tests/test_serialize.py` (renamed and modified)
* `rebound/tests/test_server.py` → `parachron/tests/test_server.py` (renamed and modified)
* `rebound/tests/test_shearingsheet.py` → `parachron/tests/test_shearingsheet.py` (renamed and modified)
* `rebound/tests/test_simulation.py` → `parachron/tests/test_simulation.py` (renamed and modified)
* `rebound/tests/test_simulationarchive.py` → `parachron/tests/test_simulationarchive.py` (renamed and modified)
* `rebound/tests/test_simulationarchive_matrix.py` → `parachron/tests/test_simulationarchive_matrix.py` (renamed and modified)
* `rebound/tests/test_size_of_simulation.py` → `parachron/tests/test_size_of_simulation.py` (renamed and modified)
* `rebound/tests/test_tponly_encounter.py` → `parachron/tests/test_tponly_encounter.py` (renamed and modified)
* `rebound/tests/test_trace.py` → `parachron/tests/test_trace.py` (renamed and modified)
* `rebound/tests/test_trace_full_bs.py` → `parachron/tests/test_trace_full_bs.py` (renamed and modified)
* `rebound/tests/test_trace_full_ias15.py` → `parachron/tests/test_trace_full_ias15.py` (renamed and modified)
* `rebound/tests/test_transformations.py` → `parachron/tests/test_transformations.py` (renamed and modified)
* `rebound/tests/test_units.py` → `parachron/tests/test_units.py` (renamed and modified)
* `rebound/tests/test_variational.py` → `parachron/tests/test_variational.py` (renamed and modified)
* `rebound/tests/test_whfast.py` → `parachron/tests/test_whfast.py` (renamed and modified)
* `rebound/tests/test_whfast_advanced.py` → `parachron/tests/test_whfast_advanced.py` (renamed and modified)
* `rebound/tests/test_whfast_testparticles.py` → `parachron/tests/test_whfast_testparticles.py` (renamed and modified)
* `rebound/tools.py` → `parachron/tools.py` (renamed and modified)
* `rebound/units.py` → `parachron/units.py` (renamed and modified)
* `rebound/variation.py` → `parachron/variation.py` (renamed and modified)
* `rebound/vectors.py` → `parachron/vectors.py` (renamed only)
* `rebound/widget.py` → `parachron/widget.py` (renamed and modified)

## New files (34)
* `CHANGES_FROM_REBOUND.md`
* `RELEASING.md` — how to publish on PyPI
* `benchmarks/whfast512_kernels/README.md` — WHFast512 kernel and MEP benchmark
* `benchmarks/whfast512_kernels/kbench.c`
* `benchmarks/whfast512_kernels/mepbench.c`
* `benchmarks/whfast512_kernels/run.sh`
* `docs/img/parachron_icon.svg`
* `docs/img/parachron_logo.svg`
* `docs/openmp.md`
* `examples/mep_solar_system/Makefile`
* `examples/mep_solar_system/problem.c`
* `ipython_examples/MEP.ipynb`
* `parachron/mep.py` — Python interface of MEP
* `parachron/tests/test_ensemble.py`
* `parachron/tests/test_ensemble_pairing.py`
* `parachron/tests/test_gravity_block_order.py`
* `parachron/tests/test_ias15_parallel_stages.py`
* `parachron/tests/test_mep.py`
* `parachron/tests/test_openmp_determinism.py`
* `parachron/tests/test_random.py`
* `parachron/tests/test_tree_group.py`
* `parachron/tests/test_whfast512.py`
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

## Removed files (15)
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
* `rebound/tests/test_random.py`
* `rebound/tests/test_whfast512.py`
