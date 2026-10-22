# Part 3 — Implementation notes (brief)

* **Code.** C (`../pint/`): `pint_coarse.c` (coarse step, helper pool, staged schedule), `pint_pipe.c` (pipelined
  driver, production runs with `PINT_OUT`), `pint_long.c` (non-pipelined validation driver), `pint_dual.h` (dual-number
  derivatives), `pint_disp6.c` (all-component transform: reference, kernel and task-parallel versions). Python
  prototypes and analysis in `../parareal/` and `../pint/*.py` (NAFF frequency analysis, ensemble statistics).
  Fine propagator: REBOUND's WHFast512 (assembly or C kernels) of this fork.
* **Bitwise reproducibility.** All parallel reductions run in a fixed order, so serial and parallel coarse steps are
  bitwise identical for any number of helpers (`test_par`); the C kernels of WHFast512 are bitwise identical across
  instruction sets, and the coarse step falls back to the AVX2 C kernel on CPUs without AVX-512.
* **Concurrency details.** Thread-local scratch arrays in the coarse step (several slice threads call it at the same
  time); slice threads pinned to physical cores; helpers on a staged pool with per-stage atomic task counters and
  stage barriers, sleeping between jobs (spinning was slower on the 7950X); persistent slice threads with double-buffered
  window state; trajectory output written by the main thread from a ring buffer filled by the deciding slice thread.
* **WHFast512 change.** Opt-in field `jacobi_recursive` (default 0 = bitwise identical to upstream): the inertial ↔
  Jacobi conversion at every synchronisation uses REBOUND's recursive transformations instead of the precomputed 8 × 8
  matrices (Part 1 §6.4). REBOUND's WHFast512 test suite passes (75 tests). A standalone reproduction of the upstream
  energy offset is `../sync_bias_repro.c`; draft report `../WHFAST512_SYNC_BIAS.md`.
* **Low-level performance work** contributed little to the multicore speed-up: the fine solver (≈ 77 % of the time)
  was already a hand-optimised kernel; the gains came from the coarse model (Part 1) and the schedule (Part 2).
* **Reproducing the numbers.** `../pint/README.md` (build and run commands, including the 5 Gyr run),
  `../MACHINES.md` (hardware), `../RESEARCH_LOG.md` (every measurement with its command), `../pint_history.csv`.
