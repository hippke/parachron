# Parallel-in-time WHFast512 (parareal with a mean-element coarse propagator)

Research code of whfaster round 2 (see `../RESEARCH_LOG.md`, iterations 1–13). One long Solar System trajectory
(Sun + 8 planets, dt = 6 d, GR potential, order-17 correctors) is integrated in windows of N slices of 986 yr
(360 000 d); N cores run the fine WHFast512 slices concurrently, a cheap coarse propagator (`pint_coarse.c`) links
them, K = 3 parareal iterations per window (N ≤ 10; K = 4 for N = 20–40).

Coarse propagator (config B): inner planets — first-order mean semi-major axes from a numerical short-period theory
(Jacobi pair accelerations on 32² / 16² grids of mean longitudes, 2-D FFT, rebuilt from the input state), Kepler +
first-order secular drift with the trapezoidal rule; giants — WHFast512 at 360 d in a Sun(+inner masses) + 4-planet
subsystem; energy projection of the corrected boundary states.

## Build (needs the fork's `src/librebound.so`)
    cd ../../src && make -j OPENGL=0 SERVER=0 librebound.so && cd -
    # edit REB= in Makefile if the repository is not at /home/claude/rebound
    make pint_pipe pint_run pint_long pint_test

## Measure the wall-clock speedup on a machine with >= 10 physical cores
    ./pint_pipe 10 360000 3 20          # 10 threads (one per core), K=3, 20 windows = 197 kyr, + serial run
    ./pint_pipe 20 360000 4 10          # 20 cores, K=4
    # AVX2-only CPUs (no AVX-512, e.g. Intel client parts): add the kernel argument 3
    ./pint_pipe 10 360000 3 20 3
Output: parareal wall time, serial wall time of the same span, speedup, core-hours ratio, and the final position
difference to the serial run (expected ~1e-6 AU after 200 kyr: boundary jumps at the round-off-noise level,
amplified by shear; statistically equivalent to a serial run — see iteration 12). Use physical cores (no SMT
siblings), a quiet machine, and fixed frequency if possible.

Model (measured here on one core of a 2.1 GHz Xeon: fine slice 11.5–12.5 ms, coarse step 1.0 ms, r = 0.083):
S = N / (N r + K (1 + r)) → 2.45 on 10 cores (K = 3), 3.3 on 20 cores (K = 4), 1.4 on 5 cores.

## Other programs
* `pint_test` — coarse-step check and timing (r = G/F).
* `pint_run N days threads kmax dt_out start heun` — non-pipelined parareal with per-iteration errors against a
  serial reference and the asm-vs-C noise.
* `pint_long prefix W N days threads K serial xtol diag serial_kernel proj heun pert` — long runs writing the state
  at every boundary (`*_pr.bin`, `*_ser<k>.bin`), boundary-jump diagnostics; `energy_ana.py`, `freq_ana.py`,
  `drift_ana.py` analyse them.
