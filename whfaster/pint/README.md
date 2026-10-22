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
    make pint_pipe pint_long test_par     # REB defaults to ../../src of this checkout

## Run on a machine with 16 physical cores (measured on an AMD Ryzen 9 7950X, see ../RESEARCH_LOG.md it. 28-38)
    # arguments: N slice_days K_min windows kernel serial pin tol_AU K_max helpers parmode
    ./pint_pipe 16 360000 3 100 1 1 1 3e-6 6 8 2    # reference accuracy: ~3.0x, frequencies within 1e-7 "/yr
    ./pint_pipe 16 360000 3 100 1 1 1 1e-5 6 8 2    # validated by a 7-run ensemble: ~3.2x
    PINT_OUT=run ./pint_pipe 16 360000 3 5070 1 0 1 3e-6 6 8 2   # 80 Myr, trajectory in run_pr.bin, no serial
    # AVX2-only CPUs (no AVX-512): kernel argument 3; on the 7950X kernel 1 (assembly) is 8 % faster than 2
    # the 5 Gyr run of arXiv 2609.25339 (Solar System ICs of tools/ss_ics.h): 316900 windows x 16 x 986 yr,
    # about 3.1 h on the 7950X at reference accuracy (serial: 9.1 h); output: one record per 986 yr (~2 GB)
    PINT_OUT=ss5gyr ./pint_pipe 16 360000 3 316900 1 0 1 3e-5 6 8 2
Environment: `PINT_OUT` (trajectory), `PINT_PERT` (IC shift), `PINT_SERONLY`/`PINT_DT` (serial reference only, its
step), `PINT_JR` (default 1: WHFast512 `jacobi_recursive`, removes an energy drift of +7e-17 per kernel call),
`PINT_CHUNKS`, `PINT_TRACE`, `PINT_UPDDIAG`, coarse-step experiments `PINT_MEANN`, `PINT_NGHI`, `PINT_ALPHAHI`,
`PINT_DTOUT`. Use physical cores only (slice n is pinned to CPU n) and an otherwise idle machine.

## Old (round 2) usage on 10 cores
    ./pint_pipe 10 360000 3 20 2 1 1 3e-6 6   # recommended: 10 cores, adaptive K (>=3, stop when update <= 3e-6 AU)
    ./pint_pipe 10 360000 3 20          # fixed K=3 (faster, but secular frequencies only to ~3e-6 arcsec/yr)
    ./pint_pipe 20 360000 4 10          # 20 cores, K=4
    # AVX2-only CPUs (no AVX-512, e.g. Intel client parts): add the kernel argument 3
    ./pint_pipe 10 360000 3 20 3
Output: parareal wall time, serial wall time of the same span, speedup, core-hours ratio, and the final position
difference to the serial run (expected ~1e-6 AU after 200 kyr: boundary jumps at the round-off-noise level,
amplified by shear; statistically equivalent to a serial run — see iteration 12). Use physical cores (no SMT
siblings), a quiet machine, and fixed frequency if possible.

Model (measured here on one core of a 2.1 GHz Xeon: fine slice 11.5–12.5 ms, coarse step 1.0 ms, r = 0.083):
S = N / (N r + K (1 + r)) → with the adaptive rule (mean K = 3.3, secular frequencies within ~1e-7 arcsec/yr of the
serial run over 10 Myr): 2.3 on 10 cores; fixed K = 3: 2.45 (frequencies only to ~3e-6); fixed K = 4: 1.9.

## Other programs
* `pint_test` — coarse-step check and timing (r = G/F).
* `pint_run N days threads kmax dt_out start heun` — non-pipelined parareal with per-iteration errors against a
  serial reference and the asm-vs-C noise.
* `pint_long prefix W N days threads K serial xtol diag serial_kernel proj heun pert` — long runs writing the state
  at every boundary (`*_pr.bin`, `*_ser<k>.bin`), boundary-jump diagnostics; `energy_ana.py`, `freq_ana.py`,
  `drift_ana.py` analyse them.
