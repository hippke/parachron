# REBOUND use-case performance benchmarks

Fourteen short benchmarks, each reproducing the hot loop of a compute-heavy
research use of REBOUND (see `../../../reports/02_use_cases.md` for the papers).
A full run takes about 10 s, so it can be repeated after every change.

```
python3 run.py run                 # build librebound + bench, run all, compare with baseline
python3 run.py set-baseline ../../../bench_results/<file>.json
python3 run.py run --only trace_tp,eos_pmlf6_tp   # iterate on a subset
python3 run.py run --repeat 3      # ~30 s; best-of-3 per chunk, tighter numbers
python3 run.py compare A.json B.json
python3 run.py probe [names] --time 8             # check how fast rates settle
```

Results are written to `bench_results/` next to this checkout (outside the git
tree), so they survive checkouts of other commits. To benchmark another commit
without touching this tree, use a worktree:

```
git worktree add ../rebound-base <commit>
python3 run.py run --src ../rebound-base --label base
```

Requirements: a C compiler, `make`, Python 3.8+ (standard library only).
The library is built with the repository's own `src/Makefile` (`OPENGL=0 SERVER=0`),
i.e. the same flags users get from `make`.

## Benchmarks

| Name | Unit | Use case (paper) | What it exercises |
|---|---|---|---|
| `ss_wh_gr` | steps/s | Solar System 5 Gyr ensembles, plain WH + GR, dt = 3.16 d (Abbot et al. 2021, 2023) | WHFast Jacobi step, Kepler solver, GR force callback |
| `ss_whckl_gr` | steps/s | Solar System, WHCKL (17th-order corrector, lazy kernel) + GR, dt = 8.06 d (Brown & Rein 2020, 2022) | Corrector, lazy-implementer kernel |
| `ss_whfast512` | steps/s | Solar System with WHFast512 + GR (Javaheri, Rein & Tamayo 2023); skipped without AVX-512 | AVX-512 kernel |
| `obertas_5p` | steps/s | Five Earth-mass planets, WHFast dt = 0.05 yr; initial conditions from the paper's data (Obertas et al. 2017) | WHFast with 5 planets |
| `spock_3p_line` | steps/s | Random compact 3-planet systems, `safe_mode=0`, line collision search (Tamayo et al. 2016, 2020; SPOCK) | WHFast + collision search every step |
| `ias15_hltau` | orbits/s | IAS15, five Saturn-mass planets in the HL Tau gaps (Rein & Tamayo 2017) | IAS15 predictor–corrector and step control |
| `ias15_embryos` | yr/s | IAS15, 60 embryos with direct collision search and merging (giant-impact phase) | IAS15 at moderate N, collisions |
| `megno_grid` | points/s | MEGNO stability map, 2 planets × 500 yr per point (Hadden & Lithwick 2018; `ipython_examples/Megno.ipynb`) | Variational equations, per-simulation setup |
| `short_posterior` | sims/s | Many 300-step WHFast runs of posterior draws (Wang et al. 2018, HR 8799) | Simulation create/add/free overhead |
| `rings_sei_tree` | steps/s | Saturn's rings, SEI + tree collisions, 200×5×20 root boxes (Rein & Latter 2013; `examples/overstability`) | Collision tree search, shear boundary |
| `selfgrav_tree` | steps/s | Self-gravitating disc, N = 10,000, tree gravity (`examples/selfgravity_disc`) | Tree build and walk |
| `mercurius_plts` | steps/s | MERCURIUS, 2 planets + 2,500 planetesimals (`examples/planetesimal_disk_migration`) | Hybrid encounters, IAS15 sub-steps, merging |
| `trace_tp` | steps/s | TRACE, 2 planets + 2,000 test particles | TRACE encounter checks (N×N bookkeeping) |
| `eos_pmlf6_tp` | steps/s | EOS PMLF6 (jerk kernel), 2 planets + 300 test particles | Jerk loop (O(N_tp²) at the time of writing) |

Adaptive integrators report simulated time per second rather than steps per
second, so a change that lets IAS15 take larger steps counts as a speedup.

## Methodology

**Fixed work, not fixed time.** Every benchmark runs a fixed number of
"chunks" (a chunk is a batch of steps, one MEGNO grid point, or one short
simulation). The counts are set in `bench.c` (`warm_chunks`, `meas_chunks`) so
that the measured part takes about 0.5 s on the reference machine. Fixed work
means a faster build integrates exactly the same simulated interval as the
baseline. Several systems evolve over time (IAS15 step sizes, encounters,
disc evolution); with a fixed time budget a faster build would sample a
different part of that evolution and bias the comparison.

**Is half a second enough?** `run.py probe` runs each benchmark for 8 s and
prints the rate in 0.5 s bins relative to the final rate. On the reference
machine all fixed-step benchmarks are within ±2% of their final rate in the
first 0.5 s bin. There is no warm-up effect beyond the first few chunks.
The ones that drift (`ias15_hltau` by ~25% over 8 s, `mercurius_plts`
and `selfgrav_tree` by ~5%) drift because the physics evolves, which fixed
work handles. Example probe output (rate / final rate per 0.5 s bin):

```
ss_wh_gr       1.00 0.99 0.99 0.99 0.97 1.00 1.00 1.00 1.00 1.00 ...
obertas_5p     1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 ...
ias15_hltau    1.13 1.11 1.11 1.11 1.10 1.09 1.09 1.09 1.07 1.04 1.01 0.99 1.00 0.95 0.95 0.87
selfgrav_tree  1.05 0.97 0.96 0.95 0.94 0.94 0.95 0.95 0.97 0.98 0.99 0.99 1.01 1.02 1.02 1.02
```

**Comparison.** Because both runs execute the same deterministic chunk
sequence, chunk *k* does the same work in both. The speedup of a benchmark is
the median over chunks of (new chunk rate / baseline chunk rate). The median
ignores chunks hit by a slow episode (another process, VM steal), which on a
shared machine were the main source of noise. With `--repeat N`, each chunk's
fastest time over the repeats is used. The single summary number is the
geometric mean of the 14 per-benchmark speedups.

**Noise.** Six identical runs on the reference machine (cloud VM, Xeon with
AVX-512, gcc 13.3) compared pairwise (15 pairs):

- geomean speedup within ±0.7% of 1.000;
- per benchmark, the worst deviation was 0.5–2.9% for 13 of the 14 benchmarks and 3.8% for `eos_pmlf6_tp`.

The comparison flags a benchmark as changed (`+`/`-`) only beyond
max(3%, 3 standard errors); adjust with `--threshold`. A dedicated workstation
is usually quieter.

**Result check.** Each run also records a signature of the final particle
state (a weighted sum of positions, velocities and masses over the fixed
work). It is bit-for-bit reproducible between runs. The comparison prints
`identical` when an optimization leaves results unchanged bit for bit, or the
relative difference otherwise. Many valid optimizations (reordered floating-point
sums) change results at the 1e-16 level per step; chaotic benchmarks then
diverge further over the fixed work, so a `CHANGED` result is a prompt to
check, not necessarily an error.

## Files

- `bench.c` — benchmark definitions and harness (public REBOUND API only, so it builds against other source trees).
- `run.py` — build, run, store, compare.
- `obertas_ics.h` — four of the widest-spaced, longest-lived systems from Obertas et al. 2017, generated by `tools/extract_obertas_ics.py` from the data files in github.com/dtamayo/spock.
