# Multicore (OpenMP)

REBOUND can use several cores of one computer through OpenMP. This page describes what runs in
parallel, which settings help, and what to expect.

## Building with OpenMP

- C: `make OPENMP=1` (in `src/` or in an example directory).
- Python: `REBOUND_OPENMP=1 pip install .` (or `pip install -e .` in a source checkout). Compiling
  with `-fopenmp` in `CFLAGS` also works.
- The number of threads is set with the environment variable `OMP_NUM_THREADS`, or with
  `rebound.omp_set_num_threads(n)` from Python. For reproducible timings, bind the threads:
  `OMP_PROC_BIND=close OMP_PLACES=cores`.

## Results do not depend on the number of threads

Every result of an OpenMP build is bitwise identical for every number of threads. With the default
settings it is also bitwise identical to the serial build (the same library compiled without
OpenMP). The opt-in modes below sum some terms in a different, fixed order. They are then still
identical between serial and OpenMP builds and for any number of threads, but differ from the
default mode at the round-off level.

The one exception is the *order in which collisions are resolved*, which is shuffled with
`sim.rand_seed`. That seed is random by default, so set it when you need reproducible runs with
collisions.

## What runs in parallel

Small loops always run serially. A parallel region costs about a microsecond, as much as a whole
time step of a small system, so it would make small systems slower.

| Part | Runs in parallel from | Notes |
|---|---|---|
| Basic gravity, all particles active | N ≈ 256 | uses about N/64 threads |
| Basic gravity with test particles | N_active · N_test ≳ 2,000 | |
| Tree gravity (also `tree_group_size > 0`) | N ≈ 32 | tree construction is serial |
| Tree collision search | N ≈ 128 | |
| WHFast Kepler drifts | N ≈ 64 | |
| MERCURIUS and TRACE pair forces | N ≈ 256 and 3 or more threads | they do twice the work of the serial loop |
| IAS15 with `parallel_stages = 1` | N ≈ 2 | see below |

Systems of a few planets (for example the Solar System with WHFast or WHFast512) do not get
faster on several cores. A time step takes about 100 ns, less than it takes to coordinate two
cores. To use many cores for such systems, run several simulations at once (ensembles, below).

## Hybrid CPUs (performance and efficiency cores)

On CPUs with performance and efficiency cores (Intel since Alder Lake, Apple silicon), REBOUND
uses only the performance cores for small and medium systems. Parallel loops wait for their
slowest thread, so the slower efficiency cores would hold them back. All cores are used from
N_active ≥ 1024 (from 512 in IAS15's stage mode), where the work is split into many small pieces.
Override this with `REBOUND_OMP_CORES=all` or `REBOUND_OMP_CORES=performance`. The limit only
applies to REBOUND's own loops.

## Opt-in: IAS15 stage mode

```python
sim.integrator = "ias15"
sim.integrator.parallel_stages = 1
```

IAS15 normally computes the forces of its 7 substeps one after another. In stage mode it computes
all 7 at once, and each of them is also split into blocks of particles, so many threads can work at
the same time. This is the setting that helps IAS15 with N ≈ 60–2000 particles. Measured on an
Intel i5-13500 (6 performance + 8 efficiency cores), compared with IAS15 on one core in the serial
build:

| N | 4 threads | 8 threads | 20 threads |
|---|---|---|---|
| 61 | 1.9x | 1.9x | 1.9x |
| 128 | 2.4x | 2.9x | 2.9x |
| 512 | 3.5x | 4.4x | 6.0x |
| 2048 | 4.1x | 6.3x | 8.6x |

Stage mode needs about 20 % more predictor-corrector iterations than the default mode. On one or
two threads it is therefore slower, and the default mode is the better choice there. It works with
basic and compensated gravity, test particles (both types), ghost boxes and `gravity_ignore_terms`.
It is ignored with additional forces, variational particles and particle maps. The setting is saved
in Simulationarchives.

## Opt-in: basic gravity in block order

```python
sim.gravity_block_order = 1
```

By default, the parallel basic gravity reproduces the serial loop bit for bit. For that it needs
about 1.5 times the work of the serial loop and a large table of pair terms. In block order, the
particles are instead split into fixed row blocks whose partial sums are added in a fixed order.
That is the work of the serial loop with nothing extra, and it splits well across threads. Results
are identical in serial and OpenMP builds and for any number of threads, and differ from the
default only at the round-off level (≈ 1e-15 relative in the accelerations).

Measured on a 2-core VM, compared with the default mode of the serial build:

| case | default, 2 threads | block order, 1 thread | block order, 2 threads |
|---|---|---|---|
| WHFast, 256 bodies | 0.90x | 1.18x | 1.90x |
| IAS15, N = 512 | 1.25x | 1.21x | 1.76x |
| IAS15, N = 2048 | 1.04x | 1.24x | 2.33x |
| IAS15, 20 + 2,000 test particles | 1.36x | 1.31x | 1.53x |

It applies to basic gravity in all integrators that use it. The setting is saved in
Simulationarchives.

## Ensembles: many simulations at once

```python
rebound.integrate_simulations(sims, tmax)   # each sim integrated to tmax
rebound.steps_simulations(sims, N_steps)    # each sim advanced by N_steps steps
```

These functions (C: `reb_simulations_integrate`, `reb_simulations_steps`) run independent
simulations on a pool of threads. They work in serial builds too, and use all logical CPUs by
default (`n_threads` sets the number). Results are bitwise identical to running the simulations
one after another. In an OpenMP build the workers are an OpenMP team, so they follow
`OMP_PROC_BIND` / `OMP_PLACES`, and each simulation runs single-threaded, so the cores are not
oversubscribed.

In `steps_simulations`, WHFast512 simulations are advanced two at a time with the paired kernel
when both can use it, with the same results. That gives about 1.2–1.5x more throughput per core:
- with the C AVX-512 kernel;
- with the AVX2 kernel for systems of up to 4 planets;
- with the generic kernel.

Eight-planet systems on CPUs without AVX-512 are not paired, because there is no gain.
