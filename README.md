<p align="center"><img src="https://raw.githubusercontent.com/hippke/parachron/main/docs/img/parachron_logo.svg" alt="Parachron" width="560"></p>

<p align="center">
<a href="https://github.com/hippke/parachron/blob/main/LICENSE"><img src="https://img.shields.io/badge/license-GPL--3.0-green.svg?style=flat" alt="GPL-3.0"></a>
<a href="https://pypi.org/project/parachron/"><img src="https://img.shields.io/pypi/v/parachron.svg?style=flat" alt="PyPI"></a>
<a href="https://github.com/hannorein/rebound"><img src="https://img.shields.io/badge/based%20on-REBOUND%205.2-blue.svg?style=flat" alt="based on REBOUND"></a>
</p>

# Parachron

Parachron is a fork of [**REBOUND**](https://github.com/hannorein/rebound) by Hanno Rein and contributors. Everything
REBOUND does, Parachron does too, with the same API and results. On top of it, it adds mean-element parareal and a
set of single- and multi-core speed-ups. If you do not need these, please use upstream REBOUND.

## Mean-element Parareal: Parallel-in-Time Integration of the Solar System

*[The accompanying paper](http://arxiv.org/abs/2611.01237) has been submitted.*

Long-term integrations of planetary systems have long been considered inherently sequential: each time step must
start from the *exact* state at the end of the previous one, because any error in that state, however small, is
amplified by orbital shear and chaos and alters the trajectory. Because a single step of a modern Wisdom–Holman
integrator takes only ~100 ns—too short to split across cores—additional cores have previously only been used to run
multiple simulations in parallel. We introduce mean-element parareal (MEP), a parallel-in-time scheme that
accelerates a single Solar System integration (eight planets, general relativity) at the same accuracy as the serial
integrator. Combined with a latency-optimized C kernel that integrates 5 Gyr on one core in 5.8 h, MEP completes a
5 Gyr integration in 52 min on 64 cores (a parallel speed-up of 6.7×, and 10.4× faster than the previous 9.0 h
single-core record) or in 1.8 h on 16 cores (3.3× parallel speed-up, 5.0× faster than the previous record). MEP
succeeds because its coarse propagator uses a numerical first-order mean-element theory to remove the short-period
oscillations from the semi-major axes, preventing Keplerian shear from amplifying them into growing phase errors.
Standard parareal with a coarse Wisdom–Holman integrator fails for planetary orbits because no inexpensive coarse
integrator can resolve these oscillations. To our knowledge, this is the first measured multicore speed-up of a single
long-term planetary integration and opens a pathway towards larger-scale parallelization. MEP is freely available in
Parachron, a modified version of REBOUND, at https://github.com/hippke/parachron.

## Example

The notebook [`MEP_SolarSystem.ipynb`](https://github.com/hippke/parachron/blob/main/ipython_examples/MEP_SolarSystem.ipynb)
shows how to integrate the Solar System with MEP, in the setup of the paper, and compares the result with serial
WHFast512. Stopping and resuming long runs from a checkpoint is shown in
[`MEP.ipynb`](https://github.com/hippke/parachron/blob/main/ipython_examples/MEP.ipynb), and a C version is in
[`examples/mep_solar_system`](https://github.com/hippke/parachron/blob/main/examples/mep_solar_system/problem.c). MEP
runs on Linux and macOS and needs a central star with exactly 8 planets. It was designed and validated for the Solar
System; other systems are on the to-do list.

## Installation

Parachron is on [PyPI](https://pypi.org/project/parachron/):

```bash
pip install parachron
```

Parachron is distributed as source only: pip compiles the C library on your machine during the installation, so a C
compiler is needed (on Debian/Ubuntu `sudo apt install build-essential`, on macOS `xcode-select --install`). MEP runs
on Linux and macOS (threads are pinned to cores on Linux only).

**The compiler sets the speed.** The WHFast512 kernels are tuned for GCC:

- Use **GCC 13 or newer**. GCC 14 is about 1 % faster than GCC 13; clang 18 is 9–13 % slower (Ryzen 9 7950X). Older
  GCC versions have not been measured.
- pip uses the compiler in `CC`, or else the default `cc` (usually GCC on Linux, clang on macOS); conda environments
  may set `CC` to their own compiler. Check with `echo $CC` and `cc --version`, and select GCC explicitly if needed:
  ```bash
  CC=gcc-14 pip install parachron
  ```
- pip caches the library it built. After changing the compiler, rebuild it with
  `CC=gcc-14 pip install --force-reinstall --no-cache-dir parachron`.
- The build does not depend on the CPU it is compiled on: kernels for AVX-512, AVX2 and plain x86-64 are all compiled
  in, and the fastest one for the CPU is chosen at run time. Their results are bitwise identical to each other, and
  between GCC 13 and GCC 14 builds.

To install the latest development version from GitHub instead:

```bash
pip install git+https://github.com/hippke/parachron
```

The Python module is `parachron` and can be installed next to REBOUND. Its API is REBOUND's: scripts written for
REBOUND run after replacing `import rebound` with `import parachron as rebound`. The C library keeps REBOUND's C API
(`reb_*` functions, `rebound.h`) and is built as before (`make` in `src/`, or see the examples).

## Contributing

All my repositories embrace collaboration between humans and AI. Contributions are welcome, whether they are written by you, by an AI assistant, or by both together, and they will be handled promptly. Please open an [issue](https://github.com/hippke/parachron/issues) or a pull request.

## Credit and citation

Parachron is built on REBOUND. Please cite REBOUND as described in its [documentation](https://rebound.hanno-rein.de/)
(`sim.cite()` lists the papers for the methods you use), and in addition the MEP paper if you use `parachron.MEP`.
REBOUND is developed by Hanno Rein and many contributors: see the
[upstream repository](https://github.com/hannorein/rebound) and its git history. All upstream copyright notices are
retained in the source files.

## License

Parachron is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
version, like REBOUND. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the [GNU General Public License](https://github.com/hippke/parachron/blob/main/LICENSE)
for more details.

Copyright © 2011–2026 Hanno Rein and REBOUND contributors (REBOUND).
Modifications and additions © 2026 Michael Hippke. The changes relative to REBOUND are listed in
[`CHANGES_FROM_REBOUND.md`](https://github.com/hippke/parachron/blob/main/CHANGES_FROM_REBOUND.md).
