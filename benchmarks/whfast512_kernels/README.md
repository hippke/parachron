# WHFast512 kernel and MEP benchmark

`run.sh` builds several variants of the Parachron library from this checkout and measures, on the machine it runs on:

1. **Single core**: speed of each WHFast512 kernel and variant relative to the assembly kernel in the paper setup
   (Sun + 8 planets, dt = 6 d, GR, 17th-order corrector), plus a hash of the final state (bitwise checks). If GitHub is
   reachable, unmodified upstream REBOUND (Parachron's base commit `262a9cc`) is built and measured too.
2. **Hardware counters** per step (cycles, instructions, FP operations by type), if `perf` is available.
3. **Multi-core load**: ns/step per core with 16 copies running at the same time (clock under load).
4. **MEP**: wall clock and speed-up against the serial reference (1600 slices, 1.58 Myr), for 16 cores with
   several kernels, and 32 and 64 cores with the assembly kernel.

The machine does not need to be idle. Every comparison runs in ABBA blocks (reference, variant, variant, reference;
MEP: serial, MEP cases, the same reversed, serial), so drifts of the background load cancel in each block's ratio, and
the runs are pinned to the physical cores that were the most idle (both hardware threads) just before. The output
gives the spread over blocks, a null test (the assembly against itself) for the size of the remaining noise, and the
idle fraction of the cores used.

| variant | what it tests |
|---|---|
| `default` kernel 1 | the original WHFast512 assembly (reference) |
| `default` kernels 2, 3 | Parachron's C kernels (AVX-512; AVX2 = 8 lanes as two 4-lane halves, built with AVX-512VL where available), 1/sqrt as sqrt + division |
| `novl` | kernel 3 without its AVX-512VL build (plain AVX2 encoding, 16 registers), for comparison |
| `horner` | C kernels with the Stiefel series by plain Horner instead of Estrin's scheme (`-DREB_W512_STIEFEL_HORNER`; the default before Oct 29, 2026) |
| `hwrsqrt` | Newton 1/sqrt with the hardware estimate as seed (`-DREB_W512_HW_RSQRT`) |
| `newton` | C kernels with the earlier Newton 1/sqrt from a bit-trick seed (`-DREB_W512_NEWTON_RSQRT`; faster on Intel Sapphire Rapids) |
| `native` | the default compiled with `-march=native` (compiler tuning for this CPU) |
| `clang`, `clang_newton` | built with clang, if installed |
| `upstream` | unmodified REBOUND, assembly kernel (speed and bitwise check of Parachron's assembly path) |

Run time: about 15 minutes. Uses one core for most of the time, and 16 to 64 cores briefly (sections 3 and 4).

## Debian 13 (trixie)

```bash
sudo apt update
sudo apt install -y build-essential wget linux-perf clang
```

`linux-perf` and `clang` are optional (sections 2 and the clang variants are skipped without them). Allow `perf` for
your user until the next reboot (optional; skip if you prefer not to change it):

```bash
sudo sysctl -w kernel.perf_event_paranoid=1
```

Download Parachron and run (no git needed; `git clone https://github.com/hippke/parachron.git`
works too where git can reach GitHub without credentials):

```bash
wget https://github.com/hippke/parachron/archive/refs/heads/main.tar.gz
tar xzf main.tar.gz && cd parachron-main
benchmarks/whfast512_kernels/run.sh
```

The script downloads upstream REBOUND the same way (curl or wget); it never asks for credentials.

The results are written to `w512bench_<host>_<date>/results.txt` in the current directory; send that file back. Optional
settings: `CORE=<cpu>` (fixed CPU for the single-core runs instead of the most idle one), `BLOCKS=<n>` (ABBA blocks per
single-core comparison, default 5; more blocks give tighter ranges on a busy machine),
`SKIP_UPSTREAM=1`, `SKIP_MEP=1`.
