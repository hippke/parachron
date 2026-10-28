# WHFast512 kernel and MEP benchmark

`run.sh` builds several variants of the Parachron library from this checkout and measures, on the machine it runs on:

1. **Single core**: ns/step of each WHFast512 kernel and variant in the paper setup (Sun + 8 planets, dt = 6 d, GR,
   17th-order corrector), three interleaved repetitions, plus a hash of the final state (bitwise checks). If GitHub is
   reachable, unmodified upstream REBOUND (Parachron's base commit `262a9cc`) is built and measured too.
2. **Hardware counters** per step (cycles, instructions, FP operations by type), if `perf` is available.
3. **All-core load**: ns/step per core with 16 and with all physical cores busy (clock under load).
4. **MEP**: wall clock and speed-up against the serial reference (1600 slices, 1.58 Myr), for 16 cores with
   several kernels, and 32 and 64 cores with the assembly kernel.

| variant | what it tests |
|---|---|
| `default` kernel 1 | the original WHFast512 assembly (reference) |
| `default` kernels 2, 3 | Parachron's C kernels (AVX-512; AVX2 = 8 lanes as two 4-lane halves) |
| `hwrsqrt` | C kernels with the hardware 1/sqrt estimate as Newton seed (`-DREB_W512_HW_RSQRT`) |
| `sqrtdiv` | C kernels with 1/sqrt as sqrt + division (`-DREB_W512_SQRTDIV`; faster on Zen 4) |
| `native`, `native_sqrtdiv` | the same compiled with `-march=native` (compiler tuning for this CPU) |
| `clang`, `clang_sqrtdiv` | built with clang, if installed |
| `upstream` | unmodified REBOUND, assembly kernel (speed and bitwise check of Parachron's assembly path) |

Run time: about 10 minutes. Uses one core for most of the time, all physical cores briefly (sections 3 and 4), so run it
on an otherwise idle machine.

## Debian 13 (trixie)

```bash
sudo apt update
sudo apt install -y build-essential git linux-perf clang
```

`linux-perf` and `clang` are optional (sections 2 and the clang variants are skipped without them). Allow `perf` for
your user until the next reboot (optional; skip if you prefer not to change it):

```bash
sudo sysctl -w kernel.perf_event_paranoid=1
```

Clone the benchmark branch and run:

```bash
git clone -b zen5-bench https://github.com/hippke/parachron.git
cd parachron
benchmarks/whfast512_kernels/run.sh
```

The results are written to `w512bench_<host>_<date>/results.txt` in the current directory; send that file back. Optional
settings: `CORE=<cpu>` (CPU for the single-core runs, default the third physical core), `REPS=<n>` (default 3),
`SKIP_UPSTREAM=1`, `SKIP_MEP=1`.
