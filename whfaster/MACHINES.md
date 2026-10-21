# Test machines

Written Oct 21, 2026. Use these machines for all whfaster / parareal timing runs. Record the machine in every
`RESEARCH_LOG.md` entry and every `pint_history.csv` row.

## host2: 16-core server (main machine for multicore runs)
| | |
|---|---|
| access | `ssh root@23.239.118.238` (key-based, no password; hostname `host2`) |
| CPU | AMD Ryzen 9 7950X, 16 cores / 32 threads, 1 socket, AVX-512 (f, dq, bw, vl, cd, ifma, vbmi) |
| memory / OS | 93 GB, Ubuntu 24.04.4 LTS, gcc 13.3.0, `perf` installed |
| Python | python3 present, **no numpy** (Python analysis such as `freq_ana.py` and `ph1eval.py` runs on the local machine) |
| repo | `/home/ubuntu/rebound`, branch `whfaster` |
| load | idle (load average ≈ 0); dedicated to these tests |

Use it for anything that needs more than a few cores, and for all speedup numbers that are compared with the
earlier 7950X results in `RESEARCH_LOG.md` (Phase 8). WHFast512 uses the AVX-512 assembly kernel here.
Fine slice ≈ 7.0 ms, serial reference for `pint_pipe ... 100` windows ≈ 11.3 s.

## Local workstation (development, Python, small runs)
| | |
|---|---|
| CPU | Intel Core i5-13500: 6 P-cores (2 threads each, up to 4.8 GHz) + 8 E-cores (1 thread, up to 3.5 GHz) = 20 threads |
| ISA | AVX2 + FMA, **no AVX-512**: WHFast512 runs on the portable C kernels (AVX2 or generic), not the assembly kernel |
| memory / OS | 125 GB, Debian 12, gcc 12.2.0, Python 3.12.2 + numpy 2.3.5 |
| repo | git clone at `/home/michael/Downloads/rebound/rebound` (branch `whfaster`); the non-git folder `/home/michael/Downloads/rebound/whfaster/` is an older copy of the documents only |
| load | interactive desktop, not idle |

Use it for building and testing, Python prototypes and analysis, and runs with few threads. Do not compare its
timings with host2: the cores are mixed (P vs E have different speeds, which skews any per-slice load balance),
the fine kernel is a different code path, and the cost ratio r = coarse/fine is therefore different. If you do time
here, pin to P-cores and record it: CPUs 0–11 are the P-core threads (siblings 0/1, 2/3, …), 12–19 the E-cores;
one thread per physical P-core = `taskset -c 0,2,4,6,8,10`.
