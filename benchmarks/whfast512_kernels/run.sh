#!/usr/bin/env bash
# WHFast512 kernel and MEP benchmark (Parachron). Builds several variants of the library from this checkout, then
# measures on this machine:
#   1. single-core ns/step of every kernel variant (paper setup), with state hashes (bitwise checks);
#   2. hardware counters per step (if `perf` works);
#   3. per-core speed with 16 (and all) physical cores busy (clock under all-core load);
#   4. MEP wall clock and speed-up against the serial reference, for several kernels and core counts.
# Usage (from anywhere):  benchmarks/whfast512_kernels/run.sh
# Environment: OUT (output directory), CORE (logical CPU for single-core runs), REPS (repetitions, default 3),
#              SKIP_UPSTREAM=1 (do not clone upstream REBOUND), SKIP_MEP=1.
# Results: $OUT/results.txt (send this file back), raw data in $OUT/raw/.
set -u
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$HERE/../.." && pwd)
OUT=${OUT:-$PWD/w512bench_$(hostname -s)_$(date +%Y%m%d_%H%M)}
REPS=${REPS:-3}
CALLS=20                                  # 20 calls x 1e6 steps per single-core run
UPSTREAM_COMMIT=262a9cce102c69ca2cf7f3c6c1264e5e60002b7d   # Parachron's upstream base
mkdir -p "$OUT/raw" "$OUT/build"
exec > >(tee "$OUT/results.txt") 2>&1
JOBS=$(nproc)

say(){ printf '\n=== %s\n' "$*"; }

# ------------------------------------------------------------------------------------------------ machine
say "machine"
date -u +"date %Y-%m-%d %H:%M UTC"
grep -m1 'model name' /proc/cpuinfo
lscpu | grep -E '^(Socket|Core|Thread|CPU\(s\)|NUMA node\(s\)|CPU max MHz|L3)' || true
echo "AVX-512 flags: $(grep -m1 -o -w -E 'avx512[a-z0-9_]*' /proc/cpuinfo | sort -u | tr '\n' ' ')"
grep -m1 microcode /proc/cpuinfo || true
echo "kernel $(uname -r); $(. /etc/os-release 2>/dev/null; echo "$PRETTY_NAME")"
echo "governor $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo n/a); boost $(cat /sys/devices/system/cpu/cpufreq/boost 2>/dev/null || echo n/a); SMT $(cat /sys/devices/system/cpu/smt/active 2>/dev/null || echo n/a)"
echo "load average: $(cut -d' ' -f1-3 /proc/loadavg)"
gcc --version | head -1
command -v clang >/dev/null && clang --version | head -1
echo "parachron $(git -C "$REPO" rev-parse --short HEAD 2>/dev/null) ($(git -C "$REPO" rev-parse --abbrev-ref HEAD 2>/dev/null))"

# physical cores (first hardware thread of each core), in CPU order
PHYS=()
for d in /sys/devices/system/cpu/cpu[0-9]*; do
    c=${d##*cpu}; f=$(cut -d, -f1 "$d/topology/thread_siblings_list" 2>/dev/null | cut -d- -f1)
    [ "$f" = "$c" ] && PHYS+=("$c")
done
IFS=$'\n' PHYS=($(sort -n <<<"${PHYS[*]}")); unset IFS
NPHYS=${#PHYS[@]}
CORE=${CORE:-${PHYS[2]:-0}}
echo "physical cores: $NPHYS; single-core runs on CPU $CORE"

# ------------------------------------------------------------------------------------------------ builds
say "building library variants (a few minutes)"
build(){   # name, extra OPT flags, compiler
    local b="$OUT/build/$1"
    rm -rf "$b"; mkdir -p "$b"; cp -r "$REPO/src" "$b/src"
    rm -f "$b"/src/*.o "$b"/src/*.asm_o "$b"/src/*.so
    [ -n "$2" ] && echo "OPT+= $2" >> "$b/src/Makefile.defs"
    if ! make -C "$b/src" -j"$JOBS" CC="$3" > "$b/build.log" 2>&1; then echo "  $1: BUILD FAILED (see $b/build.log)"; return 1; fi
    for t in kbench mepbench; do
        "$3" -O2 -DPARACHRON "$HERE/$t.c" -I"$b/src" -L"$b/src" -lrebound -lm -lpthread -Wl,-rpath,"$b/src" -o "$b/$t" || return 1
    done
    echo "  $1: ok ($2, $3)"
}
VARIANTS=()
build default ""                    gcc && VARIANTS+=(default)
build hwrsqrt "-DREB_W512_HW_RSQRT" gcc && VARIANTS+=(hwrsqrt)
build sqrtdiv "-DREB_W512_SQRTDIV"  gcc && VARIANTS+=(sqrtdiv)
build native  "-march=native"       gcc && VARIANTS+=(native)
build native_sqrtdiv "-march=native -DREB_W512_SQRTDIV" gcc && VARIANTS+=(native_sqrtdiv)
if command -v clang >/dev/null; then
    build clang "" clang && VARIANTS+=(clang)
    build clang_sqrtdiv "-DREB_W512_SQRTDIV" clang && VARIANTS+=(clang_sqrtdiv)
fi
HAVE_UP=0
if [ "${SKIP_UPSTREAM:-0}" != 1 ]; then
    u="$OUT/build/upstream"
    if git clone -q https://github.com/hannorein/rebound "$u" 2>/dev/null && git -C "$u" checkout -q "$UPSTREAM_COMMIT" 2>/dev/null \
       && make -C "$u/src" -j"$JOBS" > "$u/build.log" 2>&1 \
       && gcc -O2 "$HERE/kbench.c" -I"$u/src" -L"$u/src" -lrebound -lm -lpthread -Wl,-rpath,"$u/src" -o "$u/kbench"; then
        HAVE_UP=1; echo "  upstream REBOUND ${UPSTREAM_COMMIT:0:8}: ok"
    else
        echo "  upstream REBOUND: not available (clone or build failed); skipped"
    fi
fi

# ------------------------------------------------------------------------------------------------ 1. single core
# Runs: "variant kernel". Kernels: 1 assembly (unchanged from upstream), 2 C AVX-512, 3 C AVX2 (8 lanes as 2x4).
RUNS=()
[ $HAVE_UP = 1 ] && RUNS+=("upstream 1")
for v in "${VARIANTS[@]}"; do
    case $v in
        default) RUNS+=("default 1" "default 2" "default 3") ;;
        *)       RUNS+=("$v 2" "$v 3") ;;
    esac
done
say "1. single core: ${#RUNS[@]} runs x $REPS repetitions x $CALLS million steps, interleaved, CPU $CORE"
for rep in $(seq "$REPS"); do
    for run in "${RUNS[@]}"; do
        set -- $run
        line=$(taskset -c "$CORE" "$OUT/build/$1/kbench" "$2" "$CALLS" 2>&1 | tail -1)
        echo "$1 $2 $rep $line" >> "$OUT/raw/single.txt"
    done
    echo "  repetition $rep done"
done
# table: median over repetitions of the per-run median ns/step
awk '
{ key=$1" "$2; n[key]++; v[key,n[key]]=$4; used[key]=$7; hash[key]=$8; if(!(key in seen)){seen[key]=1; order[++m]=key} }
function med(key,   i,j,t,a,c){ c=n[key]; for(i=1;i<=c;i++)a[i]=v[key,i]; for(i=1;i<=c;i++)for(j=i+1;j<=c;j++)if(a[j]<a[i]){t=a[i];a[i]=a[j];a[j]=t}; return (c%2)?a[(c+1)/2]:(a[c/2]+a[c/2+1])/2 }
END{
  ref=med("default 1");
  printf "%-16s %-7s %-10s %-9s %-6s %s\n","variant","kernel","ns/step","vs asm","used","state hash";
  for(i=1;i<=m;i++){k=order[i]; split(k,p," "); x=med(k);
    printf "%-16s %-7s %-10.2f %-9.3f %-6s %s\n",p[1],p[2],x,ref/x,used[k],hash[k]}
}' "$OUT/raw/single.txt" | tee "$OUT/raw/single_table.txt"
echo "(vs asm > 1: faster than the assembly kernel of this build; kernel 1 = assembly, 2 = C AVX-512, 3 = C AVX2)"
h(){ awk -v k="$1 $2" '$1" "$2==k{print $8; exit}' "$OUT/raw/single.txt"; }
chk(){ [ -n "$(h $1 $2)" ] && [ -n "$(h $3 $4)" ] && { [ "$(h $1 $2)" = "$(h $3 $4)" ] && r=identical || r=DIFFERENT; echo "  $1/$2 vs $3/$4: $r"; }; }
echo "bitwise checks:"
[ $HAVE_UP = 1 ] && chk upstream 1 default 1
chk default 2 default 3
chk sqrtdiv 2 sqrtdiv 3
chk default 2 native 2
chk sqrtdiv 2 native_sqrtdiv 3
chk sqrtdiv 2 clang_sqrtdiv 2

# ------------------------------------------------------------------------------------------------ 2. counters
say "2. hardware counters per step (perf)"
PERF_RUNS=("default 1" "default 2" "default 3" "sqrtdiv 2" "sqrtdiv 3" "native_sqrtdiv 3")
if command -v perf >/dev/null && perf stat -e cycles true >/dev/null 2>&1; then
    P=fp_ops_retired_by_type
    GROUPS_=("cycles,instructions")
    if perf list 2>/dev/null | grep -q "$P.vector_mac"; then
        GROUPS_+=("$P.vector_mac,$P.vector_mul,$P.vector_add,$P.vector_sub,$P.vector_div,$P.vector_sqrt"
                  "$P.vector_shuffle,$P.vector_blend,$P.vector_cmp,$P.vector_logical,$P.vector_other")
    fi
    for run in "${PERF_RUNS[@]}"; do
        set -- $run; [ -x "$OUT/build/$1/kbench" ] || continue
        for g in "${GROUPS_[@]}"; do
            taskset -c "$CORE" perf stat -x, -e "$g" "$OUT/build/$1/kbench" "$2" 10 2>&1 >/dev/null \
                | awk -F, -v r="$1 $2" '/^[0-9]/{split($3,e,"."); printf "%-12s %-24s %9.1f\n", r, (e[2]?e[2]:e[1]), $1/1e7}'
        done
    done | tee "$OUT/raw/perf.txt"
else
    echo "perf not available or not permitted (install linux-perf; sudo sysctl -w kernel.perf_event_paranoid=1); skipped"
fi

# ------------------------------------------------------------------------------------------------ 3. all cores
say "3. per-core speed with many cores busy (one copy per physical core, 10 million steps each)"
LOADED=("default 1" "default 2" "sqrtdiv 2" "sqrtdiv 3" "native_sqrtdiv 2" "native_sqrtdiv 3")
for P in 16 $NPHYS; do
    [ "$P" -le "$NPHYS" ] || continue
    for run in "${LOADED[@]}"; do
        set -- $run; [ -x "$OUT/build/$1/kbench" ] || continue
        rm -f "$OUT/raw/load_$$_"*
        for i in $(seq 0 $((P-1))); do
            taskset -c "${PHYS[$i]}" "$OUT/build/$1/kbench" "$2" 10 > "$OUT/raw/load_$$_$i" 2>&1 &
        done
        wait
        cut -d' ' -f1 "$OUT/raw/load_$$_"* | sort -n | awk -v r="$1 $2" -v P="$P" '{a[NR]=$1} END{n=NR; printf "  %3d cores busy  %-12s median %.2f ns/step (fastest %.2f, slowest %.2f)\n", P, r, a[int((n+1)/2)], a[1], a[n]}'
        rm -f "$OUT/raw/load_$$_"*
    done
    [ "$P" = "$NPHYS" ] && break
done | tee "$OUT/raw/allcore.txt"

# ------------------------------------------------------------------------------------------------ 4. MEP
if [ "${SKIP_MEP:-0}" != 1 ]; then
    say "4. MEP: 1600 slices of 60000 steps (1.58 Myr); serial reference = one WHFast512 call per slice on CPU $CORE"
    printf "%-12s %-6s %-9s %-8s %-9s %-7s\n" "variant/k" "cores" "tol" "wall [s]" "speed-up" "K mean" | tee "$OUT/raw/mep.txt"
    mep_case(){ # variant kernel serial_time cores tol
        local line; line=$("$OUT/build/$1/mepbench" mep "$2" "$4" "$5" 2>&1 | tail -1)
        set -- "$1" "$2" "$3" "$4" "$5" $line
        printf "%-12s %-6s %-9s %-8.2f %-9.2f %-7s\n" "$1/$2" "$4" "$5" "$6" "$(awk -v a="$3" -v b="$6" 'BEGIN{print a/b}')" "$7"
    }
    for run in "default 1" "default 2" "sqrtdiv 2" "sqrtdiv 3" "native_sqrtdiv 3"; do
        set -- $run; [ -x "$OUT/build/$1/mepbench" ] || continue
        ser=$(taskset -c "$CORE" "$OUT/build/$1/mepbench" serial "$2" | cut -d' ' -f1)
        printf "%-12s %-6s %-9s %-8.2f %-9s %-7s\n" "$1/$2" 1 serial "$ser" 1.00 - | tee -a "$OUT/raw/mep.txt"
        mep_case "$1" "$2" "$ser" 16 3e-6 | tee -a "$OUT/raw/mep.txt"
        if [ "$run" = "default 1" ]; then
            mep_case "$1" "$2" "$ser" 16 3e-5 | tee -a "$OUT/raw/mep.txt"
            for N in 32 64; do [ "$N" -le "$NPHYS" ] && mep_case "$1" "$2" "$ser" "$N" 3e-6 | tee -a "$OUT/raw/mep.txt"; done
        fi
    done
fi

say "done: send $OUT/results.txt"
