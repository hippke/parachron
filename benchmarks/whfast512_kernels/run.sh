#!/usr/bin/env bash
# WHFast512 kernel and MEP benchmark (Parachron). Builds several variants of the library from this checkout, then
# measures on this machine:
#   1. single-core speed of every kernel variant relative to the assembly kernel (paper setup), with state hashes;
#   2. hardware counters per step (if `perf` works);
#   3. per-core speed with 16 copies running at the same time (clock under multi-core load);
#   4. MEP wall clock and speed-up against the serial reference, for several kernels and core counts.
# Made for machines with other load: every comparison runs in ABBA blocks (reference, variant, variant, reference),
# so that drifts of the background load cancel in each block's ratio, and the runs are pinned to the physical cores
# that were the most idle (both hardware threads) just before. The spread over blocks shows the remaining noise;
# a null test (the assembly against itself) gives its size.
# Usage (from anywhere):  benchmarks/whfast512_kernels/run.sh
# Environment: OUT (output directory), BLOCKS (ABBA blocks per single-core comparison, default 5), CORE (fixed CPU for
#              single-core runs instead of the most idle one), SKIP_UPSTREAM=1 (do not clone upstream REBOUND), SKIP_MEP=1.
# Results: $OUT/results.txt (send this file back), raw data in $OUT/raw/.
set -u
export LC_ALL=C                           # decimal points in printf/awk/sort, whatever the system locale
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$HERE/../.." && pwd)
OUT=${OUT:-$PWD/w512bench_$(hostname -s)_$(date +%Y%m%d_%H%M)}
BLOCKS=${BLOCKS:-5}
CALLS=5                                   # 5 calls x 1e6 steps per single-core run (~0.6 s)
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
gcc --version | head -1
command -v clang >/dev/null && clang --version | head -1
echo "parachron $(git -C "$REPO" rev-parse --short HEAD 2>/dev/null) ($(git -C "$REPO" rev-parse --abbrev-ref HEAD 2>/dev/null))"

# physical cores (first hardware thread of each core) and their hardware threads
PHYS=()
declare -A SIBS
for d in /sys/devices/system/cpu/cpu[0-9]*; do
    c=${d##*cpu}; l=$(cat "$d/topology/thread_siblings_list" 2>/dev/null || echo "$c")
    f=${l%%[,-]*}
    if [ "$f" = "$c" ]; then
        PHYS+=("$c"); s=""
        for part in ${l//,/ }; do
            if [[ $part == *-* ]]; then s+=" $(seq -s ' ' "${part%-*}" "${part#*-}")"; else s+=" $part"; fi
        done
        SIBS[$c]=${s# }
    fi
done
IFS=$'\n' PHYS=($(sort -n <<<"${PHYS[*]}")); unset IFS
NPHYS=${#PHYS[@]}
echo "physical cores: $NPHYS"

# Idle fraction of every physical core over 2 s (the smaller of its hardware threads), most idle first:
# lines "first_cpu idle".
rank_cores(){
    local a b; a=$(grep '^cpu[0-9]' /proc/stat); sleep 2; b=$(grep '^cpu[0-9]' /proc/stat)
    paste <(echo "$a") <(echo "$b") | awk '{n=NF/2; t1=0; t2=0; for(i=2;i<=n;i++){t1+=$i; t2+=$(n+i)};
        cpu=substr($1,4); idle=(t2>t1)?(($(n+5)+$(n+6))-($5+$6))/(t2-t1):0; print cpu, idle}' > "$OUT/raw/idle_now"
    for p in "${PHYS[@]}"; do
        awk -v s=" ${SIBS[$p]} " -v p="$p" 'BEGIN{m=1} index(s," "$1" ")>0{if($2<m)m=$2} END{printf "%s %.3f\n", p, m}' "$OUT/raw/idle_now"
    done | sort -k2,2nr -k1,1n
}
# pick N: the N most idle physical cores. Sets PICK (first threads, comma separated), PICK_ALL (with SMT siblings)
# and PICK_IDLE (smallest idle fraction among them).
pick(){
    local r; r=$(rank_cores | head -n "$1")
    PICK=$(cut -d' ' -f1 <<<"$r" | paste -sd,)
    PICK_ALL=$(for p in $(cut -d' ' -f1 <<<"$r"); do echo ${SIBS[$p]}; done | tr ' ' '\n' | paste -sd,)
    PICK_IDLE=$(tail -1 <<<"$r" | cut -d' ' -f2)
}
echo "load average $(cut -d' ' -f1-3 /proc/loadavg); idle physical cores (>= 95 % idle, both threads): $(rank_cores | awk '$2>=0.95' | wc -l) of $NPHYS"

# median of the numbers on stdin
med(){ sort -g | awk '{a[NR]=$1} END{if(NR==0){print "nan"} else if(NR%2){print a[(NR+1)/2]} else {print (a[NR/2]+a[NR/2+1])/2}}'; }

# ------------------------------------------------------------------------------------------------ builds
say "building library variants (a few minutes)"
build(){   # name, extra OPT flags, compiler
    local b="$OUT/build/$1"
    rm -rf "$b"; mkdir -p "$b"; cp -r "$REPO/src" "$b/src"
    rm -f "$b"/src/*.o "$b"/src/*.asm_o "$b"/src/*.so
    [ -n "$2" ] && echo "OPT+= $2" >> "$b/src/Makefile.defs"
    if ! make -C "$b/src" -j"$JOBS" CC="$3" > "$b/build.log" 2>&1; then echo "  $1: BUILD FAILED (see $b/build.log)"; return 1; fi
    for t in kbench mepbench pbench; do
        "$3" -O2 -DPARACHRON "$HERE/$t.c" -I"$b/src" -L"$b/src" -lrebound -lm -lpthread -Wl,-rpath,"$b/src" -o "$b/$t" || return 1
    done
    echo "  $1: ok ($2, $3)"
}
VARIANTS=()
build default ""                    gcc && VARIANTS+=(default)
build hwrsqrt "-DREB_W512_HW_RSQRT" gcc && VARIANTS+=(hwrsqrt)
build horner  "-DREB_W512_STIEFEL_HORNER" gcc && VARIANTS+=(horner)
build knewton "-DREB_W512_KEPLER_NEWTON" gcc && VARIANTS+=(knewton)
build mat8    "-DREB_W512_MAT8" gcc && VARIANTS+=(mat8)
build mix0    "-DREB_W512_RSQRT_MIX=0" gcc && VARIANTS+=(mix0)
build feval   "-DREB_W512_FINAL_EVAL" gcc && VARIANTS+=(feval)
build novl    "-DREB_W512_NO_AVX2VL" gcc && VARIANTS+=(novl)
build newton  "-DREB_W512_NEWTON_RSQRT" gcc && VARIANTS+=(newton)
build native  "-march=native"       gcc && VARIANTS+=(native)
if command -v clang >/dev/null; then
    build clang "" clang && VARIANTS+=(clang)
    build clang_newton "-DREB_W512_NEWTON_RSQRT" clang && VARIANTS+=(clang_newton)
fi
HAVE_UP=0
if [ "${SKIP_UPSTREAM:-0}" != 1 ]; then
    u="$OUT/build/upstream"
    fetch_upstream(){   # source archive of the commit (curl or wget); git clone without credential prompts as fallback
        local url="https://github.com/hannorein/rebound/archive/$UPSTREAM_COMMIT.tar.gz"
        rm -rf "$u"; mkdir -p "$u"
        { curl -fsSL "$url" 2>/dev/null || wget -qO- "$url" 2>/dev/null; } | tar xz -C "$u" --strip-components=1 2>/dev/null \
            && [ -f "$u/src/Makefile" ] && return 0
        rm -rf "$u"
        GIT_TERMINAL_PROMPT=0 GIT_ASKPASS=/bin/false git -c credential.helper= clone -q https://github.com/hannorein/rebound "$u" < /dev/null 2>/dev/null \
            && git -C "$u" checkout -q "$UPSTREAM_COMMIT" 2>/dev/null
    }
    if fetch_upstream \
       && make -C "$u/src" -j"$JOBS" > "$u/build.log" 2>&1 \
       && gcc -O2 "$HERE/kbench.c" -I"$u/src" -L"$u/src" -lrebound -lm -lpthread -Wl,-rpath,"$u/src" -o "$u/kbench"; then
        HAVE_UP=1; echo "  upstream REBOUND ${UPSTREAM_COMMIT:0:8}: ok"
    else
        echo "  upstream REBOUND: not available (clone or build failed); skipped"
    fi
fi
[ -x "$OUT/build/default/kbench" ] || { echo "the default build failed; stopping"; exit 1; }

# ------------------------------------------------------------------------------------------------ 1. single core
# Tests "variant kernel". Kernels: 1 assembly (unchanged from upstream), 2 C AVX-512, 3 C AVX2 (8 lanes as 2x4).
# Reference A = default build, assembly kernel. "default 1" against itself is the null test.
TESTS=("default 1")
[ $HAVE_UP = 1 ] && TESTS+=("upstream 1")
for v in "${VARIANTS[@]}"; do
    case $v in
        default) TESTS+=("default 2" "default 3") ;;
        *)       TESTS+=("$v 2" "$v 3") ;;
    esac
done
say "1. single core, ABBA against the assembly kernel: ${#TESTS[@]} tests x $BLOCKS blocks x 4 runs of $CALLS million steps"
kb(){ taskset -c "$1" "$OUT/build/$2/kbench" "$3" "$CALLS" 2>&1 | tail -1; }   # cpu variant kernel
for b in $(seq "$BLOCKS"); do
    if [ -n "${CORE:-}" ]; then C1=$CORE; I1=fixed; else pick 1; C1=$PICK; I1=$PICK_IDLE; fi
    for t in "${TESTS[@]}"; do
        set -- $t
        for role in A X X A; do
            if [ $role = A ]; then line=$(kb "$C1" default 1); else line=$(kb "$C1" "$1" "$2"); fi
            echo "$b $1 $2 $role $line" >> "$OUT/raw/single.txt"
        done
    done
    echo "  block $b done (CPU $C1, idle before $I1, load average $(cut -d' ' -f1 /proc/loadavg))"
done
# per test: ratio per block = (A1 + A2) / (X1 + X2) of the per-run median ns/step; median and range over blocks
{
printf "%-16s %-7s %-9s %-9s %-17s %-6s %s\n" "variant" "kernel" "ns/step" "vs asm" "range (blocks)" "used" "state hash"
for t in "${TESTS[@]}"; do
    set -- $t
    awk -v v="$1" -v k="$2" '$2==v && $3==k {if($4=="A")a[$1]+=$5; else x[$1]+=$5} END{for(b in a) if(x[b]>0) print a[b]/x[b]}' "$OUT/raw/single.txt" > "$OUT/raw/r"
    r=$(med < "$OUT/raw/r"); lo=$(sort -g "$OUT/raw/r" | head -1); hi=$(sort -g "$OUT/raw/r" | tail -1)
    ns=$(awk -v v="$1" -v k="$2" '$2==v && $3==k && $4=="X"{print $5}' "$OUT/raw/single.txt" | med)
    last=$(awk -v v="$1" -v k="$2" '$2==v && $3==k && $4=="X"{u=$8; h=$9} END{print u, h}' "$OUT/raw/single.txt")
    name=$1; [ "$t" = "default 1" ] && name="default(null)"
    printf "%-16s %-7s %-9.2f %-9.3f %-17s %-6s %s\n" "$name" "$2" "$ns" "$r" "$(printf '%.3f-%.3f' "$lo" "$hi")" $last
done
echo "(vs asm: speed relative to the assembly kernel, median over blocks; > 1 is faster. ns/step: median of all runs of the"
echo " variant, affected by the background load. kernel 1 = assembly, 2 = C AVX-512, 3 = C AVX2. null: assembly vs itself.)"
} | tee "$OUT/raw/single_table.txt"
h(){ awk -v v="$1" -v k="$2" '$2==v && $3==k && $4=="X"{print $9; exit}' "$OUT/raw/single.txt"; }
chk(){ [ -n "$(h $1 $2)" ] && [ -n "$(h $3 $4)" ] && { [ "$(h $1 $2)" = "$(h $3 $4)" ] && r=identical || r=DIFFERENT; echo "  $1/$2 vs $3/$4: $r"; }; }
echo "bitwise checks:"
[ $HAVE_UP = 1 ] && chk upstream 1 default 1
chk default 2 default 3
chk newton 2 newton 3
chk default 2 native 2
chk default 2 clang 2
chk default 3 novl 3
chk newton 2 clang_newton 2

# ------------------------------------------------------------------------------------------------ 2. counters
say "2. hardware counters per step (perf; counts do not depend on the background load, cycles a little)"
PERF_RUNS=("default 1" "default 2" "default 3" "novl 3" "newton 2" "native 3")
if command -v perf >/dev/null && perf stat -e cycles true >/dev/null 2>&1; then
    P=fp_ops_retired_by_type
    GROUPS_=("cycles,instructions")
    if perf list 2>/dev/null | grep -q "$P.vector_mac"; then
        GROUPS_+=("$P.vector_mac,$P.vector_mul,$P.vector_add,$P.vector_sub,$P.vector_div,$P.vector_sqrt"
                  "$P.vector_shuffle,$P.vector_blend,$P.vector_cmp,$P.vector_logical,$P.vector_other")
    fi
    if [ -n "${CORE:-}" ]; then C1=$CORE; else pick 1; C1=$PICK; fi
    for run in "${PERF_RUNS[@]}"; do
        set -- $run; [ -x "$OUT/build/$1/kbench" ] || continue
        for g in "${GROUPS_[@]}"; do
            taskset -c "$C1" perf stat -x, -e "$g" "$OUT/build/$1/kbench" "$2" 10 2>&1 >/dev/null \
                | awk -F, -v r="$1 $2" '/^[0-9]/{split($3,e,"."); printf "%-18s %-24s %9.1f\n", r, (e[2]?e[2]:e[1]), $1/1e7}'
        done
    done | tee "$OUT/raw/perf.txt"
else
    echo "perf not available or not permitted (install linux-perf; sudo sysctl -w kernel.perf_event_paranoid=1); skipped"
fi

# ------------------------------------------------------------------------------------------------ 3. 16 cores busy
P16=16; [ "$NPHYS" -lt 16 ] && P16=$NPHYS
say "3. per-core speed with $P16 copies at the same time (most idle cores), ABBA against the assembly kernel, 3 blocks"
LOADED=("default 1" "default 2" "default 3" "novl 3" "newton 2" "native 3")
multi(){ # variant kernel -> median ns/step over the copies (cores in $PICK)
    local i=0
    rm -f "$OUT/raw/load_$$_"*
    for c in ${PICK//,/ }; do taskset -c "$c" "$OUT/build/$1/kbench" "$2" "$CALLS" > "$OUT/raw/load_$$_$i" 2>&1 & i=$((i+1)); done
    wait
    cut -d' ' -f1 "$OUT/raw/load_$$_"* | med
    rm -f "$OUT/raw/load_$$_"*
}
for b in 1 2 3; do
    pick "$P16"
    for t in "${LOADED[@]}"; do
        set -- $t; [ -x "$OUT/build/$1/kbench" ] || continue
        for role in A X X A; do
            if [ $role = A ]; then x=$(multi default 1); else x=$(multi "$1" "$2"); fi
            echo "$b $1 $2 $role $x" >> "$OUT/raw/multi.txt"
        done
    done
    echo "  block $b done (cores $PICK, least idle before ${PICK_IDLE})"
done
{
printf "%-18s %-9s %-9s %s\n" "variant/kernel" "ns/step" "vs asm" "range (blocks)"
for t in "${LOADED[@]}"; do
    set -- $t
    awk -v v="$1" -v k="$2" '$2==v && $3==k {if($4=="A")a[$1]+=$5; else x[$1]+=$5} END{for(b in a) if(x[b]>0) print a[b]/x[b]}' "$OUT/raw/multi.txt" > "$OUT/raw/r"
    [ -s "$OUT/raw/r" ] || continue
    r=$(med < "$OUT/raw/r"); lo=$(sort -g "$OUT/raw/r" | head -1); hi=$(sort -g "$OUT/raw/r" | tail -1)
    ns=$(awk -v v="$1" -v k="$2" '$2==v && $3==k && $4=="X"{print $5}' "$OUT/raw/multi.txt" | med)
    name="$1/$2"; [ "$t" = "default 1" ] && name="default/1 (null)"
    printf "%-18s %-9.2f %-9.3f %s\n" "$name" "$ns" "$r" "$(printf '%.3f-%.3f' "$lo" "$hi")"
done
} | tee "$OUT/raw/multi_table.txt"

# ------------------------------------------------------------------------------------------------ 4. MEP
if [ "${SKIP_MEP:-0}" != 1 ]; then
    say "4. MEP: 1600 slices of 60000 steps (1.58 Myr), serial reference = one WHFast512 call per slice on one core"
    echo "Palindromic order per kernel (serial, MEP cases, MEP cases reversed, serial); speed-up = mean serial / mean MEP."
    ser(){ pick 1; taskset -c "$PICK" "$OUT/build/$1/mepbench" serial "$2" | cut -d' ' -f1; }
    mep(){ # variant kernel cores tol -> "wall kmean idle" (or "skip")
        [ "$3" -le "$NPHYS" ] || { echo skip; return; }
        pick "$3"
        echo "$(taskset -c "$PICK_ALL" "$OUT/build/$1/mepbench" mep "$2" "$3" "$4" 2>&1 | tail -1 | cut -d' ' -f1,2) $PICK_IDLE"
    }
    printf "%-20s %-6s %-7s %-9s %-9s %-7s %s\n" "variant/kernel" "cores" "tol" "wall [s]" "speed-up" "K mean" "least idle core" | tee "$OUT/raw/mep.txt"
    for run in "default 1" "default 2" "default 3" "newton 2" "native 3"; do
        set -- $run; [ -x "$OUT/build/$1/mepbench" ] || continue
        if [ "$run" = "default 1" ]; then CASES=("16 3e-6" "16 3e-5" "32 3e-6" "64 3e-6"); else CASES=("16 3e-6"); fi
        unset W K I; declare -A W K I
        s1=$(ser "$1" "$2")
        for c in "${CASES[@]}"; do set -- $run $c; r=$(mep "$1" "$2" "$3" "$4"); read -r w k i <<<"$r"; W[$c]=$w; K[$c]=${k:-}; I[$c]=${i:-}; done
        for (( j=${#CASES[@]}-1; j>=0; j-- )); do c=${CASES[$j]}; set -- $run $c; r=$(mep "$1" "$2" "$3" "$4"); W[$c]+=" ${r%% *}"; done
        set -- $run
        s2=$(ser "$1" "$2")
        sm=$(awk -v a="$s1" -v b="$s2" 'BEGIN{print (a+b)/2}')
        printf "%-20s %-6s %-7s %-9.2f %-9s %-7s %s\n" "$1/$2" 1 serial "$sm" "1.00" "-" "($s1, $s2 s)" | tee -a "$OUT/raw/mep.txt"
        for c in "${CASES[@]}"; do
            set -- $c
            if [[ ${W[$c]} == skip* ]]; then printf "%-20s %-6s %-7s skipped (fewer physical cores)\n" "${run/ //}" "$1" "$2" | tee -a "$OUT/raw/mep.txt"; continue; fi
            wm=$(awk '{print ($1+$2)/2}' <<<"${W[$c]}")
            printf "%-20s %-6s %-7s %-9.2f %-9.2f %-7s %s\n" "${run/ //}" "$1" "$2" "$wm" "$(awk -v a="$sm" -v b="$wm" 'BEGIN{print a/b}')" "${K[$c]}" "${I[$c]} (runs ${W[$c]} s)" | tee -a "$OUT/raw/mep.txt"
        done
    done
    echo "(least idle core: idle fraction, before the first MEP run, of the least idle core used; < 0.9 means MEP shared"
    echo " cores with other jobs and its speed-up is underestimated)"
fi

# ------------------------------------------------------------------------------------------------ 5. paired kernel
# Two independent simulations in one instruction stream (reb_integrator_whfast512_steps_pair; C AVX-512 kernel only)
# against the same two stepped one after the other: throughput for ensembles, or two MEP slices per thread.
say "5. paired kernel (two simulations in one instruction stream), C AVX-512 kernel, ABBA rounds"
if [ -n "${CORE:-}" ]; then C1=$CORE; else pick 1; C1=$PICK; fi
for v in default horner; do
    [ -x "$OUT/build/$v/pbench" ] || continue
    r=$(taskset -c "$C1" "$OUT/build/$v/pbench" 2 7 2>&1 | tail -1); read -r s p q k same <<<"$r"
    if [ "${k:-}" != 2 ]; then echo "$v: kernel 2 not available on this CPU (used $k); skipped"; continue; fi
    printf "%-8s single %s ns/step, paired %s ns/step per simulation, ratio %s (< 1: pairing pays), results %s\n" "$v" "$s" "$p" "$q" "$same" | tee -a "$OUT/raw/pair.txt"
done

echo "load average at the end: $(cut -d' ' -f1-3 /proc/loadavg)"
say "done: send $OUT/results.txt"
