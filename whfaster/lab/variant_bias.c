// One lab variant of the C kernel, AVX-512 ops. Compile with -DSUFFIX=name and variant macros.
#include <stdint.h>
#include <string.h>
#include <signal.h>
#include "integrator_whfast512_data.h"
#pragma GCC target("avx512f,avx512dq,avx2,fma")
#include <immintrin.h>
#define REB_W512_TARGET
#define W_INL static inline __attribute__((always_inline))
#define REB_W512_NOINLINE_FALLBACK static __attribute__((noinline))
#include "ops_avx512.h"
#define REB_W512_SUFFIX SUFFIX
#include "kernel_bias.h"
