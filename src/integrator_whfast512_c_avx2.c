/**
 * @file    integrator_whfast512_c_avx2.c
 * @brief   WHFast512 portable C kernel, AVX2+FMA instantiation.
 * @details Defines the vector operations of integrator_whfast512_kernel.h with AVX2 and FMA
 *          intrinsics. An 8-lane vector is a pair of 256-bit registers (lanes 0-3 and 4-7). The
 *          functions are compiled for AVX2+FMA through a target pragma, so the file needs no
 *          special compiler flags; the driver calls them only on CPUs that support AVX2 and FMA.
 *
 * @section LICENSE
 * Copyright (c) 2026 Hanno Rein
 * Copyright (c) 2026 Michael Hippke (portable C kernels, derived from the WHFast512 assembly)
 * Part of Parachron (https://github.com/hippke/parachron), a modified version of REBOUND.
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include <stdint.h>
#include <string.h>
#include <signal.h>
#include "integrator_whfast512_data.h"

#if (defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)) && (defined(__clang__) || (defined(__GNUC__) && __GNUC__ >= 5) || defined(_MSC_VER))

#if defined(__clang__)
#pragma clang attribute push (__attribute__((target("avx2,fma"))), apply_to = function)
#elif defined(__GNUC__)
#pragma GCC push_options
#pragma GCC target("avx2,fma")
#endif
#define REB_W512_TARGET
#include <immintrin.h>

typedef struct { __m256d l, h; } V;   // lanes 0-3, 4-7
typedef struct { __m256d l, h; } K;   // all bits set in true lanes

#if defined(_MSC_VER) && !defined(__clang__)
#define W_INL static __forceinline
#define REB_W512_NOINLINE_FALLBACK static __declspec(noinline)
#else
#define W_INL static inline __attribute__((always_inline))
#define REB_W512_NOINLINE_FALLBACK static __attribute__((noinline))
#endif

W_INL V w_make(__m256d l, __m256d h){ V r; r.l = l; r.h = h; return r; }
W_INL K w_kmake(__m256d l, __m256d h){ K r; r.l = l; r.h = h; return r; }
W_INL V w_set1(double a){ const __m256d t = _mm256_set1_pd(a); return w_make(t, t); }
W_INL V w_zero(void){ const __m256d t = _mm256_setzero_pd(); return w_make(t, t); }
W_INL V w_load(const double* p){ return w_make(_mm256_load_pd(p), _mm256_load_pd(p+4)); }
W_INL void w_store(double* p, V a){ _mm256_store_pd(p, a.l); _mm256_store_pd(p+4, a.h); }
W_INL V w_add(V a, V b){ return w_make(_mm256_add_pd(a.l, b.l), _mm256_add_pd(a.h, b.h)); }
W_INL V w_sub(V a, V b){ return w_make(_mm256_sub_pd(a.l, b.l), _mm256_sub_pd(a.h, b.h)); }
W_INL V w_mul(V a, V b){ return w_make(_mm256_mul_pd(a.l, b.l), _mm256_mul_pd(a.h, b.h)); }
W_INL V w_div(V a, V b){ return w_make(_mm256_div_pd(a.l, b.l), _mm256_div_pd(a.h, b.h)); }
W_INL V w_sqrt(V a){ return w_make(_mm256_sqrt_pd(a.l), _mm256_sqrt_pd(a.h)); }
W_INL V w_abs(V a){
    const __m256d m = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFFLL));
    return w_make(_mm256_and_pd(a.l, m), _mm256_and_pd(a.h, m));
}
W_INL V w_floor(V a){ return w_make(_mm256_round_pd(a.l, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC), _mm256_round_pd(a.h, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC)); }
W_INL V w_fmadd(V a, V b, V c){ return w_make(_mm256_fmadd_pd(a.l, b.l, c.l), _mm256_fmadd_pd(a.h, b.h, c.h)); }
W_INL V w_fnmadd(V a, V b, V c){ return w_make(_mm256_fnmadd_pd(a.l, b.l, c.l), _mm256_fnmadd_pd(a.h, b.h, c.h)); }
W_INL V w_fmsub(V a, V b, V c){ return w_make(_mm256_fmsub_pd(a.l, b.l, c.l), _mm256_fmsub_pd(a.h, b.h, c.h)); }
W_INL K w_klt(V a, V b){ return w_kmake(_mm256_cmp_pd(a.l, b.l, _CMP_LT_OQ), _mm256_cmp_pd(a.h, b.h, _CMP_LT_OQ)); }
W_INL K w_kgt(V a, V b){ return w_kmake(_mm256_cmp_pd(a.l, b.l, _CMP_GT_OQ), _mm256_cmp_pd(a.h, b.h, _CMP_GT_OQ)); }
W_INL K w_kle(V a, V b){ return w_kmake(_mm256_cmp_pd(a.l, b.l, _CMP_LE_OQ), _mm256_cmp_pd(a.h, b.h, _CMP_LE_OQ)); }
W_INL K w_kand(K a, K b){ return w_kmake(_mm256_and_pd(a.l, b.l), _mm256_and_pd(a.h, b.h)); }
W_INL K w_kandnot(K a, K b){ return w_kmake(_mm256_andnot_pd(a.l, b.l), _mm256_andnot_pd(a.h, b.h)); }
W_INL K w_knot(K a){
    const __m256d ones = _mm256_castsi256_pd(_mm256_set1_epi64x(-1));
    return w_kmake(_mm256_xor_pd(a.l, ones), _mm256_xor_pd(a.h, ones));
}
W_INL int w_kany(K a){ return (_mm256_movemask_pd(a.l) | _mm256_movemask_pd(a.h)) != 0; }
W_INL K w_kall(void){ const __m256d ones = _mm256_castsi256_pd(_mm256_set1_epi64x(-1)); return w_kmake(ones, ones); }
W_INL K w_kbits(unsigned int bits){
    const __m256i b = _mm256_set1_epi64x((long long)bits);
    const __m256i sl = _mm256_set_epi64x(8, 4, 2, 1);
    const __m256i sh = _mm256_set_epi64x(128, 64, 32, 16);
    return w_kmake(_mm256_castsi256_pd(_mm256_cmpeq_epi64(_mm256_and_si256(b, sl), sl)),
                   _mm256_castsi256_pd(_mm256_cmpeq_epi64(_mm256_and_si256(b, sh), sh)));
}
W_INL V w_sel(K k, V a, V b){ return w_make(_mm256_blendv_pd(b.l, a.l, k.l), _mm256_blendv_pd(b.h, a.h, k.h)); }
W_INL V w_maskz(K k, V a){ return w_make(_mm256_and_pd(a.l, k.l), _mm256_and_pd(a.h, k.h)); }
// Bit-trick seed for 1/sqrt(a): integer 0x5FE6EB50C7B537A9 - (bits(a) >> 1) (deterministic, same on every machine).
W_INL __m256d w_rsqrt_seed4(__m256d a){ return _mm256_castsi256_pd(_mm256_sub_epi64(_mm256_set1_epi64x(0x5FE6EB50C7B537A9LL), _mm256_srli_epi64(_mm256_castpd_si256(a), 1))); }
W_INL V w_rsqrt_seed(V a){ return w_make(w_rsqrt_seed4(a.l), w_rsqrt_seed4(a.h)); }
// Within each group of four lanes (same immediates as vpermpd in the assembly).
W_INL V w_perm_3201(V a){ return w_make(_mm256_permute4x64_pd(a.l, 0x4B), _mm256_permute4x64_pd(a.h, 0x4B)); }
W_INL V w_perm_2310(V a){ return w_make(_mm256_permute4x64_pd(a.l, 0x1E), _mm256_permute4x64_pd(a.h, 0x1E)); }
W_INL V w_perm_3012(V a){ return w_make(_mm256_permute4x64_pd(a.l, 0x93), _mm256_permute4x64_pd(a.h, 0x93)); }
W_INL V w_swap_pairs(V a){ return w_make(_mm256_permute_pd(a.l, 0x5), _mm256_permute_pd(a.h, 0x5)); }
// Across the two halves.
W_INL V w_perm_b3(V a){ return w_make(a.h, _mm256_permute4x64_pd(a.l, 0x39)); }                          // 4567 1230
W_INL V w_perm_b4(V a){ return w_make(_mm256_permute4x64_pd(a.h, 0x39), _mm256_permute4x64_pd(a.l, 0x4E)); } // 5674 2301
W_INL V w_perm_b34(V a){ return w_make(_mm256_permute4x64_pd(a.h, 0x93), a.l); }                         // 7456 0123

#define REB_W512_SUFFIX avx2
#include "integrator_whfast512_kernel.h"

#if defined(__clang__)
#pragma clang attribute pop
#elif defined(__GNUC__)
#pragma GCC pop_options
#endif

int reb_whfast512c_avx2_compiled(void){ return 1; }

#else // not x86
int reb_whfast512c_full_steps_avx2(struct simd_data* data, uint64_t* N_steps, int skip_first_kepler_step, volatile sig_atomic_t* sigint, int gr, int nsys, int encounter, int escape){
    (void)data; (void)N_steps; (void)skip_first_kepler_step; (void)sigint; (void)gr; (void)nsys; (void)encounter; (void)escape;
    return 1; // never called (see reb_whfast512c_avx2_compiled)
}
void reb_whfast512c_kepler_step_avx2(struct simd_data* data){ (void)data; }
void reb_whfast512c_corrector_step_avx2(struct simd_data* data, double inv, int gr, int nsys){ (void)data; (void)inv; (void)gr; (void)nsys; }
int reb_whfast512c_avx2_compiled(void){ return 0; }
int reb_whfast512c_full_steps_pair_avx2(struct simd_data* da, struct simd_data* db, uint64_t N_steps, uint64_t* N_done, int skip_a, int skip_b, volatile sig_atomic_t* sigint, int* status_a, int* status_b, int gr, int nsys, int encounter, int escape){
    (void)da; (void)db; (void)N_steps; (void)N_done; (void)skip_a; (void)skip_b; (void)sigint; (void)status_a; (void)status_b; (void)gr; (void)nsys; (void)encounter; (void)escape;
    return -1; // never called
}

#endif
