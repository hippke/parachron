/**
 * @file    integrator_whfast512_c_avx2x4.c
 * @brief   WHFast512 portable C kernel, AVX2+FMA instantiation with four lanes.
 * @details Same kernel as integrator_whfast512_c_avx2.c, but an 8-lane vector is replaced by one
 *          256-bit register holding lanes 0-3. The driver uses it automatically instead of the
 *          8-lane AVX2 kernel for a single system (N_systems = 1) with at most four planets. It
 *          gives the same results (the skipped pair terms only involve empty lanes), at about
 *          half the work, and the state fits in the 16 AVX2 registers.
 *
 * @section LICENSE
 * Copyright (c) 2026 Hanno Rein
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

typedef __m256d V;   // lanes 0-3
typedef __m256d K;   // all bits set in true lanes

#if defined(_MSC_VER) && !defined(__clang__)
#define W_INL static __forceinline
#define REB_W512_NOINLINE_FALLBACK static __declspec(noinline)
#else
#define W_INL static inline __attribute__((always_inline))
#define REB_W512_NOINLINE_FALLBACK static __attribute__((noinline))
#endif

W_INL V w_set1(double a){ return _mm256_set1_pd(a); }
W_INL V w_zero(void){ return _mm256_setzero_pd(); }
W_INL V w_load(const double* p){ return _mm256_load_pd(p); }
W_INL void w_store(double* p, V a){ _mm256_store_pd(p, a); }
W_INL V w_add(V a, V b){ return _mm256_add_pd(a, b); }
W_INL V w_sub(V a, V b){ return _mm256_sub_pd(a, b); }
W_INL V w_mul(V a, V b){ return _mm256_mul_pd(a, b); }
W_INL V w_div(V a, V b){ return _mm256_div_pd(a, b); }
W_INL V w_sqrt(V a){ return _mm256_sqrt_pd(a); }
W_INL V w_abs(V a){ return _mm256_and_pd(a, _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFFLL))); }
W_INL V w_floor(V a){ return _mm256_round_pd(a, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC); }
W_INL V w_fmadd(V a, V b, V c){ return _mm256_fmadd_pd(a, b, c); }
W_INL V w_fnmadd(V a, V b, V c){ return _mm256_fnmadd_pd(a, b, c); }
W_INL V w_fmsub(V a, V b, V c){ return _mm256_fmsub_pd(a, b, c); }
W_INL K w_klt(V a, V b){ return _mm256_cmp_pd(a, b, _CMP_LT_OQ); }
W_INL K w_kgt(V a, V b){ return _mm256_cmp_pd(a, b, _CMP_GT_OQ); }
W_INL K w_kle(V a, V b){ return _mm256_cmp_pd(a, b, _CMP_LE_OQ); }
W_INL K w_kand(K a, K b){ return _mm256_and_pd(a, b); }
W_INL K w_kandnot(K a, K b){ return _mm256_andnot_pd(a, b); }
W_INL K w_knot(K a){ return _mm256_xor_pd(a, _mm256_castsi256_pd(_mm256_set1_epi64x(-1))); }
W_INL int w_kany(K a){ return _mm256_movemask_pd(a) != 0; }
W_INL K w_kall(void){ return _mm256_castsi256_pd(_mm256_set1_epi64x(-1)); }
W_INL K w_kbits(unsigned int bits){
    const __m256i b = _mm256_set1_epi64x((long long)bits);
    const __m256i sl = _mm256_set_epi64x(8, 4, 2, 1);
    return _mm256_castsi256_pd(_mm256_cmpeq_epi64(_mm256_and_si256(b, sl), sl));
}
W_INL V w_sel(K k, V a, V b){ return _mm256_blendv_pd(b, a, k); }
W_INL V w_maskz(K k, V a){ return _mm256_and_pd(a, k); }
W_INL V w_rsqrt_seed(V a){ return _mm256_castsi256_pd(_mm256_sub_epi64(_mm256_set1_epi64x(0x5FE6EB50C7B537A9LL), _mm256_srli_epi64(_mm256_castpd_si256(a), 1))); }
W_INL V w_perm_3201(V a){ return _mm256_permute4x64_pd(a, 0x4B); }
W_INL V w_perm_2310(V a){ return _mm256_permute4x64_pd(a, 0x1E); }
W_INL V w_perm_3012(V a){ return _mm256_permute4x64_pd(a, 0x93); }
W_INL V w_swap_pairs(V a){ return _mm256_permute_pd(a, 0x5); }
// Only used by the N_systems = 1 code path with more than four planets, which this
// instantiation never runs (W_ONLY_NSYS2).
W_INL V w_perm_b3(V a){ return a; }
W_INL V w_perm_b4(V a){ return a; }
W_INL V w_perm_b34(V a){ return a; }

#define W_LANES 4
#define W_ONLY_NSYS2 1
#define REB_W512_SUFFIX avx2x4
#include "integrator_whfast512_kernel.h"

#if defined(__clang__)
#pragma clang attribute pop
#elif defined(__GNUC__)
#pragma GCC pop_options
#endif

#else // not x86
int reb_whfast512c_full_steps_avx2x4(struct simd_data* data, uint64_t* N_steps, int skip_first_kepler_step, volatile sig_atomic_t* sigint, int gr, int nsys, int encounter, int escape){
    (void)data; (void)N_steps; (void)skip_first_kepler_step; (void)sigint; (void)gr; (void)nsys; (void)encounter; (void)escape;
    return 1; // never called (the 8-lane AVX2 kernel is not available either)
}
void reb_whfast512c_kepler_step_avx2x4(struct simd_data* data){ (void)data; }
void reb_whfast512c_corrector_step_avx2x4(struct simd_data* data, double inv, int gr, int nsys){ (void)data; (void)inv; (void)gr; (void)nsys; }
int reb_whfast512c_full_steps_pair_avx2x4(struct simd_data* da, struct simd_data* db, uint64_t N_steps, uint64_t* N_done, int skip_a, int skip_b, volatile sig_atomic_t* sigint, int* status_a, int* status_b, int gr, int nsys, int encounter, int escape){
    (void)da; (void)db; (void)N_steps; (void)N_done; (void)skip_a; (void)skip_b; (void)sigint; (void)status_a; (void)status_b; (void)gr; (void)nsys; (void)encounter; (void)escape;
    return -1; // never called
}

#endif
