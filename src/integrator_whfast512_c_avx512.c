/**
 * @file    integrator_whfast512_c_avx512.c
 * @brief   WHFast512 portable C kernel, AVX-512 instantiation.
 * @details Defines the vector operations of integrator_whfast512_kernel.h with AVX-512F/DQ
 *          intrinsics. The functions are compiled for AVX-512 through a target pragma, so the file
 *          needs no special compiler flags; the driver calls them only on CPUs that support
 *          AVX-512F and AVX-512DQ.
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

#if (defined(__x86_64__) || defined(_M_X64)) && (defined(__clang__) || (defined(__GNUC__) && __GNUC__ >= 5) || defined(_MSC_VER))
#define REB_W512_HAVE_AVX512 1

#if defined(__clang__)
#pragma clang attribute push (__attribute__((target("avx512f,avx512dq,avx2,fma"))), apply_to = function)
#define REB_W512_TARGET
#elif defined(__GNUC__)
#pragma GCC push_options
#pragma GCC target("avx512f,avx512dq,avx2,fma")
#define REB_W512_TARGET
#else
#define REB_W512_TARGET
#endif
#include <immintrin.h>

typedef __m512d V;
typedef __mmask8 K;

#if defined(_MSC_VER) && !defined(__clang__)
#define W_INL static __forceinline
#define REB_W512_NOINLINE_FALLBACK static __declspec(noinline)
#else
#define W_INL static inline __attribute__((always_inline))
#define REB_W512_NOINLINE_FALLBACK static __attribute__((noinline))
#endif

W_INL V w_set1(double a){ return _mm512_set1_pd(a); }
W_INL V w_zero(void){ return _mm512_setzero_pd(); }
W_INL V w_load(const double* p){ return _mm512_load_pd(p); }
W_INL void w_store(double* p, V a){ _mm512_store_pd(p, a); }
W_INL V w_add(V a, V b){ return _mm512_add_pd(a, b); }
W_INL V w_sub(V a, V b){ return _mm512_sub_pd(a, b); }
W_INL V w_mul(V a, V b){ return _mm512_mul_pd(a, b); }
W_INL V w_div(V a, V b){ return _mm512_div_pd(a, b); }
W_INL V w_sqrt(V a){ return _mm512_sqrt_pd(a); }
W_INL V w_abs(V a){ return _mm512_castsi512_pd(_mm512_and_si512(_mm512_castpd_si512(a), _mm512_set1_epi64(0x7FFFFFFFFFFFFFFFLL))); }
W_INL V w_floor(V a){ return _mm512_roundscale_pd(a, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC); }
W_INL V w_fmadd(V a, V b, V c){ return _mm512_fmadd_pd(a, b, c); }
W_INL V w_fnmadd(V a, V b, V c){ return _mm512_fnmadd_pd(a, b, c); }
W_INL V w_fmsub(V a, V b, V c){ return _mm512_fmsub_pd(a, b, c); }
W_INL K w_klt(V a, V b){ return _mm512_cmp_pd_mask(a, b, _CMP_LT_OQ); }
W_INL K w_kgt(V a, V b){ return _mm512_cmp_pd_mask(a, b, _CMP_GT_OQ); }
W_INL K w_kle(V a, V b){ return _mm512_cmp_pd_mask(a, b, _CMP_LE_OQ); }
W_INL K w_kand(K a, K b){ return (K)(a & b); }
W_INL K w_kandnot(K a, K b){ return (K)(~a & b); }
W_INL K w_knot(K a){ return (K)(~a); }
W_INL int w_kany(K a){ return a != 0; }
W_INL K w_kall(void){ return (K)0xFF; }
W_INL K w_kbits(unsigned int bits){ return (K)bits; }
W_INL V w_sel(K k, V a, V b){ return _mm512_mask_blend_pd(k, b, a); }
W_INL V w_maskz(K k, V a){ return _mm512_maskz_mov_pd(k, a); }
// Bit-trick seed for 1/sqrt(a): integer 0x5FE6EB50C7B537A9 - (bits(a) >> 1) (deterministic, same on every machine).
#ifdef REB_W512_HW_RSQRT
// Experimental (compile-time option, off by default): vrsqrt14pd (14 bits) as the seed of 1/sqrt, one Newton step
// fewer (at most 1 ulp, as the default). Results differ from the other kernels at round-off level.
W_INL V w_rsqrt_seed(V a){ return _mm512_rsqrt14_pd(a); }
#define W512_RSQRT_NR 2
#else
W_INL V w_rsqrt_seed(V a){ return _mm512_castsi512_pd(_mm512_sub_epi64(_mm512_set1_epi64(0x5FE6EB50C7B537A9LL), _mm512_srli_epi64(_mm512_castpd_si512(a), 1))); }
#endif
W_INL V w_perm_3201(V a){ return _mm512_permutex_pd(a, 0x4B); }
W_INL V w_perm_2310(V a){ return _mm512_permutex_pd(a, 0x1E); }
W_INL V w_perm_3012(V a){ return _mm512_permutex_pd(a, 0x93); }
W_INL V w_swap_pairs(V a){ return _mm512_permute_pd(a, 0x55); }
// Lanes shifted up (lane i gets lane i-s) or down (lane i gets lane i+s), zeros shifted in; s = 1, 2 or 4.
W_INL V w_shup(V a, int s){ const __m512i z = _mm512_setzero_si512(), b = _mm512_castpd_si512(a);
    return _mm512_castsi512_pd(s == 1 ? _mm512_alignr_epi64(b, z, 7) : s == 2 ? _mm512_alignr_epi64(b, z, 6) : _mm512_alignr_epi64(b, z, 4)); }
W_INL V w_shdn(V a, int s){ const __m512i z = _mm512_setzero_si512(), b = _mm512_castpd_si512(a);
    return _mm512_castsi512_pd(s == 1 ? _mm512_alignr_epi64(z, b, 1) : s == 2 ? _mm512_alignr_epi64(z, b, 2) : _mm512_alignr_epi64(z, b, 4)); }
W_INL V w_perm_b3(V a){ return _mm512_permutexvar_pd(_mm512_set_epi64(0,3,2,1,7,6,5,4), a); }
W_INL V w_perm_b4(V a){ return _mm512_permutexvar_pd(_mm512_set_epi64(1,0,3,2,4,7,6,5), a); }
W_INL V w_perm_b34(V a){ return _mm512_permutexvar_pd(_mm512_set_epi64(3,2,1,0,6,5,4,7), a); }

#define REB_W512_SUFFIX avx512
#include "integrator_whfast512_kernel.h"

#if defined(__clang__)
#pragma clang attribute pop
#elif defined(__GNUC__)
#pragma GCC pop_options
#endif

int reb_whfast512c_avx512_compiled(void){ return 1; }

#else // no AVX-512 support in this compiler/architecture
int reb_whfast512c_full_steps_avx512(struct simd_data* data, uint64_t* N_steps, int skip_first_kepler_step, volatile sig_atomic_t* sigint, int gr, int nsys, int encounter, int escape){
    (void)data; (void)N_steps; (void)skip_first_kepler_step; (void)sigint; (void)gr; (void)nsys; (void)encounter; (void)escape;
    return 1; // REB_STATUS_GENERIC_ERROR, never called (see reb_whfast512c_avx512_compiled)
}
void reb_whfast512c_kepler_step_avx512(struct simd_data* data){ (void)data; }
int reb_whfast512c_full_steps_comp_avx512(struct simd_data* data, uint64_t* N_steps, int skip_first_kepler_step, volatile sig_atomic_t* sigint, int gr, int nsys, int encounter, int escape){
    (void)data; (void)N_steps; (void)skip_first_kepler_step; (void)sigint; (void)gr; (void)nsys; (void)encounter; (void)escape;
    return 1; // never called
}
void reb_whfast512c_corrector_step_avx512(struct simd_data* data, double inv, int gr, int nsys){ (void)data; (void)inv; (void)gr; (void)nsys; }
int reb_whfast512c_avx512_compiled(void){ return 0; }
int reb_whfast512c_full_steps_pair_avx512(struct simd_data* da, struct simd_data* db, uint64_t N_steps, uint64_t* N_done, int skip_a, int skip_b, volatile sig_atomic_t* sigint, int* status_a, int* status_b, int gr, int nsys, int encounter, int escape){
    (void)da; (void)db; (void)N_steps; (void)N_done; (void)skip_a; (void)skip_b; (void)sigint; (void)status_a; (void)status_b; (void)gr; (void)nsys; (void)encounter; (void)escape;
    return -1; // never called
}

#endif
