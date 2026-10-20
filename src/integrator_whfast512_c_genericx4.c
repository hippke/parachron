/**
 * @file    integrator_whfast512_c_genericx4.c
 * @brief   WHFast512 portable C kernel, plain C instantiation with four lanes (any CPU).
 * @details Same as integrator_whfast512_c_generic.c with four lanes: used automatically for a single
 *          system (N_systems = 1) with at most four planets (same results, about half the work).
 *          Defines the vector operations of integrator_whfast512_kernel.h as loops over four
 *          doubles. Fused multiply-adds use fma() from the C library, so the results are bitwise
 *          identical to the AVX2 and AVX-512 instantiations. On CPUs with FMA instructions (all
 *          64-bit ARM CPUs, x86 since 2013) compilers turn fma() into one instruction; elsewhere
 *          it is computed in software (correct, but slow).
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
#include <math.h>
#include "integrator_whfast512_data.h"

#define REB_W512_TARGET

#if defined(__aarch64__) && defined(__ARM_NEON) && !defined(REB_W512_NO_NEON)
#define REB_W512_NEON 1
#endif

#if defined(_MSC_VER) && !defined(__clang__)
#define W_INL static __forceinline
#define REB_W512_NOINLINE_FALLBACK static __declspec(noinline)
#else
#define W_INL static inline __attribute__((always_inline))
#define REB_W512_NOINLINE_FALLBACK static __attribute__((noinline))
#endif
// One out-of-line copy of the Kepler and interaction steps (instead of one per option
// combination): keeps the code small and compilation fast.
#if defined(_MSC_VER) && !defined(__clang__)
#define REB_W512_STEP_FUNC static __declspec(noinline)
#else
#define REB_W512_STEP_FUNC static __attribute__((noinline))
#endif

#ifdef REB_W512_NEON
#define W_Q 2
#include "integrator_whfast512_neon_ops.h"
#else // plain C
typedef struct { double v[4]; } V;
typedef struct { unsigned char b; } K;   // bit i set: lane i true
#define W_LOOP for (int i=0; i<4; i++)
W_INL V w_set1(double a){ V r; W_LOOP r.v[i] = a; return r; }
W_INL V w_zero(void){ return w_set1(0.0); }
W_INL V w_load(const double* p){ V r; W_LOOP r.v[i] = p[i]; return r; }
W_INL void w_store(double* p, V a){ W_LOOP p[i] = a.v[i]; }
W_INL V w_add(V a, V b){ V r; W_LOOP r.v[i] = a.v[i] + b.v[i]; return r; }
W_INL V w_sub(V a, V b){ V r; W_LOOP r.v[i] = a.v[i] - b.v[i]; return r; }
W_INL V w_mul(V a, V b){ V r; W_LOOP r.v[i] = a.v[i] * b.v[i]; return r; }
W_INL V w_div(V a, V b){ V r; W_LOOP r.v[i] = a.v[i] / b.v[i]; return r; }
W_INL V w_sqrt(V a){ V r; W_LOOP r.v[i] = sqrt(a.v[i]); return r; }
W_INL V w_abs(V a){ V r; W_LOOP r.v[i] = fabs(a.v[i]); return r; }
W_INL V w_floor(V a){ V r; W_LOOP r.v[i] = floor(a.v[i]); return r; }
W_INL V w_fmadd(V a, V b, V c){ V r; W_LOOP r.v[i] = fma(a.v[i], b.v[i], c.v[i]); return r; }
W_INL V w_fnmadd(V a, V b, V c){ V r; W_LOOP r.v[i] = fma(-a.v[i], b.v[i], c.v[i]); return r; }
W_INL V w_fmsub(V a, V b, V c){ V r; W_LOOP r.v[i] = fma(a.v[i], b.v[i], -c.v[i]); return r; }
// Ordered, quiet comparisons: false if either operand is NaN (as _CMP_LT_OQ etc.).
W_INL K w_klt(V a, V b){ K k = {0}; W_LOOP if (a.v[i] < b.v[i]) k.b |= (unsigned char)(1u<<i); return k; }
W_INL K w_kgt(V a, V b){ K k = {0}; W_LOOP if (a.v[i] > b.v[i]) k.b |= (unsigned char)(1u<<i); return k; }
W_INL K w_kle(V a, V b){ K k = {0}; W_LOOP if (a.v[i] <= b.v[i]) k.b |= (unsigned char)(1u<<i); return k; }
W_INL K w_kand(K a, K b){ K k; k.b = (unsigned char)(a.b & b.b); return k; }
W_INL K w_kandnot(K a, K b){ K k; k.b = (unsigned char)(~a.b & b.b); return k; }
W_INL K w_knot(K a){ K k; k.b = (unsigned char)(~a.b & 0x0F); return k; }
W_INL int w_kany(K a){ return a.b != 0; }
W_INL K w_kall(void){ K k; k.b = 0x0F; return k; }
W_INL K w_kbits(unsigned int bits){ K k; k.b = (unsigned char)(bits & 0x0F); return k; }
W_INL V w_sel(K k, V a, V b){ V r; W_LOOP r.v[i] = ((k.b >> i) & 1) ? a.v[i] : b.v[i]; return r; }
W_INL V w_maskz(K k, V a){ V r; W_LOOP r.v[i] = ((k.b >> i) & 1) ? a.v[i] : 0.0; return r; }
// Bit-trick seed for 1/sqrt(a): integer 0x5FE6EB50C7B537A9 - (bits(a) >> 1) (deterministic, same on every machine).
W_INL V w_rsqrt_seed(V a){ V r; W_LOOP { uint64_t b; memcpy(&b, &a.v[i], 8); b = 0x5FE6EB50C7B537A9ULL - (b >> 1); memcpy(&r.v[i], &b, 8); } return r; }
W_INL V w_permute(V a, const int* idx){ V r; W_LOOP r.v[i] = a.v[idx[i]]; return r; }
static const int w_idx_3201[4] = {3,2,0,1};
static const int w_idx_2310[4] = {2,3,1,0};
static const int w_idx_3012[4] = {3,0,1,2};
static const int w_idx_swap[4] = {1,0,3,2};
W_INL V w_perm_3201(V a){ return w_permute(a, w_idx_3201); }
W_INL V w_perm_2310(V a){ return w_permute(a, w_idx_2310); }
W_INL V w_perm_3012(V a){ return w_permute(a, w_idx_3012); }
W_INL V w_swap_pairs(V a){ return w_permute(a, w_idx_swap); }
// Only used by the N_systems = 1 code path with more than four planets (never run here).
W_INL V w_perm_b3(V a){ return a; }
W_INL V w_perm_b4(V a){ return a; }
W_INL V w_perm_b34(V a){ return a; }

#endif // REB_W512_NEON

#define W_LANES 4
#define W_ONLY_NSYS2 1
#define W_PAIR_INLINE 1
#define REB_W512_SUFFIX genericx4
#include "integrator_whfast512_kernel.h"
