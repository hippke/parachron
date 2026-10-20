/**
 * @file    integrator_whfast512_neon_ops.h
 * @brief   Vector operations of the WHFast512 C kernel for 64-bit ARM (NEON / Advanced SIMD).
 * @details Included by integrator_whfast512_c_generic.c (W_Q = 4: eight lanes) and
 *          integrator_whfast512_c_genericx4.c (W_Q = 2: four lanes) on AArch64. A vector is W_Q
 *          float64x2_t registers (lanes 0-1, 2-3, ...). All operations are IEEE-754 correctly
 *          rounded, and the fused multiply-adds are single-rounding (FMLA/FMLS), so the results
 *          are bitwise identical to the plain C (fma()) and x86 instantiations.
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
#include <arm_neon.h>

typedef struct { float64x2_t q[W_Q]; } V;
typedef struct { uint64x2_t q[W_Q]; } K;   // all bits set in true lanes

#define W_QLOOP for (int m=0; m<W_Q; m++)
W_INL V w_set1(double a){ V r; W_QLOOP r.q[m] = vdupq_n_f64(a); return r; }
W_INL V w_zero(void){ return w_set1(0.0); }
W_INL V w_load(const double* p){ V r; W_QLOOP r.q[m] = vld1q_f64(p+2*m); return r; }
W_INL void w_store(double* p, V a){ W_QLOOP vst1q_f64(p+2*m, a.q[m]); }
W_INL V w_add(V a, V b){ V r; W_QLOOP r.q[m] = vaddq_f64(a.q[m], b.q[m]); return r; }
W_INL V w_sub(V a, V b){ V r; W_QLOOP r.q[m] = vsubq_f64(a.q[m], b.q[m]); return r; }
W_INL V w_mul(V a, V b){ V r; W_QLOOP r.q[m] = vmulq_f64(a.q[m], b.q[m]); return r; }
W_INL V w_div(V a, V b){ V r; W_QLOOP r.q[m] = vdivq_f64(a.q[m], b.q[m]); return r; }
W_INL V w_sqrt(V a){ V r; W_QLOOP r.q[m] = vsqrtq_f64(a.q[m]); return r; }
W_INL V w_abs(V a){ V r; W_QLOOP r.q[m] = vabsq_f64(a.q[m]); return r; }
W_INL V w_floor(V a){ V r; W_QLOOP r.q[m] = vrndmq_f64(a.q[m]); return r; }
// a*b + c, c - a*b, a*b - c with one rounding each (= fma(a,b,c), fma(-a,b,c), fma(a,b,-c)).
W_INL V w_fmadd(V a, V b, V c){ V r; W_QLOOP r.q[m] = vfmaq_f64(c.q[m], a.q[m], b.q[m]); return r; }
W_INL V w_fnmadd(V a, V b, V c){ V r; W_QLOOP r.q[m] = vfmsq_f64(c.q[m], a.q[m], b.q[m]); return r; }
W_INL V w_fmsub(V a, V b, V c){ V r; W_QLOOP r.q[m] = vfmaq_f64(vnegq_f64(c.q[m]), a.q[m], b.q[m]); return r; }
// Ordered comparisons: false if either operand is NaN (as _CMP_LT_OQ etc.).
W_INL K w_klt(V a, V b){ K k; W_QLOOP k.q[m] = vcltq_f64(a.q[m], b.q[m]); return k; }
W_INL K w_kgt(V a, V b){ K k; W_QLOOP k.q[m] = vcgtq_f64(a.q[m], b.q[m]); return k; }
W_INL K w_kle(V a, V b){ K k; W_QLOOP k.q[m] = vcleq_f64(a.q[m], b.q[m]); return k; }
W_INL K w_kand(K a, K b){ K k; W_QLOOP k.q[m] = vandq_u64(a.q[m], b.q[m]); return k; }
W_INL K w_kandnot(K a, K b){ K k; W_QLOOP k.q[m] = vbicq_u64(b.q[m], a.q[m]); return k; }  // ~a & b
W_INL K w_knot(K a){ K k; W_QLOOP k.q[m] = veorq_u64(a.q[m], vdupq_n_u64(~0ULL)); return k; }
W_INL int w_kany(K a){
    uint64x2_t o = a.q[0];
    for (int m=1; m<W_Q; m++) o = vorrq_u64(o, a.q[m]);
    return (vgetq_lane_u64(o, 0) | vgetq_lane_u64(o, 1)) != 0;
}
W_INL K w_kall(void){ K k; W_QLOOP k.q[m] = vdupq_n_u64(~0ULL); return k; }
W_INL K w_kbits(unsigned int bits){
    K k;
    W_QLOOP {
        const uint64_t t[2] = { ((bits >> (2*m)) & 1u) ? ~0ULL : 0ULL, ((bits >> (2*m+1)) & 1u) ? ~0ULL : 0ULL };
        k.q[m] = vld1q_u64(t);
    }
    return k;
}
W_INL V w_sel(K k, V a, V b){ V r; W_QLOOP r.q[m] = vbslq_f64(k.q[m], a.q[m], b.q[m]); return r; }
W_INL V w_maskz(K k, V a){ V r; W_QLOOP r.q[m] = vreinterpretq_f64_u64(vandq_u64(k.q[m], vreinterpretq_u64_f64(a.q[m]))); return r; }
// Bit-trick seed for 1/sqrt(a): integer 0x5FE6EB50C7B537A9 - (bits(a) >> 1) (same as the other kernels).
W_INL V w_rsqrt_seed(V a){
    V r;
    W_QLOOP r.q[m] = vreinterpretq_f64_u64(vsubq_u64(vdupq_n_u64(0x5FE6EB50C7B537A9ULL), vshrq_n_u64(vreinterpretq_u64_f64(a.q[m]), 1)));
    return r;
}
// (x0,x1) -> (x1,x0), and the pairs (a1,b0) of two registers.
W_INL float64x2_t w_q_rev(float64x2_t a){ return vextq_f64(a, a, 1); }
W_INL float64x2_t w_q_ext(float64x2_t a, float64x2_t b){ return vextq_f64(a, b, 1); }
// Permutations within each group of four lanes (registers 2g, 2g+1), same lane maps as the
// other kernels: 3201 -> (a3,a2,a0,a1), 2310 -> (a2,a3,a1,a0), 3012 -> (a3,a0,a1,a2), swap -> (a1,a0,a3,a2).
W_INL V w_perm_3201(V a){ V r; for (int g=0; g<W_Q; g+=2){ r.q[g] = w_q_rev(a.q[g+1]); r.q[g+1] = a.q[g]; } return r; }
W_INL V w_perm_2310(V a){ V r; for (int g=0; g<W_Q; g+=2){ r.q[g] = a.q[g+1]; r.q[g+1] = w_q_rev(a.q[g]); } return r; }
W_INL V w_perm_3012(V a){ V r; for (int g=0; g<W_Q; g+=2){ r.q[g] = w_q_ext(a.q[g+1], a.q[g]); r.q[g+1] = w_q_ext(a.q[g], a.q[g+1]); } return r; }
W_INL V w_swap_pairs(V a){ V r; W_QLOOP r.q[m] = w_q_rev(a.q[m]); return r; }
#if W_Q == 4
// Across the two halves (lane maps as in the AVX-512 kernel):
// b3 -> (a4,a5,a6,a7,a1,a2,a3,a0), b4 -> (a5,a6,a7,a4,a2,a3,a0,a1), b34 -> (a7,a4,a5,a6,a0,a1,a2,a3).
W_INL V w_perm_b3(V a){ V r; r.q[0] = a.q[2]; r.q[1] = a.q[3]; r.q[2] = w_q_ext(a.q[0], a.q[1]); r.q[3] = w_q_ext(a.q[1], a.q[0]); return r; }
W_INL V w_perm_b4(V a){ V r; r.q[0] = w_q_ext(a.q[2], a.q[3]); r.q[1] = w_q_ext(a.q[3], a.q[2]); r.q[2] = a.q[1]; r.q[3] = a.q[0]; return r; }
W_INL V w_perm_b34(V a){ V r; r.q[0] = w_q_ext(a.q[3], a.q[2]); r.q[1] = w_q_ext(a.q[2], a.q[3]); r.q[2] = a.q[0]; r.q[3] = a.q[1]; return r; }
#else
// Only used by the N_systems = 1 code path with more than four planets (never run with four lanes).
W_INL V w_perm_b3(V a){ return a; }
W_INL V w_perm_b4(V a){ return a; }
W_INL V w_perm_b34(V a){ return a; }
#endif
