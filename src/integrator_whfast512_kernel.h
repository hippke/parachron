/**
 * @file    integrator_whfast512_kernel.h
 * @brief   Portable C kernel of WHFast512 (Kepler step, interaction step, corrector, main loop).
 * @details This file is included by integrator_whfast512_c_avx512.c, integrator_whfast512_c_avx2.c
 *          and integrator_whfast512_c_generic.c. Each of them first defines an 8-lane vector type
 *          `V`, a lane mask type `K`, the operations below (w_add, w_fmadd, ...) for its instruction
 *          set, and the macro REB_W512_SUFFIX. The algorithm is a line-by-line port of the AVX-512
 *          assembly kernel (integrator_whfast512.s, WHFast512 version 2, Dagli & Rein): same
 *          operations in the same order, and the same operands in each fused multiply-add, so that
 *          results agree with the assembly to round-off. Differences to the assembly:
 *          - In the interaction step, 1/sqrt(r^2) (planet-planet forces, Jacobi term, central body and
 *            GR terms) is a square root followed by a division (W512_RSQRT). Both are correctly
 *            rounded, so every instruction set gives the same bits. The assembly uses vrsqrt14pd plus
 *            Newton iterations for the planet-planet terms (vrsqrt14pd does not exist for AVX2, and its
 *            result is not specified bit for bit), and square roots and divisions for the other terms.
 *            Compile-time alternative REB_W512_NEWTON_RSQRT: an integer bit-trick seed, three Newton
 *            iterations and one fused-multiply-add correction (at most 1 ulp error), faster where
 *            512-bit divisions are slow (Intel Sapphire Rapids), slower on AMD Zen 4 and Zen 5.
 *            Distance checks for encounters and escapes compare squared distances.
 *          - 0.5*f*f'' in the Halley step is a multiplication by 0.5 (the assembly subtracts one from
 *            the exponent, which gives the same result for normal numbers, but not for zero).
 *          - Starting value of the Kepler solver: sixth-order series reversion of Kepler's equation in
 *            universal variables (the assembly uses the first-order guess dt/r), so
 *            that a single Halley step is followed by a single Newton step for the Solar System at
 *            dt <= 6 days. The first Newton step uses an approximation of 1/f' from the Halley step
 *            (no division on the critical path). The Stiefel functions are evaluated again at the
 *            final X as in the assembly, with plain Horner plus a correction for the rounding of the
 *            leading constants instead of compensated Horner (see w512_stiefel_Gs13_ for the
 *            round-off bias that the plain shortcuts had).
 *          - MXCSR (flush-to-zero / denormals-are-zero) is not changed.
 *          All three instantiations use explicit fused multiply-adds, so they give bitwise
 *          identical results on every machine (the generic one uses fma() from libm).
 *
 *          Operations each backend provides (a, b, c: V; k: K):
 *            w_set1(double), w_zero(), w_load(const double*), w_store(double*, V)
 *            w_add, w_sub, w_mul, w_div, w_sqrt, w_abs, w_floor
 *            w_fmadd(a,b,c) = a*b+c, w_fnmadd(a,b,c) = c-a*b, w_fmsub(a,b,c) = a*b-c  (one rounding)
 *            w_klt(a,b) = a<b, w_kgt(a,b) = a>b, w_kle(a,b) = a<=b  (ordered, quiet: false for NaN)
 *            w_kand(k,l), w_kandnot(k,l) = (~k)&l, w_knot(k), w_kany(k), w_kall(), w_kbits(unsigned)
 *            w_sel(k,a,b) = k ? a : b, w_maskz(k,a) = k ? a : 0
 *            w_rsqrt_seed(a): bits 0x5FE6EB50C7B537A9 - (bits(a) >> 1) as an integer operation (used only with REB_W512_NEWTON_RSQRT)
 *            w_perm_3201(a), w_perm_2310(a), w_perm_3012(a)  (within each group of four lanes)
 *            w_swap_pairs(a)                                  (1032 5476)
 *            w_perm_b3(a) (45671230), w_perm_b4(a) (56742301), w_perm_b34(a) (74560123)
 *
 * @section LICENSE
 * Copyright (c) 2026 Hanno Rein, Rishit Dagli
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

#ifndef REB_W512_SUFFIX
#error "Define REB_W512_SUFFIX and the vector operations before including integrator_whfast512_kernel.h"
#endif

#define REB_W512_CAT2(a,b) a##b
#define REB_W512_CAT(a,b) REB_W512_CAT2(a,b)
#define REB_W512_NAME(name) REB_W512_CAT(name, REB_W512_SUFFIX)

#if defined(_MSC_VER) && !defined(__clang__)
#define REB_W512_INLINE static __forceinline
#else
#define REB_W512_INLINE static inline __attribute__((always_inline))
#endif
// The Kepler and interaction steps. The vector instantiations inline them into each of the 24
// main loops (options become compile-time constants). The generic instantiation can define
// REB_W512_STEP_FUNC as a plain static function to keep code size and compile time small.
#ifndef REB_W512_STEP_FUNC
#define REB_W512_STEP_FUNC REB_W512_INLINE
#endif
// Number of lanes of V. 8 for the normal kernels. A 4-lane instantiation (W_LANES 4,
// W_ONLY_NSYS2) integrates one system with up to four planets in lanes 0-3 with the
// N_systems=2 code path; the pair terms it skips involve only empty lanes (zero mass), so the
// results are the same as with 8 lanes.
#ifndef W_LANES
#define W_LANES 8
#endif

// 1/k! rounded to double, and the rounding error (true - rounded). Same bits as .IF0 / .IF0_err
// in integrator_whfast512.s.
static const uint64_t REB_W512_NAME(w512_if0_bits_)[35] = {
    0x3ff0000000000000ULL, 0x3ff0000000000000ULL, 0x3fe0000000000000ULL, 0x3fc5555555555555ULL,
    0x3fa5555555555555ULL, 0x3f81111111111111ULL, 0x3f56c16c16c16c17ULL, 0x3f2a01a01a01a01aULL,
    0x3efa01a01a01a01aULL, 0x3ec71de3a556c734ULL, 0x3e927e4fb7789f5cULL, 0x3e5ae64567f544e4ULL,
    0x3e21eed8eff8d898ULL, 0x3de6124613a86d09ULL, 0x3da93974a8c07c9dULL, 0x3d6ae7f3e733b81fULL,
    0x3d2ae7f3e733b81fULL, 0x3ce952c77030ad4aULL, 0x3ca6827863b97d97ULL, 0x3c62f49b46814157ULL,
    0x3c1e542ba4020225ULL, 0x3bd71b8ef6dcf572ULL, 0x3b90ce396db7f853ULL, 0x3b4761b41316381aULL,
    0x3aff2cf01972f578ULL, 0x3ab3f3ccdd165fa9ULL, 0x3a688e85fc6a4e5aULL, 0x3a1d1ab1c2dccea3ULL,
    0x39d0a18a2635085dULL, 0x398259f98b4358adULL, 0x3933932c5047d60eULL, 0x38e434d2e783f5bcULL,
    0x389434d2e783f5bcULL, 0x3843981254dd0d52ULL, 0x37f2710231c0fd7aULL,
};
static const uint64_t REB_W512_NAME(w512_if0_err_bits_)[35] = {
    0x0000000000000000ULL, 0x0000000000000000ULL, 0x0000000000000000ULL, 0xbc65555555555555ULL,
    0xbc45555555555555ULL, 0xbc01111111111111ULL, 0x3bef49f49f49f49fULL, 0xbb6a01a01a01a01aULL,
    0xbb3a01a01a01a01aULL, 0x3b6c154f8ddc6c00ULL, 0xbb3cbbc05b4fa99aULL, 0x3afc062e06d1f209ULL,
    0x3ac2aec959e14c06ULL, 0xba8f28e0cc748ebeULL, 0xba305d6f8a2efd1fULL, 0xb9e1d8656b0ee8cbULL,
    0xb9a1d8656b0ee8cbULL, 0xb98ac981465ddc6cULL, 0xb94eec01221a8b0bULL, 0xb8f2650f61dbdcb4ULL,
    0xb87ea72b4afe3c2fULL, 0x387d043ae40c4647ULL, 0x383aebcdbd20331cULL, 0x37d3423c7d91404fULL,
    0x3789ada5fcc1ab14ULL, 0x37458ddadf344487ULL, 0x37071c37ebd16540ULL, 0xb6a054d0c78aea14ULL,
    0xb66b9e2e28e1aa54ULL, 0xb62eaf8c39dd9bc5ULL, 0xb5d832b7b530a627ULL, 0xb580b87b91be9affULL,
    0xb530b87b91be9affULL, 0x34e2b1f4c8015a2fULL, 0xb473f8a2b4af9d6bULL,
};
// 17th order corrector: 63 coefficients, alternating Kepler and interaction (same as .CORRECTOR17_AB).
static const uint64_t REB_W512_NAME(w512_corrector17_bits_)[63] = {
    0xC00AC5EB3F7AB2F8ULL, 0xBED22E64AF0557FFULL, 0x401AC5EB3F7AB2F8ULL, 0x3ED22E64AF0557FFULL,
    0xC019198C8B8307C8ULL, 0x3F14098E956E85C7ULL, 0x40176D2DD78B5C99ULL, 0xBF14098E956E85C7ULL,
    0xC015C0CF2393B16AULL, 0xBF44D7273C9751C8ULL, 0x401414706F9C063AULL, 0x3F44D7273C9751C8ULL,
    0xC0126811BBA45B0AULL, 0x3F6B2467AFD31964ULL, 0x4010BBB307ACAFDBULL, 0xBF6B2467AFD31964ULL,
    0xC00E1EA8A76A0957ULL, 0xBF88B9144F7F2B3FULL, 0x400AC5EB3F7AB2F8ULL, 0x3F88B9144F7F2B3FULL,
    0xC0076D2DD78B5C99ULL, 0x3FA099A47793A306ULL, 0x400414706F9C063AULL, 0xBFA099A47793A306ULL,
    0xC000BBB307ACAFDBULL, 0xBFB0B07AC0FE3DF1ULL, 0x3FFAC5EB3F7AB2F8ULL, 0x3FB0B07AC0FE3DF1ULL,
    0xBFF414706F9C063AULL, 0x3FB7D2865A643682ULL, 0x3FEAC5EB3F7AB2F8ULL, 0xBFC7D2865A643682ULL,
    0xBFEAC5EB3F7AB2F8ULL, 0x3FB7D2865A643682ULL, 0x3FF414706F9C063AULL, 0x3FB0B07AC0FE3DF1ULL,
    0xBFFAC5EB3F7AB2F8ULL, 0xBFB0B07AC0FE3DF1ULL, 0x4000BBB307ACAFDBULL, 0xBFA099A47793A306ULL,
    0xC00414706F9C063AULL, 0x3FA099A47793A306ULL, 0x40076D2DD78B5C99ULL, 0x3F88B9144F7F2B3FULL,
    0xC00AC5EB3F7AB2F8ULL, 0xBF88B9144F7F2B3FULL, 0x400E1EA8A76A0957ULL, 0xBF6B2467AFD31964ULL,
    0xC010BBB307ACAFDBULL, 0x3F6B2467AFD31964ULL, 0x40126811BBA45B0AULL, 0x3F44D7273C9751C8ULL,
    0xC01414706F9C063AULL, 0xBF44D7273C9751C8ULL, 0x4015C0CF2393B16AULL, 0xBF14098E956E85C7ULL,
    0xC0176D2DD78B5C99ULL, 0x3F14098E956E85C7ULL, 0x4019198C8B8307C8ULL, 0x3ED22E64AF0557FFULL,
    0xC01AC5EB3F7AB2F8ULL, 0xBED22E64AF0557FFULL, 0x400AC5EB3F7AB2F8ULL,
};

REB_W512_INLINE double REB_W512_NAME(w512_bits_)(uint64_t b){
    double d;
    memcpy(&d, &b, sizeof(d));
    return d;
}
#define W512_IF(k)     REB_W512_NAME(w512_bits_)(REB_W512_NAME(w512_if0_bits_)[k])
#define W512_IF_ERR(k) REB_W512_NAME(w512_bits_)(REB_W512_NAME(w512_if0_err_bits_)[k])

#define W512_TWOPI          6.283185307179586231995926937088 // 0x401921fb54442d18
#define W512_EPS            1e-11   // Newton convergence: |delta X| < EPS
#define W512_BISECTION_EPS  1e-15
#define W512_NEWTON_MAX     5       // Newton iterations before falling back to bisection
#define W512_BISECTION_MIN  52      // bisection iterations (elliptic lanes); hyperbolic lanes continue to convergence

// Constants that stay fixed during a call (registers in the assembly).
struct REB_W512_NAME(w512_consts_) {
    V M;        // .LM: mass used in the Kepler solver
    V DT;       // .LDT
    V MM0_DT;   // .LMM0_DT: -dt*m0 (Jacobi term)
    V ONE, HALF, EPS;
    K k1;       // lanes that hold particles
    const struct simd_data* data;
};
#define W512_CONSTS struct REB_W512_NAME(w512_consts_)

// State (registers .LX ... .LVZ in the assembly). CX ... CVZ: the low-order parts of the coordinates in the
// compensated mode (see w512_add_comp_), zero otherwise.
struct REB_W512_NAME(w512_state_) {
    V X, Y, Z, VX, VY, VZ;
    V CX, CY, CZ, CVX, CVY, CVZ;
};
#define W512_STATE struct REB_W512_NAME(w512_state_)

///////////////////////////////////////////////////////////////////////////////
// Stiefel functions (series in z = beta*X^2), as the mm_stiefel_* macros.

// Low accuracy (Gs0..Gs3) for the Halley steps. numTerms odd.
REB_W512_INLINE void REB_W512_NAME(w512_stiefel_Gs03_)(const int numTerms, const V X, const V beta, V* G0, V* G1, V* G2, V* G3){
    const V X2 = w_mul(X, X);
    V g3 = w_set1(W512_IF(numTerms));
    V g2 = w_set1(W512_IF(numTerms-1));
    const V z = w_mul(beta, X2);
    int off = numTerms - 2;
    const int it = (off - 1)/2 - 1;
    for (int k=0; k<it; k++){
        g3 = w_fnmadd(z, g3, w_set1(W512_IF(off))); off--;
        g2 = w_fnmadd(z, g2, w_set1(W512_IF(off))); off--;
    }
    g3 = w_fnmadd(z, g3, w_set1(W512_IF(3)));
    V g0 = g2;
    g2 = w_fnmadd(z, g2, w_set1(W512_IF(2)));
    g0 = w_fnmadd(z, g0, w_set1(W512_IF(1)));
    *G2 = w_mul(X2, g2);
    g3 = w_mul(X, g3);
    *G3 = w_mul(X2, g3);
    *G1 = w_fnmadd(z, g3, X);
    *G0 = g0;
}

// High accuracy (Gs1..Gs3), plain Horner with a correction for the rounding of the leading
// constants. The constants 1/k! are rounded to double, and since the leading terms
// dominate the sums, plain Horner has a systematic relative error: 1/3! is stored 0.5 ulp low, so G3
// comes out 5.6e-17 too small in every step, which biases the Kepler step (seen as a systematic
// drift of a in a grid of single Kepler orbits). The rounding errors of the leading constants (W512_IF_ERR, rounded -
// true) are added before the final rounding (fused multiply-add), in parallel with the Horner chain.
// numTerms odd.
REB_W512_INLINE void REB_W512_NAME(w512_stiefel_Gs13_)(const int numTerms, const V X, const V beta, V* G1, V* G2, V* G3){
    const V X2 = w_mul(X, X);
    V g3 = w_set1(W512_IF(numTerms));
    V g2 = w_set1(W512_IF(numTerms-1));
    const V z = w_mul(beta, X2);
    int off = numTerms - 2;
    const int it = (off - 1)/2;
    for (int k=0; k<it; k++){
        g3 = w_fnmadd(z, g3, w_set1(W512_IF(off))); off--;
        g2 = w_fnmadd(z, g2, w_set1(W512_IF(off))); off--;
    }
    // Missing parts of the sums: -(E3 - z E5 + z^2 E7) for g3 and z E4 - z^2 E6 for g2 (E = W512_IF_ERR).
    const V e3 = w_fmadd(w_fmadd(w_set1(-W512_IF_ERR(7)), z, w_set1(W512_IF_ERR(5))), z, w_set1(-W512_IF_ERR(3)));
    const V e2 = w_mul(z, w_fmsub(w_set1(-W512_IF_ERR(6)), z, w_set1(-W512_IF_ERR(4))));
    *G2 = w_fmadd(X2, g2, w_mul(X2, e2));
    const V hc = w_mul(X, e3);
    const V h = w_mul(X, g3);
    *G3 = w_fmadd(X2, h, w_mul(X2, hc));
    *G1 = w_fnmadd(z, w_fmadd(X, g3, hc), X);
}

// One step of compensated Horner evaluation (comp_horner_step): C <- IF[k] - z*C, CLO tracks the error.
REB_W512_INLINE void REB_W512_NAME(w512_comp_horner_step_)(const V z, V* C, V* CLO, const int k, const int has_err){
    const V P = w_mul(z, *C);
    const V pe = w_fmsub(*C, z, P);          // exact error of the product
    const V IFk = w_set1(W512_IF(k));
    const V Cn = w_sub(IFk, P);
    V se = w_sub(w_sub(IFk, Cn), P);         // error of the subtraction
    if (has_err){
        se = w_sub(se, w_set1(W512_IF_ERR(k)));
    }
    const V e = w_sub(se, pe);
    *CLO = w_fnmadd(z, *CLO, e);
    *C = Cn;
}

// High accuracy (Gs1..Gs3) with the last four Horner steps compensated (mm_stiefel_Gs13_comp).
REB_W512_INLINE void REB_W512_NAME(w512_stiefel_Gs13_comp_)(const int numTerms, const V X, const V beta, V* G1, V* G2, V* G3){
    const V X2 = w_mul(X, X);
    V g3 = w_set1(W512_IF(numTerms));
    V g2 = w_set1(W512_IF(numTerms-1));
    const V z = w_mul(beta, X2);
    int off = numTerms - 2;
    const int it = (off - 1)/2 - 2;
    for (int k=0; k<it; k++){
        g3 = w_fnmadd(z, g3, w_set1(W512_IF(off))); off--;
        g2 = w_fnmadd(z, g2, w_set1(W512_IF(off))); off--;
    }
    V lo3 = w_zero();
    V lo2 = w_zero();
    REB_W512_NAME(w512_comp_horner_step_)(z, &g3, &lo3, 5, 1);
    REB_W512_NAME(w512_comp_horner_step_)(z, &g2, &lo2, 4, 1);
    REB_W512_NAME(w512_comp_horner_step_)(z, &g3, &lo3, 3, 1);
    REB_W512_NAME(w512_comp_horner_step_)(z, &g2, &lo2, 2, 0);
    g3 = w_add(g3, lo3);
    g2 = w_add(g2, lo2);
    *G2 = w_mul(X2, g2);
    g3 = w_mul(X, g3);
    *G3 = w_mul(X2, g3);
    *G1 = w_fnmadd(z, g3, X);
}

///////////////////////////////////////////////////////////////////////////////
// Kepler step

// Fallback for lanes in k4 that did not converge with Newton's method: bracket the root and bisect.
REB_W512_NOINLINE_FALLBACK void REB_W512_NAME(w512_kepler_fallback_)(const W512_CONSTS* c, const W512_STATE* s, const V R, const V eta, const V zeta, const V beta, V* XXp, K k4){
    const V zero = w_zero();
    V XX = *XXp;
    V XMIN = zero;
    const K k5 = w_kand(w_kle(beta, zero), k4);   // hyperbolic/parabolic lanes
    const K k6 = w_kandnot(k5, k4);               // elliptic lanes
    // Elliptic bounds: the root lies within the orbital period containing dt.
    const V sb = w_maskz(k6, w_sqrt(beta));
    V t = w_maskz(k6, w_mul(beta, sb));
    const V tp = w_mul(w_set1(W512_TWOPI), c->M);
    t = w_maskz(k6, w_div(t, tp));                // inverse period
    t = w_maskz(k6, w_mul(c->DT, t));
    t = w_maskz(k6, w_floor(t));
    const V xpp = w_maskz(k6, w_div(w_set1(W512_TWOPI), sb));  // X per period
    XMIN = w_sel(k6, w_mul(xpp, t), XMIN);
    V XMAX = w_sel(k6, w_add(XMIN, xpp), xpp);
    // Hyperbolic bounds.
    if (w_kany(k5)){
        const V hx = w_maskz(k5, w_fnmadd(s->VY, s->Z, w_mul(s->VZ, s->Y)));
        const V hy = w_maskz(k5, w_fnmadd(s->VZ, s->X, w_mul(s->VX, s->Z)));
        const V hz = w_maskz(k5, w_fnmadd(s->VX, s->Y, w_mul(s->VY, s->X)));
        V h2 = w_sel(k5, w_mul(hx, hx), hx);
        h2 = w_sel(k5, w_fmadd(hy, hy, h2), h2);
        h2 = w_sel(k5, w_fmadd(hz, hz, h2), h2);
        const V H2 = w_maskz(k5, h2);
        const V MM = w_maskz(k5, w_mul(c->M, c->M));
        V e = w_maskz(k5, w_mul(beta, H2));
        e = w_maskz(k5, w_div(e, MM));
        e = w_maskz(k5, w_sub(c->ONE, e));
        e = w_maskz(k5, w_sqrt(e));
        const V ope = w_maskz(k5, w_add(c->ONE, e));          // 1 + e
        V q = w_maskz(k5, w_div(H2, c->M));
        q = w_maskz(k5, w_div(q, ope));                         // pericentre distance
        V vq = w_maskz(k5, w_sqrt(H2));
        vq = w_maskz(k5, w_div(vq, q));                         // pericentre velocity
        V a = w_abs(w_maskz(k5, w_mul(c->DT, vq)));
        a = w_maskz(k5, w_add(R, a));
        XMIN = w_sel(k5, w_div(c->DT, a), XMIN);
        XMAX = w_sel(k5, w_div(c->DT, q), XMAX);
        const K k7 = w_kand(w_klt(c->DT, zero), k5);           // dt < 0: swap bounds
        const V tmp = XMIN;
        XMIN = w_sel(k7, XMAX, XMIN);
        XMAX = w_sel(k7, tmp, XMAX);
    }
    XX = w_sel(k4, w_mul(w_add(XMAX, XMIN), c->HALF), XX);
    int it = 0;
    while (1){
        V G1, G2, G3;
        REB_W512_NAME(w512_stiefel_Gs13_)(19, XX, beta, &G1, &G2, &G3);
        V f = w_mul(R, XX);
        f = w_fmadd(G2, eta, f);
        f = w_fmadd(G3, zeta, f);
        const K k2 = w_kgt(f, c->DT);
        const K k3 = w_knot(k2);
        XMAX = w_sel(k2, XX, XMAX);
        XMIN = w_sel(k3, XX, XMIN);
        XX = w_sel(k4, w_mul(w_add(XMAX, XMIN), c->HALF), XX);
        it++;
        if (it < W512_BISECTION_MIN) continue;
        // Hyperbolic lanes continue until the bracket is small enough.
        const V d = w_abs(w_sub(XMAX, XMIN));
        const V tol = w_mul(w_abs(w_add(XMAX, XMIN)), w_set1(W512_BISECTION_EPS));
        k4 = w_kand(w_klt(tol, d), k5);
        if (!w_kany(k4)) break;
    }
    *XXp = XX;
}

// x += d with compensated (Kahan-Babuska) summation: c holds the rounding error of the coordinate, which is added to
// the next increment. Needs floating-point contraction off (as all builds of this file).
REB_W512_INLINE void REB_W512_NAME(w512_add_comp_)(V* x, V* c, const V d){
    const V y = w_add(d, *c);
    const V t = w_add(*x, y);
    const V bp = w_sub(t, *x);
    *c = w_add(w_sub(*x, w_sub(t, bp)), w_sub(y, bp));     // exact: (x + y) - t
    *x = t;
}

// Kepler step; comp = 1: the update is written as increments and added with compensation (compile-time constant).
REB_W512_INLINE void REB_W512_NAME(w512_kepler_core_)(const W512_CONSTS* c, W512_STATE* s, const int comp){
    const V r2 = w_fmadd(s->Z, s->Z, w_fmadd(s->Y, s->Y, w_mul(s->X, s->X)));
    const V v2 = w_fmadd(s->VZ, s->VZ, w_fmadd(s->VY, s->VY, w_mul(s->VX, s->VX)));
    V R = w_sqrt(r2);
    const V IR = w_div(c->ONE, R);
    const V beta = w_fmsub(w_add(c->M, c->M), IR, v2);
    const V eta = w_fmadd(s->Z, s->VZ, w_fmadd(s->Y, s->VY, w_mul(s->X, s->VX)));
    const V zeta = w_fnmadd(beta, R, c->M);
    // Starting value: series reversion of the universal Kepler equation to sixth
    // order in tau = dt/r,
    //   X = tau (1 + P),  P = sum_{i,j} c_ij eta^i v2^j  (0 <= i <= 5, 0 <= j <= 2),
    // where the coefficients c_ij depend on r (and dt, M) only. r, 1/r and the c_ij are available
    // while the previous interaction step is still running; after eta and v2 (which need the new
    // velocities) only a short Estrin evaluation remains. In the dimensionless variables u = eta tau/r,
    // w = v2 tau^2, m = M tau^2/r the series is
    //   P = m/6 + m^2/15 - w (1/6 + 17/120 m) + 3/40 w^2
    //     + u [-1/2 - m/3 - 43/180 m^2 + w (3/8 + 11/20 m) - 5/16 w^2]
    //     + u^2 [1/2 + 5/8 m - 3/4 w] + u^3 [-5/8 - 7/6 m + 35/24 w] + 7/8 u^4 - 21/16 u^5.
    // It replaces the third-order start of item 54: with it, one Halley step (short series) is
    // accurate to < 1e-11 for all eight Solar System planets at dt = 5-6 days, so a single Newton
    // step (which also provides the accurate Stiefel functions) suffices (the third-order start
    // needed a second Newton step for Mercury in 46 % of the steps at dt = 6 d).
    // The start is used where it is a small correction to the first-order guess dt/r
    // (|X - X1| < 0.1 |X1|, finite), otherwise the first-order guess (and the first Halley step).
    // The decision is per lane, so packed systems stay independent. Only the starting value
    // depends on this series; the converged result does not (beyond round-off).
    const V X1 = w_mul(c->DT, IR);
    V Xs;
    {
        const V tau = X1;
        const V tau2 = w_mul(tau, tau);
        const V ku = w_mul(tau, IR);                      // u = ku*eta
        const V m = w_mul(w_mul(c->M, IR), tau2);         // M tau^2 / r
        const V ku2 = w_mul(ku, ku), ku3 = w_mul(ku2, ku);
        const V ku4 = w_mul(ku2, ku2), ku5 = w_mul(ku4, ku);
        const V t4 = w_mul(tau2, tau2);
        const V c00 = w_mul(m, w_fmadd(w_set1(1.0/15.0), m, w_set1(1.0/6.0)));
        const V c01 = w_mul(tau2, w_fmadd(w_set1(-17.0/120.0), m, w_set1(-1.0/6.0)));
        const V c02 = w_mul(t4, w_set1(3.0/40.0));
        const V c10 = w_mul(ku, w_fmadd(m, w_fmadd(w_set1(-43.0/180.0), m, w_set1(-1.0/3.0)), w_set1(-0.5)));
        const V c11 = w_mul(w_mul(ku, tau2), w_fmadd(w_set1(11.0/20.0), m, w_set1(3.0/8.0)));
        const V c12 = w_mul(w_mul(ku, t4), w_set1(-5.0/16.0));
        const V c20 = w_mul(ku2, w_fmadd(w_set1(5.0/8.0), m, c->HALF));
        const V c21 = w_mul(w_mul(ku2, tau2), w_set1(-3.0/4.0));
        const V c30 = w_mul(ku3, w_fmadd(w_set1(-7.0/6.0), m, w_set1(-5.0/8.0)));
        const V c31 = w_mul(w_mul(ku3, tau2), w_set1(35.0/24.0));
        const V c40 = w_mul(ku4, w_set1(7.0/8.0));
        const V c50 = w_mul(ku5, w_set1(-21.0/16.0));
        // Short chain after eta and v2.
        const V v22 = w_mul(v2, v2);
        const V et2 = w_mul(eta, eta);
        const V et4 = w_mul(et2, et2);
        const V row0 = w_fmadd(c02, v22, w_fmadd(c01, v2, c00));
        const V row1 = w_fmadd(c12, v22, w_fmadd(c11, v2, c10));
        const V row2 = w_fmadd(c21, v2, c20);
        const V row3 = w_fmadd(c31, v2, c30);
        const V E01 = w_fmadd(eta, row1, row0);
        const V E23 = w_fmadd(eta, row3, row2);
        const V E45 = w_fmadd(eta, c50, c40);
        const V P = w_fmadd(et4, E45, w_fmadd(et2, E23, E01));
        Xs = w_fmadd(tau, P, tau);
    }
    const K good3 = w_klt(w_abs(w_sub(Xs, X1)), w_mul(w_set1(0.1), w_abs(X1)));   // false for NaN
    V XX = w_sel(good3, Xs, X1);

    V G0, G1, G2, G3;
    // Approximation of 1/f' at the result of the last Halley step (used by the first Newton step),
    // and the lanes in which it is accurate enough to be used.
    V RI = c->ONE;
    K ri_ok = w_kbits(0);
    // Two Halley steps with short series (the first only for lanes without the series start).
    const K need_first = w_knot(good3);
    for (int h=0; h<2; h++){
        if (h==0 && !w_kany(w_kand(need_first, c->k1))) continue;   // every occupied lane has the series start (empty lanes are masked out)
        REB_W512_NAME(w512_stiefel_Gs03_)(h==0 ? 9 : 11, XX, beta, &G0, &G1, &G2, &G3);
        // f summed in the order its terms become available (as in the Newton step): r X - dt right away,
        // G2 one multiplication before G3, so a single FMA follows G3 on the critical path instead of three.
        V f = w_fmsub(R, XX, c->DT);
        f = w_fmadd(eta, G2, f);
        f = w_fmadd(zeta, G3, f);
        V fp = w_fmadd(G2, zeta, R);
        fp = w_fmadd(G1, eta, fp);
        const V fpp = w_fmadd(G1, zeta, w_mul(eta, G0));
        const V hffpp = w_mul(w_mul(f, fpp), c->HALF);
        const V den = w_fmsub(fp, fp, hffpp);
        const V Xn = w_sub(XX, w_div(w_mul(f, fp), den));
        if (h==1){
            // 1/f'(Xn) to first order in the Halley correction dXh: 1/fp - fpp dXh/fp^2, computed in
            // parallel with the Halley division (off the critical path). Relative error ~(fpp dXh/fp)^2
            // plus the truncation of the short series (~1e-10); used only where |fpp dXh/fp| < 1e-4.
            const V ri0 = w_div(c->ONE, fp);
            const V q = w_mul(w_mul(fpp, w_sub(Xn, XX)), ri0);
            RI = w_fnmadd(ri0, q, ri0);
            ri_ok = w_klt(w_abs(q), w_set1(1e-4));      // false for NaN
        }
        XX = (h==0) ? w_sel(need_first, Xn, XX) : Xn;
    }

    // Newton iterations with the accurate series; lanes stop individually (masked: a converged
    // lane keeps its X, so its result does not depend on the other lanes).
    // Correction form X <- X - f(X) RI: in the first iteration RI is the
    // approximation from the Halley step (lanes with ri_ok), so that no division is on the critical
    // path; otherwise, and in later iterations, RI = 1/f'(X) exactly (quadratic convergence). A
    // relative error delta of RI changes X by delta*dX < 1e-8 * EPS.
    // Bias: the Stiefel functions are evaluated again at the final X (below), as in
    // the assembly, and w512_stiefel_Gs13_ corrects for the rounding of the constants 1/k!. The
    // earlier shortcuts (plain Horner without the correction, first-order update of
    // G1..G3 instead of a new evaluation) biased the round-off of the Kepler step: Mercury's
    // semi-major axis drifted by +1.4e-18 per step at dt = 6 d (t = +59 over 64 phases x 1e6 steps),
    // and in a grid of single Kepler orbits the drift statistic <t^2> reached ~20 instead of ~1.2.
    K k4 = c->k1;
    int it = 0;
    const int all_ri_ok = !w_kany(w_kandnot(ri_ok, c->k1));
    while (1){
        const V old = XX;
        REB_W512_NAME(w512_stiefel_Gs13_)(19, XX, beta, &G1, &G2, &G3);
        if (it > 0 || !all_ri_ok){
            // Exact derivative (rare: lanes with a large Halley correction, or a second iteration).
            const V rie = w_div(c->ONE, w_fmadd(zeta, G2, w_fmadd(eta, G1, R)));
            RI = (it > 0) ? rie : w_sel(ri_ok, RI, rie);
        }
        V fN = w_fmadd(eta, G2, w_fmsub(R, XX, c->DT));
        fN = w_fmadd(zeta, G3, fN);
        XX = w_sel(k4, w_fnmadd(fN, RI, XX), XX);
        const V dX = w_sub(XX, old);
        // Not converged: |dX| > EPS, or NaN in an occupied lane (such a lane goes to the bisection
        // fallback instead of being accepted). Empty lanes always count as converged.
        k4 = w_kand(w_kandnot(w_kle(w_abs(dX), c->EPS), k4), c->k1);
        if (!w_kany(k4)) break;
        it++;
        if (it == W512_NEWTON_MAX){
            REB_W512_NAME(w512_kepler_fallback_)(c, s, R, eta, zeta, beta, &XX, k4);
            break;
        }
    }
    // Stiefel functions at the final X (all lanes, including those from the bisection fallback).
    REB_W512_NAME(w512_stiefel_Gs13_)(19, XX, beta, &G1, &G2, &G3);

    // f and g functions.
    R = w_fmadd(G1, eta, R);
    R = w_fmadd(G2, zeta, R);
    const V irn = w_div(c->ONE, R);
    const V mg2 = w_mul(c->M, G2);
    const V nf = w_mul(IR, mg2);                    // -f
    const V ngd = w_mul(mg2, irn);                  // -gd
    const V g = w_fnmadd(G3, c->M, c->DT);
    const V nfd = w_mul(w_mul(w_mul(c->M, G1), IR), irn);  // -fd
    if (comp){
        // x' - x = g v - nf x, v' - v = -nfd x - ngd v (f = 1 - nf, gd = 1 - ngd)
        const V zero = w_zero();
        const V dx = w_maskz(c->k1, w_fmadd(g, s->VX, w_fnmadd(nf, s->X, zero)));
        const V dy = w_maskz(c->k1, w_fmadd(g, s->VY, w_fnmadd(nf, s->Y, zero)));
        const V dz = w_maskz(c->k1, w_fmadd(g, s->VZ, w_fnmadd(nf, s->Z, zero)));
        const V dvx = w_maskz(c->k1, w_fnmadd(nfd, s->X, w_fnmadd(ngd, s->VX, zero)));
        const V dvy = w_maskz(c->k1, w_fnmadd(nfd, s->Y, w_fnmadd(ngd, s->VY, zero)));
        const V dvz = w_maskz(c->k1, w_fnmadd(nfd, s->Z, w_fnmadd(ngd, s->VZ, zero)));
        REB_W512_NAME(w512_add_comp_)(&s->X, &s->CX, dx);
        REB_W512_NAME(w512_add_comp_)(&s->Y, &s->CY, dy);
        REB_W512_NAME(w512_add_comp_)(&s->Z, &s->CZ, dz);
        REB_W512_NAME(w512_add_comp_)(&s->VX, &s->CVX, dvx);
        REB_W512_NAME(w512_add_comp_)(&s->VY, &s->CVY, dvy);
        REB_W512_NAME(w512_add_comp_)(&s->VZ, &s->CVZ, dvz);
        return;
    }
    const V nx = w_maskz(c->k1, w_fmadd(g, s->VX, w_fnmadd(nf, s->X, s->X)));
    const V ny = w_maskz(c->k1, w_fmadd(g, s->VY, w_fnmadd(nf, s->Y, s->Y)));
    const V nz = w_maskz(c->k1, w_fmadd(g, s->VZ, w_fnmadd(nf, s->Z, s->Z)));
    s->VX = w_maskz(c->k1, w_fnmadd(nfd, s->X, w_fnmadd(s->VX, ngd, s->VX)));
    s->VY = w_maskz(c->k1, w_fnmadd(nfd, s->Y, w_fnmadd(s->VY, ngd, s->VY)));
    s->VZ = w_maskz(c->k1, w_fnmadd(nfd, s->Z, w_fnmadd(s->VZ, ngd, s->VZ)));
    s->X = nx;
    s->Y = ny;
    s->Z = nz;
}
REB_W512_STEP_FUNC void REB_W512_NAME(w512_kepler_step_)(const W512_CONSTS* c, W512_STATE* s){
    REB_W512_NAME(w512_kepler_core_)(c, s, 0);
}
REB_W512_STEP_FUNC void REB_W512_NAME(w512_kepler_step_comp_)(const W512_CONSTS* c, W512_STATE* s){
    REB_W512_NAME(w512_kepler_core_)(c, s, 1);
}

///////////////////////////////////////////////////////////////////////////////
// Interaction step

// 1/sqrt(a) for a > 0. Default: a square root and a division. Both are correctly rounded, so all C kernels agree bit
// for bit; fewer operations than the Newton iteration below, and faster on CPUs with fast dividers (Zen 4 +4 %,
// Zen 5 +11 % per step), slower where 512-bit divisions are slow (Sapphire Rapids).
// REB_W512_NEWTON_RSQRT (compile-time, off): integer bit-trick seed and W512_RSQRT_NR Newton steps plus a correction
// (at most 1 ulp, deterministic). REB_W512_HW_RSQRT (experimental, implies it): the hardware estimate as the seed, one
// Newton step fewer; results then differ between instruction sets at round-off level.
#if defined(REB_W512_HW_RSQRT) && !defined(REB_W512_NEWTON_RSQRT)
#define REB_W512_NEWTON_RSQRT
#endif
#ifdef REB_W512_NEWTON_RSQRT
#ifndef W512_RSQRT_NR
#define W512_RSQRT_NR 3
#endif
REB_W512_INLINE V REB_W512_NAME(w512_rsqrt_)(const V a, const V half_a){
    const V th = w_set1(1.5);
    V y = w_rsqrt_seed(a);
    for (int q=0; q<W512_RSQRT_NR; q++){
        y = w_mul(y, w_fnmadd(w_mul(half_a, y), y, th));
    }
    const V e = w_fnmadd(w_mul(a, y), y, w_set1(1.0));   // 1 - a*y^2
    return w_fmadd(w_mul(w_set1(0.5), y), e, y);
}
#define W512_RSQRT(a) REB_W512_NAME(w512_rsqrt_)((a), w_mul(w_set1(0.5), (a)))
#else
#define W512_RSQRT(a) w_div(w_set1(1.0), w_sqrt(a))
#endif

// 8x8 matrix (stored transposed: element (i,k) at mat[i+8k]) times three vectors (mat8_mul3).
REB_W512_INLINE void REB_W512_NAME(w512_mat8_mul3_)(const double* mat, const V a, const V b, const V cc, V* oa, V* ob, V* oc){
    REB_ATTRIBUTE_ALIGNED_64 double ia[8];
    REB_ATTRIBUTE_ALIGNED_64 double ib[8];
    REB_ATTRIBUTE_ALIGNED_64 double ic[8];
    w_store(ia, a);
    w_store(ib, b);
    w_store(ic, cc);
    V row0 = w_load(mat);
    V row1 = w_load(mat+8);
    V o0 = w_mul(row0, w_set1(ia[0]));
    V o1 = w_mul(row0, w_set1(ib[0]));
    V o2 = w_mul(row0, w_set1(ic[0]));
    V t0 = w_mul(row1, w_set1(ia[1]));
    V t1 = w_mul(row1, w_set1(ib[1]));
    V t2 = w_mul(row1, w_set1(ic[1]));
    for (int k=2; k<W_LANES; k+=2){
        row0 = w_load(mat+8*k);
        row1 = w_load(mat+8*(k+1));
        o0 = w_fmadd(row0, w_set1(ia[k]), o0);
        o1 = w_fmadd(row0, w_set1(ib[k]), o1);
        o2 = w_fmadd(row0, w_set1(ic[k]), o2);
        t0 = w_fmadd(row1, w_set1(ia[k+1]), t0);
        t1 = w_fmadd(row1, w_set1(ib[k+1]), t1);
        t2 = w_fmadd(row1, w_set1(ic[k+1]), t2);
    }
    *oa = w_add(o0, t0);
    *ob = w_add(o1, t1);
    *oc = w_add(o2, t2);
}

// multiplier/r^3 for separations (d0, d1, d2); with_mult=0 means multiplier 1.
REB_W512_INLINE V REB_W512_NAME(w512_gravity_prefactor_)(const W512_CONSTS* c, const V d0, const V d1, const V d2, const int with_mult, const V mult, const int encounter, int* status){
    const V r2 = w_fmadd(d2, d2, w_fmadd(d1, d1, w_mul(d0, d0)));
    const V y = W512_RSQRT(r2);
    if (encounter){
        if (w_kany(w_kand(w_kgt(y, w_load(c->data->exit_min_distance_r.lanes)), c->k1))) *status = 3; // REB_STATUS_ENCOUNTER
    }
    const V y2 = w_mul(y, y);
    if (with_mult){
        return w_mul(y2, w_mul(y, mult));
    }else{
        return w_mul(y, y2);
    }
}

// Interaction step; comp = 1: the velocity kick is added with compensation (compile-time constant).
REB_W512_INLINE void REB_W512_NAME(w512_interaction_core_)(const W512_CONSTS* c, W512_STATE* s, const int gr, const int nsys, const int encounter, const int escape, int* status, const int comp){
    const struct simd_data* const data = c->data;
    // Jacobi term in Jacobi coordinates (undoes the Kepler-step force of the interior mass).
    V JTX = w_zero(), JTY = w_zero(), JTZ = w_zero();    // comp: the Jacobi-term kick, added with the interaction kick
    {
        const V r2 = w_fmadd(s->Z, s->Z, w_fmadd(s->Y, s->Y, w_mul(s->X, s->X)));
        const V y = W512_RSQRT(r2);
        const V a = w_mul(w_mul(c->DT, c->M), w_mul(w_mul(y, y), y));    // dt*M/r^3
        if (comp){
            JTX = w_maskz(c->k1, w_mul(a, s->X));
            JTY = w_maskz(c->k1, w_mul(a, s->Y));
            JTZ = w_maskz(c->k1, w_mul(a, s->Z));
        }else{
            s->VX = w_maskz(c->k1, w_fmadd(a, s->X, s->VX));
            s->VY = w_maskz(c->k1, w_fmadd(a, s->Y, s->VY));
            s->VZ = w_maskz(c->k1, w_fmadd(a, s->Z, s->VZ));
        }
    }
    V HX, HY, HZ;
    REB_W512_NAME(w512_mat8_mul3_)(data->mat8_jacobi_to_heliocentric, s->X, s->Y, s->Z, &HX, &HY, &HZ);

    const V hr2 = w_fmadd(HZ, HZ, w_fmadd(HY, HY, w_mul(HX, HX)));
    if (escape){
        const V dmax = w_load(data->exit_max_distance.lanes);
        if (w_kany(w_kand(w_kgt(hr2, w_mul(dmax, dmax)), c->k1))) *status = 4; // REB_STATUS_ESCAPE
    }
    if (encounter){
        const V dmin = w_load(data->exit_min_distance.lanes);
        if (w_kany(w_kand(w_klt(hr2, w_mul(dmin, dmin)), c->k1))) *status = 3; // REB_STATUS_ENCOUNTER
    }
    // Central-body term in heliocentric coordinates.
    const V yh = W512_RSQRT(hr2);
    const V yh2 = w_mul(yh, yh);
    const V j = w_maskz(c->k1, w_mul(c->MM0_DT, w_mul(yh2, yh)));   // -m0*dt/r^3
    V HVX = w_mul(HX, j);
    V HVY = w_mul(HY, j);
    V HVZ = w_mul(HZ, j);
    if (gr){
        const V g = w_mul(c->DT, w_load(data->gr_prefac.lanes));
        const V gg = w_maskz(c->k1, w_mul(g, w_mul(yh2, yh2)));     // g/r^4
        HVX = w_maskz(c->k1, w_fmadd(HX, gg, HVX));
        HVY = w_maskz(c->k1, w_fmadd(HY, gg, HVY));
        HVZ = w_maskz(c->k1, w_fmadd(HZ, gg, HVZ));
    }
    const V dm = w_mul(c->DT, w_load(data->m.lanes));   // dt*m
    const V one = c->ONE;

    // Planet-planet terms. Lanes are paired by permutations; each pair is evaluated once and
    // both directions are applied (blocks as in the assembly).
    if (nsys < 4){
        const V d0 = w_sub(HX, w_perm_3201(HX));
        const V d1 = w_sub(HY, w_perm_3201(HY));
        const V d2 = w_sub(HZ, w_perm_3201(HZ));
        const V pf = REB_W512_NAME(w512_gravity_prefactor_)(c, d0, d1, d2, 0, one, encounter, status);
        const V w = w_mul(w_perm_3201(dm), pf);
        HVX = w_fnmadd(d0, w, HVX);
        HVY = w_fnmadd(d1, w, HVY);
        HVZ = w_fnmadd(d2, w, HVZ);
        const V w2 = w_perm_2310(w_mul(dm, pf));
        HVX = w_fmadd(w_perm_2310(d0), w2, HVX);
        HVY = w_fmadd(w_perm_2310(d1), w2, HVY);
        HVZ = w_fmadd(w_perm_2310(d2), w2, HVZ);
    }
    {
        const V d0 = w_sub(HX, w_swap_pairs(HX));
        const V d1 = w_sub(HY, w_swap_pairs(HY));
        const V d2 = w_sub(HZ, w_swap_pairs(HZ));
        const V pf = REB_W512_NAME(w512_gravity_prefactor_)(c, d0, d1, d2, 1, w_swap_pairs(dm), encounter, status);
        HVX = w_fnmadd(d0, pf, HVX);
        HVY = w_fnmadd(d1, pf, HVY);
        HVZ = w_fnmadd(d2, pf, HVZ);
    }
    V A0, A1, A2;
    if (nsys == 1){
        V CX, CY, CZ;
        {
            const V d0 = w_sub(HX, w_perm_b3(HX));
            const V d1 = w_sub(HY, w_perm_b3(HY));
            const V d2 = w_sub(HZ, w_perm_b3(HZ));
            const V pf = REB_W512_NAME(w512_gravity_prefactor_)(c, d0, d1, d2, 0, one, encounter, status);
            const V w = w_mul(w_perm_b3(dm), pf);
            HVX = w_fnmadd(d0, w, HVX);
            HVY = w_fnmadd(d1, w, HVY);
            HVZ = w_fnmadd(d2, w, HVZ);
            const V w2 = w_mul(dm, pf);
            CX = w_mul(d0, w2);
            CY = w_mul(d1, w2);
            CZ = w_mul(d2, w2);
        }
        {
            const V d0 = w_sub(HX, w_perm_b4(HX));
            const V d1 = w_sub(HY, w_perm_b4(HY));
            const V d2 = w_sub(HZ, w_perm_b4(HZ));
            const V pf = REB_W512_NAME(w512_gravity_prefactor_)(c, d0, d1, d2, 0, one, encounter, status);
            const V w = w_mul(w_perm_b4(dm), pf);
            HVX = w_fnmadd(d0, w, HVX);
            HVY = w_fnmadd(d1, w, HVY);
            HVZ = w_fnmadd(d2, w, HVZ);
            const V w2 = w_perm_3012(w_mul(dm, pf));
            CX = w_fmadd(w_perm_3012(d0), w2, CX);
            CY = w_fmadd(w_perm_3012(d1), w2, CY);
            CZ = w_fmadd(w_perm_3012(d2), w2, CZ);
        }
        A0 = w_maskz(c->k1, w_add(HVX, w_perm_b34(CX)));
        A1 = w_maskz(c->k1, w_add(HVY, w_perm_b34(CY)));
        A2 = w_maskz(c->k1, w_add(HVZ, w_perm_b34(CZ)));
    }else{
        A0 = w_maskz(c->k1, HVX);
        A1 = w_maskz(c->k1, HVY);
        A2 = w_maskz(c->k1, HVZ);
    }
    // Velocity kicks from heliocentric to Jacobi.
    V J0, J1, J2;
    REB_W512_NAME(w512_mat8_mul3_)(data->mat8_inertial_to_jacobi, A0, A1, A2, &J0, &J1, &J2);
    if (comp){
        REB_W512_NAME(w512_add_comp_)(&s->VX, &s->CVX, w_add(J0, JTX));
        REB_W512_NAME(w512_add_comp_)(&s->VY, &s->CVY, w_add(J1, JTY));
        REB_W512_NAME(w512_add_comp_)(&s->VZ, &s->CVZ, w_add(J2, JTZ));
        return;
    }
    s->VX = w_add(J0, s->VX);
    s->VY = w_add(J1, s->VY);
    s->VZ = w_add(J2, s->VZ);
}
REB_W512_STEP_FUNC void REB_W512_NAME(w512_interaction_step_)(const W512_CONSTS* c, W512_STATE* s, const int gr, const int nsys, const int encounter, const int escape, int* status){
    REB_W512_NAME(w512_interaction_core_)(c, s, gr, nsys, encounter, escape, status, 0);
}
REB_W512_STEP_FUNC void REB_W512_NAME(w512_interaction_step_comp_)(const W512_CONSTS* c, W512_STATE* s, const int gr, const int nsys, const int encounter, const int escape, int* status){
    REB_W512_NAME(w512_interaction_core_)(c, s, gr, nsys, encounter, escape, status, 1);
}

///////////////////////////////////////////////////////////////////////////////
// Entry points

REB_W512_INLINE void REB_W512_NAME(w512_load_)(const struct simd_data* data, W512_CONSTS* c, W512_STATE* s){
    c->data = data;
    c->k1 = w_kbits(data->mask);
    c->DT = w_load(data->dt.lanes);
    c->M = w_load(data->M.lanes);
    c->MM0_DT = w_sub(w_zero(), w_mul(w_load(data->M0.lanes), c->DT)); // -(m0*dt): sign flip, exact
    c->ONE = w_set1(1.0);
    c->HALF = w_set1(0.5);
    c->EPS = w_set1(W512_EPS);
    s->X = w_load(data->x.lanes);
    s->Y = w_load(data->y.lanes);
    s->Z = w_load(data->z.lanes);
    s->VX = w_load(data->vx.lanes);
    s->VY = w_load(data->vy.lanes);
    s->VZ = w_load(data->vz.lanes);
    s->CX = s->CY = s->CZ = s->CVX = s->CVY = s->CVZ = w_zero();
}

REB_W512_INLINE void REB_W512_NAME(w512_store_)(struct simd_data* data, const W512_STATE* s){
    w_store(data->x.lanes, s->X);
    w_store(data->y.lanes, s->Y);
    w_store(data->z.lanes, s->Z);
    w_store(data->vx.lanes, s->VX);
    w_store(data->vy.lanes, s->VY);
    w_store(data->vz.lanes, s->VZ);
}

// Main loop for one combination of options (all arguments except data/N_steps/skip/sigint are
// compile-time constants after inlining).
REB_W512_INLINE int REB_W512_NAME(w512_full_steps_)(struct simd_data* data, uint64_t* N_steps, const int skip_first_kepler_step, volatile sig_atomic_t* sigint, const int gr, const int nsys, const int encounter, const int escape){
    W512_CONSTS c;
    W512_STATE s;
    REB_W512_NAME(w512_load_)(data, &c, &s);
    int64_t remaining = (int64_t)*N_steps;
    int status = 0;
    int ret;
    int skip = skip_first_kepler_step;
    while (1){
        if (!skip){
            REB_W512_NAME(w512_kepler_step_)(&c, &s);
        }
        skip = 0;
        REB_W512_NAME(w512_interaction_step_)(&c, &s, gr, nsys, encounter, escape, &status);
        if (*sigint){
            ret = 6; // REB_STATUS_SIGINT
            remaining--;
            *N_steps -= (uint64_t)remaining;
            break;
        }
        if ((encounter || escape) && status){
            ret = status;
            remaining--;
            *N_steps -= (uint64_t)remaining;
            break;
        }
        remaining--;
        if (remaining > 0) continue;
        ret = -1; // REB_STATUS_RUNNING
        break;
    }
    REB_W512_NAME(w512_store_)(data, &s);
    return ret;
}

// Compensated main loop (option `compensated`): as w512_full_steps_, with the coordinate updates of the Kepler and
// interaction steps added with compensated summation. The low-order parts live in registers for the duration of the
// call and are added to the coordinates (rounded) at its end. encounter/escape are run-time flags here.
REB_W512_INLINE int REB_W512_NAME(w512_full_steps_comp_)(struct simd_data* data, uint64_t* N_steps, const int skip_first_kepler_step, volatile sig_atomic_t* sigint, const int gr, const int nsys, const int encounter, const int escape){
    W512_CONSTS c;
    W512_STATE s;
    REB_W512_NAME(w512_load_)(data, &c, &s);
    int64_t remaining = (int64_t)*N_steps;
    int status = 0;
    int ret;
    int skip = skip_first_kepler_step;
    while (1){
        if (!skip){
            REB_W512_NAME(w512_kepler_step_comp_)(&c, &s);
        }
        skip = 0;
        REB_W512_NAME(w512_interaction_step_comp_)(&c, &s, gr, nsys, encounter, escape, &status);
        if (*sigint){
            ret = 6; // REB_STATUS_SIGINT
            remaining--;
            *N_steps -= (uint64_t)remaining;
            break;
        }
        if ((encounter || escape) && status){
            ret = status;
            remaining--;
            *N_steps -= (uint64_t)remaining;
            break;
        }
        remaining--;
        if (remaining > 0) continue;
        ret = -1; // REB_STATUS_RUNNING
        break;
    }
    s.X = w_add(s.X, s.CX); s.Y = w_add(s.Y, s.CY); s.Z = w_add(s.Z, s.CZ);
    s.VX = w_add(s.VX, s.CVX); s.VY = w_add(s.VY, s.CVY); s.VZ = w_add(s.VZ, s.CVZ);
    REB_W512_NAME(w512_store_)(data, &s);
    return ret;
}

// Out-of-line Kepler and interaction steps for the paired loop (eight-lane instantiations; the
// four-lane ones define W_PAIR_INLINE and inline them, because there the two systems must overlap inside one
// function) (one copy of the Kepler step and one
// interaction step per gr/nsys/encounter/escape combination used): the out-of-order core still
// overlaps the two simulations' chains across the calls, and the paired kernels stay small.
#if defined(_MSC_VER) && !defined(__clang__)
#define REB_W512_NOINLINE static __declspec(noinline)
#else
#define REB_W512_NOINLINE static __attribute__((noinline, unused))
#endif
REB_W512_NOINLINE void REB_W512_NAME(w512_kepler_step_ool_)(const W512_CONSTS* c, W512_STATE* s){
    REB_W512_NAME(w512_kepler_step_)(c, s);
}
#define W512_INTERACTION_OOL(GR, NSYS) \
REB_W512_NOINLINE void REB_W512_CAT(REB_W512_NAME(w512_interaction_ool_), GR##_##NSYS)(const W512_CONSTS* c, W512_STATE* s, const int encounter, const int escape, int* status){ \
    REB_W512_NAME(w512_interaction_step_)(c, s, GR, NSYS, encounter, escape, status); \
}
W512_INTERACTION_OOL(0, 1) W512_INTERACTION_OOL(0, 2) W512_INTERACTION_OOL(0, 4)
W512_INTERACTION_OOL(1, 1) W512_INTERACTION_OOL(1, 2) W512_INTERACTION_OOL(1, 4)
#define W512_INTERACTION_OOL_CALL(c, s, st) do { \
    if (gr==0 && nsys==1) REB_W512_CAT(REB_W512_NAME(w512_interaction_ool_), 0_1)(c, s, encounter, escape, st); \
    else if (gr==0 && nsys==2) REB_W512_CAT(REB_W512_NAME(w512_interaction_ool_), 0_2)(c, s, encounter, escape, st); \
    else if (gr==0) REB_W512_CAT(REB_W512_NAME(w512_interaction_ool_), 0_4)(c, s, encounter, escape, st); \
    else if (nsys==1) REB_W512_CAT(REB_W512_NAME(w512_interaction_ool_), 1_1)(c, s, encounter, escape, st); \
    else if (nsys==2) REB_W512_CAT(REB_W512_NAME(w512_interaction_ool_), 1_2)(c, s, encounter, escape, st); \
    else REB_W512_CAT(REB_W512_NAME(w512_interaction_ool_), 1_4)(c, s, encounter, escape, st); \
} while (0)
// Two independent simulations per call (item 49): the steps of A and B alternate in one loop,
// so the out-of-order core overlaps the two dependency chains (each step is one long chain of
// dependent operations; the vector units are idle most of the time with a single chain). Each
// simulation gets exactly the operations of a single-simulation call, so the results are
// bitwise identical to stepping them separately. Both must use the same options (gr, nsys,
// encounter, escape). Returns a bit mask: bit 0 set if A stopped (status in *status_a), bit 1
// if B stopped; *N_done is the number of steps both completed.
REB_W512_INLINE int REB_W512_NAME(w512_full_steps_pair_)(struct simd_data* da, struct simd_data* db, uint64_t N_steps, uint64_t* N_done, const int skip_a, const int skip_b, volatile sig_atomic_t* sigint, int* status_a, int* status_b, const int gr, const int nsys, const int encounter, const int escape){
    W512_CONSTS ca, cb;
    W512_STATE sa, sb;
    REB_W512_NAME(w512_load_)(da, &ca, &sa);
    REB_W512_NAME(w512_load_)(db, &cb, &sb);
    int sta = 0, stb = 0;
    uint64_t n = 0;
    int stopped = 0;
    int ska = skip_a, skb = skip_b;
    while (n < N_steps){
#ifdef W_PAIR_INLINE
        // Small (four-lane) instantiations: steps inlined, the compiler interleaves both chains.
        if (!ska) REB_W512_NAME(w512_kepler_step_)(&ca, &sa);
        if (!skb) REB_W512_NAME(w512_kepler_step_)(&cb, &sb);
        ska = 0; skb = 0;
        REB_W512_NAME(w512_interaction_step_)(&ca, &sa, gr, nsys, encounter, escape, &sta);
        REB_W512_NAME(w512_interaction_step_)(&cb, &sb, gr, nsys, encounter, escape, &stb);
#else
        if (!ska) REB_W512_NAME(w512_kepler_step_ool_)(&ca, &sa);
        if (!skb) REB_W512_NAME(w512_kepler_step_ool_)(&cb, &sb);
        ska = 0; skb = 0;
        W512_INTERACTION_OOL_CALL(&ca, &sa, &sta);
        W512_INTERACTION_OOL_CALL(&cb, &sb, &stb);
#endif
        n++;
        if (*sigint){ sta = 6; stb = 6; stopped = 3; break; }
        if ((encounter || escape) && (sta || stb)){
            stopped = (sta ? 1 : 0) | (stb ? 2 : 0);
            break;
        }
    }
    REB_W512_NAME(w512_store_)(da, &sa);
    REB_W512_NAME(w512_store_)(db, &sb);
    *N_done = n;
    *status_a = sta;
    *status_b = stb;
    return stopped;
}

REB_W512_TARGET int REB_W512_CAT(reb_whfast512c_full_steps_pair_, REB_W512_SUFFIX)(struct simd_data* da, struct simd_data* db, uint64_t N_steps, uint64_t* N_done, int skip_a, int skip_b, volatile sig_atomic_t* sigint, int* status_a, int* status_b, int gr, int nsys, int encounter, int escape){
    gr = gr ? 1 : 0;
    encounter = encounter ? 1 : 0;
    escape = escape ? 1 : 0;
#ifdef W_ONLY_NSYS2
    nsys = 2;
#else
    if (nsys != 2 && nsys != 4) nsys = 1;
#endif
    // Specialised on gr and nsys only; the encounter/escape checks are run-time flags here (a
    // compare per step), which keeps code size and compile time small (24 -> 6 copies).
#define W512_PAIR_CASE(GR, NSYS) \
    if (gr==GR && nsys==NSYS) \
        return REB_W512_NAME(w512_full_steps_pair_)(da, db, N_steps, N_done, skip_a, skip_b, sigint, status_a, status_b, GR, NSYS, encounter, escape);
#ifdef W_ONLY_NSYS2
    W512_PAIR_CASE(0, 2) W512_PAIR_CASE(1, 2)
#else
    W512_PAIR_CASE(0, 1) W512_PAIR_CASE(0, 2) W512_PAIR_CASE(0, 4)
    W512_PAIR_CASE(1, 1) W512_PAIR_CASE(1, 2) W512_PAIR_CASE(1, 4)
#endif
#undef W512_PAIR_CASE
    return -1;
}

#define W512_FULL_STEPS_CASE(GR, NSYS, ENC, ESC) \
    if (gr==GR && nsys==NSYS && encounter==ENC && escape==ESC) \
        return REB_W512_NAME(w512_full_steps_)(data, N_steps, skip_first_kepler_step, sigint, GR, NSYS, ENC, ESC);
#define W512_FULL_STEPS_CASES_NSYS(GR, NSYS) \
    W512_FULL_STEPS_CASE(GR, NSYS, 0, 0) W512_FULL_STEPS_CASE(GR, NSYS, 1, 0) \
    W512_FULL_STEPS_CASE(GR, NSYS, 0, 1) W512_FULL_STEPS_CASE(GR, NSYS, 1, 1)

REB_W512_TARGET int REB_W512_CAT(reb_whfast512c_full_steps_, REB_W512_SUFFIX)(struct simd_data* data, uint64_t* N_steps, int skip_first_kepler_step, volatile sig_atomic_t* sigint, int gr, int nsys, int encounter, int escape){
    gr = gr ? 1 : 0;
    encounter = encounter ? 1 : 0;
    escape = escape ? 1 : 0;
    if (nsys != 2 && nsys != 4) nsys = 1;
#ifdef W_ONLY_NSYS2
    nsys = 2;
    W512_FULL_STEPS_CASES_NSYS(0, 2)
    W512_FULL_STEPS_CASES_NSYS(1, 2)
#else
    W512_FULL_STEPS_CASES_NSYS(0, 1)
    W512_FULL_STEPS_CASES_NSYS(0, 2)
    W512_FULL_STEPS_CASES_NSYS(0, 4)
    W512_FULL_STEPS_CASES_NSYS(1, 1)
    W512_FULL_STEPS_CASES_NSYS(1, 2)
    W512_FULL_STEPS_CASES_NSYS(1, 4)
#endif
    return -1;
}

REB_W512_TARGET int REB_W512_CAT(reb_whfast512c_full_steps_comp_, REB_W512_SUFFIX)(struct simd_data* data, uint64_t* N_steps, int skip_first_kepler_step, volatile sig_atomic_t* sigint, int gr, int nsys, int encounter, int escape){
    gr = gr ? 1 : 0;
    encounter = encounter ? 1 : 0;
    escape = escape ? 1 : 0;
    if (nsys != 2 && nsys != 4) nsys = 1;
#define W512_COMP_CASE(GR, NSYS) \
    if (gr==GR && nsys==NSYS) \
        return REB_W512_NAME(w512_full_steps_comp_)(data, N_steps, skip_first_kepler_step, sigint, GR, NSYS, encounter, escape);
#ifdef W_ONLY_NSYS2
    nsys = 2;
    W512_COMP_CASE(0, 2) W512_COMP_CASE(1, 2)
#else
    W512_COMP_CASE(0, 1) W512_COMP_CASE(0, 2) W512_COMP_CASE(0, 4)
    W512_COMP_CASE(1, 1) W512_COMP_CASE(1, 2) W512_COMP_CASE(1, 4)
#endif
#undef W512_COMP_CASE
    return -1;
}

REB_W512_TARGET void REB_W512_CAT(reb_whfast512c_kepler_step_, REB_W512_SUFFIX)(struct simd_data* data){
    W512_CONSTS c;
    W512_STATE s;
    REB_W512_NAME(w512_load_)(data, &c, &s);
    REB_W512_NAME(w512_kepler_step_)(&c, &s);
    REB_W512_NAME(w512_store_)(data, &s);
}

// Symplectic corrector (17th order): alternating Kepler (four sub-steps) and interaction steps.
REB_W512_INLINE void REB_W512_NAME(w512_corrector_)(struct simd_data* data, const double inv, const int gr, const int nsys){
    W512_CONSTS c;
    W512_STATE s;
    REB_W512_NAME(w512_load_)(data, &c, &s);
    const V dt = w_load(data->dt.lanes);
    const V M0 = w_load(data->M0.lanes);
    int status = 0;
    for (int i=0; i<63; i++){
        const double coef = REB_W512_NAME(w512_bits_)(REB_W512_NAME(w512_corrector17_bits_)[i]);
        if (i%2 == 0){ // Kepler
            c.DT = w_mul(w_mul(w_mul(w_set1(coef), dt), c.HALF), c.HALF);
            for (int k=0; k<4; k++){
                REB_W512_NAME(w512_kepler_step_)(&c, &s);
            }
        }else{ // Interaction
            c.DT = w_mul(w_mul(w_set1(coef), w_set1(inv)), dt);
            c.MM0_DT = w_sub(w_zero(), w_mul(M0, c.DT));
            REB_W512_NAME(w512_interaction_step_)(&c, &s, gr, nsys, 0, 0, &status);
        }
    }
    REB_W512_NAME(w512_store_)(data, &s);
}

REB_W512_TARGET void REB_W512_CAT(reb_whfast512c_corrector_step_, REB_W512_SUFFIX)(struct simd_data* data, double inv, int gr, int nsys){
    gr = gr ? 1 : 0;
    if (nsys != 2 && nsys != 4) nsys = 1;
#ifdef W_ONLY_NSYS2
    nsys = 2;
    if (gr){
        REB_W512_NAME(w512_corrector_)(data, inv, 1, 2);
    }else{
        REB_W512_NAME(w512_corrector_)(data, inv, 0, 2);
    }
    return;
#endif
    if (gr){
        if (nsys==1) REB_W512_NAME(w512_corrector_)(data, inv, 1, 1);
        if (nsys==2) REB_W512_NAME(w512_corrector_)(data, inv, 1, 2);
        if (nsys==4) REB_W512_NAME(w512_corrector_)(data, inv, 1, 4);
    }else{
        if (nsys==1) REB_W512_NAME(w512_corrector_)(data, inv, 0, 1);
        if (nsys==2) REB_W512_NAME(w512_corrector_)(data, inv, 0, 2);
        if (nsys==4) REB_W512_NAME(w512_corrector_)(data, inv, 0, 4);
    }
}
