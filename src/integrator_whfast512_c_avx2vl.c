/**
 * @file    integrator_whfast512_c_avx2vl.c
 * @brief   WHFast512 portable C kernel, AVX2 instantiation compiled with AVX-512VL.
 * @details The kernel of integrator_whfast512_c_avx2.c (an 8-lane vector as two 256-bit registers), compiled with
 *          AVX-512VL: the EVEX encoding gives it 32 vector registers instead of 16. Same operations in the same
 *          order, so the results are bitwise identical to the other C kernels. The driver uses it for kernel 3
 *          (C AVX2) on CPUs with AVX-512VL, where the 256-bit kernel can be faster than the 512-bit one (AMD Zen 4).
 *
 * @section LICENSE
 * Copyright (c) 2026 Michael Hippke. Part of Parachron (https://github.com/hippke/parachron), a modified version
 * of REBOUND. GNU General Public License v3 or later (see LICENSE).
 */
#define REB_W512_AVX2_VL
#include "integrator_whfast512_c_avx2.c"
