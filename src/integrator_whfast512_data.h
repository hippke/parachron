/**
 * @file    integrator_whfast512_data.h
 * @brief   Data layout shared by the WHFast512 driver and its kernels.
 * @details The driver (integrator_whfast512.c) fills a `struct simd_data`; the kernels
 *          (hand-written AVX-512 assembly in integrator_whfast512.s, and the portable C
 *          kernels in integrator_whfast512_c_*.c) read and update it. The layout must not
 *          change: the assembly uses fixed offsets, and the structure is written to
 *          Simulationarchive files.
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
#ifndef _INTEGRATOR_WHFAST512_DATA_H
#define _INTEGRATOR_WHFAST512_DATA_H

#include <stdint.h>
#include <signal.h>

#if !defined(_WIN32)
#define REB_ATTRIBUTE_ALIGNED_64 __attribute__((aligned(64)))
#else
#define REB_ATTRIBUTE_ALIGNED_64 __declspec(align(64))
#endif

// Eight doubles, 64-byte aligned (one AVX-512 register). Plain struct so that the driver
// does not need to be compiled with any vector extensions.
typedef struct {
    double lanes[8];
} REB_ATTRIBUTE_ALIGNED_64 reb_simd8d;

// The main data structure. Passed as a pointer to the kernels.
struct simd_data{
    // Various constants
    REB_ATTRIBUTE_ALIGNED_64 reb_simd8d M;             //  Masses used in Kepler-Solver
    REB_ATTRIBUTE_ALIGNED_64 reb_simd8d dt;            //  Timestep
    REB_ATTRIBUTE_ALIGNED_64 reb_simd8d gr_prefac;     //  Prefactor for GR
    REB_ATTRIBUTE_ALIGNED_64 reb_simd8d m;
    REB_ATTRIBUTE_ALIGNED_64 reb_simd8d x;
    REB_ATTRIBUTE_ALIGNED_64 reb_simd8d y;
    REB_ATTRIBUTE_ALIGNED_64 reb_simd8d z;
    REB_ATTRIBUTE_ALIGNED_64 reb_simd8d vx;
    REB_ATTRIBUTE_ALIGNED_64 reb_simd8d vy;
    REB_ATTRIBUTE_ALIGNED_64 reb_simd8d vz;
    REB_ATTRIBUTE_ALIGNED_64 double mat8_inertial_to_jacobi[64]; // Coordinate transformation matricies. Can be recalculated from particle masses.
    REB_ATTRIBUTE_ALIGNED_64 double mat8_jacobi_to_heliocentric[64];
    REB_ATTRIBUTE_ALIGNED_64 reb_simd8d M0;            //  Masses used in Jacobi Term
    REB_ATTRIBUTE_ALIGNED_64 unsigned char mask;       // Mask for cases with less than 8 planets
    REB_ATTRIBUTE_ALIGNED_64 reb_simd8d exit_max_distance;
    REB_ATTRIBUTE_ALIGNED_64 reb_simd8d exit_min_distance;
    REB_ATTRIBUTE_ALIGNED_64 reb_simd8d exit_min_distance_r; // Inverse of exit_min_distance
    REB_ATTRIBUTE_ALIGNED_64 double mat8_jacobi_to_inertial[64];
#ifdef DEBUG_AVX512
    REB_ATTRIBUTE_ALIGNED_64 uint64_t counter[8];
#endif // DEBUG_AVX512
};

// Kernels. Values of reb_integrator_whfast512_state.kernel.
enum REB_WHFAST512_KERNEL {
    REB_WHFAST512_KERNEL_AUTO = 0,      // Fastest available on this machine (default)
    REB_WHFAST512_KERNEL_ASM = 1,       // Hand-written AVX-512 assembly (x86-64, AVX-512F/DQ)
    REB_WHFAST512_KERNEL_C_AVX512 = 2,  // Portable C kernel, AVX-512 instantiation
    REB_WHFAST512_KERNEL_C_AVX2 = 3,    // Portable C kernel, AVX2+FMA instantiation
    REB_WHFAST512_KERNEL_C_GENERIC = 4, // Portable C kernel, plain C (any CPU)
};

// Portable C kernels (integrator_whfast512_c_avx512.c, _c_avx2.c, _c_generic.c).
// All three instantiate the same source (integrator_whfast512_kernel.h) with explicit
// fused multiply-adds, so they give bitwise identical results on any machine.
// Return values: REB_STATUS_RUNNING (-1), REB_STATUS_ENCOUNTER (3), REB_STATUS_ESCAPE (4), REB_STATUS_SIGINT (6).
#define REB_WHFAST512_C_KERNEL_DECLARE(suffix) \
    int reb_whfast512c_full_steps_##suffix(struct simd_data* data, uint64_t* N_steps, int skip_first_kepler_step, volatile sig_atomic_t* sigint, int gr, int nsys, int encounter, int escape); \
    void reb_whfast512c_kepler_step_##suffix(struct simd_data* data); \
    void reb_whfast512c_corrector_step_##suffix(struct simd_data* data, double inv, int gr, int nsys);
REB_WHFAST512_C_KERNEL_DECLARE(avx512)
REB_WHFAST512_C_KERNEL_DECLARE(avx2)
REB_WHFAST512_C_KERNEL_DECLARE(generic)

// 1 if the C kernel for this instruction set was compiled in (x86 compilers that support the target).
int reb_whfast512c_avx512_compiled(void);
int reb_whfast512c_avx2_compiled(void);

#endif // _INTEGRATOR_WHFAST512_DATA_H
