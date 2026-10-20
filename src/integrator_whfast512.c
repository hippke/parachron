/**
 * integrator_whfast512.c: ASM version of WHFast512 in Jacobi Coordinates
 * 
 * Copyright (c) 2026 Rishit Dagli, Hanno Rein, Pejvak Javaheri
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

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include "rebound.h"
#include "rebound_internal.h"
#include "integrator_whfast512_data.h"
#if defined(__GNUC__) || defined(__clang__)
#define REB_SYSV_ABI __attribute__((sysv_abi))
#else // defined(__GNUC__) || defined(__clang__)
#define REB_SYSV_ABI 
#endif // defined(__GNUC__) || defined(__clang__)
#include "particle.h"
#include "tools.h"
#include "gravity.h"
#include "boundary.h"
#include "integrator_whfast.h"
#include "integrator_whfast512.h"

//#define DEBUG_AVX512 1

void reb_integrator_whfast512_free(void* state);		
void* reb_integrator_whfast512_create();		
void reb_integrator_whfast512_step(struct reb_simulation* r, void* state);
void reb_integrator_whfast512_synchronize(struct reb_simulation* r, void* state);
const struct reb_binarydata_field_descriptor reb_integrator_whfast512_field_descriptor_list[];

// The main data structure (struct simd_data) is defined in integrator_whfast512_data.h. It is
// passed as a pointer to the kernels: the AVX-512 assembly (integrator_whfast512.s) or the
// portable C kernels (integrator_whfast512_c_*.c).

// The assembly kernel exists for x86-64 with GCC/clang-style toolchains (also clang on Windows).
#if (defined(__x86_64__) || defined(_M_AMD64)) && ((defined(__GNUC__) || defined(__clang__)) || !defined(_WIN32)) && !defined(REB_WHFAST512_NO_ASM)
#define REB_WHFAST512_HAVE_ASM 1
#endif

const struct reb_integrator reb_integrator_whfast512 = {
    .documentation =
    "WHFast512 is a highly optimized implementation of the symplectic [Wisdom & Holman (1991)] integrator. " 
    "It supports simulations with up to 9 particles (8 planets + 1 central object). " 
    "Note that by default WHFast512 combines 1e6 timesteps to improve speed. "
    "You can set `concatenate_steps` to a lower value for more fine grained output. "
    "WHFast512 uses Jacobi coordinates and support symplectic correctors as well as general relativistic corrections. "
    "\n\n"
    "WHFast512 runs on any CPU. The `kernel` setting selects the implementation: "
    "by default (0) the hand-written AVX-512 assembly is used if the CPU supports AVX-512, "
    "otherwise a C kernel compiled for AVX2+FMA (x86 CPUs since about 2013), otherwise a generic C kernel "
    "(any CPU; fast on CPUs with fused multiply-add instructions such as all 64-bit ARM CPUs, slow on old x86 CPUs without FMA). "
    "The C kernels give bitwise identical results on all machines and instruction sets. "
    "For a single system with at most four planets the AVX2 and generic C kernels automatically use a four-lane variant "
    "(same results, faster; on AVX-512 machines `kernel=3` is then also faster than the assembly). "
    "The assembly kernel agrees with them to round-off. "
    "`kernel_used` reports which kernel ran. "
    "\n\n"
    "The algorithm is described in two papers [Javaheri et al. (2023)] and [Dagli & Rein (in prep)]. "
    "Note that in July 2026 significant changes have been made. "
    "See [Dagli & Rein (in prep)] for details on WHFast512 version 2. "
    "\n\n"
    "[Wisdom & Holman (1991)]: https://ui.adsabs.harvard.edu/abs/1991AJ....102.1528W/abstract\n"
    "[Javaheri et al. (2023)]: https://ui.adsabs.harvard.edu/abs/2023OJAp....6E..29J/abstract\n"
    "[Dagli & Rein (in prep)]: https://hanno-rein.de/publications.html\n"
    ,
    .step = reb_integrator_whfast512_step,
    .create = reb_integrator_whfast512_create,
    .free = reb_integrator_whfast512_free,
    .synchronize = reb_integrator_whfast512_synchronize,
    .field_descriptor_list = reb_integrator_whfast512_field_descriptor_list,
};

const struct reb_binarydata_field_descriptor reb_integrator_whfast512_field_descriptor_list[] = {
    { "If this flag is set to 1 (default is 0) then general relativistic corrections are included. "
        "The corrections are in the form of an additional potential term and reproduce the correct precession rate. "
        "The constants are hard coded for this effect and assume that the simulation is in units of G=1 and one length unit corresponds to one astronomical unit. ",
        REB_UINT,        "gr_potential",    offsetof(struct reb_integrator_whfast512_state, gr_potential), 0, 0, 0},
    { "If this flag is set to 17 (default is 0), then symplectic correctors are used.", 
        REB_UINT,        "corrector",       offsetof(struct reb_integrator_whfast512_state, corrector), 0, 0, 0},
    { "If this is set to a number other than 1 then timesteps are combined. "
        "By doing multiple timesteps in a row, WHFast512 can keep all simulation data in registers which significantly speeds up the calculation. "
        "This number should be as large as the output cadence allows. "
        "The default is 1e6. ",
        REB_UINT64,      "concatenate_steps", offsetof(struct reb_integrator_whfast512_state, concatenate_steps), 0, 0, 0},
    { "By default this value is set to 1, implying all 8 particles in the simulation correspond to one system. "
        "By setting N_systems to either 2 or 4, one can integrate multiple planetary systems with 2, 3, or 4 particles at the same time. "
        "See the example problems on how to setup the particles for this case. ",
        REB_UINT,        "N_systems",       offsetof(struct reb_integrator_whfast512_state, N_systems), 0, 0, 0},
    { "Kernel used for the integration: 0 (default) picks the fastest one available on this machine "
        "(AVX-512 assembly if the CPU supports AVX-512, otherwise the C kernel for AVX2+FMA, otherwise the generic C kernel). "
        "1: AVX-512 assembly, 2: C kernel compiled for AVX-512, 3: C kernel for AVX2+FMA, 4: generic C kernel (any CPU). "
        "The C kernels give bitwise identical results on all machines; the assembly kernel differs from them at round-off level. ",
        REB_UINT,        "kernel",          offsetof(struct reb_integrator_whfast512_state, kernel), 0, 0, 0},
    { "Kernel used by the last step (read only). Same values as `kernel`.",
        REB_UINT,        "kernel_used",     offsetof(struct reb_integrator_whfast512_state, kernel_used), 0, 0, 0},
    { "", REB_POINTER_ALIGNED, "data",        offsetof(struct reb_integrator_whfast512_state, data), SIZE_MAX, sizeof(struct simd_data), 0},
    { "", REB_DOUBLE,      "last_synchronization", offsetof(struct reb_integrator_whfast512_state, last_synchronization), 0, 0, 0},
    { 0 }, // Null terminated list
};

void* reb_integrator_whfast512_create(){
    struct reb_integrator_whfast512_state* whfast512 = calloc(sizeof(struct reb_integrator_whfast512_state),1);
    whfast512->N_systems = 1;
    whfast512->gr_potential = 0;
    whfast512->concatenate_steps = 1e6;
    return whfast512;
}

void reb_integrator_whfast512_free(void* state){
    struct reb_integrator_whfast512_state* whfast512 = state;
    reb_aligned_free(whfast512->data);
    free(whfast512);
}

// Helper macro to print out offsets in structure for assembly code
#define SIMD_DATA_MEMBERS X(M) X(dt) X(gr_prefac) X(m) X(x) X(y) X(z) X(vx) X(vy) X(vz) \
X(mat8_inertial_to_jacobi) \
X(mat8_jacobi_to_heliocentric) \
X(M0) X(mask) \
X(exit_max_distance) X(exit_min_distance_r) \
X(mat8_jacobi_to_inertial)\
X(counter) 

#ifdef REB_WHFAST512_HAVE_ASM
// External function definitions. Implemented in integrator_whfast512.s.
REB_SYSV_ABI extern enum REB_STATUS reb_whfast512_kepler_step(struct simd_data* data);
REB_SYSV_ABI extern void reb_whfast512_set1_pd(void* address, double value);
REB_SYSV_ABI extern void reb_whfast512_movu_pd(void* destination, void* source);
// _n2 = two systems of up to 4 planets, _n4 = four systems of 2 planets.
REB_SYSV_ABI extern enum REB_STATUS reb_whfast512_full_steps_gr0_n1_encounter0_escape0(struct simd_data* data, uint64_t* N_steps, int skip_first_kepler_step, volatile sig_atomic_t* sigint);
REB_SYSV_ABI extern enum REB_STATUS reb_whfast512_full_steps_gr0_n1_encounter1_escape0(struct simd_data* data, uint64_t* N_steps, int skip_first_kepler_step, volatile sig_atomic_t* sigint);
REB_SYSV_ABI extern enum REB_STATUS reb_whfast512_full_steps_gr0_n1_encounter0_escape1(struct simd_data* data, uint64_t* N_steps, int skip_first_kepler_step, volatile sig_atomic_t* sigint);
REB_SYSV_ABI extern enum REB_STATUS reb_whfast512_full_steps_gr0_n1_encounter1_escape1(struct simd_data* data, uint64_t* N_steps, int skip_first_kepler_step, volatile sig_atomic_t* sigint);
REB_SYSV_ABI extern enum REB_STATUS reb_whfast512_full_steps_gr0_n2_encounter0_escape0(struct simd_data* data, uint64_t* N_steps, int skip_first_kepler_step, volatile sig_atomic_t* sigint);
REB_SYSV_ABI extern enum REB_STATUS reb_whfast512_full_steps_gr0_n2_encounter1_escape0(struct simd_data* data, uint64_t* N_steps, int skip_first_kepler_step, volatile sig_atomic_t* sigint);
REB_SYSV_ABI extern enum REB_STATUS reb_whfast512_full_steps_gr0_n2_encounter0_escape1(struct simd_data* data, uint64_t* N_steps, int skip_first_kepler_step, volatile sig_atomic_t* sigint);
REB_SYSV_ABI extern enum REB_STATUS reb_whfast512_full_steps_gr0_n2_encounter1_escape1(struct simd_data* data, uint64_t* N_steps, int skip_first_kepler_step, volatile sig_atomic_t* sigint);
REB_SYSV_ABI extern enum REB_STATUS reb_whfast512_full_steps_gr0_n4_encounter0_escape0(struct simd_data* data, uint64_t* N_steps, int skip_first_kepler_step, volatile sig_atomic_t* sigint);
REB_SYSV_ABI extern enum REB_STATUS reb_whfast512_full_steps_gr0_n4_encounter1_escape0(struct simd_data* data, uint64_t* N_steps, int skip_first_kepler_step, volatile sig_atomic_t* sigint);
REB_SYSV_ABI extern enum REB_STATUS reb_whfast512_full_steps_gr0_n4_encounter0_escape1(struct simd_data* data, uint64_t* N_steps, int skip_first_kepler_step, volatile sig_atomic_t* sigint);
REB_SYSV_ABI extern enum REB_STATUS reb_whfast512_full_steps_gr0_n4_encounter1_escape1(struct simd_data* data, uint64_t* N_steps, int skip_first_kepler_step, volatile sig_atomic_t* sigint);
REB_SYSV_ABI extern enum REB_STATUS reb_whfast512_full_steps_gr1_n1_encounter0_escape0(struct simd_data* data, uint64_t* N_steps, int skip_first_kepler_step, volatile sig_atomic_t* sigint);
REB_SYSV_ABI extern enum REB_STATUS reb_whfast512_full_steps_gr1_n1_encounter1_escape0(struct simd_data* data, uint64_t* N_steps, int skip_first_kepler_step, volatile sig_atomic_t* sigint);
REB_SYSV_ABI extern enum REB_STATUS reb_whfast512_full_steps_gr1_n1_encounter0_escape1(struct simd_data* data, uint64_t* N_steps, int skip_first_kepler_step, volatile sig_atomic_t* sigint);
REB_SYSV_ABI extern enum REB_STATUS reb_whfast512_full_steps_gr1_n1_encounter1_escape1(struct simd_data* data, uint64_t* N_steps, int skip_first_kepler_step, volatile sig_atomic_t* sigint);
REB_SYSV_ABI extern enum REB_STATUS reb_whfast512_full_steps_gr1_n2_encounter0_escape0(struct simd_data* data, uint64_t* N_steps, int skip_first_kepler_step, volatile sig_atomic_t* sigint);
REB_SYSV_ABI extern enum REB_STATUS reb_whfast512_full_steps_gr1_n2_encounter1_escape0(struct simd_data* data, uint64_t* N_steps, int skip_first_kepler_step, volatile sig_atomic_t* sigint);
REB_SYSV_ABI extern enum REB_STATUS reb_whfast512_full_steps_gr1_n2_encounter0_escape1(struct simd_data* data, uint64_t* N_steps, int skip_first_kepler_step, volatile sig_atomic_t* sigint);
REB_SYSV_ABI extern enum REB_STATUS reb_whfast512_full_steps_gr1_n2_encounter1_escape1(struct simd_data* data, uint64_t* N_steps, int skip_first_kepler_step, volatile sig_atomic_t* sigint);
REB_SYSV_ABI extern enum REB_STATUS reb_whfast512_full_steps_gr1_n4_encounter0_escape0(struct simd_data* data, uint64_t* N_steps, int skip_first_kepler_step, volatile sig_atomic_t* sigint);
REB_SYSV_ABI extern enum REB_STATUS reb_whfast512_full_steps_gr1_n4_encounter1_escape0(struct simd_data* data, uint64_t* N_steps, int skip_first_kepler_step, volatile sig_atomic_t* sigint);
REB_SYSV_ABI extern enum REB_STATUS reb_whfast512_full_steps_gr1_n4_encounter0_escape1(struct simd_data* data, uint64_t* N_steps, int skip_first_kepler_step, volatile sig_atomic_t* sigint);
REB_SYSV_ABI extern enum REB_STATUS reb_whfast512_full_steps_gr1_n4_encounter1_escape1(struct simd_data* data, uint64_t* N_steps, int skip_first_kepler_step, volatile sig_atomic_t* sigint);
REB_SYSV_ABI extern void reb_whfast512_corrector_step_gr0_n1(struct simd_data* data, double inv);
REB_SYSV_ABI extern void reb_whfast512_corrector_step_gr0_n2(struct simd_data* data, double inv);
REB_SYSV_ABI extern void reb_whfast512_corrector_step_gr0_n4(struct simd_data* data, double inv);
REB_SYSV_ABI extern void reb_whfast512_corrector_step_gr1_n1(struct simd_data* data, double inv);
REB_SYSV_ABI extern void reb_whfast512_corrector_step_gr1_n2(struct simd_data* data, double inv);
REB_SYSV_ABI extern void reb_whfast512_corrector_step_gr1_n4(struct simd_data* data, double inv);
#endif // REB_WHFAST512_HAVE_ASM

// Small helpers to fill the 64-byte vectors in struct simd_data (plain C, no vector instructions).
static void w512_set1_pd(void* address, double value){
    double tmp[8];
    for (int i=0; i<8; i++) tmp[i] = value;
    memcpy(address, tmp, sizeof(tmp));
}
static void w512_movu_pd(void* destination, const void* source){
    memcpy(destination, source, 8*sizeof(double));
}

// Returns 1 if this CPU (and OS) supports AVX2 and FMA.
static int whfast512_avx2_fma_available(void){
#if (defined(__i386__) || defined(__x86_64__)) && (defined(__GNUC__) || defined(__clang__)) && !defined(_WIN32)
    __builtin_cpu_init();
    return __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
#elif defined(_WIN32) && (defined(_M_X64) || defined(_M_IX86) || defined(__x86_64__))
    return reb_avx2_fma_available();
#else
    return 0;
#endif
}

static const char* whfast512_kernel_name(unsigned int k){
    switch (k){
        case REB_WHFAST512_KERNEL_ASM: return "AVX-512 assembly";
        case REB_WHFAST512_KERNEL_C_AVX512: return "C, AVX-512";
        case REB_WHFAST512_KERNEL_C_AVX2: return "C, AVX2+FMA";
        case REB_WHFAST512_KERNEL_C_GENERIC: return "C, generic";
        default: return "unknown";
    }
}

// Returns 1 if kernel k can run on this machine.
static int whfast512_kernel_available(unsigned int k){
    switch (k){
        case REB_WHFAST512_KERNEL_ASM:
#ifdef REB_WHFAST512_HAVE_ASM
            return reb_avx512_available();
#else
            return 0;
#endif
        case REB_WHFAST512_KERNEL_C_AVX512: return reb_whfast512c_avx512_compiled() && reb_avx512_available() && whfast512_avx2_fma_available();
        case REB_WHFAST512_KERNEL_C_AVX2: return reb_whfast512c_avx2_compiled() && whfast512_avx2_fma_available();
        case REB_WHFAST512_KERNEL_C_GENERIC: return 1;
        default: return 0;
    }
}

// Picks the kernel: the requested one, or the fastest available one for REB_WHFAST512_KERNEL_AUTO.
// Returns 0 (and prints an error) if the requested kernel is not available.
static unsigned int whfast512_select_kernel(struct reb_simulation* r, unsigned int requested){
    if (requested == REB_WHFAST512_KERNEL_AUTO){
        static const unsigned int order[] = {REB_WHFAST512_KERNEL_ASM, REB_WHFAST512_KERNEL_C_AVX2, REB_WHFAST512_KERNEL_C_GENERIC};
        for (unsigned int i=0; i<sizeof(order)/sizeof(order[0]); i++){
            if (whfast512_kernel_available(order[i])) return order[i];
        }
        return REB_WHFAST512_KERNEL_C_GENERIC;
    }
    if (requested > REB_WHFAST512_KERNEL_C_GENERIC){
        reb_simulation_error(r, "WHFast512: unknown kernel. Use 0 (auto), 1 (assembly), 2 (C AVX-512), 3 (C AVX2), or 4 (C generic).");
        return 0;
    }
    if (!whfast512_kernel_available(requested)){
        char msg[256];
        snprintf(msg, sizeof(msg), "WHFast512: the kernel '%s' is not available on this machine or in this build.", whfast512_kernel_name(requested));
        reb_simulation_error(r, msg);
        return 0;
    }
    return requested;
}


// Hepler function to load particle data into avx512 registers
static void load_into_m512d(void* destination, struct reb_simulation* r, size_t offset, const double* transformation, int N_systems){
    struct reb_particle* particles = r->particles;
    const unsigned int p_per_system = 8/N_systems;
    const unsigned int N_per_system = r->N/N_systems;
    double tmp[8] = {0}; 
    for (int s=0; s<N_systems; s++){
        for (unsigned int i=1; i<N_per_system; i++){
            tmp[s*p_per_system+i-1] = *(double*)((char*)(&particles[s*N_per_system+i])+offset);
        }
    }
    if (transformation != NULL){
        double tmp2[8] = {0}; 
        for (int i=0; i<8; i++) {
            for (int j=0; j<8; j++) {
                tmp2[i] += tmp[j]*transformation[8*j+i];
            }
        }
        w512_movu_pd(destination, tmp2);
    }else{
        w512_movu_pd(destination, tmp);
    }
}

// Hepler function to load particle data from avx512 registers
static void load_from_m512d(struct reb_simulation* r, size_t offset, const double* transformation, int N_systems, void* vector){
    struct reb_particle* particles = r->particles;
    const unsigned int p_per_system = 8/N_systems;
    const unsigned int N_per_system = r->N/N_systems;
    double tmp[8];
    w512_movu_pd(tmp, vector);
    double tmp2[8] = {0}; 
    for (int i=0; i<8; i++) {
        for (int j=0; j<8; j++) {
            tmp2[i] += tmp[j]*transformation[8*j+i];
        }
    }
    for (int s=0; s<N_systems; s++){
        for (unsigned int i=1; i<N_per_system; i++){
            *(double*)((char*)(&particles[s*N_per_system+i])+offset) = tmp2[s*p_per_system+i-1];
        }
    }
}

// Convert jacobi coordinates to inertial coordinates
// Also performs com step (assume original particles are unmodified)
// Note: Speed is not a concern here 
static void jacobi_to_inertial_posvel_and_com(struct reb_simulation* r, struct simd_data* data, double dt_com, unsigned int N_systems){
    const unsigned int N_per_system = r->N/N_systems;
    struct reb_particle com[4];
    for (unsigned s=0;s<N_systems;s++){
        com[s] = reb_simulation_com_range(r,s*N_per_system, (s+1)*N_per_system); // original com
    }
    struct reb_particle* particles = r->particles;
    load_from_m512d(r, offsetof(struct reb_particle, x), data->mat8_jacobi_to_inertial, N_systems, &data->x);
    load_from_m512d(r, offsetof(struct reb_particle, y), data->mat8_jacobi_to_inertial, N_systems, &data->y);
    load_from_m512d(r, offsetof(struct reb_particle, z), data->mat8_jacobi_to_inertial, N_systems, &data->z);
    load_from_m512d(r, offsetof(struct reb_particle, vx), data->mat8_jacobi_to_inertial, N_systems, &data->vx);
    load_from_m512d(r, offsetof(struct reb_particle, vy), data->mat8_jacobi_to_inertial, N_systems, &data->vy);
    load_from_m512d(r, offsetof(struct reb_particle, vz), data->mat8_jacobi_to_inertial, N_systems, &data->vz);
    for (unsigned s=0;s<N_systems;s++){
        particles[s*N_per_system+0].x  = 0.0;
        particles[s*N_per_system+0].y  = 0.0;
        particles[s*N_per_system+0].z  = 0.0;
        particles[s*N_per_system+0].vx = 0.0;
        particles[s*N_per_system+0].vy = 0.0;
        particles[s*N_per_system+0].vz = 0.0;
        for (unsigned int i=1;i<N_per_system;i++){
            particles[s*N_per_system+0].x  -= particles[s*N_per_system+i].m/particles[s*N_per_system+0].m*particles[s*N_per_system+i].x;
            particles[s*N_per_system+0].y  -= particles[s*N_per_system+i].m/particles[s*N_per_system+0].m*particles[s*N_per_system+i].y;
            particles[s*N_per_system+0].z  -= particles[s*N_per_system+i].m/particles[s*N_per_system+0].m*particles[s*N_per_system+i].z;
            particles[s*N_per_system+0].vx -= particles[s*N_per_system+i].m/particles[s*N_per_system+0].m*particles[s*N_per_system+i].vx;
            particles[s*N_per_system+0].vy -= particles[s*N_per_system+i].m/particles[s*N_per_system+0].m*particles[s*N_per_system+i].vy;
            particles[s*N_per_system+0].vz -= particles[s*N_per_system+i].m/particles[s*N_per_system+0].m*particles[s*N_per_system+i].vz;
        }
        for (unsigned int i=0;i<N_per_system;i++){
            particles[s*N_per_system+i].x  += com[s].x + com[s].vx*dt_com;
            particles[s*N_per_system+i].y  += com[s].y + com[s].vy*dt_com;
            particles[s*N_per_system+i].z  += com[s].z + com[s].vz*dt_com;
            particles[s*N_per_system+i].vx += com[s].vx;
            particles[s*N_per_system+i].vy += com[s].vy;
            particles[s*N_per_system+i].vz += com[s].vz;
        }
    }
}


static enum REB_STATUS whfast512_full_steps(struct reb_simulation* r, struct reb_integrator_whfast512_state* whfast512, uint64_t* N_steps, int skip){
    struct simd_data* data = whfast512->data;
    const int encounter = r->exit_min_distance ? 1 : 0;
    const int escape = r->exit_max_distance ? 1 : 0;
    switch (whfast512->kernel_used){
        case REB_WHFAST512_KERNEL_C_AVX512:
            return reb_whfast512c_full_steps_avx512(data, N_steps, skip, &reb_sigint, whfast512->gr_potential, whfast512->N_systems, encounter, escape);
        case REB_WHFAST512_KERNEL_C_AVX2:
            if (whfast512->four_lanes){
                return reb_whfast512c_full_steps_avx2x4(data, N_steps, skip, &reb_sigint, whfast512->gr_potential, 2, encounter, escape);
            }
            return reb_whfast512c_full_steps_avx2(data, N_steps, skip, &reb_sigint, whfast512->gr_potential, whfast512->N_systems, encounter, escape);
        case REB_WHFAST512_KERNEL_C_GENERIC:
            if (whfast512->four_lanes){
                return reb_whfast512c_full_steps_genericx4(data, N_steps, skip, &reb_sigint, whfast512->gr_potential, 2, encounter, escape);
            }
            return reb_whfast512c_full_steps_generic(data, N_steps, skip, &reb_sigint, whfast512->gr_potential, whfast512->N_systems, encounter, escape);
        default:
            break;
    }
#ifdef REB_WHFAST512_HAVE_ASM
    // Assembly: all the 24 possible options are separate functions.
#define MACRO_5(encounter,escape) \
    if (whfast512->gr_potential){\
        switch (whfast512->N_systems){\
            case 2: return reb_whfast512_full_steps_gr1_n2_encounter##encounter##_escape##escape(data, N_steps, skip, &reb_sigint); break;\
            case 4: return reb_whfast512_full_steps_gr1_n4_encounter##encounter##_escape##escape(data, N_steps, skip, &reb_sigint); break;\
            default: return reb_whfast512_full_steps_gr1_n1_encounter##encounter##_escape##escape(data, N_steps, skip, &reb_sigint); break;\
        }\
    }else{\
        switch (whfast512->N_systems){\
            case 2: return reb_whfast512_full_steps_gr0_n2_encounter##encounter##_escape##escape(data, N_steps, skip, &reb_sigint); break;\
            case 4: return reb_whfast512_full_steps_gr0_n4_encounter##encounter##_escape##escape(data, N_steps, skip, &reb_sigint); break;\
            default: return reb_whfast512_full_steps_gr0_n1_encounter##encounter##_escape##escape(data, N_steps, skip, &reb_sigint); break;\
        }\
    }
    if (encounter){
        if (escape){
            MACRO_5(1,1)
        }else{
            MACRO_5(1,0)
        }
    }else{
        if (escape){
            MACRO_5(0,1)
        }else{
            MACRO_5(0,0)
        }
    }
#undef MACRO_5
#endif // REB_WHFAST512_HAVE_ASM
    reb_simulation_error(r, "WHFast512: no kernel selected.");
    return REB_STATUS_GENERIC_ERROR;
}

static void whfast512_corrector_step(struct reb_integrator_whfast512_state* whfast512, double inv){
    struct simd_data* data = whfast512->data;
    switch (whfast512->kernel_used){
        case REB_WHFAST512_KERNEL_C_AVX512:
            reb_whfast512c_corrector_step_avx512(data, inv, whfast512->gr_potential, whfast512->N_systems); return;
        case REB_WHFAST512_KERNEL_C_AVX2:
            if (whfast512->four_lanes){
                reb_whfast512c_corrector_step_avx2x4(data, inv, whfast512->gr_potential, 2); return;
            }
            reb_whfast512c_corrector_step_avx2(data, inv, whfast512->gr_potential, whfast512->N_systems); return;
        case REB_WHFAST512_KERNEL_C_GENERIC:
            if (whfast512->four_lanes){
                reb_whfast512c_corrector_step_genericx4(data, inv, whfast512->gr_potential, 2); return;
            }
            reb_whfast512c_corrector_step_generic(data, inv, whfast512->gr_potential, whfast512->N_systems); return;
        default:
            break;
    }
#ifdef REB_WHFAST512_HAVE_ASM

    if (whfast512->gr_potential){
        switch (whfast512->N_systems){
            case 2: reb_whfast512_corrector_step_gr1_n2(data, inv); break;
            case 4: reb_whfast512_corrector_step_gr1_n4(data, inv); break;
            default: reb_whfast512_corrector_step_gr1_n1(data, inv); break;
        }
    }else{
        switch (whfast512->N_systems){
            case 2: reb_whfast512_corrector_step_gr0_n2(data, inv); break;
            case 4: reb_whfast512_corrector_step_gr0_n4(data, inv); break;
            default: reb_whfast512_corrector_step_gr0_n1(data, inv); break;
        }
    }
#endif // REB_WHFAST512_HAVE_ASM
}

// One Kepler step (drift) with the timestep stored in data->dt.
static void whfast512_kepler_step(struct reb_integrator_whfast512_state* whfast512){
    struct simd_data* data = whfast512->data;
    switch (whfast512->kernel_used){
        case REB_WHFAST512_KERNEL_C_AVX512: reb_whfast512c_kepler_step_avx512(data); return;
        case REB_WHFAST512_KERNEL_C_AVX2:
            if (whfast512->four_lanes){
                reb_whfast512c_kepler_step_avx2x4(data); return;
            }
            reb_whfast512c_kepler_step_avx2(data); return;
        case REB_WHFAST512_KERNEL_C_GENERIC:
            if (whfast512->four_lanes){
                reb_whfast512c_kepler_step_genericx4(data); return;
            }
            reb_whfast512c_kepler_step_generic(data); return;
        default: break;
    }
#ifdef REB_WHFAST512_HAVE_ASM
    reb_whfast512_kepler_step(data);
#endif // REB_WHFAST512_HAVE_ASM
}

static void inertial_to_jacobi_posvel(struct reb_simulation* r, struct simd_data* data, unsigned int N_systems){
    const unsigned int N_per_system = r->N/N_systems;
    // Transformations assume system is in COM frame.
    struct reb_particle com[4];
    struct reb_particle* p_tmp = malloc(sizeof(struct reb_particle)*r->N);
    memcpy(p_tmp, r->particles, sizeof(struct reb_particle)*r->N);
    for (unsigned int s=0;s<N_systems;s++){
        com[s] = reb_simulation_com_range(r,s*N_per_system, (s+1)*N_per_system); // original com
    }
    for (unsigned int s=0;s<N_systems;s++){
        for (unsigned int i=0;i<N_per_system;i++){
            r->particles[s*N_per_system+i].x -= com[s].x;
            r->particles[s*N_per_system+i].y -= com[s].y;
            r->particles[s*N_per_system+i].z -= com[s].z;
            r->particles[s*N_per_system+i].vx -= com[s].vx;
            r->particles[s*N_per_system+i].vy -= com[s].vy;
            r->particles[s*N_per_system+i].vz -= com[s].vz;
        }
    }
    reb_simulation_move_to_com(r);
    // Same layout as for democratic heliocentric
    load_into_m512d(&data->x,  r, offsetof(struct reb_particle,x),data->mat8_inertial_to_jacobi, N_systems);
    load_into_m512d(&data->y,  r, offsetof(struct reb_particle,y),data->mat8_inertial_to_jacobi, N_systems);
    load_into_m512d(&data->z,  r, offsetof(struct reb_particle,z),data->mat8_inertial_to_jacobi, N_systems);
    load_into_m512d(&data->vx, r, offsetof(struct reb_particle,vx),data->mat8_inertial_to_jacobi, N_systems);
    load_into_m512d(&data->vy, r, offsetof(struct reb_particle,vy),data->mat8_inertial_to_jacobi, N_systems);
    load_into_m512d(&data->vz, r, offsetof(struct reb_particle,vz),data->mat8_inertial_to_jacobi, N_systems);
    load_into_m512d(&data->m,  r, offsetof(struct reb_particle,m),NULL, N_systems);
    // Undo COM transformation. COM will be applied in jacobi_to_inertial_posvel_and_com().
    memcpy(r->particles, p_tmp, sizeof(struct reb_particle)*r->N);
    free(p_tmp);
}


// Precalculate various constants and put them in 512 bit vectors.
static void recalculate_constants(struct reb_simulation* r, unsigned int N_systems){
    struct reb_integrator_whfast512_state* whfast512 = r->integrator.state;
    reb_aligned_free(whfast512->data); // free in case previously allocated
    whfast512->data = reb_aligned_alloc(64,sizeof(struct simd_data));
    memset(whfast512->data, 0, sizeof(struct simd_data));
    if (!whfast512->data){
        reb_simulation_error(r, "WHFast512 was not able to allocate memory.");
        return;
    }
    struct simd_data* data = whfast512->data;
    struct reb_particle* particles = r->particles;
    const unsigned int p_per_system = 8/N_systems;
    const unsigned int N_per_system = r->N/N_systems;
    double mat8_inertial_to_heliocentric[64];
    double M[8] = {0.0};
    double M0[8] = {0.0};
    switch (N_systems){
        case 1:
            data->mask = (1 << (r->N -1)) - 1;
            break;
        case 2:
            if (N_per_system==5) data->mask = 0xFF;
            if (N_per_system==4) data->mask = 0x77;
            if (N_per_system==3) data->mask = 0x33;
            if (N_per_system==2) data->mask = 0x11;
            break;
        case 4:
            if (N_per_system==3) data->mask = 0xFF;
            if (N_per_system==2) data->mask = 0x55;
            break;
        default:
            reb_simulation_error(r,"Invalid value for N_systems.");
    }
    // Zeroing.
    for (unsigned int i=1; i<9; i++){
        for (unsigned int j=1; j<9; j++){
            data->mat8_inertial_to_jacobi[(i-1)+8*(j-1)] = 0.0;
            data->mat8_jacobi_to_inertial[(i-1)+8*(j-1)] = 0.0;
            data->mat8_jacobi_to_heliocentric[(i-1)+8*(j-1)] = 0.0;
            mat8_inertial_to_heliocentric[(i-1)+8*(j-1)] = 0.0;
        }
    }

    // Filling vector
    for (unsigned int s=0; s<N_systems; s++){
        for (unsigned int i=0;i<N_per_system-1;i++){
            for (unsigned int j=0;j<i+2;j++){
                M[s*p_per_system+i] += r->particles[s*N_per_system+j].m;
            }
        }
    }

    // Fill matricies
    for (unsigned int s=0; s<N_systems; s++){
        double ms = particles[s*N_per_system+0].m;
        for (unsigned int i=1; i<N_per_system; i++){
            for (unsigned int j=i; j<N_per_system; j++){
                data->mat8_inertial_to_jacobi[(s*p_per_system+i-1)+8*(s*p_per_system+j-1)] += particles[s*N_per_system+j].m/ms;
            }
            for (unsigned int j=1; j<N_per_system; j++){
                mat8_inertial_to_heliocentric[(s*p_per_system+i-1)+8*(s*p_per_system+j-1)] += particles[s*N_per_system+j].m/particles[s*N_per_system+0].m;
            }
            mat8_inertial_to_heliocentric[(s*p_per_system+i-1)+8*(s*p_per_system+i-1)] += 1.0;
            data->mat8_inertial_to_jacobi[(s*p_per_system+i-1)+8*(s*p_per_system+i-1)] += 1.0;
            data->mat8_jacobi_to_inertial[(s*p_per_system+i-1)+8*(s*p_per_system+i-1)] += ms/(ms + particles[s*N_per_system+i].m);
            ms += particles[s*N_per_system+i].m;
            if (i<N_per_system-1){ 
                for (unsigned int ii=i; ii>0; ii--){
                    int jj = i+1;
                    data->mat8_jacobi_to_inertial[(s*p_per_system+ii-1)+8*(s*p_per_system+jj-1)] -= particles[s*N_per_system+jj].m/(ms+particles[s*N_per_system+jj].m);
                }
            }
            M0[(s*p_per_system+i-1)] = particles[s*N_per_system+0].m;
        }
    }

    // Might be numerically more stable to calculate this manually rather than do a matrix multiplication.
    for (unsigned int i=1; i<9; i++){
        for (unsigned int j=1; j<9; j++){
            for (unsigned int k=1; k<9; k++){
                data->mat8_jacobi_to_heliocentric[(i-1)+8*(j-1)] += mat8_inertial_to_heliocentric[(i-1)+8*(k-1)] * data->mat8_jacobi_to_inertial[(k-1)+8*(j-1)];
            }
        }
    }

    w512_movu_pd(&data->M, &M);
    w512_movu_pd(&data->M0, &M0); //  = particles[0].m 

    // GR prefactors. Note: assumes units of AU, year/2pi.
    double c = 10065.32;
    double _gr_prefac[8];
    for(unsigned int i=0;i<8;i++){
        _gr_prefac[i] = 0; // for when N<8
    }
    for (unsigned int s=0; s<N_systems; s++){
        double m0 = r->particles[s*N_per_system].m;
        for (unsigned int p=1; p<N_per_system; p++){
            _gr_prefac[s*p_per_system+(p-1)] = -6.*m0*m0/(c*c);
        }
    }
    w512_movu_pd(&data->gr_prefac, &_gr_prefac);
    w512_set1_pd(&data->dt,r->dt); 
    w512_set1_pd(&data->exit_max_distance, r->exit_max_distance);
    w512_set1_pd(&data->exit_min_distance, r->exit_min_distance);
    w512_set1_pd(&data->exit_min_distance_r, 1.0/r->exit_min_distance);
#define X(name) printf(".set P512_" #name ", %zu\n", offsetof(struct simd_data, name));
    //    SIMD_DATA_MEMBERS
#undef X

}

static int reb_integrator_whfast512_verify_setup(struct reb_simulation* const r){
    struct reb_integrator_whfast512_state* whfast512 = r->integrator.state;
    // Check if all assumptions are satisfied.
    // Note: These are not checked every timestep. 
    // So it is possible for the user to screw things up.
    whfast512->kernel_used = whfast512_select_kernel(r, whfast512->kernel);
    if (!whfast512->kernel_used){
        return 1;
    }
    // A single system with at most four planets fits into lanes 0-3: the AVX2 kernel then uses
    // its four-lane variant (same results, about half the work).
    whfast512->four_lanes = (whfast512->N_systems == 1 && r->N <= 5) ? 1 : 0;
    if (r->N_var!=0){
        reb_simulation_error(r, "WHFast512 does not support variational particles.");
        return 1;
    }
    if (r->exact_finish_time!=0){
        reb_simulation_error(r, "WHFast512 requires exact_finish_time=0.");
        return 1;
    }
    if (r->N>9 && whfast512->N_systems == 1) {
        reb_simulation_error(r, "WHFast512 supports a maximum of 9 particles when N_systems is set to 1.");
        return 1;
    }
    if (r->N>10 && whfast512->N_systems == 2) {
        reb_simulation_error(r, "WHFast512 supports a maximum of 10 particles when N_systems is set to 2.");
        return 1;
    }
    if (r->N>12 && whfast512->N_systems == 4) {
        reb_simulation_error(r, "WHFast512 supports a maximum of 12 particles when N_systems is set to 4.");
        return 1;
    }
    if (whfast512->N_systems != 1 && whfast512->N_systems !=2 && whfast512->N_systems != 4){
        reb_simulation_error(r, "WHFast512 supports 1, 2, or 4 systems only.");
        return 1;
    }
    if (r->N % whfast512->N_systems != 0){
        reb_simulation_error(r, "Number of particles must be a multiple of whfast512.N_systems.");
        return 1;
    }
    if (r->G!=1.0){
        reb_simulation_error(r, "WHFast512 requires units in which G=1. Please rescale your system.");
        return 1;
    }
    if (r->N_active!=SIZE_MAX && r->N_active!=r->N){
        reb_simulation_error(r, "WHFast512 does not support test particles.");
        return 1;
    }
    r->gravity = REB_GRAVITY_NONE; // WHFast512 uses its own gravity routine.
    return 0; // success
}

// If the simulation is synchronized: set up the data, apply the corrector and do the first half
// drift. Returns 1 if the first Kepler step of the main loop has to be skipped.
static int whfast512_prepare(struct reb_simulation* const r, struct reb_integrator_whfast512_state* whfast512){
    if (!r->is_synchronized){
        return 0;
    }
    const double dt = r->dt;
    recalculate_constants(r, whfast512->N_systems);
    struct simd_data* data = whfast512->data;
    // Use WHFast to apply the correctors.
    inertial_to_jacobi_posvel(r, data, whfast512->N_systems);
    whfast512->last_synchronization = r->t;
    if (whfast512->corrector){
        whfast512_corrector_step(whfast512, 1.0);
    }
    // First half DRIFT step.
    w512_set1_pd(&data->dt, dt/2.0); 
    whfast512_kepler_step(whfast512);
    w512_set1_pd(&data->dt, dt); // Reset
    return 1;
}

// Optimized main loops allowing for concatenate_steps
void reb_integrator_whfast512_step(struct reb_simulation* const r, void* state){
    struct reb_integrator_whfast512_state* whfast512 = state;
    const double dt = r->dt;
    uint64_t N_steps = whfast512->concatenate_steps;

    if (reb_integrator_whfast512_verify_setup(r)){
        r->status = REB_STATUS_GENERIC_ERROR;
        return; // Error occured
    }

    const int skip_first_kepler_step = whfast512_prepare(r, whfast512);
    r->status = whfast512_full_steps(r, whfast512, &N_steps, skip_first_kepler_step);

    r->is_synchronized = 0;
    r->t += dt*N_steps;     // Note: N_steps might have been changed by whfast512_full_steps.
    r->dt_last_done = dt;
}

// Synchronization routine. Called every time an output is needed.
void reb_integrator_whfast512_synchronize(struct reb_simulation* const r, void* state){
    struct reb_integrator_whfast512_state* const whfast512 = state;
    if (!r->is_synchronized){
        struct simd_data * data = whfast512->data;
        if (!data){
            reb_simulation_error(r, "WHFast512 is unable to synchronize. data is NULL.");
            r->status = REB_STATUS_GENERIC_ERROR;
            return;
        }
        // Select again: the state may come from a Simulationarchive written on another machine.
        whfast512->kernel_used = whfast512_select_kernel(r, whfast512->kernel);
        if (!whfast512->kernel_used){
            r->status = REB_STATUS_GENERIC_ERROR;
            return;
        }
        whfast512->four_lanes = (whfast512->N_systems == 1 && r->N <= 5) ? 1 : 0;
        w512_set1_pd(&data->dt, r->dt/2.0); 
        whfast512_kepler_step(whfast512);
        w512_set1_pd(&data->dt, r->dt); // Reset
        if (whfast512->corrector){
            whfast512_corrector_step(whfast512, -1.0);
        }
        double dt_com = r->t - whfast512->last_synchronization;
        jacobi_to_inertial_posvel_and_com(r, data, dt_com, whfast512->N_systems);
        whfast512->last_synchronization = r->t;
        r->is_synchronized = 1;
        reb_aligned_free(whfast512->data);
        whfast512->data = NULL;
    }
}


// Two simulations per call. See rebound.h (reb_integrator_whfast512_steps_pair).
typedef int (*whfast512_pair_kernel)(struct simd_data* da, struct simd_data* db, uint64_t N_steps, uint64_t* N_done, int skip_a, int skip_b, volatile sig_atomic_t* sigint, int* status_a, int* status_b, int gr, int nsys, int encounter, int escape);

static void whfast512_steps_single(struct reb_simulation* const r, uint64_t N_steps){
    struct reb_integrator_whfast512_state* whfast512 = r->integrator.state;
    const uint64_t saved = whfast512->concatenate_steps;
    whfast512->concatenate_steps = N_steps;
    reb_simulation_steps(r, 1);
    whfast512->concatenate_steps = saved;
}

void reb_integrator_whfast512_steps_pair(struct reb_simulation* const r1, struct reb_simulation* const r2, uint64_t N_steps){
    struct reb_simulation* rs[2] = {r1, r2};
    for (int i=0; i<2; i++){
        if (strcmp(rs[i]->integrator.name, "whfast512") != 0){
            reb_simulation_error(rs[i], "reb_integrator_whfast512_steps_pair() requires the WHFast512 integrator.");
            rs[i]->status = REB_STATUS_GENERIC_ERROR;
            return;
        }
        if (reb_integrator_whfast512_verify_setup(rs[i])){
            rs[i]->status = REB_STATUS_GENERIC_ERROR;
            return;
        }
    }
    struct reb_integrator_whfast512_state* w1 = r1->integrator.state;
    struct reb_integrator_whfast512_state* w2 = r2->integrator.state;
    const int enc1 = r1->exit_min_distance ? 1 : 0, esc1 = r1->exit_max_distance ? 1 : 0;
    const int enc2 = r2->exit_min_distance ? 1 : 0, esc2 = r2->exit_max_distance ? 1 : 0;
    whfast512_pair_kernel kernel = NULL;
    if (r1 != r2 && w1->kernel_used == w2->kernel_used && w1->four_lanes == w2->four_lanes
            && w1->gr_potential == w2->gr_potential && w1->N_systems == w2->N_systems
            && enc1 == enc2 && esc1 == esc2){
        switch (w1->kernel_used){
            case REB_WHFAST512_KERNEL_C_AVX512: kernel = reb_whfast512c_full_steps_pair_avx512; break;
            // Eight-lane AVX2: one simulation already fills the 16 registers, pairing gains nothing
            // (0.97-1.00x measured), so the two are stepped separately.
            case REB_WHFAST512_KERNEL_C_AVX2: kernel = w1->four_lanes ? reb_whfast512c_full_steps_pair_avx2x4 : NULL; break;
            case REB_WHFAST512_KERNEL_C_GENERIC: kernel = w1->four_lanes ? reb_whfast512c_full_steps_pair_genericx4 : reb_whfast512c_full_steps_pair_generic; break;
            default: break; // the assembly kernel has no paired version
        }
    }
    if (!kernel){
        // Cannot be paired: integrate one after the other (same results).
        whfast512_steps_single(r1, N_steps);
        if (r1 != r2) whfast512_steps_single(r2, N_steps);
        return;
    }
    r1->status = REB_STATUS_RUNNING;
    r2->status = REB_STATUS_RUNNING;
    const int skip1 = whfast512_prepare(r1, w1);
    const int skip2 = whfast512_prepare(r2, w2);
    uint64_t done = 0;
    int st1 = 0, st2 = 0;
    const int nsys = w1->four_lanes ? 2 : (int)w1->N_systems;
    const int stopped = kernel(w1->data, w2->data, N_steps, &done, skip1, skip2, &reb_sigint, &st1, &st2, w1->gr_potential, nsys, enc1, esc1);
    for (int i=0; i<2; i++){
        struct reb_simulation* r = rs[i];
        r->t += r->dt*done;
        r->dt_last_done = r->dt;
        r->is_synchronized = 0;
    }
    if (stopped & 1) r1->status = st1;
    if (stopped & 2) r2->status = st2;
    // If only one of them stopped (encounter or escape), the other one continues alone.
    for (int i=0; i<2; i++){
        struct reb_simulation* r = rs[i];
        if (r->status < 0 && done < N_steps && !reb_sigint){
            uint64_t n = N_steps - done;
            r->status = whfast512_full_steps(r, r->integrator.state, &n, 0);
            r->t += r->dt*n;
        }
    }
    for (int i=0; i<2; i++){
        reb_simulation_synchronize(rs[i]);
        if (rs[i]->status <= 0) rs[i]->status = REB_STATUS_SUCCESS;
    }
}

#ifdef DEBUG_AVX512
uint64_t reb_integrator_whfast512_debug_counter(struct reb_simulation* r, int test_p){
    struct reb_integrator_whfast512_state* whfast512 = r->integrator.state;
    struct simd_data* data = whfast512->data;
    uint64_t i[8];
    w512_movu_pd(&i[0], &data->counter);
    return i[test_p];
}

void reb_integrator_whfast512_debug_kepler_step(struct reb_simulation* r, long N_steps){
    recalculate_constants(r, 1);
    struct reb_integrator_whfast512_state* whfast512 = r->integrator.state;
    struct simd_data* data = whfast512->data;
    inertial_to_jacobi_posvel(r, data, 1);
    for(long i=0;i<N_steps;i++){
        whfast512_kepler_step(whfast512);
    }
    jacobi_to_inertial_posvel_and_com(r, data,0.0, 1);
}

// Debug function to print vectors
static inline void printavx512(void* a) {
    double _nax[8];
    w512_movu_pd(&_nax[0], a);
    printf("avx = {%.17g, %.17g, %.17g, %.17g, %.17g, %.17g, %.17g, %.17g}\n", _nax[0], _nax[1], _nax[2], _nax[3], _nax[4], _nax[5], _nax[6], _nax[7]);
}

// Print mask in binary format for debuggin
static inline void printmask8(__mmask8 mask) {
    for (int i = 7; i >= 0; i--) {
        printf("%d", (mask >> i) & 1);
        if (i == 4) printf(" "); 
    }
    printf("\n");
}

// Debug function to print 8x8 matrix
static inline void printmat8(double* a) {
    for (int i=0; i<8;i++){
        for (int j=0; j<8;j++){
            printf("%.16f ", a[i*8+j]);
        }
        printf("\n");
    }
}
#endif // DEBUG_AVX512

