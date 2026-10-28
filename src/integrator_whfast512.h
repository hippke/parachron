/**
 * integrator_whfast512.h: The AVX512 accelerated symplectic integrator WHFast512 in ASM
 * 
 * Copyright (c) 2026 Hanno Rein, Rishit Dagli, Pejvak Javaheri
 * Modified in 2026 by Michael Hippke for Parachron (https://github.com/hippke/parachron),
 * a modified version of REBOUND; see CHANGES_FROM_REBOUND.md for the list of changes.
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
#ifndef _INTEGRATOR_WHFAST512_H
#define _INTEGRATOR_WHFAST512_H

#include "rebound.h"

extern const struct reb_integrator reb_integrator_whfast512;

struct reb_integrator_whfast512_state {
    unsigned int gr_potential;
    unsigned int N_systems;
    unsigned int corrector;             
    uint64_t concatenate_steps;
    unsigned int kernel;                // 0 (default): fastest available. 1: AVX-512 assembly, 2: C AVX-512, 3: C AVX2, 4: C generic. See enum REB_WHFAST512_KERNEL.
    unsigned int jacobi_recursive;      // 0 (default): inertial <-> Jacobi by the precomputed 8x8 matrices (as upstream). 1: by the recursive transformations (no systematic energy change per synchronization).
    unsigned int compensated;           // 0 (default): plain updates. 1: coordinate updates with compensated summation (C kernels only; lower round-off).

    // Internal use
    double last_synchronization;
    void* data;
    unsigned int kernel_used;           // Kernel used by the last step (set by the integrator, for information).
    unsigned int four_lanes;            // Internal: 1 if the AVX2 kernel runs its four-lane variant (one system, at most four planets).
    // Internal: cache of the constants in data (see recalculate_constants), not saved.
    void* cache;
    double* cache_m;
    size_t cache_N, cache_N_allocated;
    unsigned int cache_N_systems;
    double cache_dt, cache_exit_max, cache_exit_min;
};

#endif
