/**
 * @file 	gravity.h
 * @brief 	Default gravity modules. 
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 *
 * @section LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
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
#ifndef _GRAVITY_H
#define _GRAVITY_H
struct reb_simulation;

// Basic O(N^2) gravity routine
void reb_gravity_basic_calculate_acceleration(struct reb_simulation* r);

// Tree based gravity routine
void reb_gravity_tree_calculate_acceleration(struct reb_simulation* r);

// O(N^2) gravity routine with compensated summation
void reb_gravity_compensated_calculate_acceleration(struct reb_simulation* r);

// Gravity routine for WH-style integrators in Jacobi coordinates
void reb_gravity_jacobi_calculate_acceleration(struct reb_simulation* r);

// Calculates the acceleration for the variational equations.
void reb_gravity_basic_calculate_acceleration_var(struct reb_simulation* r);

// Calculates the jerk (derivative of the acceleration) and applies it to the particles' velocity.
void reb_gravity_basic_calculate_and_apply_jerk(struct reb_simulation* r, const double v);

// Serial half-pair gravity on plain arrays (see gravity.c; used by IAS15's stage-parallel mode).
void reb_gravity_basic_arrays(const size_t N, const double* x, const double* y, const double* z, const double* m, double* ax, double* ay, double* az, const double G, const double softening);
void reb_gravity_basic_rows_arrays(const double* x, const double* y, const double* z, const double* m, double* ax, double* ay, double* az, const size_t i0, const size_t i1, const double G, const double softening);
void reb_gravity_basic_rows_arrays_ex(const double* x, const double* y, const double* z, const double* m, double* ax, double* ay, double* az, const size_t i0, const size_t i1, const size_t starti, const size_t startj, const size_t jhi, const int partner, const struct reb_vec6d gb, const double G, const double softening);
void reb_gravity_compensated_range_arrays(const double* x, const double* y, const double* z, const double* m, double* ax, double* ay, double* az, const size_t k0, const size_t k1, const size_t N, const size_t N_active, const int testparticle_type, const unsigned int gravity_ignore_terms, const double G, const double softening);

#endif
