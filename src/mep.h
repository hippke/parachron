/**
 * @file    mep.h
 * @brief   Internal header of mean-element parareal (MEP): coarse propagator and element conversions.
 * @details Planet state: Jacobi elements per planet [ln a, lambda, k, h, q, p] (k = e cos pomega, h = e sin pomega,
 *          q = sin(i/2) cos Omega, p = sin(i/2) sin Omega); units G = 1, AU, yr/2pi (as WHFast512).
 *
 * @section LICENSE
 * Copyright (c) 2026 Michael Hippke
 *
 * This file is part of Parachron (https://github.com/hippke/parachron), a modified version of rebound.
 *
 * Parachron is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Parachron is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Parachron.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#ifndef _MEP_H
#define _MEP_H
#include "rebound.h"
#define MEP_NP 8    // planets
#define MEP_NI 4    // inner planets treated by the mean-element theory
struct mep_cfg {
    double masses[MEP_NP+1];
    double slice_days;      // slice length [d]
    double dt_out_days;     // WH step of the giants' coarse run [d]
    int ng_ii;              // grid for inner-inner pairs (power of two)
    int ng_io;              // grid for inner-giant pairs (power of two, <= ng_ii, divides it)
    int maxj;               // last planet index paired with the inner planets (5 = Saturn)
    int kernel;             // WHFast512 C kernel of the giants' coarse run (2 = AVX-512, 3 = AVX2, 4 = generic)
    int heun;               // 1: secular drift with the trapezoidal rule
    double pcut_yr;         // drop short-period terms with period > pcut_yr from the transform (0: keep all)
    int ng_hi;              // grid for close inner-inner pairs (0: off)
    double alpha_hi;
    int meann;              // mean motions of the mean ln a in the small divisors of the start tables: 1: tables rebuilt
                            // at the mean ln a; 2 (default): only the divisor sums re-evaluated (same iteration counts,
                            // nearly free)
    int dual;               // 1: exact dE/dv by dual numbers
    int gcorr;              // symplectic corrector order of the giants' coarse run
    int gr;                 // 1: GR potential of WHFast512 included in the coarse step
    // Synodic grids for inner-inner pairs with a_i/a_j >= alpha_an (default 0.7: Venus-Earth): the pair's rates are
    // tabulated on ng_th x ng_ph points of (lambda_i, phi = lambda_i - lambda_j) instead of ng_ii x ng_ii points of
    // (lambda_i, lambda_j). A close pair's spectrum is wide in phi (synodic harmonics) and narrow in the other direction
    // (eccentricity order), so 16 x 64 resolves it like a 64 x 64 grid at the cost of 32 x 32. ng_th = 0: off.
    int ng_th, ng_ph; double alpha_an;
};
void reb_mep_default_cfg(struct mep_cfg* c, const double* masses, double slice_days, int gr);
void reb_mep_coarse(const struct mep_cfg* c, const double E_in[MEP_NP][6], double E_out[MEP_NP][6]);
void reb_mep_coarse_par(const struct mep_cfg* c, const double E_in[MEP_NP][6], double E_out[MEP_NP][6], int par);
void reb_mep_coarse_giants(const struct mep_cfg* c, const double E_in[MEP_NP][6], double Eo[MEP_NP][6]);
void reb_mep_coarse_pre(const struct mep_cfg* c, const double E_in[MEP_NP][6], double E_out[MEP_NP][6], int par, const double (*Mpre)[6]);
void reb_mep_sim_to_el(struct reb_simulation* r, double E[MEP_NP][6]);
struct reb_simulation* reb_mep_el_to_sim(const double* masses, const double E[MEP_NP][6]);
void reb_mep_el_to_xv(const double* masses, const double E[MEP_NP][6], struct reb_particle p[MEP_NP+1]);
void reb_mep_xv_to_el(const struct reb_particle p[MEP_NP+1], double E[MEP_NP][6]);
void reb_mep_set_threads(int nhelpers, int pin);
int reb_mep_pool_guest(void);
#endif
