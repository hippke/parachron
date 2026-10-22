// Parallel-in-time (parareal) support for WHFast512 Solar System runs: coarse propagator G and element
// conversions. Planet state: Jacobi elements per planet [ln a, lambda, k, h, q, p] (k = e cos pomega,
// h = e sin pomega, q = sin(i/2) cos Omega, p = sin(i/2) sin Omega), units G = 1, AU, yr/2pi.
// G ("config A", whfaster RESEARCH_LOG round 2): inner planets (Mercury..Mars) first-order mean semi-major axes
// (numerical short-period theory on 2-D grids, rebuilt from the input state) + Kepler/secular Euler drift;
// giants by WHFast512 (C kernel) at a large step in a Sun(+inner masses) + 4-planet subsystem.
#ifndef PINT_COARSE_H
#define PINT_COARSE_H
#include "rebound.h"
#define PINT_NP 8
#define PINT_NI 4
struct pint_cfg {
    double masses[PINT_NP+1];
    double slice_days;      // slice length
    double dt_out_days;     // WH step of the giants
    int ng_ii;              // grid for inner-inner pairs (power of two)
    int ng_io;              // grid for inner-giant pairs (power of two, <= ng_ii, divides it)
    int maxj;               // last planet index paired with the inner planets (5 = Saturn)
    int kernel;             // WHFast512 kernel for the giants (2 = C AVX-512, 3 = C AVX2, 4 = generic)
    int heun;               // 1: secular drift with the trapezoidal rule (rates at start and at the predicted end)
    double pcut_yr;         // drop short-period terms with period > pcut_yr from the transform (0: keep all)
    int ng_hi;              // grid for inner-inner pairs with a_i/a_j >= alpha_hi (power of two, >= ng_ii; 0: off)
    double alpha_hi;
    int meann;              // 1: start tables rebuilt at the inner planets' mean ln a (mean motions in the small divisors)
    int dual;               // 1 (default): dE/dv by dual numbers (exact) instead of finite differences in the tables
    int gcorr;              // symplectic corrector order of the giants' coarse WH run (default 17)
};
void pint_default_cfg(struct pint_cfg* c, const double* masses, double slice_days);
// G: E_in -> E_out over one slice
void pint_coarse(const struct pint_cfg* c, const double E_in[PINT_NP][6], double E_out[PINT_NP][6]);
// giants' part of G only (rows PINT_NI.. of Eo), and G with those rows precomputed (bitwise identical to G)
void pint_coarse_giants(const struct pint_cfg* c, const double E_in[PINT_NP][6], double Eo[PINT_NP][6]);
void pint_coarse_pre(const struct pint_cfg* c, const double E_in[PINT_NP][6], double E_out[PINT_NP][6], int par, const double (*Mpre)[6]);
// conversions (Jacobi, as REBOUND's default orbits)
void pint_sim_to_el(struct reb_simulation* r, double E[PINT_NP][6]);
struct reb_simulation* pint_el_to_sim(const double* masses, const double E[PINT_NP][6]);
// diagnostics: time spent in parts of the last pint_coarse call (seconds)
void pint_lna_disp(const struct pint_cfg* c, const double E[PINT_NP][6], double disp[PINT_NI], double R[PINT_NI][6]);
// parallel variant: par != 0 runs the tables on the helper pool (pint_set_threads) when it is free
void pint_coarse_par(const struct pint_cfg* c, const double E_in[PINT_NP][6], double E_out[PINT_NP][6], int par);
void pint_set_threads(int nhelpers);
void pint_disp6(const struct pint_cfg* c, const double E[PINT_NP][6], const double* n_div, double disp[PINT_NI][6]);
extern __thread double pint_t_tables, pint_t_outer, pint_t_total;
#endif
