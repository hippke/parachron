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
};
void pint_default_cfg(struct pint_cfg* c, const double* masses, double slice_days);
// G: E_in -> E_out over one slice
void pint_coarse(const struct pint_cfg* c, const double E_in[PINT_NP][6], double E_out[PINT_NP][6]);
// conversions (Jacobi, as REBOUND's default orbits)
void pint_sim_to_el(struct reb_simulation* r, double E[PINT_NP][6]);
struct reb_simulation* pint_el_to_sim(const double* masses, const double E[PINT_NP][6]);
// diagnostics: time spent in parts of the last pint_coarse call (seconds)
extern __thread double pint_t_tables, pint_t_outer, pint_t_total;
#endif
