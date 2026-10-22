// Round 5, step 1: pint_disp6k (one kernel per pair) vs pint_disp6 (one FFT + evaluation per element): agreement, time
#include <math.h>
#include <stdio.h>
#include <time.h>
#include "pint_coarse.h"
#include "../tools/ss_ics.h"
void pint_disp6(const struct pint_cfg* c, const double E[PINT_NP][6], const double* n_div, double disp[PINT_NI][6]);
void pint_disp6k(const struct pint_cfg* c, const double E[PINT_NP][6], const double* n_div, double disp[PINT_NI][6]);
extern double pint_disp6k_t[3]; extern int pint_disp6k_dual;
static double now(){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec; }
int main(){
    struct reb_simulation* r=ss_create(6.,0); double E[PINT_NP][6]; pint_sim_to_el(r,E);
    for (int dual=0; dual<=1; dual++) for (int hi=0; hi<=1; hi++){ pint_disp6k_dual=dual;
        struct pint_cfg c; pint_default_cfg(&c,ss_mass,360000.); if (hi){ c.ng_hi=64; c.alpha_hi=0.7; }
        double A[PINT_NI][6], B[PINT_NI][6]; pint_disp6(&c,E,NULL,A); pint_disp6k(&c,E,NULL,B);
        double md=0, mr=0; for (int i=0;i<PINT_NI;i++) for (int e=0;e<6;e++){ const double d=fabs(A[i][e]-B[i][e]); if (d>md) md=d; if (fabs(A[i][e])>mr) mr=fabs(A[i][e]); }
        const int R=200; double t0=now(); for (int q=0;q<R;q++) pint_disp6(&c,E,NULL,A); double ta=(now()-t0)/R;
        double tg=0,tf=0,tk=0; t0=now(); for (int q=0;q<R;q++){ pint_disp6k(&c,E,NULL,B); tg+=pint_disp6k_t[0]; tf+=pint_disp6k_t[1]; tk+=pint_disp6k_t[2]; } double tb=(now()-t0)/R;
        printf("dual %d ng_hi %2d: max |diff| %.2e (max |disp| %.2e); pint_disp6 %.3f ms, pint_disp6k %.3f ms (grids %.3f, pair fills %.3f, kernels+sums %.3f)\n",
               dual, hi?64:0, md, mr, ta*1e3, tb*1e3, tg/R*1e3, tf/R*1e3, tk/R*1e3);
    }
}
