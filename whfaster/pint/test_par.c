#define _GNU_SOURCE
#include <math.h>
// bitwise check and timing of the parallel coarse step
#include "pint_coarse.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "../tools/ss_ics.h"
static double now(){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec; }
int main(int argc, char** argv){
    int nh=argc>1?atoi(argv[1]):3;
    struct pint_cfg cfg; pint_default_cfg(&cfg,ss_mass,360000.);
    struct reb_simulation* r0=ss_create(6.,0); double E[PINT_NP][6]; pint_sim_to_el(r0,E); reb_simulation_free(r0);
    double A[PINT_NP][6], B[PINT_NP][6];
    pint_coarse(&cfg,E,A);
    pint_set_threads(nh);
    pint_coarse_par(&cfg,E,B,1);
    printf("bitwise identical: %s\n", memcmp(A,B,sizeof A)==0 ? "yes":"NO");
    int R=200; double t0=now(); for (int q=0;q<R;q++) pint_coarse(&cfg,E,A); double ts=(now()-t0)/R;
    t0=now(); for (int q=0;q<R;q++) pint_coarse_par(&cfg,E,B,1); double tp=(now()-t0)/R;
    printf("coarse step: serial %.3f ms, parallel (%d helpers + caller) %.3f ms\n",ts*1e3,nh,tp*1e3);
    double d[PINT_NI], Rr[PINT_NI][6];
    t0=now(); for (int q=0;q<R;q++) pint_lna_disp(&cfg,E,d,Rr); double tt=(now()-t0)/R;
    printf("one table pass (serial, with rates) %.3f ms -> giants WH + rest ~ %.3f ms\n",tt*1e3,(ts-2*tt)*1e3);
    return 0;
}
