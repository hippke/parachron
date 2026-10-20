// ns/step for WHFast512 kernels, paper setup: dt=6d, GR, corrector 17. Interleaved repeats.
#include "rebound.h"
#include "integrator_whfast512.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "ss_ics.h"
static double now(){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec; }
int main(int argc, char** argv){
    long N = argc>1 ? atol(argv[1]) : 20000000;
    int rep = argc>2 ? atoi(argv[2]) : 5;
    double dtd = argc>3 ? atof(argv[3]) : 6.0;
    int kernels[3] = {1,2,3};
    double best[3]={1e9,1e9,1e9}; static double all[3][4096];
    for (int k=0;k<rep;k++) for (int ki=0;ki<3;ki++){
        struct reb_simulation* r = ss_create(dtd, 0);
        reb_simulation_set_integrator(r, "whfast512");
        struct reb_integrator_whfast512_state* w = r->integrator.state;
        w->gr_potential=1; w->corrector=17; w->kernel=kernels[ki]; w->concatenate_steps=N;
        double E0 = reb_simulation_energy(r);
        double t0=now(); reb_simulation_steps(r,1); double t1=now();
        reb_simulation_synchronize(r);
        double ns = (t1-t0)/N*1e9; all[ki][k]=ns; if (ns<best[ki]) best[ki]=ns;
        reb_simulation_free(r);
    }
    for (int ki=0;ki<3;ki++){
        double rr[4096]; for(int k=0;k<rep;k++) rr[k]=all[ki][k]/all[0][k];
        for(int i=0;i<rep;i++)for(int j=i+1;j<rep;j++) if(rr[j]<rr[i]){double t=rr[i];rr[i]=rr[j];rr[j]=t;}
        printf("kernel %d: paired time ratio vs kernel 1: median %.4f [q25 %.4f q75 %.4f] -> speedup %.3fx\n", kernels[ki], rr[rep/2], rr[rep/4], rr[3*rep/4], 1/rr[rep/2]);
        double v[4096]; for(int k=0;k<rep;k++) v[k]=all[ki][k];
        for(int i=0;i<rep;i++)for(int j=i+1;j<rep;j++) if(v[j]<v[i]){double t=v[i];v[i]=v[j];v[j]=t;}
        printf("kernel %d: min %.1f median %.1f ns/step -> 5 Gyr (dt=%gd) = %.2f h\n", kernels[ki], best[ki], v[rep/2], dtd, best[ki]*5e9*365.25/dtd/3.6e12);
    }
}
