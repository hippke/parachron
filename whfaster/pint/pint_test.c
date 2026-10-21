// Checks pint_coarse against the Python prototype and times G against the fine slice F.
#include "pint_coarse.h"
#include "integrator_whfast512.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "../tools/ss_ics.h"
static double now(){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec; }
int main(int argc, char** argv){
    double slice = argc>1 ? atof(argv[1]) : 360000.;
    int reps = argc>2 ? atoi(argv[2]) : 20;
    struct reb_simulation* r=ss_create(6.,0);
    double E[PINT_NP][6], Eo[PINT_NP][6];
    pint_sim_to_el(r,E);
    struct pint_cfg c; pint_default_cfg(&c, ss_mass, slice);
    if (argc>3) c.dt_out_days=atof(argv[3]);
    pint_coarse(&c,E,Eo);
    printf("E_in/E_out:\n");
    for (int i=0;i<PINT_NP;i++){ for (int e=0;e<6;e++) printf("%.17g ",Eo[i][e]); printf("\n"); }
    double tt=0, ta=0, to=0; double best=1e9;
    for (int k=0;k<reps;k++){ pint_coarse(&c,E,Eo); tt+=pint_t_total; ta+=pint_t_tables; to+=pint_t_outer; if (pint_t_total<best) best=pint_t_total; }
    extern __thread double pint_prof[4];
    printf("profile per G call: grid r,v,dE/dv %.3f ms, pair loops %.3f ms, FFT %.3f ms, eval %.3f ms\n",1e3*pint_prof[0]/(reps+1),1e3*pint_prof[1]/(reps+1),1e3*pint_prof[2]/(reps+1),1e3*pint_prof[3]/(reps+1));
    printf("G: mean %.3f ms (tables %.3f ms, giants WH %.3f ms), min %.3f ms\n", 1e3*tt/reps, 1e3*ta/reps, 1e3*to/reps, 1e3*best);
    // fine slice
    double bf=1e9;
    for (int k=0;k<3;k++){
        struct reb_simulation* s=pint_el_to_sim(ss_mass,E);
        s->dt=6.*DAY; s->exact_finish_time=0;
        reb_simulation_set_integrator(s,"whfast512");
        struct reb_integrator_whfast512_state* w=s->integrator.state;
        w->gr_potential=1; w->corrector=17; w->kernel=2; w->concatenate_steps=lround(slice/6.);
        double t0=now(); reb_simulation_steps(s,1); double t1=now();
        if (t1-t0<bf) bf=t1-t0;
        reb_simulation_free(s);
    }
    printf("F: min %.3f ms (%.1f ns/step)  ->  r = G/F = %.4f (mean G) / %.4f (min G)\n", 1e3*bf, bf/(slice/6.)*1e9, (tt/reps)/bf, best/bf);
    return 0;
}
