// serial WHFast512 run (dt = 6 d, GR, corr17) with all planets' x perturbed by `pert` AU; records like pint_long
#include "pint_coarse.h"
#include "integrator_whfast512.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../tools/ss_ics.h"
static double energy_sim(struct reb_simulation* s){
    double E=reb_simulation_energy(s); const double c=10065.32, m0=s->particles[0].m;
    for (unsigned i=1;i<s->N;i++){ const double dx=s->particles[i].x-s->particles[0].x, dy=s->particles[i].y-s->particles[0].y, dz=s->particles[i].z-s->particles[0].z;
        E += -3.*m0*m0*s->particles[i].m/(c*c*(dx*dx+dy*dy+dz*dz)); }
    return E; }
int main(int argc, char** argv){
    const char* out=argv[1]; const long nsl=atol(argv[2]); const double days=atof(argv[3]); const double pert=atof(argv[4]);
    const int kern = argc>5 ? atoi(argv[5]) : 2; const double dtd = argc>6 ? atof(argv[6]) : 6.;
    struct reb_simulation* s=ss_create(dtd,pert);
    reb_simulation_set_integrator(s,"whfast512");
    struct reb_integrator_whfast512_state* w=s->integrator.state;
    w->gr_potential=1; w->corrector=17; w->kernel=kern; w->concatenate_steps=lround(days/dtd);
    FILE* f=fopen(out,"wb"); double E[PINT_NP][6];
    for (long q=0;q<=nsl;q++){
        if (q) reb_simulation_steps(s,1);
        pint_sim_to_el(s,E); double t=s->t, e=energy_sim(s); fwrite(&t,8,1,f); fwrite(&e,8,1,f); fwrite(E,8,48,f);
    }
    fclose(f); return 0;
}
