// Energy ensemble / bias test (paper Appendix C style).
// usage: ens kernel N_steps n_runs run_offset perturb out.txt [dt_days]
// Each run: Solar System, dt (default 6 d), GR, corrector 17; x of every planet perturbed by
// perturb*x*(uniform -1..1) (deterministic per run index). Records dE/E (incl. GR potential) at
// ~log-spaced checkpoints. Output: run step dE/E
#include "rebound.h"
#include "integrator_whfast512.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ss_ics.h"
static const double C_LIGHT=10065.32;
static double energy(struct reb_simulation* r){
    double E = reb_simulation_energy(r); struct reb_particle* p=r->particles; double m0=p[0].m;
    for (int i=1;i<r->N;i++){ double dx=p[i].x-p[0].x, dy=p[i].y-p[0].y, dz=p[i].z-p[0].z;
        E += p[i].m*(-3.*m0*m0/(C_LIGHT*C_LIGHT))/(dx*dx+dy*dy+dz*dz); }
    return E;
}
static uint64_t rng_state;
static double urand(){ rng_state = rng_state*6364136223846793005ULL + 1442695040888963407ULL; return ((rng_state>>11)*(1.0/9007199254740992.0))*2-1; }
int main(int argc, char** argv){
    int kernel=atoi(argv[1]); double N=atof(argv[2]); int nruns=atoi(argv[3]); int off=atoi(argv[4]); double pert=atof(argv[5]);
    FILE* f=fopen(argv[6],"w"); double dtd = argc>7 ? atof(argv[7]) : 6.0;
    for (int run=off; run<off+nruns; run++){
        rng_state = 0x9E3779B97F4A7C15ULL*(uint64_t)(run+1);
        struct reb_simulation* r = ss_create(dtd, 0);
        for (int i=1;i<9;i++){ r->particles[i].x *= 1+pert*urand(); r->particles[i].y *= 1+pert*urand(); }
        reb_simulation_move_to_com(r);
        reb_simulation_set_integrator(r, "whfast512");
        struct reb_integrator_whfast512_state* w = r->integrator.state;
        w->gr_potential=1; w->corrector=17; w->kernel=kernel;
        double E0=energy(r); double done=0;
        for (double cp=1e3; done<N; cp*=pow(10,0.25)){
            double target = cp<N ? floor(cp) : N;
            uint64_t n = (uint64_t)(target-done); if (n==0) continue;
            w->concatenate_steps = n; reb_simulation_steps(r,1); done=target;
            fprintf(f,"%d %.0f %.6e\n", run, done, (energy(r)-E0)/E0); fflush(f);
        }
        reb_simulation_free(r);
    }
    fclose(f);
}
