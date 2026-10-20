// Accuracy driver for lead S1. Usage: s1_acc method dt_days T_years sample_years gr out.txt
// methods: whc (WHFast Jacobi, corr 17), whckl, whckm, whckc (kernels + corr17 + corrector2),
//          saba:<type int>, w512 (WHFast512 C kernel 2, corr 17), ias15
// Output per sample: t, E_rel_err, then for planets 1..8: a e inc Omega pomega lambda (heliocentric)
#include "rebound.h"
#include "integrator_whfast.h"
#include "integrator_saba.h"
#include "integrator_whfast512.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ss_ics.h"
static int use_gr=0;
static const double C_LIGHT=10065.32;
static void gr_force(struct reb_simulation* r){
    struct reb_particle* p = r->particles; const double m0=p[0].m;
    const double pf = -6.*m0*m0/(C_LIGHT*C_LIGHT);
    for (int i=1;i<r->N;i++){
        double dx=p[i].x-p[0].x, dy=p[i].y-p[0].y, dz=p[i].z-p[0].z; double r2=dx*dx+dy*dy+dz*dz;
        double f=pf/(r2*r2); p[i].ax+=f*dx; p[i].ay+=f*dy; p[i].az+=f*dz;
    }
}
static double energy(struct reb_simulation* r){
    double E = reb_simulation_energy(r);
    if (use_gr){ struct reb_particle* p=r->particles; double m0=p[0].m;
        for (int i=1;i<r->N;i++){ double dx=p[i].x-p[0].x, dy=p[i].y-p[0].y, dz=p[i].z-p[0].z;
            E += p[i].m*(-3.*m0*m0/(C_LIGHT*C_LIGHT))/(dx*dx+dy*dy+dz*dz); } }
    return E;
}
int main(int argc, char** argv){
    if (argc<7){ fprintf(stderr,"usage\n"); return 1; }
    const char* method=argv[1]; double dtd=atof(argv[2]); double T=atof(argv[3])*2*M_PI; double ds=atof(argv[4])*DAY;
    use_gr=atoi(argv[5]); double perturb = argc>7? atof(argv[7]) : 0;
    struct reb_simulation* r = ss_create(dtd, perturb);
    if (!strcmp(method,"w512")){
        reb_simulation_set_integrator(r,"whfast512"); struct reb_integrator_whfast512_state* w=r->integrator.state;
        w->gr_potential=use_gr; w->corrector=17; w->kernel=2; w->concatenate_steps=1;
    } else if (!strncmp(method,"saba:",5)){
        reb_simulation_set_integrator(r,"saba"); struct reb_integrator_saba_state* s=r->integrator.state;
        s->type=atoi(method+5); s->safe_mode=0; s->keep_unsynchronized=1;
        if (use_gr) r->additional_forces=gr_force;
    } else if (!strcmp(method,"ias15")){
        reb_simulation_set_integrator(r,"ias15"); if (use_gr) r->additional_forces=gr_force;
        r->exact_finish_time=1;
    } else {
        reb_simulation_set_integrator(r,"whfast"); struct reb_integrator_whfast_state* w=r->integrator.state;
        w->coordinates=REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI; w->safe_mode=0; w->keep_unsynchronized=1; w->corrector=17;
        if (!strcmp(method,"whckl")){ w->kernel=REB_INTEGRATOR_WHFAST_KERNEL_LAZY; w->corrector2=1; }
        if (!strcmp(method,"whckm")){ w->kernel=REB_INTEGRATOR_WHFAST_KERNEL_MODIFIEDKICK; w->corrector2=1; }
        if (!strcmp(method,"whckc")){ w->kernel=REB_INTEGRATOR_WHFAST_KERNEL_COMPOSITION; w->corrector2=1; }
        if (!strcmp(method,"wh")){ w->corrector=0; }
        if (use_gr) r->additional_forces=gr_force;
    }
    FILE* f=fopen(argv[6],"w");
    double E0=energy(r);
    long nps = lround(ds/r->dt); if (fabs(nps*r->dt-ds)>1e-9*ds){ fprintf(stderr,"dt does not divide sample\n"); return 1; }
    long nsamp = lround(T/ds);
    for (long k=0; k<=nsamp; k++){
        if (k) { if (!strcmp(method,"ias15")) reb_simulation_integrate(r, k*ds); else reb_simulation_steps(r, nps); }
        reb_simulation_synchronize(r);
        fprintf(f,"%.10e %.6e", r->t/(2*M_PI), (energy(r)-E0)/E0);
        for (int i=1;i<9;i++){ struct reb_orbit o=reb_orbit_from_particle(1.,r->particles[i],r->particles[0]);
            fprintf(f," %.15e %.15e %.15e %.15e %.15e %.15e",o.a,o.e,o.inc,o.Omega,o.pomega,o.l); }
        fprintf(f,"\n");
    }
    fclose(f);
    return 0;
}
