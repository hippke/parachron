// Parareal driver around WHFast512 (fine: C kernel, dt = 6 d, GR, corrector 17) with the mean-element coarse
// propagator of pint_coarse.c. Threads (pthreads) run the fine slices of each iteration; the coarse sweep is serial.
// Reports per-iteration errors against a continuous serial fine run of the same window, wall-clock times,
// and the measured costs needed by the speedup model.
//   usage: pint_run [N slices] [slice days] [threads] [kmax] [dt_out days] [start days]
#include "pint_coarse.h"
#include "integrator_whfast512.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <time.h>
#include "../tools/ss_ics.h"
#define NMAX 1024
static double now(){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec; }
typedef double el[PINT_NP][6];
static double slice_days; static int kernel_F=2;
static void fine(const el Ein, el Eout, int kernel){
    struct reb_simulation* s=pint_el_to_sim(ss_mass,Ein);
    s->dt=6.*DAY; s->exact_finish_time=0;
    reb_simulation_set_integrator(s,"whfast512");
    struct reb_integrator_whfast512_state* w=s->integrator.state;
    w->gr_potential=1; w->corrector=17; w->kernel=kernel; w->concatenate_steps=lround(slice_days/6.);
    reb_simulation_steps(s,1);
    pint_sim_to_el(s,Eout);
    reb_simulation_free(s);
}
static void el_pos(const el E, double X[PINT_NP+1][3]){
    struct reb_simulation* s=pint_el_to_sim(ss_mass,E);
    for (int i=0;i<=PINT_NP;i++){ X[i][0]=s->particles[i].x; X[i][1]=s->particles[i].y; X[i][2]=s->particles[i].z; }
    reb_simulation_free(s);
}
static double poserr(const el A, const el B){
    double XA[PINT_NP+1][3], XB[PINT_NP+1][3]; el_pos(A,XA); el_pos(B,XB);
    double m=0; for (int i=1;i<=PINT_NP;i++) for (int d=0;d<3;d++){ double e=fabs(XA[i][d]-XB[i][d]); if (e>m) m=e; }
    return m;
}
static double wrapd(double x){ return x - 2.*M_PI*floor((x+M_PI)/(2.*M_PI)); }
// thread work
static el *gU, *gF; static int gN, gP, gfirst; static double tF[NMAX];
static void* worker(void* arg){
    const int t=(int)(long)arg;
    for (int n=gfirst+t; n<gN; n+=gP){ double t0=now(); fine(gU[n],gF[n],kernel_F); tF[n]=now()-t0; }
    return NULL;
}
int main(int argc, char** argv){
    const int N = argc>1 ? atoi(argv[1]) : 10;
    slice_days = argc>2 ? atof(argv[2]) : 360000.;
    const int P = argc>3 ? atoi(argv[3]) : 2;
    const int kmax = argc>4 ? atoi(argv[4]) : 5;
    const double dt_out = argc>5 ? atof(argv[5]) : 360.;
    const double start = argc>6 ? atof(argv[6]) : 0.;
    struct pint_cfg c; pint_default_cfg(&c, ss_mass, slice_days); c.dt_out_days=dt_out;
    if (argc>7) c.heun=atoi(argv[7]);
    // initial state (optionally after a fine pre-run)
    struct reb_simulation* r0=ss_create(6.,0);
    reb_simulation_set_integrator(r0,"whfast512");
    { struct reb_integrator_whfast512_state* w=r0->integrator.state; w->gr_potential=1; w->corrector=17; w->kernel=2;
      if (start>0){ w->concatenate_steps=lround(start/6.); reb_simulation_steps(r0,1); } }
    el E0; pint_sim_to_el(r0,E0);
    reb_simulation_free(r0);
    // serial references: continuous fine run (kernel 2) and the noise run (kernel 1, asm)
    static el ref[NMAX+1], ref1[NMAX+1];
    double tser=0;
    for (int kk=0; kk<2; kk++){
        el* R = kk ? ref1 : ref;
        struct reb_simulation* s=pint_el_to_sim(ss_mass,E0);
        s->dt=6.*DAY; s->exact_finish_time=0;
        reb_simulation_set_integrator(s,"whfast512");
        struct reb_integrator_whfast512_state* w=s->integrator.state;
        w->gr_potential=1; w->corrector=17; w->kernel= kk ? 1 : 2; w->concatenate_steps=lround(slice_days/6.);
        memcpy(R[0],E0,sizeof(el));
        double t0=now();
        for (int n=0;n<N;n++){ reb_simulation_steps(s,1); pint_sim_to_el(s,R[n+1]); }
        if (!kk) tser=now()-t0;
        reb_simulation_free(s);
    }
    double noise=0; for (int n=1;n<=N;n++){ double e=poserr(ref[n],ref1[n]); if (e>noise) noise=e; }
    printf("N=%d slices x %.0f d (%.0f yr), threads=%d, giants dt=%.0f d\n", N, slice_days, slice_days/365.25, P, dt_out);
    printf("serial fine run: %.3f s (%.3f ms/slice); noise (asm vs C) max pos err %.2e AU\n", tser, 1e3*tser/N, noise);
    // parareal
    static el U[NMAX+1], Gold[NMAX+1], Fv[NMAX+1], Un[NMAX+1];
    double twall0=now(), tG=0;
    memcpy(U[0],E0,sizeof(el));
    for (int n=0;n<N;n++){ double t0=now(); pint_coarse(&c,U[n],U[n+1]); tG+=now()-t0; memcpy(Gold[n+1],U[n+1],sizeof(el)); }
    double e0=0; for (int n=1;n<=N;n++){ double e=poserr(U[n],ref[n]); if (e>e0) e0=e; }
    printf("k=0  max pos err %.2e AU   (coarse sweep %.3f ms)\n", e0, 1e3*tG);
    double tGsum=tG, tFwall=0, tFtask=0; int kconv=-1; int nG=N;
    for (int k=1;k<=kmax;k++){
        // fine sweep, slices k-1..N-1 (the first k-1 boundaries are exact already)
        gU=U; gF=Fv; gN=N; gP=P; gfirst=k-1;
        double t0=now();
        pthread_t th[256];
        for (int t=0;t<P;t++) pthread_create(&th[t],NULL,worker,(void*)(long)t);
        for (int t=0;t<P;t++) pthread_join(th[t],NULL);
        double t1=now(); tFwall+=t1-t0;
        for (int n=k-1;n<N;n++) tFtask+=tF[n];
        // coarse sweep + correction
        memcpy(Un[0],E0,sizeof(el));
        for (int n=0;n<k;n++) memcpy(Un[n],U[n],sizeof(el));
        memcpy(Un[k],Fv[k-1],sizeof(el));     // exact
        double t2=now();
        for (int n=k;n<N;n++){
            el g; pint_coarse(&c,Un[n],g); nG++;
            for (int i=0;i<PINT_NP;i++) for (int e=0;e<6;e++){
                double d=Fv[n][i][e]-Gold[n+1][i][e]; if (e==1) d=wrapd(d);
                Un[n+1][i][e]=g[i][e]+d;
            }
            Un[n+1][0][1]=Un[n+1][0][1]; for (int i=0;i<PINT_NP;i++) Un[n+1][i][1]=wrapd(Un[n+1][i][1]);
            memcpy(Gold[n+1],g,sizeof(el));
        }
        double t3=now(); tGsum+=t3-t2;
        double dmax=0; for (int n=1;n<=N;n++){ double e=poserr(Un[n],U[n]); if (e>dmax) dmax=e; }
        memcpy(U,Un,sizeof(el)*(N+1));
        double em=0; for (int n=1;n<=N;n++){ double e=poserr(U[n],ref[n]); if (e>em) em=e; }
        printf("k=%d  max pos err %.2e AU  (update %.1e)  fine sweep %.3f s, coarse sweep %.3f ms%s\n", k, em, dmax, t1-t0, 1e3*(t3-t2),
               em<=3*noise ? "  CONVERGED (<= 3x noise)" : "");
        if (em<=3*noise && kconv<0){ kconv=k; break; }
    }
    const double twall=now()-twall0;
    const double TF=tser/N, TGm=tGsum/nG, r=TGm/TF;
    printf("\nmeasured: wall %.3f s on %d threads vs serial %.3f s -> speedup %.2f (this machine)\n", twall, P, tser, tser/twall);
    printf("costs: fine slice %.3f ms, coarse step %.3f ms -> r = %.4f; k = %d; fine core-time %.3f s\n", 1e3*TF, 1e3*TGm, r, kconv, tFtask);
    if (kconv>0){
        printf("model on N=%d cores (one slice per core, pipelined): S = N/(N r + k(1+r)) = %.2f, efficiency %.0f%%\n",
               N, N/(N*r+kconv*(1+r)), 100.*N/(N*r+kconv*(1+r))/N);
        printf("model, non-pipelined (this driver's schedule): S = N/((k+1) N r + k) = %.2f\n", N/((kconv+1)*N*r+kconv));
    }
    return 0;
}
