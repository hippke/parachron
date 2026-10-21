// Pipelined parareal (one thread per slice) around WHFast512 with the mean-element coarse propagator.
// Thread n owns slice n: for each iteration k it waits for its start state U^k_n, computes the coarse step,
// forms U^k_{n+1} = G(U^k_n) + F(U^{k-1}_n) - G(U^{k-1}_n) (+ energy projection), hands it to thread n+1, and then
// runs the fine slice F(U^k_n) for the next iteration. Wall time ~ N T_G + K (T_F + T_G) per window.
// Windows are chained (the next window starts from the converged end state).
//   usage: pint_pipe N slice_days K windows [fine_kernel=2] [serial=1] [pin=1] [tol_AU=0] [Kmax=6] [helpers=0] [parmode=1]
// helpers > 0: the coarse step's tables (and the giants' run) are spread over a pool of helper threads, in the
// k = 0 sweep only (parmode 1) or in every iteration (parmode 2); results are bitwise identical.
// tol > 0: adaptive iteration count per window: at least K iterations, stop at the first iteration whose largest
// boundary update (position change) is <= tol (at most Kmax). The next iteration's fine slices run speculatively
// (in 10 chunks) and are abandoned when the window stops.
// Prints wall-clock of the parareal run and of the serial run of the same span, and the final-state difference.
#define _GNU_SOURCE
#include "pint_coarse.h"
#include "integrator_whfast512.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <semaphore.h>
#include <sched.h>
#include <time.h>
#include "../tools/ss_ics.h"
#define NMAX 128
#define KMAX 12
typedef double el[PINT_NP][6];
static double now(){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec; }
#include <stdatomic.h>
static double slice_days; static int fkern=2, K, N, pin=1, Kmax=6, npar=0, parmode=1; static double tol=0.;
static atomic_int stop_at;               // iteration at which the current window stopped (-1: running)
static double upd_max[KMAX+1]; static pthread_mutex_t upd_mx=PTHREAD_MUTEX_INITIALIZER;
static double poserr_el(const double A[PINT_NP][6], const double B[PINT_NP][6]);
static struct pint_cfg cfg;
static double energy_sim(struct reb_simulation* s){
    double E=reb_simulation_energy(s); const double c=10065.32, m0=s->particles[0].m;
    for (unsigned i=1;i<s->N;i++){
        const double dx=s->particles[i].x-s->particles[0].x, dy=s->particles[i].y-s->particles[0].y, dz=s->particles[i].z-s->particles[0].z;
        E += -3.*m0*m0*s->particles[i].m/(c*c*(dx*dx+dy*dy+dz*dz)); }
    return E;
}
static double energy_el(const el E){ struct reb_simulation* s=pint_el_to_sim(ss_mass,E); double e=energy_sim(s); reb_simulation_free(s); return e; }
static void project_energy(el E, double Et){
    struct reb_simulation* s=pint_el_to_sim(ss_mass,E);
    for (int it=0; it<2; it++){
        const double e=energy_sim(s);
        double Kin=0; for (unsigned i=0;i<s->N;i++){ struct reb_particle p=s->particles[i]; Kin+=0.5*p.m*(p.vx*p.vx+p.vy*p.vy+p.vz*p.vz); }
        const double f=sqrt(1.+(Et-e)/Kin);
        for (unsigned i=0;i<s->N;i++){ s->particles[i].vx*=f; s->particles[i].vy*=f; s->particles[i].vz*=f; }
    }
    pint_sim_to_el(s,E); reb_simulation_free(s);
}
// fine slice; returns 0 if abandoned because the window stopped at an iteration < kneed
static int fine(const el Ein, el Eout, int kneed){
    struct reb_simulation* s=pint_el_to_sim(ss_mass,Ein);
    s->dt=6.*DAY; s->exact_finish_time=0;
    reb_simulation_set_integrator(s,"whfast512");
    struct reb_integrator_whfast512_state* w=s->integrator.state;
    w->gr_potential=1; w->corrector=17; w->kernel=fkern;
    const long nst=lround(slice_days/6.);
    if (tol>0. && nst%10==0){
        w->concatenate_steps=nst/10;
        for (int q=0;q<10;q++){
            const int st=atomic_load(&stop_at);
            if (st>=0 && st<kneed){ reb_simulation_free(s); return 0; }
            reb_simulation_steps(s,1);
        }
    } else { w->concatenate_steps=nst; reb_simulation_steps(s,1); }
    pint_sim_to_el(s,Eout); reb_simulation_free(s);
    return 1;
}
static double poserr_el(const double A[PINT_NP][6], const double B[PINT_NP][6]){
    struct reb_simulation* a=pint_el_to_sim(ss_mass,A); struct reb_simulation* b=pint_el_to_sim(ss_mass,B);
    double m=0; for (int i=1;i<=PINT_NP;i++){ double d[3]={a->particles[i].x-b->particles[i].x,a->particles[i].y-b->particles[i].y,a->particles[i].z-b->particles[i].z}; for(int x=0;x<3;x++) if (fabs(d[x])>m) m=fabs(d[x]); }
    reb_simulation_free(a); reb_simulation_free(b); return m;
}
static double wrapd(double x){ return x - 2.*M_PI*floor((x+M_PI)/(2.*M_PI)); }
static el U[KMAX+1][NMAX+1];
static sem_t ready[KMAX+1][NMAX+1];
static double tbusyF[NMAX], tbusyG[NMAX];
// optional event trace (env PINT_TRACE=file): per thread (window, k, n, type, t_start, t_end); type 0 = coarse+update,
// 1 = fine completed, 2 = fine abandoned, 3 = waiting for the start state
typedef struct { int w,k,type; double t0,t1; } tev;
static tev* trace[NMAX]; static int ntr[NMAX], trace_on=0, cur_w=0; static double trace_t0;
#define TRACE_MAX 200000
static void tr(int n,int k,int type,double t0,double t1){ if (trace_on && ntr[n]<TRACE_MAX) trace[n][ntr[n]++]=(tev){cur_w,k,type,t0-trace_t0,t1-trace_t0}; }
static void* slice_thread(void* arg){
    const int n=(int)(long)arg;
    if (pin){ cpu_set_t cs; CPU_ZERO(&cs); CPU_SET(n % CPU_SETSIZE, &cs); pthread_setaffinity_np(pthread_self(), sizeof cs, &cs); }
    el Fprev, Gold; double EFprev=0, EUprev=0; int have_F=0;
    const int kend = tol>0. ? Kmax : K;
    for (int k=0;k<=kend;k++){
        { const int st=atomic_load(&stop_at); if (st>=0 && k>st) break; }
        { double tw=now(); if (n>0) sem_wait(&ready[k][n]); tr(n,k,3,tw,now()); }
        { const int st=atomic_load(&stop_at); if (st>=0 && k>st) break; }
        double t0=now();
        el g;
        // G is not needed for boundaries that are already exact (n <= k-1): skip it (saves coarse time on the chain)
        if (k==0 || n>k-1) pint_coarse_par(&cfg,U[k][n],g, npar>0 && (parmode==2 || k==0));
        if (k==0 || n<k-1){
            if (k==0) memcpy(U[0][n+1],g,sizeof(el)); else memcpy(U[k][n+1],U[k-1][n+1],sizeof(el));
        } else if (n==k-1){
            memcpy(U[k][n+1],Fprev,sizeof(el));    // exact
        } else {
            if (!have_F) break;                   // cannot happen unless the window stopped
            for (int i=0;i<PINT_NP;i++) for (int e=0;e<6;e++){
                double d=Fprev[i][e]-Gold[i][e]; if (e==1) d=wrapd(d);
                U[k][n+1][i][e]=g[i][e]+d;
            }
            for (int i=0;i<PINT_NP;i++) U[k][n+1][i][1]=wrapd(U[k][n+1][i][1]);
            project_energy(U[k][n+1], energy_el(U[k][n]) + (EFprev-EUprev));
        }
        memcpy(Gold,g,sizeof(el));
        if (tol>0. && k>0){
            const double u = n>=k-1 ? poserr_el(U[k][n+1],U[k-1][n+1]) : 0.;
            pthread_mutex_lock(&upd_mx); if (u>upd_max[k]) upd_max[k]=u; pthread_mutex_unlock(&upd_mx);
        }
        double t1=now(); tbusyG[n]+=t1-t0; tr(n,k,0,t0,t1);
        if (n<N-1) sem_post(&ready[k][n+1]);
        else {
            // last slice: all updates of iteration k are known -> decide
            int stop = (tol>0.) ? ((k>=K && upd_max[k]<=tol) || k==Kmax) : (k==K);
            if (stop){
                atomic_store(&stop_at,k);
                for (int kk=k+1; kk<=KMAX; kk++) for (int q=1;q<=N;q++) sem_post(&ready[kk][q]);   // release waiting threads
            }
            sem_post(&ready[k][N]);
        }
        have_F=0;
        if (k<kend && n>=k){     // fine slice for the next iteration (speculative when adaptive)
            if (fine(U[k][n],Fprev,k+1)){ EUprev=energy_el(U[k][n]); EFprev=energy_el(Fprev); have_F=1; tr(n,k,1,t1,now()); }
            else tr(n,k,2,t1,now());
        }
        tbusyF[n]+=now()-t1;
    }
    return NULL;
}
int main(int argc, char** argv){
    if (argc<5){ fprintf(stderr,"usage: pint_pipe N slice_days K windows [fine_kernel=2] [serial=1] [pin=1]\n"); return 1; }
    N=atoi(argv[1]); slice_days=atof(argv[2]);
    if (fabs(slice_days/6.-lround(slice_days/6.))>1e-9){ fprintf(stderr,"slice length must be a multiple of the 6 d step\n"); return 1; } K=atoi(argv[3]); const int W=atoi(argv[4]);
    if (argc>5) fkern=atoi(argv[5]);
    const int ser = argc>6 ? atoi(argv[6]) : 1;
    if (argc>7) pin=atoi(argv[7]);
    if (argc>8) tol=atof(argv[8]);
    if (argc>9) Kmax=atoi(argv[9]);
    if (argc>10) npar=atoi(argv[10]);
    if (argc>11) parmode=atoi(argv[11]);
    if (npar>0) pint_set_threads(npar);
    if (Kmax>KMAX-1) Kmax=KMAX-1;
    int khist[KMAX+1]={0};
    pint_default_cfg(&cfg, ss_mass, slice_days);
    if (fkern!=2) cfg.kernel = fkern;
    struct reb_simulation* r0=ss_create(6.,0); el E0; pint_sim_to_el(r0,E0); reb_simulation_free(r0);
    el Eser; double tser=0;
    if (ser){
        struct reb_simulation* s=pint_el_to_sim(ss_mass,E0);
        s->dt=6.*DAY; s->exact_finish_time=0; reb_simulation_set_integrator(s,"whfast512");
        struct reb_integrator_whfast512_state* w=s->integrator.state;
        w->gr_potential=1; w->corrector=17; w->kernel=fkern; w->concatenate_steps=lround(slice_days/6.);
        double t0=now(); for (long q=0;q<(long)W*N;q++) reb_simulation_steps(s,1); tser=now()-t0;
        pint_sim_to_el(s,Eser); reb_simulation_free(s);
    }
    el Ecur; memcpy(Ecur,E0,sizeof(el));
    const char* trf=getenv("PINT_TRACE");
    if (trf){ trace_on=1; for (int n=0;n<N;n++) trace[n]=malloc(TRACE_MAX*sizeof(tev)); }
    const double t0=now(); trace_t0=t0;
    for (int w=0;w<W;w++){
        cur_w=w;
        for (int k=0;k<=KMAX;k++){ memcpy(U[k][0],Ecur,sizeof(el)); upd_max[k]=0; for (int n=0;n<=N;n++) sem_init(&ready[k][n],0,0); }
        atomic_store(&stop_at,-1);
        pthread_t th[NMAX];
        for (int n=0;n<N;n++) pthread_create(&th[n],NULL,slice_thread,(void*)(long)n);
        for (int n=0;n<N;n++) pthread_join(th[n],NULL);
        const int ks=atomic_load(&stop_at);
        khist[ks]++;
        memcpy(Ecur,U[ks][N],sizeof(el));
    }
    const double tpr=now()-t0;
    if (trf){ FILE* f=fopen(trf,"w"); for (int n=0;n<N;n++) for (int i=0;i<ntr[n];i++) fprintf(f,"%d %d %d %d %.7f %.7f\n",trace[n][i].w,n,trace[n][i].k,trace[n][i].type,trace[n][i].t0,trace[n][i].t1); fclose(f); }
    double sF=0,sG=0; for (int n=0;n<N;n++){ sF+=tbusyF[n]; sG+=tbusyG[n]; }
    printf("iterations per window:"); double km=0; for (int k=0;k<=KMAX;k++) if (khist[k]){ printf(" K=%d: %d",k,khist[k]); km+=k*khist[k]; } printf(" (mean %.2f)\n",km/W);
    printf("pipelined parareal: N=%d threads (+%d coarse helpers, mode %d), K=%d, %d windows of %d x %.0f d: wall %.3f s, busy fine %.3f s, coarse %.3f s (core-s)\n",
           N,npar,parmode,K,W,N,slice_days,tpr,sF,sG);
    if (ser){
        struct reb_simulation* a=pint_el_to_sim(ss_mass,Ecur); struct reb_simulation* b=pint_el_to_sim(ss_mass,Eser);
        double m=0; for (int i=1;i<=PINT_NP;i++){ double d[3]={a->particles[i].x-b->particles[i].x,a->particles[i].y-b->particles[i].y,a->particles[i].z-b->particles[i].z}; for(int x=0;x<3;x++) if (fabs(d[x])>m) m=fabs(d[x]); }
        printf("serial: wall %.3f s -> speedup %.2f; core-hours ratio %.2f; final max position difference %.2e AU\n", tser, tser/tpr, (tpr*N)/tser, m);
    }
    return 0;
}
