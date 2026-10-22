// Pipelined parareal (one thread per slice) around WHFast512 with the mean-element coarse propagator.
// Thread n owns slice n: for each iteration k it waits for its start state U^k_n, computes the coarse step,
// forms U^k_{n+1} = G(U^k_n) + F(U^{k-1}_n) - G(U^{k-1}_n) (+ energy projection), hands it to thread n+1, and then
// runs the fine slice F(U^k_n) for the next iteration. Wall time ~ N T_G + K (T_F + T_G) per window.
// Windows are chained (the next window starts from the converged end state).
//   usage: pint_pipe N slice_days K windows [fine_kernel=2] [serial=1] [pin=1] [tol_AU=0] [Kmax=6] [helpers=0] [parmode=1] [M=1]
// M = 2: each thread owns two consecutive slices and runs their fine solves in one paired kernel call
// (reb_integrator_whfast512_steps_pair, bitwise identical to separate calls; N/M threads).
// helpers > 0: the coarse step's tables (and the giants' run) are spread over a pool of helper threads, in the
// k = 0 sweep only (parmode 1) or in every iteration (parmode 2); results are bitwise identical.
// tol > 0: adaptive iteration count per window: at least K iterations, stop at the first iteration whose largest
// boundary update (position change) is <= tol (at most Kmax). The next iteration's fine slices run speculatively
// (in 10 chunks) and are abandoned when the window stops.
// Prints wall-clock of the parareal run and of the serial run of the same span, and the final-state difference.
// Environment: PINT_OUT=prefix writes the accepted boundary states (<prefix>_pr.bin) and, with serial=1, the serial
// run (<prefix>_ser<kernel>.bin), records (t, E_tot incl. GR, E[8][6]) as in pint_long; PINT_PERT=d shifts the x of
// all planets by d AU; PINT_SERONLY=1 runs only the serial reference, PINT_DT its time step; PINT_JR=0/1 WHFast512 jacobi_recursive (default 1); PINT_TRACE, PINT_UPDDIAG: diagnostics.
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
// PINT_JR: WHFast512 jacobi_recursive for fine and serial runs (default 1, no energy offset per call);
// PINT_CHUNKS: kernel calls per speculative fine slice (abandonment granularity; default 2: as fast as 10, 5x fewer
// kernel calls, each of which carries a small energy offset, see RESEARCH_LOG it. 39)
static int jr=1, nchunk=2;
// PINT_BAL=d: load balance between the two CCDs of the 7950X — slices n < N/2 get (1 + d) x the steps, slices
// n >= N/2 (cores 8-15, measured 3-5 % slower under load) (1 - d) x; same window length (rounded to 10 steps)
#define NSL 128
static long nstep_s[NSL]; static double tcum[NSL+1]; static struct pint_cfg cfgs[NSL];
static double slice_days; static int fkern=2, K, N, pin=1, Kmax=6, npar=0, parmode=1, M=1; static double tol=0.;   // M: slices per thread (1 or 2)
// window state is double-buffered by window parity, so that the persistent slice threads can start window w+1 while
// others still finish (abandon) speculative work of window w
static atomic_int stop_at[2];            // iteration at which the window stopped (-1: running)
static double upd_max[2][KMAX+1];
static double upd_pl[2][KMAX+1][PINT_NP]; static int upd_diag=0;
static double upd_el[2][KMAX+1][PINT_NI][4];   // PINT_UPDDIAG=2: inner planets' element updates x a: ln a, lambda, |dz|, |dzeta|   // PINT_UPDDIAG=1: per-planet maxima printed per window
static double poserr_pl(const double A[PINT_NP][6], const double B[PINT_NP][6], double* pl); static pthread_mutex_t upd_mx=PTHREAD_MUTEX_INITIALIZER;
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
// PINT_EDIAG=1: energy bookkeeping of the parareal chain — mean and spread of (a) the element round trip of every
// completed fine result (energy after sim -> elements -> sim minus the energy of the integrated state) and (b) the
// residual of every energy projection (energy of the projected boundary minus its target); relative to |E0|.
static int ediag=0; static double ediag_E0=1.;
typedef struct { double s, s2; long n; } eacc;
static eacc acc_rt, acc_r1; static pthread_mutex_t acc_mx=PTHREAD_MUTEX_INITIALIZER;
static void eadd(eacc* a, double x){ pthread_mutex_lock(&acc_mx); a->s+=x; a->s2+=x*x; a->n++; pthread_mutex_unlock(&acc_mx); }
static void eprint(const char* name, const eacc* a){ const double m=a->s/a->n, sd=sqrt(a->s2/a->n-m*m);
    printf("ediag %s: n %ld, mean %+.3e +- %.1e, rms %.2e (relative to |E0|)\n", name, a->n, m, sd/sqrt(a->n), sd); }
static double energy_sim(struct reb_simulation* s);
// PINT_PROJ: 2 (default) scales the velocities as v += v*d with d = (Et - E)/(2 Kin), which keeps the full relative
// precision of the tiny correction; 1 = the original v *= sqrt(1 + (Et - E)/Kin), whose factor is quantised to the
// spacing of doubles near 1 and left a systematic residual of -1e-16 per projection (iteration 49).
static int projmode=2;
static void project_energy(el E, double Et){
    struct reb_simulation* s=pint_el_to_sim(ss_mass,E);
    for (int it=0; it<2; it++){
        const double e=energy_sim(s);
        double Kin=0; for (unsigned i=0;i<s->N;i++){ struct reb_particle p=s->particles[i]; Kin+=0.5*p.m*(p.vx*p.vx+p.vy*p.vy+p.vz*p.vz); }
        if (projmode==1){
            const double f=sqrt(1.+(Et-e)/Kin);
            for (unsigned i=0;i<s->N;i++){ s->particles[i].vx*=f; s->particles[i].vy*=f; s->particles[i].vz*=f; }
        } else {
            const double d=(Et-e)/(2.*Kin);
            for (unsigned i=0;i<s->N;i++){ s->particles[i].vx+=s->particles[i].vx*d; s->particles[i].vy+=s->particles[i].vy*d; s->particles[i].vz+=s->particles[i].vz*d; }
        }
    }
    pint_sim_to_el(s,E); reb_simulation_free(s);
}
// fine slice; returns 0 if abandoned because the window stopped at an iteration < kneed
static int fine(const el Ein, el Eout, int kneed, atomic_int* stp, int n){
    struct reb_simulation* s=pint_el_to_sim(ss_mass,Ein);
    s->dt=6.*DAY; s->exact_finish_time=0;
    reb_simulation_set_integrator(s,"whfast512");
    struct reb_integrator_whfast512_state* w=s->integrator.state;
    w->gr_potential=1; w->corrector=17; w->jacobi_recursive=jr; w->kernel=fkern;
    const long nst=nstep_s[n];
    if (tol>0. && nchunk>1 && nst%nchunk==0){
        w->concatenate_steps=nst/nchunk;
        for (int q=0;q<nchunk;q++){
            const int st=atomic_load(stp);
            if (st>=0 && st<kneed){ reb_simulation_free(s); return 0; }
            reb_simulation_steps(s,1);
        }
    } else { w->concatenate_steps=nst; reb_simulation_steps(s,1); }
    if (ediag){ const double ed=energy_sim(s); pint_sim_to_el(s,Eout); eadd(&acc_rt,(energy_el(Eout)-ed)/ediag_E0); }
    else pint_sim_to_el(s,Eout);
    reb_simulation_free(s);
    return 1;
}
// two fine slices in one paired kernel call (bitwise identical to two separate calls); same abandonment rule
static struct reb_simulation* fine_sim(const el Ein){
    struct reb_simulation* s=pint_el_to_sim(ss_mass,Ein);
    s->dt=6.*DAY; s->exact_finish_time=0;
    reb_simulation_set_integrator(s,"whfast512");
    struct reb_integrator_whfast512_state* w=s->integrator.state;
    w->gr_potential=1; w->corrector=17; w->jacobi_recursive=jr; w->kernel=fkern;
    return s;
}
static int fine2(const el Ea, const el Eb, el Oa, el Ob, int kneed, atomic_int* stp){
    struct reb_simulation* a=fine_sim(Ea); struct reb_simulation* b=fine_sim(Eb);
    const long nst=lround(slice_days/6.);
    const int nch = (tol>0. && nchunk>1 && nst%nchunk==0) ? nchunk : 1;
    for (int q=0;q<nch;q++){
        const int st=atomic_load(stp);
        if (nch>1 && st>=0 && st<kneed){ reb_simulation_free(a); reb_simulation_free(b); return 0; }
        reb_integrator_whfast512_steps_pair(a,b,nst/nch);
    }
    pint_sim_to_el(a,Oa); pint_sim_to_el(b,Ob); reb_simulation_free(a); reb_simulation_free(b);
    return 1;
}
static double poserr_el(const double A[PINT_NP][6], const double B[PINT_NP][6]){
    struct reb_simulation* a=pint_el_to_sim(ss_mass,A); struct reb_simulation* b=pint_el_to_sim(ss_mass,B);
    double m=0; for (int i=1;i<=PINT_NP;i++){ double d[3]={a->particles[i].x-b->particles[i].x,a->particles[i].y-b->particles[i].y,a->particles[i].z-b->particles[i].z}; for(int x=0;x<3;x++) if (fabs(d[x])>m) m=fabs(d[x]); }
    reb_simulation_free(a); reb_simulation_free(b); return m;
}
static double poserr_pl(const double A[PINT_NP][6], const double B[PINT_NP][6], double* pl){
    struct reb_simulation* a=pint_el_to_sim(ss_mass,A); struct reb_simulation* b=pint_el_to_sim(ss_mass,B);
    double m=0; for (int i=1;i<=PINT_NP;i++){ double q=0; double d[3]={a->particles[i].x-b->particles[i].x,a->particles[i].y-b->particles[i].y,a->particles[i].z-b->particles[i].z}; for(int x=0;x<3;x++) if (fabs(d[x])>q) q=fabs(d[x]); pl[i-1]=q; if (q>m) m=q; }
    reb_simulation_free(a); reb_simulation_free(b); return m;
}
static void write_rec(FILE* f, double t, const el E){ double e=energy_el(E); fwrite(&t,8,1,f); fwrite(&e,8,1,f); fwrite(E,8,48,f); }
static double wrapd(double x){ return x - 2.*M_PI*floor((x+M_PI)/(2.*M_PI)); }
static el U[2][KMAX+1][NMAX+1];
static sem_t ready[2][KMAX+1][NMAX+1], win_start[2], win_done[2];
static int W, khist[KMAX+1]; static el Ecur;
// trajectory output (env PINT_OUT=prefix): the accepted boundary states of every window are copied into a ring by the
// deciding thread and written by the main thread as records (t, E_tot incl. GR, E[8][6]) to <prefix>_pr.bin
#define RING 64
static el* ring[RING]; static sem_t out_ready, out_free; static int out_on=0;
// prepare the state of a window (parity p) that starts from Ecur
static void window_init(int p){
    for (int k=0;k<=KMAX;k++){ memcpy(U[p][k][0],Ecur,sizeof(el)); upd_max[p][k]=0; for (int i=0;i<PINT_NP;i++) upd_pl[p][k][i]=0; memset(upd_el[p][k],0,sizeof upd_el[p][k]); for (int n=0;n<=N;n++){ sem_destroy(&ready[p][k][n]); sem_init(&ready[p][k][n],0,0); } }
    sem_destroy(&win_done[p]); sem_init(&win_done[p],0,0);
    atomic_store(&stop_at[p],-1);
}
static double tbusyF[NMAX], tbusyG[NMAX];
// optional event trace (env PINT_TRACE=file): per thread (window, k, n, type, t_start, t_end); type 0 = coarse+update,
// 1 = fine completed, 2 = fine abandoned, 3 = waiting for the start state
typedef struct { int w,k,type; double t0,t1; } tev;
static tev* trace[NMAX]; static int ntr[NMAX], trace_on=0; static double trace_t0;
#define TRACE_MAX 200000
static void tr(int w,int n,int k,int type,double t0,double t1){ if (trace_on && ntr[n]<TRACE_MAX) trace[n][ntr[n]++]=(tev){w,k,type,t0-trace_t0,t1-trace_t0}; }
// thread c owns slices c*M ... c*M+M-1 (M = 2: both coarse hops back to back, then both fine slices in one paired call)
static void* slice_thread(void* arg){
    const int c=(int)(long)arg, n0=c*M;
    if (pin){ cpu_set_t cs; CPU_ZERO(&cs); CPU_SET(c % CPU_SETSIZE, &cs); pthread_setaffinity_np(pthread_self(), sizeof cs, &cs); }
    el Fprev[2], Gold[2];
    const int kend = tol>0. ? Kmax : K;
    for (int w=0; w<W; w++){
    const int p=w&1; el (*Up)[NMAX+1]=U[p]; sem_t (*rdy)[NMAX+1]=ready[p]; atomic_int* stp=&stop_at[p];
    if (c==0 && w>0) sem_wait(&win_start[p]);
    double EFprev[2]={0,0}, EUprev[2]={0,0}; int have_F[2]={0,0};
    for (int k=0;k<=kend;k++){
        { const int st=atomic_load(stp); if (st>=0 && k>st) break; }
        { double tw=now(); if (n0>0) sem_wait(&rdy[k][n0]); tr(w,c,k,3,tw,now()); }
        { const int st=atomic_load(stp); if (st>=0 && k>st) break; }
        int brk=0;
        double t1=now();
        for (int j=0;j<M;j++){
        const int n=n0+j;
        double t0=now();
        el g;
        // G is not needed for boundaries that are already exact (n <= k-1): skip it (saves coarse time on the chain)
        if (k==0 || n>k-1) pint_coarse_par(&cfgs[n],Up[k][n],g, npar>0 && (parmode==2 || k==0));
        if (k==0 || n<k-1){
            if (k==0) memcpy(Up[0][n+1],g,sizeof(el)); else memcpy(Up[k][n+1],Up[k-1][n+1],sizeof(el));
        } else if (n==k-1){
            memcpy(Up[k][n+1],Fprev[j],sizeof(el));    // exact
        } else {
            if (!have_F[j]){ brk=1; break; }      // cannot happen unless the window stopped
            for (int i=0;i<PINT_NP;i++) for (int e=0;e<6;e++){
                double d=Fprev[j][i][e]-Gold[j][i][e]; if (e==1) d=wrapd(d);
                Up[k][n+1][i][e]=g[i][e]+d;
            }
            for (int i=0;i<PINT_NP;i++) Up[k][n+1][i][1]=wrapd(Up[k][n+1][i][1]);
            { const double Et=energy_el(Up[k][n]) + (EFprev[j]-EUprev[j]);
              project_energy(Up[k][n+1], Et);
              if (ediag) eadd(&acc_r1,(energy_el(Up[k][n+1])-Et)/ediag_E0); }
        }
        memcpy(Gold[j],g,sizeof(el));
        if (tol>0. && k>0){
            double pl[PINT_NP]={0};
            const double u = n>=k-1 ? (upd_diag ? poserr_pl(Up[k][n+1],Up[k-1][n+1],pl) : poserr_el(Up[k][n+1],Up[k-1][n+1])) : 0.;
            pthread_mutex_lock(&upd_mx); if (u>upd_max[p][k]) upd_max[p][k]=u;
            for (int i=0;i<PINT_NP;i++) if (pl[i]>upd_pl[p][k][i]) upd_pl[p][k][i]=pl[i];
            if (upd_diag>1 && n>=k-1) for (int i=0;i<PINT_NI;i++){
                const double* A=Up[k][n+1][i]; const double* B=Up[k-1][n+1][i]; const double a=exp(A[0]);
                const double v[4]={a*fabs(A[0]-B[0]), a*fabs(wrapd(A[1]-B[1])), a*hypot(A[2]-B[2],A[3]-B[3]), a*hypot(A[4]-B[4],A[5]-B[5])};
                for (int e=0;e<4;e++) if (v[e]>upd_el[p][k][i][e]) upd_el[p][k][i][e]=v[e];
            }
            pthread_mutex_unlock(&upd_mx);
        }
        t1=now(); tbusyG[c]+=t1-t0; tr(w,c,k,0,t0,t1);
        if (j<M-1) continue;
        if (n<N-1) sem_post(&rdy[k][n+1]);
        else {
            // last slice: all updates of iteration k are known -> decide
            int stop = (tol>0.) ? ((k>=K && upd_max[p][k]<=tol) || k==Kmax) : (k==K);
            if (stop){
                // every thread has entered this window, so none uses the other parity's state any more: prepare the
                // next window there before releasing anybody
                khist[k]++; memcpy(Ecur,Up[k][N],sizeof(el));
                if (out_on){ sem_wait(&out_free); memcpy(ring[w%RING],&Up[k][1],sizeof(el)*N); sem_post(&out_ready); }
                if (upd_diag) for (int kk=1;kk<=k;kk++){ fprintf(stderr,"UPD %d %d",w,kk); for (int i=0;i<PINT_NP;i++) fprintf(stderr," %.3e",upd_pl[p][kk][i]); if (upd_diag>1) for (int i=0;i<PINT_NI;i++) for (int e=0;e<4;e++) fprintf(stderr," %.3e",upd_el[p][kk][i][e]); fprintf(stderr,"\n"); }
                if (w+1<W) window_init(p^1);
                atomic_store(stp,k);
                for (int q=0;q<N/M;q++) sem_post(&win_done[p]);
                for (int kk=k+1; kk<=KMAX; kk++) for (int q=1;q<=N;q++) sem_post(&rdy[kk][q]);   // release waiting threads
                if (w+1<W) sem_post(&win_start[p^1]);
            }
        }
        }
        if (brk) break;
        // fine slices for the next iteration (speculative when adaptive)
        int need[2]={0,0}, nneed=0;
        for (int j=0;j<M;j++){ have_F[j]=0; if (k<kend && n0+j>=k){ need[j]=1; nneed++; } }
        if (nneed==2){
            if (fine2(Up[k][n0],Up[k][n0+1],Fprev[0],Fprev[1],k+1,stp)){
                for (int j=0;j<2;j++){ EUprev[j]=energy_el(Up[k][n0+j]); EFprev[j]=energy_el(Fprev[j]); have_F[j]=1; }
                tr(w,c,k,1,t1,now());
            } else tr(w,c,k,2,t1,now());
        } else if (nneed==1){
            const int j=need[0]?0:1;
            if (fine(Up[k][n0+j],Fprev[j],k+1,stp,n0+j)){ EUprev[j]=energy_el(Up[k][n0+j]); EFprev[j]=energy_el(Fprev[j]); have_F[j]=1; tr(w,c,k,1,t1,now()); }
            else tr(w,c,k,2,t1,now());
        }
        tbusyF[c]+=now()-t1;
    }
    sem_wait(&win_done[p]);    // a thread can finish k = Kmax before the window is decided
    }
    return NULL;
}
int main(int argc, char** argv){
    if (argc<5){ fprintf(stderr,"usage: pint_pipe N slice_days K windows [fine_kernel=2] [serial=1] [pin=1]\n"); return 1; }
    N=atoi(argv[1]); slice_days=atof(argv[2]);
    if (fabs(slice_days/6.-lround(slice_days/6.))>1e-9){ fprintf(stderr,"slice length must be a multiple of the 6 d step\n"); return 1; } K=atoi(argv[3]); W=atoi(argv[4]);
    if (argc>5) fkern=atoi(argv[5]);
    const int ser = argc>6 ? atoi(argv[6]) : 1;
    if (argc>7) pin=atoi(argv[7]);
    if (argc>8) tol=atof(argv[8]);
    if (argc>9) Kmax=atoi(argv[9]);
    if (argc>10) npar=atoi(argv[10]);
    if (argc>11) parmode=atoi(argv[11]);
    if (argc>12) M=atoi(argv[12]);
    if (M<1 || M>2 || N%M){ fprintf(stderr,"slices per thread must be 1 or 2 and divide N\n"); return 1; }
    if (npar>0) pint_set_threads(npar);
    if (getenv("PINT_UPDDIAG")) upd_diag=atoi(getenv("PINT_UPDDIAG"));
    if (getenv("PINT_JR")) jr=atoi(getenv("PINT_JR"));
    if (getenv("PINT_EDIAG")) ediag=atoi(getenv("PINT_EDIAG"));
    if (getenv("PINT_PROJ")) projmode=atoi(getenv("PINT_PROJ"));
    if (getenv("PINT_CHUNKS")) nchunk=atoi(getenv("PINT_CHUNKS"));
    if (Kmax>KMAX-1) Kmax=KMAX-1;
    pint_default_cfg(&cfg, ss_mass, slice_days);
    {   const double bal = getenv("PINT_BAL") ? atof(getenv("PINT_BAL")) : 0.;
        const long nst=lround(slice_days/6.), d=10*lround(nst*bal/10.);
        if (bal!=0. && (N%2 || M!=1)){ fprintf(stderr,"PINT_BAL needs an even N and M = 1\n"); return 1; }
        tcum[0]=0.;
        for (int n=0;n<N;n++){
            nstep_s[n] = bal==0. ? nst : (n<N/2 ? nst+d : nst-d);
            if (bal==0.) cfgs[n]=cfg; else pint_default_cfg(&cfgs[n], ss_mass, 6.*nstep_s[n]);
            if (fkern!=2) cfgs[n].kernel=fkern;
            tcum[n+1]=tcum[n]+6.*nstep_s[n]*DAY;
        }
        if (bal!=0.) fprintf(stderr,"load balance: %ld / %ld steps per slice\n",nstep_s[0],nstep_s[N-1]);
    }
    if (fkern!=2) cfg.kernel = fkern;
    if (getenv("PINT_GKERNEL")){ cfg.kernel=atoi(getenv("PINT_GKERNEL")); for (int n=0;n<N;n++) cfgs[n].kernel=cfg.kernel; }   // giants' coarse kernel
    const char* outp=getenv("PINT_OUT"); const double pert = getenv("PINT_PERT") ? atof(getenv("PINT_PERT")) : 0.;
    const double Tsl=slice_days*DAY;
    struct reb_simulation* r0=ss_create(6.,pert); el E0; pint_sim_to_el(r0,E0); reb_simulation_free(r0);
    el Eser; double tser=0;
    if (ser){
        struct reb_simulation* s=pint_el_to_sim(ss_mass,E0);
        // PINT_DT: time step of the serial reference [d] (default 6; must divide the slice)
        const double sdt = getenv("PINT_DT") ? atof(getenv("PINT_DT")) : 6.;
        if (fabs(slice_days/sdt-lround(slice_days/sdt))>1e-9){ fprintf(stderr,"PINT_DT must divide the slice\n"); return 1; }
        s->dt=sdt*DAY; s->exact_finish_time=0; reb_simulation_set_integrator(s,"whfast512");
        struct reb_integrator_whfast512_state* w=s->integrator.state;
        w->gr_potential=1; w->corrector=17; w->jacobi_recursive=jr; w->kernel=fkern; w->concatenate_steps=lround(slice_days/sdt);
        FILE* fs=NULL;
        if (outp){ char fn[512]; snprintf(fn,sizeof fn,"%s_ser%d.bin",outp,fkern); fs=fopen(fn,"wb"); write_rec(fs,0.,E0); }
        double t0=now();
        for (long q=0;q<(long)W*N;q++){ reb_simulation_steps(s,1); if (fs){ el E; pint_sim_to_el(s,E); write_rec(fs,(q+1)*Tsl,E); } }
        tser=now()-t0; if (fs) fclose(fs);
        pint_sim_to_el(s,Eser); reb_simulation_free(s);
    }
    if (getenv("PINT_SERONLY")){ printf("serial: wall %.3f s (%d x %d slices)\n",tser,W,N); return 0; }
    memcpy(Ecur,E0,sizeof(el)); ediag_E0=fabs(energy_el(E0));
    const char* trf=getenv("PINT_TRACE");
    if (trf){ trace_on=1; for (int n=0;n<N;n++) trace[n]=malloc(TRACE_MAX*sizeof(tev)); }
    const double t0=now(); trace_t0=t0;
    for (int p=0;p<2;p++){ sem_init(&win_start[p],0,0); sem_init(&win_done[p],0,0); for (int k=0;k<=KMAX;k++) for (int n=0;n<=N;n++) sem_init(&ready[p][k][n],0,0); }
    window_init(0);
    FILE* fo=NULL;
    if (outp){
        char fn[512]; snprintf(fn,sizeof fn,"%s_pr.bin",outp); fo=fopen(fn,"wb"); write_rec(fo,0.,E0);
        out_on=1; sem_init(&out_ready,0,0); sem_init(&out_free,0,RING); for (int q=0;q<RING;q++) ring[q]=malloc(sizeof(el)*N);
    }
    pthread_t th[NMAX];
    for (int c=0;c<N/M;c++) pthread_create(&th[c],NULL,slice_thread,(void*)(long)c);
    if (fo) for (int w=0;w<W;w++){
        sem_wait(&out_ready);
        for (int n=0;n<N;n++) write_rec(fo,(double)w*N*Tsl+tcum[n+1],(const double (*)[6])ring[w%RING][n]);
        sem_post(&out_free);
        if (w%1000==999) fflush(fo);
    }
    for (int c=0;c<N/M;c++) pthread_join(th[c],NULL);
    if (fo) fclose(fo);
    const double tpr=now()-t0;
    if (trf){ FILE* f=fopen(trf,"w"); for (int n=0;n<N;n++) for (int i=0;i<ntr[n];i++) fprintf(f,"%d %d %d %d %.7f %.7f\n",trace[n][i].w,n,trace[n][i].k,trace[n][i].type,trace[n][i].t0,trace[n][i].t1); fclose(f); }
    double sF=0,sG=0; for (int n=0;n<N;n++){ sF+=tbusyF[n]; sG+=tbusyG[n]; }
    printf("iterations per window:"); double km=0; for (int k=0;k<=KMAX;k++) if (khist[k]){ printf(" K=%d: %d",k,khist[k]); km+=k*khist[k]; } printf(" (mean %.2f)\n",km/W);
    printf("pipelined parareal: N=%d threads (+%d coarse helpers, mode %d), K=%d, %d windows of %d x %.0f d: wall %.3f s, busy fine %.3f s, coarse %.3f s (core-s)\n",
           N,npar,parmode,K,W,N,slice_days,tpr,sF,sG);
    if (ediag){ eprint("round trip of fine results",&acc_rt); eprint("projection residual",&acc_r1); }
    if (ser){
        struct reb_simulation* a=pint_el_to_sim(ss_mass,Ecur); struct reb_simulation* b=pint_el_to_sim(ss_mass,Eser);
        double m=0; for (int i=1;i<=PINT_NP;i++){ double d[3]={a->particles[i].x-b->particles[i].x,a->particles[i].y-b->particles[i].y,a->particles[i].z-b->particles[i].z}; for(int x=0;x<3;x++) if (fabs(d[x])>m) m=fabs(d[x]); }
        printf("serial: wall %.3f s -> speedup %.2f; core-hours ratio %.2f; final max position difference %.2e AU\n", tser, tser/tpr, (tpr*N)/tser, m);
    }
    return 0;
}
