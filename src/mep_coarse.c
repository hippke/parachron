/**
 * @file    mep_coarse.c
 * @brief   Coarse propagator of mean-element parareal (MEP) for WHFast512 planetary systems.
 * @details Inner four planets: numerical first-order mean-element theory (Jacobi pair accelerations on 2-D grids of
 *          the mean longitudes, 2-D FFT, rebuilt from the input state at every call), Kepler + first-order secular
 *          drift with the trapezoidal rule. Outer four planets: WHFast512 at a large step in a reduced system.
 *          Derived from the research code of the MEP paper (results bitwise identical).
 *          Linux only (pthreads, CPU affinity); on other platforms the file compiles to nothing.
 *
 * @section LICENSE
 * Copyright (c) 2026 Michael Hippke
 *
 * This file is part of Parachron (https://github.com/hippke/parachron), a modified version of rebound.
 *
 * Parachron is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Parachron is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Parachron.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#if defined(__linux__) && !defined(__EMSCRIPTEN__)
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "mep.h"
#include <stdio.h>
#include "integrator_whfast512.h"
// The coarse step is compiled for the baseline instruction set (portable builds, e.g. wheels). Its hot functions get an
// x86-64-v3 (AVX2) clone that is selected at load time by CPU features (ifunc), which makes MEP ~5% faster. Results are
// bitwise identical either way (the builds use -ffp-contract=off, so the clone uses no fused multiply-adds).
#if defined(__x86_64__) && !defined(__AVX2__) && defined(__GNUC__) && !defined(__clang__) && __GNUC__>=12
#define MEP_CLONES __attribute__((target_clones("arch=x86-64-v3","default")))
#else
#define MEP_CLONES
#endif
#include "mep_dual.h"
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>
#define DAY (2.*M_PI/365.25)
#define NGMAX 128
#define CEL(x) ((const double (*)[6])(x))     // non-const -> const 2-D array (ISO C before C2X needs the cast)
static const double CGR = 10065.32;
#if defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
#define MEP_RELAX() _mm_pause()
#elif defined(__aarch64__)
#define MEP_RELAX() __asm__ __volatile__("yield")
#else
#define MEP_RELAX() do {} while (0)
#endif

void reb_mep_default_cfg(struct mep_cfg* c, const double* masses, double slice_days, int gr){
    memcpy(c->masses, masses, sizeof(double)*(MEP_NP+1));
    c->slice_days = slice_days; c->dt_out_days = 360.; c->ng_ii = 32; c->ng_io = 16; c->maxj = 5; c->heun = 1; c->pcut_yr = 0.; c->ng_hi = 0; c->alpha_hi = 0.6;
    c->meann = 0; c->gcorr = 17; c->dual = 1; c->gr = gr;
    // giants' kernel: the C kernels are bitwise identical across instruction sets; use the fastest available one
    c->kernel = 4;
#if (defined(__x86_64__) || defined(__i386__)) && (defined(__GNUC__) || defined(__clang__))
    if (__builtin_cpu_supports("avx512f")) c->kernel = 2;
    else if (__builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma")) c->kernel = 3;
#endif
}
static double wrap(double x){ return x - 2.*M_PI*floor((x+M_PI)/(2.*M_PI)); }

MEP_CLONES static void el_to_rv(double mu, const double* E, double* r, double* v){
    const double a=exp(E[0]), lam=E[1], k=E[2], h=E[3], q=E[4], p=E[5];
    const double e=hypot(k,h), pom=atan2(h,k), si=hypot(q,p), Om=atan2(p,q);
    const double inc=2.*asin(si<1.?si:1.);
    const double M=lam-pom;
    double Ea=M;
    for (int it=0; it<30; it++){ const double dE=(Ea-e*sin(Ea)-M)/(1.-e*cos(Ea)); Ea-=dE; if (fabs(dE)<1e-15) break; }
    const double cE=cos(Ea), sE=sin(Ea), b=sqrt(1.-e*e), n=sqrt(mu/(a*a*a));
    const double xp=a*(cE-e), yp=a*b*sE;
    const double vxp=-a*n*sE/(1.-e*cE), vyp=a*n*b*cE/(1.-e*cE);
    const double om=pom-Om, co=cos(om), so=sin(om), cO=cos(Om), sO=sin(Om), ci=cos(inc), s_i=sin(inc);
    const double P[3]={cO*co-sO*so*ci, sO*co+cO*so*ci, so*s_i};
    const double Q[3]={-cO*so-sO*co*ci, -sO*so+cO*co*ci, co*s_i};
    for (int d=0; d<3; d++){ r[d]=xp*P[d]+yp*Q[d]; v[d]=vxp*P[d]+vyp*Q[d]; }
}
MEP_CLONES static void rv_to_el(double mu, const double* r, const double* v, double* E){
    const double rn=sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]), v2=v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
    const double a=1./(2./rn-v2/mu);
    const double hv[3]={r[1]*v[2]-r[2]*v[1], r[2]*v[0]-r[0]*v[2], r[0]*v[1]-r[1]*v[0]};
    const double hn=sqrt(hv[0]*hv[0]+hv[1]*hv[1]+hv[2]*hv[2]);
    double cz=hv[2]/hn; if (cz>1.) cz=1.; if (cz<-1.) cz=-1.;
    const double inc=acos(cz);
    const double Om=atan2(hv[0],-hv[1]);
    const double vxh[3]={v[1]*hv[2]-v[2]*hv[1], v[2]*hv[0]-v[0]*hv[2], v[0]*hv[1]-v[1]*hv[0]};
    double ev[3]; for(int d=0;d<3;d++) ev[d]=vxh[d]/mu - r[d]/rn;
    const double e=sqrt(ev[0]*ev[0]+ev[1]*ev[1]+ev[2]*ev[2]);
    const double nv[3]={cos(Om),sin(Om),0.};
    const double hu[3]={hv[0]/hn,hv[1]/hn,hv[2]/hn};
    const double nxe[3]={nv[1]*ev[2]-nv[2]*ev[1], nv[2]*ev[0]-nv[0]*ev[2], nv[0]*ev[1]-nv[1]*ev[0]};
    const double om=atan2(nxe[0]*hu[0]+nxe[1]*hu[1]+nxe[2]*hu[2], nv[0]*ev[0]+nv[1]*ev[1]+nv[2]*ev[2]);
    const double pom=Om+om;
    const double rv=r[0]*v[0]+r[1]*v[1]+r[2]*v[2];
    const double sE=rv/sqrt(mu*a), cE=1.-rn/a;
    const double Ea=atan2(sE,cE);
    const double M=Ea-e*sin(Ea);
    const double si=sin(inc/2.);
    E[0]=log(a); E[1]=M+pom; E[2]=e*cos(pom); E[3]=e*sin(pom); E[4]=si*cos(Om); E[5]=si*sin(Om);
}
// dE/dv: ln a analytically (d ln a = 2a v.dv/mu); the other components (if ncomp > 1) by forward differences
static __thread int dual_t=0;
MEP_CLONES static void dEdv(double mu, const double* r, const double* v, double D[6][3], int ncomp){
    if (dual_t && ncomp>1){ mep_dedv_dual(mu,r,v,D); return; }
    const double rn=sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]), v2=v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
    const double a=1./(2./rn-v2/mu);
    for (int c=0;c<3;c++) D[0][c]=2.*a*v[c]/mu;
    if (ncomp<=1) return;
    const double vs=sqrt(v2)*1e-7;
    double E0[6]; rv_to_el(mu,r,v,E0);
    for (int c=0;c<3;c++){
        double vp[3]={v[0],v[1],v[2]}; vp[c]+=vs;
        double Ep[6]; rv_to_el(mu,r,vp,Ep);
        for (int e=1;e<ncomp;e++){
            double d=Ep[e]-E0[e];
            if (e==1) d=wrap(d);
            D[e][c]=d/vs;
        }
    }
}
// in-place radix-2 complex FFT (forward, e^{-i}), n power of two (<= NGMAX), stride s; cached twiddles
static double complex tw[NGMAX+1][NGMAX/2];
static pthread_once_t tw_once = PTHREAD_ONCE_INIT;
static void tw_init(void){ for (int n=2;n<=NGMAX;n<<=1) for (int j=0;j<n/2;j++) tw[n][j]=cos(-2.*M_PI*j/n)+I*sin(-2.*M_PI*j/n); }
MEP_CLONES static void fft1(double complex* a, int n, int s){
    pthread_once(&tw_once, tw_init);
    double complex t[NGMAX];
    for (int i=0,j=0;i<n;i++){
        t[j]=a[i*s];
        int bit=n>>1; for (; j&bit; bit>>=1) j^=bit; j^=bit;
    }
    for (int len=2; len<=n; len<<=1){
        const int h=len/2, step=n/len;
        for (int i=0;i<n;i+=len){
            for (int j=0;j<h;j++){
                const double complex u=t[i+j], v=t[i+j+h]*tw[n][j*step];
                t[i+j]=u+v; t[i+j+h]=u-v;
            }
        }
    }
    for (int i=0;i<n;i++) a[i*s]=t[i];
}
MEP_CLONES static void fft2(double complex* a, int n){
    for (int r=0;r<n;r++) fft1(a+r*n,n,1);
    for (int c=0;c<n;c++) fft1(a+c,n,n);
}
static inline int kf(int a, int n){ return a<n/2 ? a : a-n; }   // numpy fftfreq(n, 1/n)
                                                                // short-period displacements of ln a for the two real grids packed in Z = fft2(fa + i fb) (unnormalised):
                                                                // disp = sum_{k != 0} Im(F_k e^{i k.theta}) / (k.nu), F = fft2(f)/n^2, with F_a(k) = (Z(k) + conj Z(-k))/2,
                                                                // F_b(k) = (Z(k) - conj Z(-k))/(2i)
static __thread double nu_cut_t=0.;
MEP_CLONES static void sp_eval_pair(const double complex* Z, int n, double ni, double nj, double li, double lj, double* da, double* db){
    double complex ei[NGMAX], ej[NGMAX];
    for (int a=0;a<n;a++){ const int K=kf(a,n); ei[a]=cos(K*li)+I*sin(K*li); ej[a]=cos(K*lj)+I*sin(K*lj); }
    double sa=0, sb=0; const double nrm=1./(n*n);
    for (int a=0;a<n;a++){
        const int Ka=kf(a,n); const int am=(n-a)%n;
        for (int b=0;b<n;b++){
            const double kk=Ka*ni+kf(b,n)*nj;
            if (kk==0. || fabs(kk)<nu_cut_t) continue;
            const int bm=(n-b)%n;
            const double complex z=Z[a*n+b], zc=conj(Z[am*n+bm]);
            const double complex Fa=0.5*(z+zc), Fb=-0.5*I*(z-zc);
            const double complex e=ei[a]*ej[b];
            sa+=cimag(Fa*e)/kk; sb+=cimag(Fb*e)/kk;
        }
    }
    *da=sa*nrm; *db=sb*nrm;
}
// ---------------------------------------------------------------------------------------------------------------
// table(): first-order ln a short-period displacements (disp) and secular rates (R) of the inner planets.
// Organised as independent tasks (grids per planet; one task per inner-inner pair / per packed pair of giants; GR per
// planet) whose results are reduced in the fixed serial order, so that the result is bitwise independent of the
// number of threads. With a helper pool (reb_mep_set_threads) the tasks of one call run on several cores.
#define NTASKMAX 32
struct tctx {
    const struct mep_cfg* c; const double (*E)[6]; int nc;
    double Mint[MEP_NP], mu[MEP_NP], n[MEP_NP];
    int N1, N2, Nlo;
    double ri[MEP_NI][NGMAX][3], vi[MEP_NI][NGMAX][3], Di[MEP_NI][NGMAX][6][3], rj[MEP_NP][NGMAX][3];
    // tasks
    int ntask; int ti[NTASKMAX], tj0[NTASKMAX], tj1[NTASKMAX], tkind[NTASKMAX];   // kind 0 inner pair, 1 giants, 2 GR, 3 grid
    double tda[NTASKMAX], tdb[NTASKMAX], tRi[NTASKMAX][2][6], tRj[NTASKMAX][6], tNg[NTASKMAX];
};
MEP_CLONES static void task_grid(struct tctx* x, int p){
    const struct mep_cfg* c=x->c; dual_t=c->dual;
    if (p<MEP_NI){
        double Eg[6]; memcpy(Eg,x->E[p],sizeof Eg);
        for (int g=0; g<x->N1; g++){ Eg[1]=2.*M_PI*g/x->N1; el_to_rv(x->mu[p],Eg,x->ri[p][g],x->vi[p][g]); dEdv(x->mu[p],x->ri[p][g],x->vi[p][g],x->Di[p][g],x->nc); }
    } else if (p<=c->maxj){
        double Eg[6]; memcpy(Eg,x->E[p],sizeof Eg); double vdum[3];
        for (int g=0; g<x->N2; g++){ Eg[1]=2.*M_PI*g/x->N2; el_to_rv(x->mu[p],Eg,x->rj[p][g],vdum); }
    }
}
MEP_CLONES static void task_pair(struct tctx* x, int t){
    const struct mep_cfg* c=x->c; const double* m=c->masses; const double (*E)[6]=x->E; const int nc=x->nc;
    static __thread double complex buf[NGMAX*NGMAX];
    nu_cut_t = c->pcut_yr>0. ? 1./c->pcut_yr : 0.;
    const int i=x->ti[t];
    const int jlo=x->tj0[t], jhi=x->tj1[t];
    const int inner_j = x->tkind[t]==0;
    for (int s=0;s<2;s++) for (int e=0;e<6;e++) x->tRi[t][s][e]=0;
    for (int e=0;e<6;e++) x->tRj[t][e]=0;
    int Ng=0;
    for (int j=jlo;j<=jhi;j++){
        const int hi = inner_j && c->ng_hi>0 && exp(E[i][0]-E[j][0])>=c->alpha_hi;
        Ng = inner_j ? (hi ? x->N1 : x->Nlo) : x->N2; const int si = inner_j ? x->N1/Ng : x->N1/x->N2;
        const int slot=j-jlo;
        double Ri[6]={0}, Rj[6]={0};
        const double pf=inner_j ? -(x->Mint[j]+m[j+1])/x->Mint[j]*m[i+1] : 0.;
        for (int a=0;a<Ng;a++){
            const double* r1=x->ri[i][a*si]; const double* D1=x->Di[i][a*si][0];
            for (int b=0;b<Ng;b++){
                const double* r2 = inner_j ? x->ri[j][b*(x->N1/Ng)] : x->rj[j][b];
                const double d[3]={r2[0]-r1[0],r2[1]-r1[1],r2[2]-r1[2]};
                const double d2=d[0]*d[0]+d[1]*d[1]+d[2]*d[2], id3=1./(d2*sqrt(d2));
                const double rj2=r2[0]*r2[0]+r2[1]*r2[1]+r2[2]*r2[2], irj3=1./(rj2*sqrt(rj2));
                double ai[3];
                for (int q=0;q<3;q++) ai[q]=m[j+1]*(d[q]*id3-r2[q]*irj3);
                const double f=D1[0]*ai[0]+D1[1]*ai[1]+D1[2]*ai[2];
                if (nc>1) for (int e=0;e<6;e++){ const double* De=x->Di[i][a*si][e]; Ri[e]+=De[0]*ai[0]+De[1]*ai[1]+De[2]*ai[2]; }
                if (inner_j){
                    const double rdot=r2[0]*r1[0]+r2[1]*r1[1]+r2[2]*r1[2];
                    double aj[3];
                    for (int q=0;q<3;q++) aj[q]=pf*(d[q]*id3-r2[q]*irj3+r1[q]*irj3-3.*r2[q]*(rdot/rj2*irj3));
                    const double* D2=x->Di[j][b*si][0];
                    const double g=D2[0]*aj[0]+D2[1]*aj[1]+D2[2]*aj[2];
                    buf[a*Ng+b]=f+I*g;
                    if (nc>1) for (int e=0;e<6;e++){ const double* De=x->Di[j][b*si][e]; Rj[e]+=De[0]*aj[0]+De[1]*aj[1]+De[2]*aj[2]; }
                } else {
                    if (slot==0) buf[a*Ng+b]=f; else buf[a*Ng+b]+=I*f;
                }
            }
        }
        for (int e=0;e<6;e++){ x->tRi[t][slot][e]=Ri[e]; if (inner_j) x->tRj[t][e]=Rj[e]; }
    }
    x->tNg[t]=Ng;
    fft2(buf,Ng);
    double da,db;
    if (inner_j){ sp_eval_pair(buf,Ng,x->n[i],x->n[jlo],E[i][1],E[jlo][1],&da,&db); x->tda[t]=da; x->tdb[t]=db; }
    else {
        sp_eval_pair(buf,Ng,x->n[i],x->n[jlo],E[i][1],E[jlo][1],&da,&db); x->tda[t]=da; x->tdb[t]=0;
        if (jhi>jlo){ sp_eval_pair(buf,Ng,x->n[i],x->n[jhi],E[i][1],E[jhi][1],&da,&db); x->tdb[t]=db; }
    }
}
MEP_CLONES static void task_gr(struct tctx* x, int t){
    const double* m=x->c->masses; const int i=x->ti[t]; const int N1=x->N1;
    double complex f1[NGMAX]; double Rg[6]={0};
    const double pre= x->c->gr ? -6.*m[0]*m[0]/(CGR*CGR) : 0.;
    for (int a=0;a<N1;a++){
        const double* r1=x->ri[i][a];
        const double rr=sqrt(r1[0]*r1[0]+r1[1]*r1[1]+r1[2]*r1[2]), r4=rr*rr*rr*rr;
        double ag[3]; for (int q=0;q<3;q++) ag[q]=pre*r1[q]/r4;
        double f=0; for (int q=0;q<3;q++) f+=x->Di[i][a][0][q]*ag[q];
        f1[a]=f;
        if (x->nc>1) for (int e=0;e<6;e++){ double s=0; for (int q=0;q<3;q++) s+=x->Di[i][a][e][q]*ag[q]; Rg[e]+=s; }
    }
    for (int e=0;e<6;e++) x->tRj[t][e]=Rg[e];
    fft1(f1,N1,1);
    double complex s=0.;
    for (int a=0;a<N1;a++){
        const int K=kf(a,N1); if (K==0) continue;
        s += f1[a]/N1*(cos(K*x->E[i][1])+I*sin(K*x->E[i][1]))/(I*K*x->n[i]);
    }
    x->tda[t]=creal(s);
}
// ---- helper pool (fork-join over task indices); workers spin briefly before sleeping (wake-up latency of sleeping
// threads is ~0.1 ms, comparable to a table pass) ----
#include <stdatomic.h>
// A job is a list of stages; each stage is a set of independent tasks, and a stage starts only when all tasks of the
// previous one are finished. Helpers (and the caller, in pool_join) take tasks from per-stage counters.
static pthread_mutex_t pool_mx=PTHREAD_MUTEX_INITIALIZER, pool_busy=PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t pool_cv=PTHREAD_COND_INITIALIZER;
#define MAXST 8
struct pjob { int nst; void (*fn[MAXST])(void*,int); void* arg[MAXST]; int ntask[MAXST]; };
static struct pjob job;
static int pool_n=0; static atomic_int pool_gen=0, pool_active=0, st_next[MAXST], st_done[MAXST];
static void pool_run_tasks(void){
    for (int s=0;s<job.nst;s++){
        for (;;){
            const int t=atomic_fetch_add(&st_next[s],1); if (t>=job.ntask[s]) break;
            job.fn[s](job.arg[s],t); atomic_fetch_add_explicit(&st_done[s],1,memory_order_release);
        }
        while (atomic_load_explicit(&st_done[s],memory_order_acquire)<job.ntask[s]) MEP_RELAX();
    }
}
static void* pool_worker(void* a){
    int gen=(int)(long)a;          // generation at creation (no job can be in flight then)
    for (;;){
        if (atomic_load(&pool_gen)==gen){
            pthread_mutex_lock(&pool_mx);
            while (atomic_load(&pool_gen)==gen) pthread_cond_wait(&pool_cv,&pool_mx);
            pthread_mutex_unlock(&pool_mx);
        }
        gen=atomic_load(&pool_gen);
        pool_run_tasks();
        atomic_fetch_sub_explicit(&pool_active,1,memory_order_release);
    }
    return NULL;
}
// the other hardware thread of logical CPU q's core (from sysfs), -1 if none
static int smt_sibling(int q){
    char fn[128]; snprintf(fn,sizeof fn,"/sys/devices/system/cpu/cpu%d/topology/thread_siblings_list",q);
    FILE* f=fopen(fn,"r"); if (!f) return -1;
    char buf[64]={0}; const int ok = fgets(buf,sizeof buf,f)!=NULL; fclose(f); if (!ok) return -1;
    int a=-1,b=-1; char sep=0;
    if (sscanf(buf,"%d%c%d",&a,&sep,&b)!=3) return -1;      // "0,16" or "0-1"
    if (a==q) return b;
    if (b==q) return a;
    return -1;
}
// Helper threads of the coarse step (created once, kept sleeping between jobs). pin != 0: helper q is pinned to the
// SMT sibling of logical CPU q, i.e. of the core of slice q (in iteration k the cores of slices 0..k-1 have no fine
// work, so these hardware threads are the least contended ones; +5 % on a Ryzen 9 7950X).
void reb_mep_set_threads(int nhelpers, int pin){
    pthread_mutex_lock(&pool_busy);
    for (int q=pool_n;q<nhelpers;q++){ pthread_t th; pthread_create(&th,NULL,pool_worker,(void*)(long)atomic_load(&pool_gen));
        const int cpu = pin ? smt_sibling(q) : -1;
        if (cpu>=0){ cpu_set_t cs; CPU_ZERO(&cs); CPU_SET(cpu,&cs); pthread_setaffinity_np(th,sizeof cs,&cs); }
        pthread_detach(th); }
    if (nhelpers>pool_n) pool_n=nhelpers;
    pthread_mutex_unlock(&pool_busy);
}
// the caller must hold pool_busy; pool_launch starts the helpers on j, pool_join lets the caller help and waits
static void pool_launch(const struct pjob* j){
    job=*j;
    for (int s=0;s<j->nst;s++){ atomic_store(&st_next[s],0); atomic_store(&st_done[s],0); }
    atomic_store(&pool_active,pool_n);
    pthread_mutex_lock(&pool_mx); atomic_fetch_add_explicit(&pool_gen,1,memory_order_release); pthread_cond_broadcast(&pool_cv); pthread_mutex_unlock(&pool_mx);
}
static void pool_join(void){
    pool_run_tasks();
    while (atomic_load_explicit(&pool_active,memory_order_acquire)>0) MEP_RELAX();
}
static int pool_acquire(int par){ return par && pool_n>0 && pthread_mutex_trylock(&pool_busy)==0; }
static void pool_release(void){ pthread_mutex_unlock(&pool_busy); }
// run fn(arg, 0..ntask-1); on the pool if par and the pool is free, else serially in the caller
static void run_tasks(void (*fn)(void*,int), void* arg, int ntask, int par){
    if (!pool_acquire(par)){ for (int t=0;t<ntask;t++) fn(arg,t); return; }
    struct pjob j={1,{fn},{arg},{ntask}}; pool_launch(&j); pool_join(); pool_release();
}
static void grid_task(void* a, int t){ task_grid((struct tctx*)a, t); }
static void pairgr_task(void* a, int t){ struct tctx* x=a; if (x->tkind[t]==2) task_gr(x,t); else task_pair(x,t); }
// table setup: masses, mean motions (planets < nn only; the others are set by table_nmotion) and the task list in the
// serial loop order
static void table_setup(struct tctx* x, const struct mep_cfg* c, const double E[MEP_NP][6], int rates, int nn){
    x->c=c; x->E=E; x->nc= rates ? 6 : 1;
    const double* m=c->masses; double acc=0;
    for (int i=0;i<MEP_NP;i++){ acc+=m[i]; x->Mint[i]=acc; x->mu[i]=x->Mint[i]+m[i+1]; if (i<nn) x->n[i]=sqrt(x->mu[i]/exp(3.*E[i][0])); }
    x->N1=(c->ng_hi>c->ng_ii ? c->ng_hi : c->ng_ii); x->N2=c->ng_io; x->Nlo=c->ng_ii;
    nu_cut_t = c->pcut_yr>0. ? 1./c->pcut_yr : 0.;
    int nt=0;
    for (int i=0;i<MEP_NI;i++){
        for (int j=i+1;j<=c->maxj;j++){
            if (j<MEP_NI){ x->ti[nt]=i; x->tj0[nt]=j; x->tj1[nt]=j; x->tkind[nt]=0; nt++; }
            else if (((j-MEP_NI)&1)==0){ x->ti[nt]=i; x->tj0[nt]=j; x->tj1[nt]= j+1<=c->maxj ? j+1 : j; x->tkind[nt]=1; nt++; }
        }
        x->ti[nt]=i; x->tkind[nt]=2; nt++;
    }
    x->ntask=nt;
}
static void table_nmotion(struct tctx* x, int from){ for (int i=from;i<MEP_NP;i++) x->n[i]=sqrt(x->mu[i]/exp(3.*x->E[i][0])); }
// reduction in the serial order
static void table_reduce(struct tctx* x, double disp[MEP_NI], double R[MEP_NI][6]){
    const int nt=x->ntask;
    for (int i=0;i<MEP_NI;i++){ disp[i]=0; if (R) for (int e=0;e<6;e++) R[i][e]=0; }
    for (int t=0;t<nt;t++){
        const int i=x->ti[t];
        if (x->tkind[t]==0){
            const int j=x->tj0[t]; const double Ng=x->tNg[t];
            if (R) for (int e=0;e<6;e++){ R[i][e]+=x->tRi[t][0][e]/(Ng*Ng); R[j][e]+=x->tRj[t][e]/(Ng*Ng); }
            disp[i]+=x->tda[t]; disp[j]+=x->tdb[t];
        } else if (x->tkind[t]==1){
            const double Ng=x->tNg[t];
            for (int j=x->tj0[t]; j<=x->tj1[t]; j++) if (R) for (int e=0;e<6;e++) R[i][e]+=x->tRi[t][j-x->tj0[t]][e]/(Ng*Ng);
            disp[i]+=x->tda[t]; if (x->tj1[t]>x->tj0[t]) disp[i]+=x->tdb[t];
        } else {
            if (R) for (int e=0;e<6;e++) R[i][e]+=x->tRj[t][e]/x->N1;
            disp[i]+=x->tda[t];
        }
    }
}
static void table_par(const struct mep_cfg* c, const double E[MEP_NP][6], double disp[MEP_NI], double R[MEP_NI][6], int par){
    static __thread struct tctx* x=NULL; if (!x) x=malloc(sizeof *x);
    table_setup(x,c,E,R!=NULL,MEP_NP);
    run_tasks(grid_task,x,c->maxj+1,par);
    run_tasks(pairgr_task,x,x->ntask,par);
    table_reduce(x,disp,R);
}
static void outer_wh(const struct mep_cfg* c, const double E[MEP_NP][6], double Eo[MEP_NP][6]){
    const double* m=c->masses;
    double Mi=0; for (int i=0;i<=MEP_NI;i++) Mi+=m[i];
    struct reb_simulation* s=reb_simulation_create();
    reb_simulation_add_fmt(s,"m",Mi);
    for (int j=MEP_NI;j<MEP_NP;j++){
        const double* e=E[j];
        const double ecc=hypot(e[2],e[3]), pom=atan2(e[3],e[2]), si=hypot(e[4],e[5]), Om=atan2(e[5],e[4]);
        reb_simulation_add_fmt(s,"m a e inc Omega pomega l", m[j+1], exp(e[0]), ecc, 2.*asin(si), Om, pom, e[1]);
    }
    reb_simulation_move_to_com(s);
    // number of giant steps: the step is adjusted (<= dt_out_days) so that the steps cover the slice exactly
    const long nst=(long)ceil(c->slice_days/c->dt_out_days-1e-9);
    s->dt=c->slice_days/nst*DAY; s->exact_finish_time=0;
    reb_simulation_set_integrator(s,"whfast512");
    struct reb_integrator_whfast512_state* w=s->integrator.state;
    w->gr_potential=c->gr; w->corrector=c->gcorr; w->kernel=c->kernel;
    w->concatenate_steps=nst;
    reb_simulation_steps(s,1);
    for (int j=MEP_NI;j<MEP_NP;j++){
        const int idx=j-MEP_NI+1;
        struct reb_particle prim=reb_simulation_com_range(s,0,idx);
        struct reb_orbit o=reb_orbit_from_particle(s->G,s->particles[idx],prim);
        const double si=sin(o.inc/2.);
        Eo[j][0]=log(o.a); Eo[j][1]=o.l; Eo[j][2]=o.e*cos(o.pomega); Eo[j][3]=o.e*sin(o.pomega);
        Eo[j][4]=si*cos(o.Omega); Eo[j][5]=si*sin(o.Omega);
    }
    reb_simulation_free(s);
}
// first (predictor) part of the inner planets' step: mean ln a, Kepler + secular drift with the start rates
static void inner_predict(const struct mep_cfg* c, const double* mu, const double E_in[MEP_NP][6], const double disp[MEP_NI], const double R[MEP_NI][6], double M2[MEP_NP][6]){
    const double T=c->slice_days*DAY;
    for (int i=0;i<MEP_NI;i++){
        const double lna=E_in[i][0]-disp[i];
        const double nbar=sqrt(mu[i]/exp(3.*lna));
        M2[i][0]=lna;
        M2[i][1]=wrap(E_in[i][1]+(nbar+R[i][1])*T);
        for (int e=2;e<6;e++) M2[i][e]=E_in[i][e]+R[i][e]*T;
    }
}
// Heun corrector with the end rates R1, and the end displacement
static void inner_correct(const struct mep_cfg* c, const double* mu, const double E_in[MEP_NP][6], const double R[MEP_NI][6], const double R1[MEP_NI][6], const double disp1[MEP_NI], double M2[MEP_NP][6], double E_out[MEP_NP][6]){
    const double T=c->slice_days*DAY;
    for (int i=0;i<MEP_NI;i++){
        const double lna=M2[i][0], nbar=sqrt(mu[i]/exp(3.*lna));
        M2[i][1]=wrap(E_in[i][1]+(nbar+0.5*(R[i][1]+R1[i][1]))*T);
        for (int e=2;e<6;e++) M2[i][e]=E_in[i][e]+0.5*(R[i][e]+R1[i][e])*T;
    }
    memcpy(E_out,M2,sizeof(double)*MEP_NP*6);
    for (int i=0;i<MEP_NI;i++) E_out[i][0]=M2[i][0]+disp1[i];
}
// parallel schedule (bitwise identical to the serial one): the caller integrates the giants (the longest single
// task) while the helpers build the start tables, form the predictor and build the parts of the end tables that need
// only the inner planets (their grids, the inner-inner pairs, GR); then the giants' grids and the inner-giant pairs.
struct pctx {
    const struct mep_cfg* c; const double (*E_in)[6]; double (*M2)[6]; const double* mu;
    struct tctx *x0, *x1; double disp[MEP_NI], R[MEP_NI][6], E1[MEP_NP][6];
    int ninner, inner_t[NTASKMAX], nouter, outer_t[NTASKMAX];
};
static void p_grid0(void* a, int t){ struct pctx* p=a; task_grid(p->x0,t); }
static void p_pair0(void* a, int t){ struct pctx* p=a; pairgr_task(p->x0,t); }
// mean-n: start tables again with the inner planets' ln a replaced by the mean ln a (the giants' grids stay valid)
static void mean_state(const double E_in[MEP_NP][6], const double disp[MEP_NI], double E1[MEP_NP][6]){
    memcpy(E1,E_in,sizeof(double)*MEP_NP*6);
    for (int i=0;i<MEP_NI;i++) E1[i][0]=E_in[i][0]-disp[i];
}
static void p_mean(void* a, int t){ (void)t; struct pctx* p=a;
    table_reduce(p->x0,p->disp,NULL);
    mean_state(p->E_in,p->disp,p->E1);
    table_setup(p->x0,p->c,(const double (*)[6])p->E1,1,MEP_NP);
}
static void p_pred(void* a, int t){ (void)t; struct pctx* p=a;
    table_reduce(p->x0,p->disp,p->R);
    inner_predict(p->c,p->mu,p->E_in,p->disp,CEL(p->R),p->M2);
    table_setup(p->x1,p->c,(const double (*)[6])p->M2,p->c->heun,MEP_NI);
    p->ninner=p->nouter=0;
    for (int t2=0;t2<p->x1->ntask;t2++){ if (p->x1->tkind[t2]==1) p->outer_t[p->nouter++]=t2; else p->inner_t[p->ninner++]=t2; }
}
static void p_grid1i(void* a, int t){ struct pctx* p=a; task_grid(p->x1,t); }
static void p_pair1i(void* a, int t){ struct pctx* p=a; pairgr_task(p->x1,p->inner_t[t]); }
static void p_grid1o(void* a, int t){ struct pctx* p=a; task_grid(p->x1,MEP_NI+t); }
static void p_pair1o(void* a, int t){ struct pctx* p=a; pairgr_task(p->x1,p->outer_t[t]); }
static int coarse_staged(const struct mep_cfg* c, const double* mu, const double E_in[MEP_NP][6], double E_out[MEP_NP][6], const double (*Mpre)[6]){
    if (!c->heun || !pool_acquire(1)) return 0;
    static __thread struct tctx *x0=NULL, *x1=NULL; if (!x0){ x0=malloc(sizeof *x0); x1=malloc(sizeof *x1); }
    static __thread struct pctx p; double M2[MEP_NP][6], Mg[MEP_NP][6];
    p.c=c; p.E_in=E_in; p.M2=M2; p.mu=mu; p.x0=x0; p.x1=x1;
    table_setup(x0,c,E_in,1,MEP_NP);
    // the number of inner tasks of the end table equals that of the start table (same task list)
    int ni=0; for (int t=0;t<x0->ntask;t++) if (x0->tkind[t]!=1) ni++;
    struct pjob a={5,{p_grid0,p_pair0,p_pred,p_grid1i,p_pair1i},{&p,&p,&p,&p,&p},{c->maxj+1,x0->ntask,1,MEP_NI,ni}};
    if (c->meann){
        struct pjob am={8,{p_grid0,p_pair0,p_mean,p_grid0,p_pair0,p_pred,p_grid1i,p_pair1i},{&p,&p,&p,&p,&p,&p,&p,&p},{c->maxj+1,x0->ntask,1,MEP_NI,x0->ntask,1,MEP_NI,ni}};
        a=am;
    }
    pool_launch(&a);
    if (Mpre) memcpy(Mg,Mpre,sizeof Mg); else outer_wh(c,E_in,Mg);
    pool_join();
    for (int j=MEP_NI;j<MEP_NP;j++) memcpy(M2[j],Mg[j],sizeof(double)*6);
    table_nmotion(x1,MEP_NI);
    struct pjob b={2,{p_grid1o,p_pair1o},{&p,&p},{c->maxj+1-MEP_NI,p.nouter}};
    pool_launch(&b); pool_join(); pool_release();
    double disp1[MEP_NI], R1[MEP_NI][6];
    table_reduce(x1,disp1,R1);
    inner_correct(c,mu,E_in,CEL(p.R),CEL(R1),disp1,M2,E_out);
    return 1;
}
void reb_mep_coarse_par(const struct mep_cfg* c, const double E_in[MEP_NP][6], double E_out[MEP_NP][6], int par);
void reb_mep_coarse_pre(const struct mep_cfg* c, const double E_in[MEP_NP][6], double E_out[MEP_NP][6], int par, const double (*Mpre)[6]);
void reb_mep_coarse(const struct mep_cfg* c, const double E_in[MEP_NP][6], double E_out[MEP_NP][6]){ reb_mep_coarse_par(c,E_in,E_out,0); }
// the giants' part of G alone (their WH run; depends only on the giants' rows of E_in): rows MEP_NI.. of Eo
void reb_mep_coarse_giants(const struct mep_cfg* c, const double E_in[MEP_NP][6], double Eo[MEP_NP][6]){ outer_wh(c,E_in,Eo); }
void reb_mep_coarse_par(const struct mep_cfg* c, const double E_in[MEP_NP][6], double E_out[MEP_NP][6], int par){ reb_mep_coarse_pre(c,E_in,E_out,par,NULL); }
// as reb_mep_coarse_par, with the giants' end rows Mpre (= reb_mep_coarse_giants(E_in)) precomputed; bitwise identical
void reb_mep_coarse_pre(const struct mep_cfg* c, const double E_in[MEP_NP][6], double E_out[MEP_NP][6], int par, const double (*Mpre)[6]){
    const double* m=c->masses;
    double Mint[MEP_NP], mu[MEP_NP]; double acc=0;
    for (int i=0;i<MEP_NP;i++){ acc+=m[i]; Mint[i]=acc; mu[i]=Mint[i]+m[i+1]; }
    if (par && coarse_staged(c,mu,E_in,E_out,Mpre)) return;
    double disp[MEP_NI], R[MEP_NI][6], M2[MEP_NP][6];
    table_par(c,E_in,disp,R,par);
    if (c->meann){ double E1[MEP_NP][6]; mean_state(E_in,disp,E1); table_par(c,(const double (*)[6])E1,disp,R,par); }
    inner_predict(c,mu,E_in,disp,CEL(R),M2);
    if (Mpre){ for (int j=MEP_NI;j<MEP_NP;j++) memcpy(M2[j],Mpre[j],sizeof(double)*6); } else outer_wh(c,E_in,M2);
    double disp1[MEP_NI];
    if (c->heun){
        double R1[MEP_NI][6];
        table_par(c,CEL(M2),disp1,R1,par);
        inner_correct(c,mu,E_in,CEL(R),CEL(R1),disp1,M2,E_out);
    } else {
        table_par(c,CEL(M2),disp1,NULL,par);
        memcpy(E_out,M2,sizeof M2);
        for (int i=0;i<MEP_NI;i++) E_out[i][0]=M2[i][0]+disp1[i];
    }
}
void reb_mep_sim_to_el(struct reb_simulation* r, double E[MEP_NP][6]){
    for (int i=1;i<=MEP_NP;i++){
        struct reb_particle prim=reb_simulation_com_range(r,0,i);
        struct reb_orbit o=reb_orbit_from_particle(r->G,r->particles[i],prim);
        const double si=sin(o.inc/2.);
        E[i-1][0]=log(o.a); E[i-1][1]=o.l; E[i-1][2]=o.e*cos(o.pomega); E[i-1][3]=o.e*sin(o.pomega);
        E[i-1][4]=si*cos(o.Omega); E[i-1][5]=si*sin(o.Omega);
    }
}
struct reb_simulation* reb_mep_el_to_sim(const double* masses, const double E[MEP_NP][6]){
    struct reb_simulation* s=reb_simulation_create();
    reb_simulation_add_fmt(s,"m",masses[0]);
    for (int j=0;j<MEP_NP;j++){
        const double* e=E[j];
        const double ecc=hypot(e[2],e[3]), pom=atan2(e[3],e[2]), si=hypot(e[4],e[5]), Om=atan2(e[5],e[4]);
        reb_simulation_add_fmt(s,"m a e inc Omega pomega l", masses[j+1], exp(e[0]), ecc, 2.*asin(si), Om, pom, e[1]);
    }
    reb_simulation_move_to_com(s);
    return s;
}
#else
typedef int mep_coarse_unavailable;    // avoids an empty translation unit on other platforms
#endif // __linux__
