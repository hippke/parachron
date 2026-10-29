/**
 * @file    mep_coarse.c
 * @brief   Coarse propagator of mean-element parareal (MEP) for WHFast512 planetary systems.
 * @details Inner four planets: numerical first-order mean-element theory (Jacobi pair accelerations on 2-D grids of
 *          the mean longitudes, 2-D FFT, rebuilt from the input state at every call; close pairs on a grid of
 *          (lambda_i, lambda_i - lambda_j) that resolves the synodic harmonics, small divisors at the mean motions of
 *          the mean semi-major axes), Kepler + first-order secular drift with the trapezoidal rule. Outer four
 *          planets: WHFast512 at a large step in a reduced system.
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
#define CDD(x) ((const double (*)[3])(x))
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
    c->gcorr = 17; c->dual = 1; c->gr = gr;
    // synodic grids for the close inner pairs (Venus-Earth in the Solar System) and the mean motions of the mean ln a
    // in the small divisors (research log it. 54): 3.1x -> 3.4x on 16 cores at the default tolerance
    c->ng_th = 16; c->ng_ph = 64; c->alpha_an = 0.7; c->meann = 2;
    // giants' kernel: the C kernels are bitwise identical across instruction sets. The AVX2 kernel runs a system of at
    // most four planets in its four-lane variant, which is faster for the four giants than the eight-lane AVX-512 kernel
    // with half of its lanes empty (Ryzen 9 7950X: giants' run 111 -> 88 us, coarse hop 0.126 -> 0.103 ms, MEP on 16
    // cores 3.4 % faster at the default tolerance, 5.5 % at k_min = 2; results bitwise unchanged).
    c->kernel = 4;
#if (defined(__x86_64__) || defined(__i386__)) && (defined(__GNUC__) || defined(__clang__))
    if (__builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma")) c->kernel = 3;
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
// m-point FFTs (as fft1) of the n columns of the m x n complex array a (row-major): the butterflies of fft1 for every
// column, vectorised across the columns (bitwise identical to fft1 on each column)
MEP_CLONES static void fft_cols(double complex* a, int m, int n){
    pthread_once(&tw_once, tw_init);
    for (int i=0,j=0;i<m;i++){
        if (i<j) for (int c=0;c<n;c++){ const double complex x=a[i*n+c]; a[i*n+c]=a[j*n+c]; a[j*n+c]=x; }
        int bit=m>>1; for (; j&bit; bit>>=1) j^=bit; j^=bit;
    }
    for (int len=2; len<=m; len<<=1){
        const int h=len/2, step=m/len;
        for (int i=0;i<m;i+=len){
            for (int j=0;j<h;j++){
                double* restrict p=(double*)(a+(i+j)*n); double* restrict q=(double*)(a+(i+j+h)*n);
                const double wr=creal(tw[m][j*step]), wi=cimag(tw[m][j*step]);
                for (int c=0;c<n;c++){
                    const double xr=q[2*c], xi=q[2*c+1], vr=xr*wr-xi*wi, vi=xr*wi+xi*wr;
                    const double ur=p[2*c], ui=p[2*c+1];
                    p[2*c]=ur+vr; p[2*c+1]=ui+vi; q[2*c]=ur-vr; q[2*c+1]=ui-vi;
                }
            }
        }
    }
}
// 2-D FFT of an n1 x n2 array (row-major, rows of length n2): the n1 row transforms, then the n2 column transforms
// (the rows as the columns of the transposed array)
MEP_CLONES static void fft2(double complex* a, int n1, int n2){
    static __thread double complex* T=NULL; static __thread int tn=0;
    if (n1*n2>tn){ free(T); tn=n1*n2; T=malloc(sizeof(double complex)*tn); }
    for (int r=0;r<n1;r++) for (int c=0;c<n2;c++) T[c*n1+r]=a[r*n2+c];
    fft_cols(T,n2,n1);
    for (int r=0;r<n1;r++) for (int c=0;c<n2;c++) a[r*n2+c]=T[c*n1+r];
    fft_cols(a,n1,n2);
}
static inline int kf(int a, int n){ return a<n/2 ? a : a-n; }   // numpy fftfreq(n, 1/n)
// short-period displacements of ln a for the two real grids packed in Z = fft2(fa + i fb) (unnormalised, n1 x n2
// points of the angles (theta, phi) with rates (ni, nj) and current values (li, lj)):
// disp = sum_{k != 0} Im(F_k e^{i k.theta}) / (k.nu), F = fft2(f)/(n1 n2), with F_a(k) = (Z(k) + conj Z(-k))/2,
// F_b(k) = (Z(k) - conj Z(-k))/(2i). For a synodic grid the caller passes phi = lambda_i - lambda_j, nj = n_i - n_j.
static __thread double nu_cut_t=0.;
// The sum is done in two passes: the terms' numerators q_a = Im(F_a(k) e^{i k.theta}), q_b = Im(F_b(k) e^{i k.theta})
// (sp_terms; want: 1 q_a, 2 q_b, 3 both), then the divisor sums (sp_sums). The numerators do not depend on the mean
// motions, so the re-evaluation with the mean motions of the mean ln a (meann 2) only redoes the sums. The operations
// are those of the complex expressions Fa = (z + conj z')/2, Fb = -i (z - conj z')/2, Im(Fa e), Im(Fb e), so the sums
// are bitwise unchanged.
MEP_CLONES static void sp_terms(const double complex* Z, int n1, int n2, double li, double lj, int want, double* qa, double* qb){
    double complex ei[NGMAX], ej[NGMAX];
    for (int a=0;a<n1;a++){ const int K=kf(a,n1); ei[a]=cos(K*li)+I*sin(K*li); }
    for (int b=0;b<n2;b++){ const int K=kf(b,n2); ej[b]=cos(K*lj)+I*sin(K*lj); }
    double ejr[NGMAX], eji[NGMAX], cr[NGMAX], ci[NGMAX];
    for (int b=0;b<n2;b++){ ejr[b]=creal(ej[b]); eji[b]=cimag(ej[b]); }
    for (int a=0;a<n1;a++){
        const int am=(n1-a)%n1; const double* z=(const double*)(Z+a*n2);
        for (int b=0;b<n2;b++){ const double complex zc=conj(Z[am*n2+(n2-b)%n2]); cr[b]=creal(zc); ci[b]=cimag(zc); }
        const double ar=creal(ei[a]), ai=cimag(ei[a]);
        double* restrict pa=qa+a*n2; double* restrict pb=qb+a*n2;
        for (int b=0;b<n2;b++){
            const double er=ar*ejr[b]-ai*eji[b], eim=ar*eji[b]+ai*ejr[b];
            const double zr=z[2*b], zi=z[2*b+1];
            if (want&1){ const double far=0.5*(zr+cr[b]), fai=0.5*(zi+ci[b]); pa[b]=far*eim+fai*er; }
            if (want&2){ const double fbr=-(-0.5*(zi-ci[b])), fbi=-0.5*(zr-cr[b]); pb[b]=fbr*eim+fbi*er; }
        }
    }
}
// disp = sum_{k != 0} q_k / (k.nu) (terms with k.nu = 0 or |k.nu| < nu_cut skipped), summed in the order of k
MEP_CLONES static void sp_sums(const double* qa, const double* qb, int n1, int n2, double ni, double nj, int want, double* da, double* db){
    double Kb[NGMAX]; for (int b=0;b<n2;b++) Kb[b]=kf(b,n2);
    const double cut=nu_cut_t;
    double sa=0, sb=0; const double nrm=1./((double)n1*n2);
    for (int a=0;a<n1;a++){
        const int Ka=kf(a,n1);
        double ra[NGMAX], rb[NGMAX], kk[NGMAX];
        const double* restrict pa=qa+a*n2; const double* restrict pb=qb+a*n2;
        for (int b=0;b<n2;b++) kk[b]=Ka*ni+Kb[b]*nj;
        if (want&1) for (int b=0;b<n2;b++) ra[b]=pa[b]/kk[b];
        if (want&2) for (int b=0;b<n2;b++) rb[b]=pb[b]/kk[b];
        for (int b=0;b<n2;b++){
            if (kk[b]==0. || fabs(kk[b])<cut) continue;
            if (want&1) sa+=ra[b];
            if (want&2) sb+=rb[b];
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
    int N1, N2, Nlo, Ng1[MEP_NI];       // N1: largest inner grid, N2: giants' grid; Ng1[p]: grid points of inner planet p
    double ri[MEP_NI][NGMAX][3], vi[MEP_NI][NGMAX][3], Di[MEP_NI][NGMAX][6][3], rj[MEP_NP][NGMAX][3];
    // tasks
    int ntask; int ti[NTASKMAX], tj0[NTASKMAX], tj1[NTASKMAX], tkind[NTASKMAX];   // kind 0 inner pair, 1 giants, 2 GR, 3 grid
    int tNa[NTASKMAX], tNb[NTASKMAX], tsyn[NTASKMAX], zoff[NTASKMAX];   // grid of the task (Na x Nb; GR: Ng1[i] x 1), synodic flag, offset of its transform in Z
    double tda[NTASKMAX], tdb[NTASKMAX], tRi[NTASKMAX][2][6], tRj[NTASKMAX][6], tNg[NTASKMAX];   // tNg: grid points of the task
    double complex* Z; int zcap;        // transforms of all tasks, kept for the re-evaluation with other mean motions (meann 2)
    double* Q; int qcap;                // numerators of the divisor sums (sp_terms), at offsets 2 zoff (q_a, then q_b)
};
// Grid tasks: the inner planets' grids are built in chunks of GCH points (the dual-number dE/dv per point is the
// expensive part, and the stage time on the helpers is that of its longest task), the giants' grids in chunks of GCHO
// points (they are on the critical path, right after the giants' integration).
// Task t < ngrid_inner: planet t / nch, chunk t % nch (empty when beyond the planet's grid), where nch = largest
// inner grid of the configuration / GCH; the giants follow (nchg chunks each). The task count depends only on the
// configuration.
#define GCH 8
#define GCHO 4
static int gmax_cfg(const struct mep_cfg* c){
    int g=c->ng_ii; if (c->ng_hi>g) g=c->ng_hi;
    if (c->ng_th>0){ if (c->ng_th>g) g=c->ng_th; if (c->ng_ph>g) g=c->ng_ph; }
    return g;
}
static int ngrid_inner(const struct mep_cfg* c){ return MEP_NI*(gmax_cfg(c)/GCH); }
static int nchg_cfg(const struct mep_cfg* c){ return c->ng_io>GCHO ? c->ng_io/GCHO : 1; }
static int ngrid_all(const struct mep_cfg* c){ return ngrid_inner(c)+(c->maxj+1-MEP_NI)*nchg_cfg(c); }
MEP_CLONES static void task_grid(struct tctx* x, int t){
    const struct mep_cfg* c=x->c; dual_t=c->dual;
    const int nch=gmax_cfg(c)/GCH, ni=MEP_NI*nch;
    if (t<ni){
        const int p=t/nch, g0=(t%nch)*GCH, ng=x->Ng1[p]; if (g0>=ng) return;
        double Eg[6]; memcpy(Eg,x->E[p],sizeof Eg);
        for (int g=g0; g<g0+GCH; g++){ Eg[1]=2.*M_PI*g/ng; el_to_rv(x->mu[p],Eg,x->ri[p][g],x->vi[p][g]); dEdv(x->mu[p],x->ri[p][g],x->vi[p][g],x->Di[p][g],x->nc); }
    } else {
        const int nchg=nchg_cfg(c), p=MEP_NI+(t-ni)/nchg, gn=x->N2/nchg, g0=((t-ni)%nchg)*gn; if (p>c->maxj) return;
        double Eg[6]; memcpy(Eg,x->E[p],sizeof Eg); double vdum[3];
        for (int g=g0; g<g0+gn; g++){ Eg[1]=2.*M_PI*g/x->N2; el_to_rv(x->mu[p],Eg,x->rj[p][g],vdum); }
    }
}
// short-period ln a of a pair task from its transform Z (x->n: mean motions of the small divisors); terms: 1 for the
// first evaluation of the table (numerators computed), 0 for a re-evaluation with other mean motions
MEP_CLONES static void pair_eval(struct tctx* x, int t, int terms){
    const double (*E)[6]=x->E; const int i=x->ti[t], jlo=x->tj0[t], jhi=x->tj1[t], Na=x->tNa[t], Nb=x->tNb[t];
    const double complex* Z=x->Z+x->zoff[t]; double* qa=x->Q+2*x->zoff[t]; double* qb=qa+Na*Nb;
    double da,db;
    if (x->tkind[t]==0){
        if (terms) sp_terms(Z,Na,Nb,E[i][1],x->tsyn[t] ? E[i][1]-E[jlo][1] : E[jlo][1],3,qa,qb);
        sp_sums(qa,qb,Na,Nb,x->n[i],x->tsyn[t] ? x->n[i]-x->n[jlo] : x->n[jlo],3,&da,&db);
        x->tda[t]=da; x->tdb[t]=db;
    } else {
        if (terms) sp_terms(Z,Na,Nb,E[i][1],E[jlo][1],1,qa,qb);
        sp_sums(qa,qb,Na,Nb,x->n[i],x->n[jlo],1,&da,&db); x->tda[t]=da; x->tdb[t]=0;
        if (jhi>jlo){
            if (terms) sp_terms(Z,Na,Nb,E[i][1],E[jhi][1],2,qa,qb);
            sp_sums(qa,qb,Na,Nb,x->n[i],x->n[jhi],2,&da,&db); x->tdb[t]=db;
        }
    }
}
// Planet j's grid points of a pair task in the order of b, as structures of arrays (so that a row of the pair loop
// vectorises): position, |r|^2, 1/|r|^3, r/|r|^3, 3 r and dE/dv. A synodic grid is stored twice over (row a starts at
// an offset). Everything is computed with the operations of the point-by-point loop, so the tables are bitwise
// unchanged.
struct jgrid { double x[2*NGMAX], y[2*NGMAX], z[2*NGMAX], rr[2*NGMAX], ir3[2*NGMAX], sx[2*NGMAX], sy[2*NGMAX], sz[2*NGMAX],
    tx[2*NGMAX], ty[2*NGMAX], tz[2*NGMAX], D[6][3][2*NGMAX]; };
static inline void jgrid_set(struct jgrid* J, int t, const double* r, const double (*D)[3], int nd){
    const double rr=r[0]*r[0]+r[1]*r[1]+r[2]*r[2], ir3=1./(rr*sqrt(rr));
    J->x[t]=r[0]; J->y[t]=r[1]; J->z[t]=r[2]; J->rr[t]=rr; J->ir3[t]=ir3;
    J->sx[t]=r[0]*ir3; J->sy[t]=r[1]*ir3; J->sz[t]=r[2]*ir3; J->tx[t]=3.*r[0]; J->ty[t]=3.*r[1]; J->tz[t]=3.*r[2];
    for (int e=0;e<nd;e++) for (int q=0;q<3;q++) J->D[e][q][t]=D[e][q];
}
// row a of a pair task (planet i at r1 with dE/dv D1, planet j's points J from offset o): f (and g for an inner j)
// per point into bufd (re, im of the transform's input), the rate contributions summed in the order of b
static inline void pair_row(const struct jgrid* restrict J, int o, int Nb, const double* r1, const double (*D1)[3], double mj,
        double pf, int inner, int rates, int slot, double* restrict bufd, double Ri[6], double Rj[6]){
    double ID3[NGMAX], F[NGMAX], G[NGMAX], CI[6][NGMAX], CJ[6][NGMAX];
    const double* restrict X=J->x+o; const double* restrict Y=J->y+o; const double* restrict Z=J->z+o;
    const double* restrict RR=J->rr+o; const double* restrict IR3=J->ir3+o;
    const double* restrict SX=J->sx+o; const double* restrict SY=J->sy+o; const double* restrict SZ=J->sz+o;
    const double* restrict TX=J->tx+o; const double* restrict TY=J->ty+o; const double* restrict TZ=J->tz+o;
    const double r1x=r1[0], r1y=r1[1], r1z=r1[2];
    // 1/|d|^3 in a loop of its own: sqrt sets errno, which keeps a loop that calls it from vectorising
    for (int b=0;b<Nb;b++){
        const double dx=X[b]-r1x, dy=Y[b]-r1y, dz=Z[b]-r1z;
        const double d2=dx*dx+dy*dy+dz*dz; ID3[b]=1./(d2*sqrt(d2));
    }
    // the accelerations of the row first, then one simple loop per element: vectorised by any compiler (GCC 14 left two
    // variants of the fused loop scalar: serial hop 0.79 -> 0.45 ms with GCC 14, 0.50 -> 0.45 ms with GCC 13), bitwise
    // unchanged
    double AX[NGMAX], AY[NGMAX], AZ[NGMAX], BX[NGMAX], BY[NGMAX], BZ[NGMAX];
    for (int b=0;b<Nb;b++){
        const double dx=X[b]-r1x, dy=Y[b]-r1y, dz=Z[b]-r1z, id3=ID3[b];
        AX[b]=mj*(dx*id3-SX[b]); AY[b]=mj*(dy*id3-SY[b]); AZ[b]=mj*(dz*id3-SZ[b]);
        if (inner){
            const double rdot=X[b]*r1x+Y[b]*r1y+Z[b]*r1z, w=rdot/RR[b]*IR3[b];
            BX[b]=pf*(dx*id3-SX[b]+r1x*IR3[b]-TX[b]*w); BY[b]=pf*(dy*id3-SY[b]+r1y*IR3[b]-TY[b]*w); BZ[b]=pf*(dz*id3-SZ[b]+r1z*IR3[b]-TZ[b]*w);
        }
    }
    for (int b=0;b<Nb;b++) F[b]=D1[0][0]*AX[b]+D1[0][1]*AY[b]+D1[0][2]*AZ[b];
    if (rates) for (int e=0;e<6;e++){ const double d0=D1[e][0], d1=D1[e][1], d2=D1[e][2]; double* restrict ce=CI[e];
        for (int b=0;b<Nb;b++) ce[b]=d0*AX[b]+d1*AY[b]+d2*AZ[b]; }
    if (inner){
        { const double* restrict q0=J->D[0][0]+o; const double* restrict q1=J->D[0][1]+o; const double* restrict q2=J->D[0][2]+o;
          for (int b=0;b<Nb;b++) G[b]=q0[b]*BX[b]+q1[b]*BY[b]+q2[b]*BZ[b]; }
        if (rates) for (int e=0;e<6;e++){ const double* restrict q0=J->D[e][0]+o; const double* restrict q1=J->D[e][1]+o; const double* restrict q2=J->D[e][2]+o;
            double* restrict ce=CJ[e]; for (int b=0;b<Nb;b++) ce[b]=q0[b]*BX[b]+q1[b]*BY[b]+q2[b]*BZ[b]; }
    }

    if (inner) for (int b=0;b<Nb;b++){ bufd[2*b]=F[b]; bufd[2*b+1]=G[b]; }
    else if (slot==0) for (int b=0;b<Nb;b++){ bufd[2*b]=F[b]; bufd[2*b+1]=0.; }
    else for (int b=0;b<Nb;b++) bufd[2*b+1]+=F[b];
    if (rates) for (int b=0;b<Nb;b++){
        for (int e=0;e<6;e++) Ri[e]+=CI[e][b];
        if (inner) for (int e=0;e<6;e++) Rj[e]+=CJ[e][b];
    }
}
MEP_CLONES static void task_pair(struct tctx* x, int t){
    const struct mep_cfg* c=x->c; const double* m=c->masses; const int nc=x->nc;
    double complex* buf=x->Z+x->zoff[t];
    nu_cut_t = c->pcut_yr>0. ? 1./c->pcut_yr : 0.;
    const int i=x->ti[t];
    const int jlo=x->tj0[t], jhi=x->tj1[t];
    const int inner_j = x->tkind[t]==0;
    for (int s=0;s<2;s++) for (int e=0;e<6;e++) x->tRi[t][s][e]=0;
    for (int e=0;e<6;e++) x->tRj[t][e]=0;
    // grid of the pair: Na x Nb points; plain: (lambda_i, lambda_j) = 2 pi (a/Na, b/Nb) with Na = Nb;
    // synodic (inner pair): (lambda_i, phi) = 2 pi (a/Na, b/Nb), lambda_j = lambda_i - phi. Planet i's grid has Ng1[i]
    // points, planet j's Ng1[j] (inner) / N2 (giants); sa, sb are the strides of a and b on them.
    const int Na=x->tNa[t], Nb=x->tNb[t], synodic=x->tsyn[t];
    const int sa=x->Ng1[i]/Na;
    struct jgrid J;
    for (int j=jlo;j<=jhi;j++){
        const int Ngj = inner_j ? x->Ng1[j] : x->N2, sb=Ngj/Nb, saj=Ngj/Na;
        const int slot=j-jlo, nd = inner_j ? (nc>1 ? 6 : 1) : 0;
        double Ri[6]={0}, Rj[6]={0};
        const double pf=inner_j ? -(x->Mint[j]+m[j+1])/x->Mint[j]*m[i+1] : 0.;
        // planet j's point of (a, b): plain b sb; synodic (a saj - b sb) mod Ngj = sb ((a saj/sb - b) mod Nb), i.e.
        // point t = b + o_a of the sequence g(t) = sb ((-t) mod Nb), t < 2 Nb (if sb divides saj; else row by row)
        const int shift = synodic && saj%sb==0;
        if (!synodic || shift){
            for (int tt=0; tt<(shift?2*Nb:Nb); tt++){
                const int g = shift ? sb*((Nb-tt%Nb)%Nb) : tt*sb;
                jgrid_set(&J,tt, inner_j ? x->ri[j][g] : x->rj[j][g], inner_j ? x->Di[j][g] : NULL, nd);
            }
        }
        for (int a=0;a<Na;a++){
            const double* r1=x->ri[i][a*sa]; const double (*D1)[3]=CDD(x->Di[i][a*sa]);
            int o=0;
            if (shift) o=(Nb-(a*saj/sb)%Nb)%Nb;
            else if (synodic) for (int b=0, gj=a*saj; b<Nb; b++, gj = gj>=sb ? gj-sb : gj-sb+Ngj) jgrid_set(&J,b,x->ri[j][gj],x->Di[j][gj],nd);
            double* bufd=(double*)(buf+a*Nb);
            if (inner_j){ if (nc>1) pair_row(&J,o,Nb,r1,D1,m[j+1],pf,1,1,slot,bufd,Ri,Rj); else pair_row(&J,o,Nb,r1,D1,m[j+1],pf,1,0,slot,bufd,Ri,Rj); }
            else { if (nc>1) pair_row(&J,o,Nb,r1,D1,m[j+1],pf,0,1,slot,bufd,Ri,Rj); else pair_row(&J,o,Nb,r1,D1,m[j+1],pf,0,0,slot,bufd,Ri,Rj); }
        }
        for (int e=0;e<6;e++){ x->tRi[t][slot][e]=Ri[e]; if (inner_j) x->tRj[t][e]=Rj[e]; }
    }
    x->tNg[t]=(double)Na*Nb;
    fft2(buf,Na,Nb);
    pair_eval(x,t,1);
}
// GR term of inner planet i: its 1-D transform f1 (in Z) and the divisor sum
MEP_CLONES static void gr_eval(struct tctx* x, int t){
    const int i=x->ti[t], N1=x->tNa[t]; const double complex* f1=x->Z+x->zoff[t];
    double complex s=0.;
    for (int a=0;a<N1;a++){
        const int K=kf(a,N1); if (K==0) continue;
        s += f1[a]/N1*(cos(K*x->E[i][1])+I*sin(K*x->E[i][1]))/(I*K*x->n[i]);
    }
    x->tda[t]=creal(s);
}
MEP_CLONES static void task_gr(struct tctx* x, int t){
    const double* m=x->c->masses; const int i=x->ti[t]; const int N1=x->tNa[t];
    double complex* f1=x->Z+x->zoff[t]; double Rg[6]={0};
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
    gr_eval(x,t);
}
// ---- helper pool (fork-join over task indices); workers spin briefly before sleeping (wake-up latency of sleeping
// threads is ~0.1 ms, comparable to a table pass) ----
#include <stdatomic.h>
// A job is a list of stages; each stage is a set of independent tasks, and a stage starts only when all tasks of the
// previous one are finished. Helpers (and the caller, in pool_join) take tasks from per-stage counters.
static pthread_mutex_t pool_mx=PTHREAD_MUTEX_INITIALIZER, pool_busy=PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t pool_cv=PTHREAD_COND_INITIALIZER;
// A stage with a gate starts only when *gate != 0 (set by the caller while the helpers work on the earlier stages).
#define MAXST 12
struct pjob { int nst; void (*fn[MAXST])(void*,int); void* arg[MAXST]; int ntask[MAXST]; atomic_int* gate[MAXST]; };
static struct pjob job;
static int pool_n=0; static atomic_int pool_gen=0, pool_active=0, st_next[MAXST], st_done[MAXST];
// Guests: threads that would otherwise wait (MEP slice threads waiting for the hop chain or the window's end) take
// tasks of the running job too. job_live is set from pool_launch until the caller's tasks and the helpers are done;
// pool_join then waits until no guest is inside the job (guests_in), so a guest never sees the next job half set up.
static atomic_int job_live=0, guests_in=0;
static void pool_run_tasks(void){
    for (int s=0;s<job.nst;s++){
        if (job.gate[s]) while (!atomic_load_explicit(job.gate[s],memory_order_acquire)) MEP_RELAX();
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
    atomic_store(&job_live,1);
    pthread_mutex_lock(&pool_mx); atomic_fetch_add_explicit(&pool_gen,1,memory_order_release); pthread_cond_broadcast(&pool_cv); pthread_mutex_unlock(&pool_mx);
}
static void pool_join(void){
    pool_run_tasks();
    while (atomic_load_explicit(&pool_active,memory_order_acquire)>0) MEP_RELAX();
    atomic_store(&job_live,0);
    while (atomic_load(&guests_in)>0) MEP_RELAX();
}
// one pass as a guest: 1 if a job was running (its remaining tasks are done when it returns), 0 otherwise
int reb_mep_pool_guest(void){
    if (!atomic_load_explicit(&job_live,memory_order_relaxed)) return 0;
    atomic_fetch_add(&guests_in,1);
    int ran=0;
    if (atomic_load(&job_live)){ pool_run_tasks(); ran=1; }
    atomic_fetch_sub(&guests_in,1);
    return ran;
}
static int pool_acquire(int par){ return par && pool_n>0 && pthread_mutex_trylock(&pool_busy)==0; }
static void pool_release(void){ pthread_mutex_unlock(&pool_busy); }
// run fn(arg, 0..ntask-1); on the pool if par and the pool is free, else serially in the caller
static void run_tasks(void (*fn)(void*,int), void* arg, int ntask, int par){
    if (!pool_acquire(par)){ for (int t=0;t<ntask;t++) fn(arg,t); return; }
    struct pjob j={1,{fn},{arg},{ntask},{0}}; pool_launch(&j); pool_join(); pool_release();
}
static void grid_task(void* a, int t){ task_grid((struct tctx*)a, t); }
static void pairgr_task(void* a, int t){ struct tctx* x=a; if (x->tkind[t]==2) task_gr(x,t); else task_pair(x,t); }
// table setup: masses, mean motions (planets < nn only; the others are set by table_nmotion) and the task list in the
// serial loop order
static void table_setup(struct tctx* x, const struct mep_cfg* c, const double E[MEP_NP][6], int rates, int nn){
    x->c=c; x->E=E; x->nc= rates ? 6 : 1;
    const double* m=c->masses; double acc=0;
    for (int i=0;i<MEP_NP;i++){ acc+=m[i]; x->Mint[i]=acc; x->mu[i]=x->Mint[i]+m[i+1]; if (i<nn) x->n[i]=sqrt(x->mu[i]/exp(3.*E[i][0])); }
    x->N2=c->ng_io; x->Nlo=c->ng_ii;
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
    // grids: per pair task (plain ng_ii, ng_hi for close pairs, synodic ng_th x ng_ph), per inner planet the largest
    // grid any of its pair tasks reads (a synodic task reads planet j's grid at lambda_i - phi: Na | Ng1[j], Nb | Ng1[j])
    const int nhi = c->ng_hi>c->ng_ii ? c->ng_hi : c->ng_ii;
    for (int p=0;p<MEP_NI;p++) x->Ng1[p]=c->ng_ii;
    for (int t=0;t<nt;t++){
        const int i=x->ti[t], j=x->tj0[t]; x->tsyn[t]=0;
        if (x->tkind[t]==0){
            const double al=exp(E[i][0]-E[j][0]);
            const int hi = c->ng_hi>0 && al>=c->alpha_hi;
            const int syn = c->ng_th>0 && al>=c->alpha_an;
            x->tsyn[t]=syn;
            if (syn){ x->tNa[t]=c->ng_th; x->tNb[t]=c->ng_ph; } else x->tNa[t]=x->tNb[t]= hi ? nhi : c->ng_ii;
            if (x->tNa[t]>x->Ng1[i]) x->Ng1[i]=x->tNa[t];
            if (x->tNa[t]>x->Ng1[j]) x->Ng1[j]=x->tNa[t];
            if (x->tNb[t]>x->Ng1[j]) x->Ng1[j]=x->tNb[t];
        } else if (x->tkind[t]==1){ x->tNa[t]=x->tNb[t]=c->ng_io; }
    }
    x->N1=0; for (int p=0;p<MEP_NI;p++) if (x->Ng1[p]>x->N1) x->N1=x->Ng1[p];
    int off=0;
    for (int t=0;t<nt;t++){
        if (x->tkind[t]==2){ x->tNa[t]=x->Ng1[x->ti[t]]; x->tNb[t]=1; }
        x->zoff[t]=off; off+=x->tNa[t]*x->tNb[t];
    }
    if (off>x->zcap){ free(x->Z); x->Z=malloc(off*sizeof(double complex)); x->zcap=off; }
    if (off>x->qcap){ free(x->Q); x->Q=malloc(2*off*sizeof(double)); x->qcap=off; }
}
// mean motions of the inner planets from their mean ln a = ln a - disp (small divisors of the mean-n variant)
static void table_mean_n(struct tctx* x, const double disp[MEP_NI]){
    for (int i=0;i<MEP_NI;i++) x->n[i]=sqrt(x->mu[i]/exp(3.*(x->E[i][0]-disp[i])));
}
static void reeval_task(void* a, int t){ struct tctx* x=a; if (x->tkind[t]==2) gr_eval(x,t); else pair_eval(x,t,0); }
static void table_nmotion(struct tctx* x, int from){ for (int i=from;i<MEP_NP;i++) x->n[i]=sqrt(x->mu[i]/exp(3.*x->E[i][0])); }
// reduction in the serial order
static void table_reduce(struct tctx* x, double disp[MEP_NI], double R[MEP_NI][6]){
    const int nt=x->ntask;
    for (int i=0;i<MEP_NI;i++){ disp[i]=0; if (R) for (int e=0;e<6;e++) R[i][e]=0; }
    for (int t=0;t<nt;t++){
        const int i=x->ti[t];
        if (x->tkind[t]==0){
            const int j=x->tj0[t]; const double Np=x->tNg[t];
            if (R) for (int e=0;e<6;e++){ R[i][e]+=x->tRi[t][0][e]/Np; R[j][e]+=x->tRj[t][e]/Np; }
            disp[i]+=x->tda[t]; disp[j]+=x->tdb[t];
        } else if (x->tkind[t]==1){
            const double Np=x->tNg[t];
            for (int j=x->tj0[t]; j<=x->tj1[t]; j++) if (R) for (int e=0;e<6;e++) R[i][e]+=x->tRi[t][j-x->tj0[t]][e]/Np;
            disp[i]+=x->tda[t]; if (x->tj1[t]>x->tj0[t]) disp[i]+=x->tdb[t];
        } else {
            if (R) for (int e=0;e<6;e++) R[i][e]+=x->tRj[t][e]/x->tNa[t];
            disp[i]+=x->tda[t];
        }
    }
}
// mn = 2: mean-n by re-evaluating the divisor sums with the mean motions of the mean ln a (tables unchanged)
static void table_par(const struct mep_cfg* c, const double E[MEP_NP][6], double disp[MEP_NI], double R[MEP_NI][6], int par, int mn){
    static __thread struct tctx* x=NULL; if (!x) x=calloc(1,sizeof *x);
    table_setup(x,c,E,R!=NULL,MEP_NP);
    run_tasks(grid_task,x,ngrid_all(c),par);
    run_tasks(pairgr_task,x,x->ntask,par);
    table_reduce(x,disp,R);
    if (mn==2){ table_mean_n(x,disp); run_tasks(reeval_task,x,x->ntask,par); table_reduce(x,disp,R); }
}
// centre of mass of p[0..n-1]: the operations of reb_simulation_com_range (positions, velocities, mass)
static struct reb_particle com_xv(const struct reb_particle* p, int n){
    double x=0., y=0., z=0., vx=0., vy=0., vz=0., m=0.;
    for (int i=0;i<n;i++){
        x=x*m+p[i].x*p[i].m; y=y*m+p[i].y*p[i].m; z=z*m+p[i].z*p[i].m;
        vx=vx*m+p[i].vx*p[i].m; vy=vy*m+p[i].vy*p[i].m; vz=vz*m+p[i].vz*p[i].m;
        m+=p[i].m;
        if (m>0.){ x/=m; y/=m; z/=m; vx/=m; vy/=m; vz/=m; }
    }
    struct reb_particle c={0}; c.x=x; c.y=y; c.z=z; c.vx=vx; c.vy=vy; c.vz=vz; c.m=m;
    return c;
}
// particle n of a system p[0..n-1] + planet of mass m with Jacobi elements e (ln a, lambda, k, h, q, p): the operations
// of reb_simulation_add_fmt(s, "m a e inc Omega pomega l", ...) (G = 1, primary = centre of mass of p[0..n-1])
static struct reb_particle orbit_particle(const struct reb_particle* p, int n, double m, const double* e){
    const double ecc=hypot(e[2],e[3]), pom=atan2(e[3],e[2]), si=hypot(e[4],e[5]), Om=atan2(e[5],e[4]);
    const double inc=2.*asin(si), l=e[1];
    const struct reb_particle prim=com_xv(p,n);
    const double omega = cos(inc)>0. ? pom-Om : Om-pom;
    const double M = cos(inc)>0. ? l-Om-omega : Om-omega-l;
    int err=0;
    struct reb_particle q=reb_particle_from_orbit_err(1.,prim,m,exp(e[0]),ecc,inc,Om,omega,reb_M_to_f(ecc,M),&err);
    q.r=0.; q.name=NULL; q.sim=NULL;
    return q;
}
static void outer_wh(const struct mep_cfg* c, const double E[MEP_NP][6], double Eo[MEP_NP][6]){
    const double* m=c->masses;
    double Mi=0; for (int i=0;i<=MEP_NI;i++) Mi+=m[i];
    // the system of add_fmt calls, built without parsing the format
    struct reb_particle ps[MEP_NP-MEP_NI+1]; memset(ps,0,sizeof ps); ps[0].m=Mi;
    for (int j=MEP_NI;j<MEP_NP;j++) ps[j-MEP_NI+1]=orbit_particle(ps,j-MEP_NI+1,m[j+1],E[j]);
    struct reb_simulation* s=reb_simulation_create();
    for (int q=0;q<=MEP_NP-MEP_NI;q++) reb_simulation_add(s,ps[q]);
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
// only the inner planets (their grids, the inner-inner pairs, GR); then (gated by the end of the giants' run, in the
// same job, so that no helper has to be woken again) the giants' grids and the inner-giant pairs.
struct pctx {
    const struct mep_cfg* c; const double (*E_in)[6]; double (*M2)[6]; const double* mu;
    struct tctx *x0, *x1; double disp[MEP_NI], R[MEP_NI][6], E1[MEP_NP][6];
    int ninner, inner_t[NTASKMAX], nouter, outer_t[NTASKMAX];
    atomic_int giants;                  // the giants' end rows are in M2
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
// mean-n by divisors only (meann 2): the start tables' divisor sums again with the mean motions of the mean ln a
static void p_mean2(void* a, int t){ (void)t; struct pctx* p=a; table_reduce(p->x0,p->disp,NULL); table_mean_n(p->x0,p->disp); }
static void p_reeval(void* a, int t){ struct pctx* p=a; reeval_task(p->x0,t); }
static void p_pred(void* a, int t){ (void)t; struct pctx* p=a;
    table_reduce(p->x0,p->disp,p->R);
    inner_predict(p->c,p->mu,p->E_in,p->disp,CEL(p->R),p->M2);
    table_setup(p->x1,p->c,(const double (*)[6])p->M2,p->c->heun,MEP_NI);
    p->ninner=p->nouter=0;
    for (int t2=0;t2<p->x1->ntask;t2++){ if (p->x1->tkind[t2]==1) p->outer_t[p->nouter++]=t2; else p->inner_t[p->ninner++]=t2; }
}
static void p_grid1i(void* a, int t){ struct pctx* p=a; task_grid(p->x1,t); }
static void p_pair1i(void* a, int t){ struct pctx* p=a; pairgr_task(p->x1,p->inner_t[t]); }
static void p_nmo1(void* a, int t){ (void)t; struct pctx* p=a; table_nmotion(p->x1,MEP_NI); }
static void p_grid1o(void* a, int t){ struct pctx* p=a; task_grid(p->x1,ngrid_inner(p->c)+t); }
static void p_pair1o(void* a, int t){ struct pctx* p=a; pairgr_task(p->x1,p->outer_t[t]); }
static int coarse_staged(const struct mep_cfg* c, const double* mu, const double E_in[MEP_NP][6], double E_out[MEP_NP][6], const double (*Mpre)[6]){
    if (!c->heun || !pool_acquire(1)) return 0;
    static __thread struct tctx *x0=NULL, *x1=NULL; if (!x0){ x0=calloc(1,sizeof *x0); x1=calloc(1,sizeof *x1); }
    static __thread struct pctx p; double M2[MEP_NP][6], Mg[MEP_NP][6];
    p.c=c; p.E_in=E_in; p.M2=M2; p.mu=mu; p.x0=x0; p.x1=x1;
    table_setup(x0,c,E_in,1,MEP_NP);
    // the numbers of inner and of inner-giant tasks of the end table equal those of the start table (same task list)
    int ni=0, no=0; for (int t=0;t<x0->ntask;t++){ if (x0->tkind[t]!=1) ni++; else no++; }
    const int ng_all=ngrid_all(c), ng_in=ngrid_inner(c);
    struct pjob a={5,{p_grid0,p_pair0,p_pred,p_grid1i,p_pair1i},{&p,&p,&p,&p,&p},{ng_all,x0->ntask,1,ng_in,ni},{0}};
    if (c->meann==1){
        struct pjob am={8,{p_grid0,p_pair0,p_mean,p_grid0,p_pair0,p_pred,p_grid1i,p_pair1i},{&p,&p,&p,&p,&p,&p,&p,&p},{ng_all,x0->ntask,1,ng_in,x0->ntask,1,ng_in,ni},{0}};
        a=am;
    } else if (c->meann==2){
        struct pjob am={7,{p_grid0,p_pair0,p_mean2,p_reeval,p_pred,p_grid1i,p_pair1i},{&p,&p,&p,&p,&p,&p,&p},{ng_all,x0->ntask,1,x0->ntask,1,ng_in,ni},{0}};
        a=am;
    }
    // giants' part of the end tables, after the gate
    const int s0=a.nst; a.nst+=3;
    a.fn[s0]=p_nmo1; a.fn[s0+1]=p_grid1o; a.fn[s0+2]=p_pair1o;
    a.ntask[s0]=1; a.ntask[s0+1]=ng_all-ng_in; a.ntask[s0+2]=no;
    for (int q=s0;q<a.nst;q++){ a.arg[q]=&p; a.gate[q]=NULL; }
    a.gate[s0]=&p.giants;
    atomic_store(&p.giants,0);
    pool_launch(&a);
    if (Mpre) memcpy(Mg,Mpre,sizeof Mg); else outer_wh(c,E_in,Mg);
    for (int j=MEP_NI;j<MEP_NP;j++) memcpy(M2[j],Mg[j],sizeof(double)*6);    // rows the earlier stages do not read
    atomic_store_explicit(&p.giants,1,memory_order_release);
    pool_join(); pool_release();
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
    table_par(c,E_in,disp,R,par,c->meann==2?2:0);
    if (c->meann==1){ double E1[MEP_NP][6]; mean_state(E_in,disp,E1); table_par(c,(const double (*)[6])E1,disp,R,par,0); }
    inner_predict(c,mu,E_in,disp,CEL(R),M2);
    if (Mpre){ for (int j=MEP_NI;j<MEP_NP;j++) memcpy(M2[j],Mpre[j],sizeof(double)*6); } else outer_wh(c,E_in,M2);
    double disp1[MEP_NI];
    if (c->heun){
        double R1[MEP_NI][6];
        table_par(c,CEL(M2),disp1,R1,par,0);
        inner_correct(c,mu,E_in,CEL(R),CEL(R1),disp1,M2,E_out);
    } else {
        table_par(c,CEL(M2),disp1,NULL,par,0);
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
// the particles of reb_mep_el_to_sim without a simulation (the operations of reb_simulation_add_fmt and
// reb_simulation_move_to_com in the same order: bitwise identical, ~5x faster)
void reb_mep_el_to_xv(const double* masses, const double E[MEP_NP][6], struct reb_particle p[MEP_NP+1]){
    memset(p,0,sizeof(struct reb_particle)*(MEP_NP+1));
    p[0].m=masses[0];
    for (int j=0;j<MEP_NP;j++) p[j+1]=orbit_particle(p,j+1,masses[j+1],E[j]);
    const struct reb_particle c=com_xv(p,MEP_NP+1);
    for (int i=0;i<=MEP_NP;i++){ p[i].x-=c.x; p[i].y-=c.y; p[i].z-=c.z; p[i].vx-=c.vx; p[i].vy-=c.vy; p[i].vz-=c.vz; }
}
// reb_mep_sim_to_el of particles p
void reb_mep_xv_to_el(const struct reb_particle p[MEP_NP+1], double E[MEP_NP][6]){
    for (int i=1;i<=MEP_NP;i++){
        const struct reb_orbit o=reb_orbit_from_particle(1.,p[i],com_xv(p,i));
        const double si=sin(o.inc/2.);
        E[i-1][0]=log(o.a); E[i-1][1]=o.l; E[i-1][2]=o.e*cos(o.pomega); E[i-1][3]=o.e*sin(o.pomega);
        E[i-1][4]=si*cos(o.Omega); E[i-1][5]=si*sin(o.Omega);
    }
}
#else
typedef int mep_coarse_unavailable;    // avoids an empty translation unit on other platforms
#endif // __linux__
