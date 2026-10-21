// Coarse propagator for parareal around WHFast512 (see pint_coarse.h). Mirrors parareal/mgc.py (config A).
#include "pint_coarse.h"
#include "integrator_whfast512.h"
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#define DAY (2.*M_PI/365.25)
#define NGMAX 64
static const double CGR = 10065.32;
double pint_t_tables=0, pint_t_outer=0, pint_t_total=0;
double pint_prof[4]={0,0,0,0};
static double now(){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec; }

void pint_default_cfg(struct pint_cfg* c, const double* masses, double slice_days){
    memcpy(c->masses, masses, sizeof(double)*(PINT_NP+1));
    c->slice_days = slice_days; c->dt_out_days = 360.; c->ng_ii = 32; c->ng_io = 16; c->maxj = 5; c->kernel = 2;
}
static double wrap(double x){ return x - 2.*M_PI*floor((x+M_PI)/(2.*M_PI)); }

static void el_to_rv(double mu, const double* E, double* r, double* v){
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
static void rv_to_el(double mu, const double* r, const double* v, double* E){
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
static void dEdv(double mu, const double* r, const double* v, double D[6][3], int ncomp){
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
static double complex tw[NGMAX+1][NGMAX/2]; static int tw_ok[NGMAX+1];
static void fft1(double complex* a, int n, int s){
    if (!tw_ok[n]){ for (int j=0;j<n/2;j++) tw[n][j]=cos(-2.*M_PI*j/n)+I*sin(-2.*M_PI*j/n); tw_ok[n]=1; }
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
static void fft2(double complex* a, int n){
    for (int r=0;r<n;r++) fft1(a+r*n,n,1);
    for (int c=0;c<n;c++) fft1(a+c,n,n);
}
static inline int kf(int a, int n){ return a<n/2 ? a : a-n; }   // numpy fftfreq(n, 1/n)
// short-period displacements of ln a for the two real grids packed in Z = fft2(fa + i fb) (unnormalised):
// disp = sum_{k != 0} Im(F_k e^{i k.theta}) / (k.nu), F = fft2(f)/n^2, with F_a(k) = (Z(k) + conj Z(-k))/2,
// F_b(k) = (Z(k) - conj Z(-k))/(2i)
static void sp_eval_pair(const double complex* Z, int n, double ni, double nj, double li, double lj, double* da, double* db){
    double complex ei[NGMAX], ej[NGMAX];
    for (int a=0;a<n;a++){ const int K=kf(a,n); ei[a]=cos(K*li)+I*sin(K*li); ej[a]=cos(K*lj)+I*sin(K*lj); }
    double sa=0, sb=0; const double nrm=1./(n*n);
    for (int a=0;a<n;a++){
        const int Ka=kf(a,n); const int am=(n-a)%n;
        for (int b=0;b<n;b++){
            const double kk=Ka*ni+kf(b,n)*nj;
            if (kk==0.) continue;
            const int bm=(n-b)%n;
            const double complex z=Z[a*n+b], zc=conj(Z[am*n+bm]);
            const double complex Fa=0.5*(z+zc), Fb=-0.5*I*(z-zc);
            const double complex e=ei[a]*ej[b];
            sa+=cimag(Fa*e)/kk; sb+=cimag(Fb*e)/kk;
        }
    }
    *da=sa*nrm; *db=sb*nrm;
}
// first-order ln a short-period displacements of the inner planets (disp[i]) and, if R != NULL, the
// secular (grid-mean) element rates R[i][0..5] of the inner planets, from the elements E
static void table(const struct pint_cfg* c, const double E[PINT_NP][6], double disp[PINT_NI], double R[PINT_NI][6]){
    const double* m=c->masses;
    double Mint[PINT_NP], mu[PINT_NP], n[PINT_NP];
    double acc=0; for (int i=0;i<PINT_NP;i++){ acc+=m[i]; Mint[i]=acc; mu[i]=Mint[i]+m[i+1]; n[i]=sqrt(mu[i]/exp(3.*E[i][0])); }
    const int N1=c->ng_ii, N2=c->ng_io, st=N1/N2, nc = R ? 6 : 1;
    double tp0=now();
    static double ri[PINT_NI][NGMAX][3], vi[PINT_NI][NGMAX][3], Di[PINT_NI][NGMAX][6][3], rj[PINT_NP][NGMAX][3];
    for (int i=0;i<PINT_NI;i++){
        double Eg[6]; memcpy(Eg,E[i],sizeof Eg);
        for (int g=0; g<N1; g++){
            Eg[1]=2.*M_PI*g/N1;
            el_to_rv(mu[i],Eg,ri[i][g],vi[i][g]);
            dEdv(mu[i],ri[i][g],vi[i][g],Di[i][g],nc);
        }
        disp[i]=0; if (R) for (int e=0;e<6;e++) R[i][e]=0;
    }
    for (int j=PINT_NI;j<=c->maxj;j++){
        double Eg[6]; memcpy(Eg,E[j],sizeof Eg); double vdum[3];
        for (int g=0; g<N2; g++){ Eg[1]=2.*M_PI*g/N2; el_to_rv(mu[j],Eg,rj[j][g],vdum); }
    }
    double tp1=now(); pint_prof[0]+=tp1-tp0; double tl=0, tf=0, te=0;
    static double complex fi[NGMAX*NGMAX], fg[NGMAX*NGMAX];
    for (int i=0;i<PINT_NI;i++){
        for (int j=i+1;j<=c->maxj;j++){
            const int inner_j = j<PINT_NI;
            const int Ng = inner_j ? N1 : N2, si = inner_j ? 1 : st;
            const int giant_slot = inner_j ? 0 : (j-PINT_NI)&1;     // giants: two grids packed per FFT
            double complex* buf = inner_j ? fi : fg;
            double Ri[6]={0}, Rj[6]={0}; double tq0=now();
            const double pf=inner_j ? -(Mint[j]+m[j+1])/Mint[j]*m[i+1] : 0.;
            for (int a=0;a<Ng;a++){
                const double* r1=ri[i][a*si]; const double* D1=Di[i][a*si][0];
                for (int b=0;b<Ng;b++){
                    const double* r2 = inner_j ? ri[j][b] : rj[j][b];
                    const double d[3]={r2[0]-r1[0],r2[1]-r1[1],r2[2]-r1[2]};
                    const double d2=d[0]*d[0]+d[1]*d[1]+d[2]*d[2], id3=1./(d2*sqrt(d2));
                    const double rj2=r2[0]*r2[0]+r2[1]*r2[1]+r2[2]*r2[2], irj3=1./(rj2*sqrt(rj2));
                    double ai[3];
                    for (int x=0;x<3;x++) ai[x]=m[j+1]*(d[x]*id3-r2[x]*irj3);
                    const double f=D1[0]*ai[0]+D1[1]*ai[1]+D1[2]*ai[2];
                    if (R) for (int e=0;e<6;e++){ const double* De=Di[i][a*si][e]; Ri[e]+=De[0]*ai[0]+De[1]*ai[1]+De[2]*ai[2]; }
                    if (inner_j){
                        const double rdot=r2[0]*r1[0]+r2[1]*r1[1]+r2[2]*r1[2];
                        double aj[3];
                        for (int x=0;x<3;x++) aj[x]=pf*(d[x]*id3-r2[x]*irj3+r1[x]*irj3-3.*r2[x]*(rdot/rj2*irj3));
                        const double* D2=Di[j][b][0];
                        const double g=D2[0]*aj[0]+D2[1]*aj[1]+D2[2]*aj[2];
                        buf[a*Ng+b]=f+I*g;
                        if (R) for (int e=0;e<6;e++){ const double* De=Di[j][b][e]; Rj[e]+=De[0]*aj[0]+De[1]*aj[1]+De[2]*aj[2]; }
                    } else {
                        if (giant_slot==0) buf[a*Ng+b]=f; else buf[a*Ng+b]+=I*f;
                    }
                }
            }
            if (R) for (int e=0;e<6;e++){ R[i][e]+=Ri[e]/(Ng*Ng); if (inner_j) R[j][e]+=Rj[e]/(Ng*Ng); }
            double tq1=now(); tl+=tq1-tq0;
            if (inner_j){
                fft2(buf,Ng);
                double tq2=now(); tf+=tq2-tq1;
                double da,db; sp_eval_pair(buf,Ng,n[i],n[j],E[i][1],E[j][1],&da,&db);
                disp[i]+=da; disp[j]+=db;
                te+=now()-tq2;
            } else if (giant_slot==1 || j==c->maxj){
                fft2(buf,Ng);
                double tq2=now(); tf+=tq2-tq1;
                const int j0 = giant_slot==1 ? j-1 : j;
                double da,db;
                sp_eval_pair(buf,Ng,n[i],n[j0],E[i][1],E[j0][1],&da,&db); disp[i]+=da;
                if (giant_slot==1){ sp_eval_pair(buf,Ng,n[i],n[j],E[i][1],E[j][1],&da,&db); disp[i]+=db; }
                te+=now()-tq2;
            }
        }
        // GR single-planet term (1-D, N1 points)
        double complex f1[NGMAX]; double Rg[6]={0};
        const double pre=-6.*m[0]*m[0]/(CGR*CGR);
        for (int a=0;a<N1;a++){
            const double* r1=ri[i][a];
            const double rr=sqrt(r1[0]*r1[0]+r1[1]*r1[1]+r1[2]*r1[2]), r4=rr*rr*rr*rr;
            double ag[3]; for (int x=0;x<3;x++) ag[x]=pre*r1[x]/r4;
            double f=0; for (int x=0;x<3;x++) f+=Di[i][a][0][x]*ag[x];
            f1[a]=f;
            if (R) for (int e=0;e<6;e++){ double s=0; for (int x=0;x<3;x++) s+=Di[i][a][e][x]*ag[x]; Rg[e]+=s; }
        }
        if (R) for (int e=0;e<6;e++) R[i][e]+=Rg[e]/N1;
        fft1(f1,N1,1);
        double complex s=0.;
        for (int a=0;a<N1;a++){
            const int K=kf(a,N1); if (K==0) continue;
            s += f1[a]/N1*(cos(K*E[i][1])+I*sin(K*E[i][1]))/(I*K*n[i]);
        }
        disp[i]+=creal(s);
    }
    pint_prof[1]+=tl; pint_prof[2]+=tf; pint_prof[3]+=te;
}
static void outer_wh(const struct pint_cfg* c, const double E[PINT_NP][6], double Eo[PINT_NP][6]){
    const double* m=c->masses;
    double Mi=0; for (int i=0;i<=PINT_NI;i++) Mi+=m[i];
    struct reb_simulation* s=reb_simulation_create();
    reb_simulation_add_fmt(s,"m",Mi);
    for (int j=PINT_NI;j<PINT_NP;j++){
        const double* e=E[j];
        const double ecc=hypot(e[2],e[3]), pom=atan2(e[3],e[2]), si=hypot(e[4],e[5]), Om=atan2(e[5],e[4]);
        reb_simulation_add_fmt(s,"m a e inc Omega pomega l", m[j+1], exp(e[0]), ecc, 2.*asin(si), Om, pom, e[1]);
    }
    reb_simulation_move_to_com(s);
    s->dt=c->dt_out_days*DAY; s->exact_finish_time=0;
    reb_simulation_set_integrator(s,"whfast512");
    struct reb_integrator_whfast512_state* w=s->integrator.state;
    w->gr_potential=1; w->corrector=17; w->kernel=c->kernel;
    const long nst=lround(c->slice_days/c->dt_out_days);
    w->concatenate_steps=nst;
    reb_simulation_steps(s,1);
    for (int j=PINT_NI;j<PINT_NP;j++){
        const int idx=j-PINT_NI+1;
        struct reb_particle prim=reb_simulation_com_range(s,0,idx);
        struct reb_orbit o=reb_orbit_from_particle(s->G,s->particles[idx],prim);
        const double si=sin(o.inc/2.);
        Eo[j][0]=log(o.a); Eo[j][1]=o.l; Eo[j][2]=o.e*cos(o.pomega); Eo[j][3]=o.e*sin(o.pomega);
        Eo[j][4]=si*cos(o.Omega); Eo[j][5]=si*sin(o.Omega);
    }
    reb_simulation_free(s);
}
void pint_coarse(const struct pint_cfg* c, const double E_in[PINT_NP][6], double E_out[PINT_NP][6]){
    const double t0=now();
    const double* m=c->masses;
    double Mint[PINT_NP], mu[PINT_NP]; double acc=0;
    for (int i=0;i<PINT_NP;i++){ acc+=m[i]; Mint[i]=acc; mu[i]=Mint[i]+m[i+1]; }
    double disp[PINT_NI], R[PINT_NI][6], M2[PINT_NP][6];
    table(c,E_in,disp,R);
    const double t1=now();
    const double T=c->slice_days*DAY;
    memcpy(M2,E_in,sizeof M2);
    for (int i=0;i<PINT_NI;i++){
        const double lna=E_in[i][0]-disp[i];
        const double nbar=sqrt(mu[i]/exp(3.*lna));
        M2[i][0]=lna;
        M2[i][1]=wrap(E_in[i][1]+(nbar+R[i][1])*T);
        for (int e=2;e<6;e++) M2[i][e]=E_in[i][e]+R[i][e]*T;
    }
    outer_wh(c,E_in,M2);
    const double t2=now();
    double disp1[PINT_NI];
    table(c,M2,disp1,NULL);
    memcpy(E_out,M2,sizeof M2);
    for (int i=0;i<PINT_NI;i++) E_out[i][0]=M2[i][0]+disp1[i];
    const double t3=now();
    pint_t_tables=(t1-t0)+(t3-t2); pint_t_outer=t2-t1; pint_t_total=t3-t0;
}
void pint_sim_to_el(struct reb_simulation* r, double E[PINT_NP][6]){
    for (int i=1;i<=PINT_NP;i++){
        struct reb_particle prim=reb_simulation_com_range(r,0,i);
        struct reb_orbit o=reb_orbit_from_particle(r->G,r->particles[i],prim);
        const double si=sin(o.inc/2.);
        E[i-1][0]=log(o.a); E[i-1][1]=o.l; E[i-1][2]=o.e*cos(o.pomega); E[i-1][3]=o.e*sin(o.pomega);
        E[i-1][4]=si*cos(o.Omega); E[i-1][5]=si*sin(o.Omega);
    }
}
struct reb_simulation* pint_el_to_sim(const double* masses, const double E[PINT_NP][6]){
    struct reb_simulation* s=reb_simulation_create();
    reb_simulation_add_fmt(s,"m",masses[0]);
    for (int j=0;j<PINT_NP;j++){
        const double* e=E[j];
        const double ecc=hypot(e[2],e[3]), pom=atan2(e[3],e[2]), si=hypot(e[4],e[5]), Om=atan2(e[5],e[4]);
        reb_simulation_add_fmt(s,"m a e inc Omega pomega l", masses[j+1], exp(e[0]), ecc, 2.*asin(si), Om, pom, e[1]);
    }
    reb_simulation_move_to_com(s);
    return s;
}
