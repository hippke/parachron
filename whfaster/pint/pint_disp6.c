// Round 3, Phase 5 (cost study): first-order short-period displacements of all six elements (ln a, lambda with the
// Kepler shear, k, h, q, p) of the inner planets, on the same grids as the config-B ln a transform; optional mean
// motions for the small divisors (n_div != NULL). Compiled into libpint.so; used to measure the C cost of the
// all-component transform of the round-3 coarse step.
#include "pint_coarse.h"
#include <math.h>
#include <complex.h>
#include <string.h>
#include <time.h>
#define NGM 128
#define CGR 10065.32
static double wrapd(double x){ return x - 2.*M_PI*floor((x+M_PI)/(2.*M_PI)); }
static void el_to_rv6(double mu, const double* E, double* r, double* v){
    const double a=exp(E[0]), lam=E[1], k=E[2], h=E[3], q=E[4], p=E[5];
    const double e=hypot(k,h), pom=atan2(h,k), si=hypot(q,p), Om=atan2(p,q), inc=2.*asin(si<1?si:1);
    const double M=lam-pom; double Ea=M;
    for (int it=0;it<30;it++){ const double d=(Ea-e*sin(Ea)-M)/(1-e*cos(Ea)); Ea-=d; if (fabs(d)<1e-15) break; }
    const double cE=cos(Ea), sE=sin(Ea), b=sqrt(1-e*e), n=sqrt(mu/(a*a*a));
    const double xp=a*(cE-e), yp=a*b*sE, vxp=-a*n*sE/(1-e*cE), vyp=a*n*b*cE/(1-e*cE);
    const double om=pom-Om, co=cos(om), so=sin(om), cO=cos(Om), sO=sin(Om), ci=cos(inc), s_i=sin(inc);
    const double P[3]={cO*co-sO*so*ci, sO*co+cO*so*ci, so*s_i}, Q[3]={-cO*so-sO*co*ci, -sO*so+cO*co*ci, co*s_i};
    for (int x=0;x<3;x++){ r[x]=xp*P[x]+yp*Q[x]; v[x]=vxp*P[x]+vyp*Q[x]; }
}
static void rv_to_el6(double mu, const double* r, const double* v, double* E){
    const double rn=sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]), v2=v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
    const double a=1./(2./rn-v2/mu);
    const double hv[3]={r[1]*v[2]-r[2]*v[1], r[2]*v[0]-r[0]*v[2], r[0]*v[1]-r[1]*v[0]};
    const double hn=sqrt(hv[0]*hv[0]+hv[1]*hv[1]+hv[2]*hv[2]);
    const double inc=acos(fmax(-1.,fmin(1.,hv[2]/hn))), Om=atan2(hv[0],-hv[1]);
    const double vxh[3]={v[1]*hv[2]-v[2]*hv[1], v[2]*hv[0]-v[0]*hv[2], v[0]*hv[1]-v[1]*hv[0]};
    double ev[3]; for (int x=0;x<3;x++) ev[x]=vxh[x]/mu-r[x]/rn;
    const double e=sqrt(ev[0]*ev[0]+ev[1]*ev[1]+ev[2]*ev[2]);
    const double nv[3]={cos(Om),sin(Om),0.}, hu[3]={hv[0]/hn,hv[1]/hn,hv[2]/hn};
    const double cx[3]={nv[1]*ev[2]-nv[2]*ev[1], nv[2]*ev[0]-nv[0]*ev[2], nv[0]*ev[1]-nv[1]*ev[0]};
    const double om=atan2(cx[0]*hu[0]+cx[1]*hu[1]+cx[2]*hu[2], nv[0]*ev[0]+nv[1]*ev[1]+nv[2]*ev[2]);
    const double pom=Om+om, rv=r[0]*v[0]+r[1]*v[1]+r[2]*v[2];
    const double sE=rv/sqrt(mu*a), cE=1-rn/a, Ea=atan2(sE,cE), M=Ea-e*sin(Ea);
    const double si=sin(inc/2);
    E[0]=log(a); E[1]=M+pom; E[2]=e*cos(pom); E[3]=e*sin(pom); E[4]=si*cos(Om); E[5]=si*sin(Om);
}
static void dEdv6(double mu, const double* r, const double* v, double D[6][3]){
    const double vs=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])*1e-6;
    for (int c=0;c<3;c++){
        double vp[3]={v[0],v[1],v[2]}, vm[3]={v[0],v[1],v[2]}; vp[c]+=vs; vm[c]-=vs;
        double Ep[6], Em[6]; rv_to_el6(mu,r,vp,Ep); rv_to_el6(mu,r,vm,Em);
        for (int e=0;e<6;e++){ double d=Ep[e]-Em[e]; if (e==1) d=wrapd(d); D[e][c]=d/(2*vs); }
    }
}
#include "pint_dual.h"
#define dEdv6d pint_dedv_dual
int pint_disp6k_dual=0;   // 1: dE/dv by dual numbers in pint_disp6k (else central differences as in pint_disp6)
static double complex twd[NGM+1][NGM/2]; static int twd_ok=0;
static void fftd(double complex* a, int n, int s){
    if (!twd_ok){ for (int m=2;m<=NGM;m<<=1) for (int j=0;j<m/2;j++) twd[m][j]=cexp(-2.*M_PI*I*j/m); twd_ok=1; }
    double complex t[NGM];
    for (int i=0,j=0;i<n;i++){ t[j]=a[i*s]; int bit=n>>1; for (; j&bit; bit>>=1) j^=bit; j^=bit; }
    for (int len=2; len<=n; len<<=1){ const int h=len/2, st=n/len;
        for (int i=0;i<n;i+=len) for (int j=0;j<h;j++){ const double complex u=t[i+j], w=t[i+j+h]*twd[n][j*st]; t[i+j]=u+w; t[i+j+h]=u-w; } }
    for (int i=0;i<n;i++) a[i*s]=t[i];
}
static void fft2d(double complex* a, int n){ for (int r=0;r<n;r++) fftd(a+r*n,n,1); for (int c=0;c<n;c++) fftd(a+c,n,n); }
static inline int kfd(int a, int n){ return a<n/2 ? a : a-n; }
// displacement contributions of the two real grids packed in Z (components: comp of planet a / of planet b)
static void eval2(const double complex* Z, int n, double na, double nb, double la, double lb, double* sa, double* sb,
                  double* ra, double* rb){
    double complex ea[NGM], eb[NGM];
    for (int a=0;a<n;a++){ const int K=kfd(a,n); ea[a]=cexp(I*K*la); eb[a]=cexp(I*K*lb); }
    double s1=0,s2=0,r1=0,r2=0; const double nrm=1./(n*n);
    for (int a=0;a<n;a++){ const int Ka=kfd(a,n), am=(n-a)%n;
        for (int b=0;b<n;b++){
            const double kk=Ka*na+kfd(b,n)*nb; if (kk==0.) continue;
            const int bm=(n-b)%n;
            const double complex z=Z[a*n+b], zc=conj(Z[am*n+bm]);
            const double complex Fa=0.5*(z+zc), Fb=-0.5*I*(z-zc), e=ea[a]*eb[b];
            const double complex pa=Fa*e, pb=Fb*e;
            s1+=cimag(pa)/kk; s2+=cimag(pb)/kk; r1+=creal(pa)/(kk*kk); r2+=creal(pb)/(kk*kk);
        }
    }
    *sa=s1*nrm; *sb=s2*nrm; if (ra){ *ra=r1*nrm; *rb=r2*nrm; }
}
// disp[i][c] (inner planets) for pairs (i<j<=maxj); grids: inner-inner ng_ii (ng_hi for a_i/a_j >= alpha_hi), giants ng_io
void pint_disp6(const struct pint_cfg* c, const double E[PINT_NP][6], const double* n_div, double disp[PINT_NI][6]){
    const double* m=c->masses; double Mint[PINT_NP], mu[PINT_NP], n[PINT_NP]; double acc=0;
    for (int i=0;i<PINT_NP;i++){ acc+=m[i]; Mint[i]=acc; mu[i]=Mint[i]+m[i+1]; n[i]= n_div ? n_div[i] : sqrt(mu[i]/exp(3.*E[i][0])); }
    const int N1=(c->ng_hi>c->ng_ii ? c->ng_hi : c->ng_ii);
    static __thread double ri[PINT_NI][NGM][3], vi[PINT_NI][NGM][3], Di[PINT_NI][NGM][6][3], rj[PINT_NP][NGM][3];
    static __thread double complex buf[6][NGM*NGM];
    for (int i=0;i<PINT_NI;i++){
        double Eg[6]; memcpy(Eg,E[i],sizeof Eg);
        for (int g=0;g<N1;g++){ Eg[1]=2.*M_PI*g/N1; el_to_rv6(mu[i],Eg,ri[i][g],vi[i][g]); dEdv6(mu[i],ri[i][g],vi[i][g],Di[i][g]); }
        for (int e=0;e<6;e++) disp[i][e]=0;
    }
    for (int j=PINT_NI;j<=c->maxj;j++){ double Eg[6]; memcpy(Eg,E[j],sizeof Eg); double vd[3];
        for (int g=0;g<c->ng_io;g++){ Eg[1]=2.*M_PI*g/c->ng_io; el_to_rv6(mu[j],Eg,rj[j][g],vd); } }
    for (int i=0;i<PINT_NI;i++){
        for (int j=i+1;j<=c->maxj;j++){
            const int inner_j=j<PINT_NI;
            const int hi=inner_j && c->ng_hi>0 && exp(E[i][0]-E[j][0])>=c->alpha_hi;
            const int Ng= inner_j ? (hi ? N1 : c->ng_ii) : c->ng_io, si= inner_j ? N1/Ng : 1;
            const int sio = inner_j ? 0 : N1/c->ng_io;
            const double pf=inner_j ? -(Mint[j]+m[j+1])/Mint[j]*m[i+1] : 0.;
            for (int a=0;a<Ng;a++){
                const int ia = inner_j ? a*si : a*sio;
                const double* r1=ri[i][ia];
                for (int b=0;b<Ng;b++){
                    const double* r2= inner_j ? ri[j][b*si] : rj[j][b];
                    const double d[3]={r2[0]-r1[0],r2[1]-r1[1],r2[2]-r1[2]};
                    const double d2=d[0]*d[0]+d[1]*d[1]+d[2]*d[2], id3=1./(d2*sqrt(d2));
                    const double rj2=r2[0]*r2[0]+r2[1]*r2[1]+r2[2]*r2[2], irj3=1./(rj2*sqrt(rj2));
                    double ai[3]; for (int x=0;x<3;x++) ai[x]=m[j+1]*(d[x]*id3-r2[x]*irj3);
                    double aj[3]={0,0,0};
                    if (inner_j){ const double rdot=r2[0]*r1[0]+r2[1]*r1[1]+r2[2]*r1[2];
                        for (int x=0;x<3;x++) aj[x]=pf*(d[x]*id3-r2[x]*irj3+r1[x]*irj3-3.*r2[x]*(rdot/rj2*irj3)); }
                    for (int e=0;e<6;e++){
                        const double* De=Di[i][ia][e]; const double f=De[0]*ai[0]+De[1]*ai[1]+De[2]*ai[2];
                        double g=0; if (inner_j){ const double* Dj=Di[j][b*si][e]; g=Dj[0]*aj[0]+Dj[1]*aj[1]+Dj[2]*aj[2]; }
                        buf[e][a*Ng+b]=f+I*g;
                    }
                }
            }
            double sa[6],sb[6],ra=0,rb=0,ra0=0,rb0=0;
            for (int e=0;e<6;e++){ fft2d(buf[e],Ng);
                eval2(buf[e],Ng,n[i],n[j],E[i][1],E[j][1],&sa[e],&sb[e], e==0?&ra0:&ra, e==0?&rb0:&rb); }
            for (int e=0;e<6;e++){ disp[i][e]+=sa[e]; if (inner_j) disp[j][e]+=sb[e]; }
            disp[i][1]+=1.5*n[i]*ra0; if (inner_j) disp[j][1]+=1.5*n[j]*rb0;
        }
        // GR single-planet term
        double complex f1[6][NGM]; const double pre=-6.*m[0]*m[0]/(CGR*CGR);
        for (int a=0;a<N1;a++){ const double* r1=ri[i][a]; const double rr=sqrt(r1[0]*r1[0]+r1[1]*r1[1]+r1[2]*r1[2]), r4=rr*rr*rr*rr;
            double ag[3]; for (int x=0;x<3;x++) ag[x]=pre*r1[x]/r4;
            for (int e=0;e<6;e++){ double s=0; for (int x=0;x<3;x++) s+=Di[i][a][e][x]*ag[x]; f1[e][a]=s; } }
        double complex sl0=0;
        for (int e=0;e<6;e++){ fftd(f1[e],N1,1); double complex s=0;
            for (int a=0;a<N1;a++){ const int K=kfd(a,N1); if (!K) continue; const double kk=K*n[i];
                const double complex ph=cexp(I*K*E[i][1])*f1[e][a]/N1; s+=ph/(I*kk); if (e==0) sl0+=ph/(kk*kk); }
            disp[i][e]+=creal(s); }
        disp[i][1]+=1.5*n[i]*creal(sl0);
    }
}

// Round 5, step 1: the same displacements with one kernel per pair instead of one FFT + evaluation per element.
// disp = sum_k F_k e^{i k.theta}/(k.nu) = (1/n^2) sum_g f(theta_g) W1(g) with W1 = FFT(e^{i k.theta}/(k.nu)), and the
// shear term sum_k F_k e^{i k.theta}/(k.nu)^2 likewise with W2; W1, W2 do not depend on the element, so per pair two
// FFTs replace six FFTs and six evaluations. Exact identity (agrees with pint_disp6 to round-off).
static void kernels(int n, double na, double nb, double la, double lb, double complex* W1, double complex* W2){
    double complex ea[NGM], eb[NGM];
    for (int a=0;a<n;a++){ const int K=kfd(a,n); ea[a]=cexp(I*K*la); eb[a]=cexp(I*K*lb); }
    for (int a=0;a<n;a++){ const int Ka=kfd(a,n);
        for (int b=0;b<n;b++){ const double kk=Ka*na+kfd(b,n)*nb;
            if (kk==0.){ W1[a*n+b]=0; W2[a*n+b]=0; continue; }
            const double complex e=ea[a]*eb[b]; W1[a*n+b]=e/kk; W2[a*n+b]=e/(kk*kk); } }
    fft2d(W1,n); fft2d(W2,n);
}
double pint_disp6k_t[3];   // timing of the last call: grids, pair fills, kernels + sums
static double nowd(){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec; }
void pint_disp6k(const struct pint_cfg* c, const double E[PINT_NP][6], const double* n_div, double disp[PINT_NI][6]){
    const double t0=nowd(); double tf=0, tk=0;
    const double* m=c->masses; double Mint[PINT_NP], mu[PINT_NP], n[PINT_NP]; double acc=0;
    for (int i=0;i<PINT_NP;i++){ acc+=m[i]; Mint[i]=acc; mu[i]=Mint[i]+m[i+1]; n[i]= n_div ? n_div[i] : sqrt(mu[i]/exp(3.*E[i][0])); }
    const int N1=(c->ng_hi>c->ng_ii ? c->ng_hi : c->ng_ii);
    static __thread double ri[PINT_NI][NGM][3], vi[PINT_NI][NGM][3], Di[PINT_NI][NGM][6][3], rj[PINT_NP][NGM][3];
    static __thread double fa[6][NGM*NGM], fb[6][NGM*NGM];
    static __thread double complex W1[NGM*NGM], W2[NGM*NGM];
    for (int i=0;i<PINT_NI;i++){
        double Eg[6]; memcpy(Eg,E[i],sizeof Eg);
        for (int g=0;g<N1;g++){ Eg[1]=2.*M_PI*g/N1; el_to_rv6(mu[i],Eg,ri[i][g],vi[i][g]); if (pint_disp6k_dual) dEdv6d(mu[i],ri[i][g],vi[i][g],Di[i][g]); else dEdv6(mu[i],ri[i][g],vi[i][g],Di[i][g]); }
        for (int e=0;e<6;e++) disp[i][e]=0;
    }
    for (int j=PINT_NI;j<=c->maxj;j++){ double Eg[6]; memcpy(Eg,E[j],sizeof Eg); double vd[3];
        for (int g=0;g<c->ng_io;g++){ Eg[1]=2.*M_PI*g/c->ng_io; el_to_rv6(mu[j],Eg,rj[j][g],vd); } }
    const double t1=nowd();
    for (int i=0;i<PINT_NI;i++){
        for (int j=i+1;j<=c->maxj;j++){
            const double ta=nowd();
            const int inner_j=j<PINT_NI;
            const int hi=inner_j && c->ng_hi>0 && exp(E[i][0]-E[j][0])>=c->alpha_hi;
            const int Ng= inner_j ? (hi ? N1 : c->ng_ii) : c->ng_io, si= inner_j ? N1/Ng : 1;
            const int sio = inner_j ? 0 : N1/c->ng_io;
            const double pf=inner_j ? -(Mint[j]+m[j+1])/Mint[j]*m[i+1] : 0.;
            for (int a=0;a<Ng;a++){
                const int ia = inner_j ? a*si : a*sio;
                const double* r1=ri[i][ia];
                for (int b=0;b<Ng;b++){
                    const double* r2= inner_j ? ri[j][b*si] : rj[j][b];
                    const double d[3]={r2[0]-r1[0],r2[1]-r1[1],r2[2]-r1[2]};
                    const double d2=d[0]*d[0]+d[1]*d[1]+d[2]*d[2], id3=1./(d2*sqrt(d2));
                    const double rj2=r2[0]*r2[0]+r2[1]*r2[1]+r2[2]*r2[2], irj3=1./(rj2*sqrt(rj2));
                    double ai[3]; for (int x=0;x<3;x++) ai[x]=m[j+1]*(d[x]*id3-r2[x]*irj3);
                    double aj[3]={0,0,0};
                    if (inner_j){ const double rdot=r2[0]*r1[0]+r2[1]*r1[1]+r2[2]*r1[2];
                        for (int x=0;x<3;x++) aj[x]=pf*(d[x]*id3-r2[x]*irj3+r1[x]*irj3-3.*r2[x]*(rdot/rj2*irj3)); }
                    for (int e=0;e<6;e++){
                        const double* De=Di[i][ia][e]; fa[e][a*Ng+b]=De[0]*ai[0]+De[1]*ai[1]+De[2]*ai[2];
                        if (inner_j){ const double* Dj=Di[j][b*si][e]; fb[e][a*Ng+b]=Dj[0]*aj[0]+Dj[1]*aj[1]+Dj[2]*aj[2]; }
                    }
                }
            }
            const double tb=nowd(); tf+=tb-ta;
            kernels(Ng,n[i],n[j],E[i][1],E[j][1],W1,W2);
            const double nrm=1./((double)Ng*Ng*Ng*Ng);    // 1/n^2 of F_k and 1/n^2 of the kernel sum
            double sa[6]={0}, sb[6]={0}, ra0=0, rb0=0;
            for (int g=0;g<Ng*Ng;g++){
                const double w1=cimag(W1[g]);
                for (int e=0;e<6;e++){ sa[e]+=fa[e][g]*w1; if (inner_j) sb[e]+=fb[e][g]*w1; }
                ra0+=fa[0][g]*creal(W2[g]); if (inner_j) rb0+=fb[0][g]*creal(W2[g]);
            }
            for (int e=0;e<6;e++){ disp[i][e]+=sa[e]*nrm*Ng*Ng; if (inner_j) disp[j][e]+=sb[e]*nrm*Ng*Ng; }
            disp[i][1]+=1.5*n[i]*ra0*nrm*Ng*Ng; if (inner_j) disp[j][1]+=1.5*n[j]*rb0*nrm*Ng*Ng;
            tk+=nowd()-tb;
        }
        // GR single-planet term (as in pint_disp6)
        double complex f1[6][NGM]; const double pre=-6.*m[0]*m[0]/(CGR*CGR);
        for (int a=0;a<N1;a++){ const double* r1=ri[i][a]; const double rr=sqrt(r1[0]*r1[0]+r1[1]*r1[1]+r1[2]*r1[2]), r4=rr*rr*rr*rr;
            double ag[3]; for (int x=0;x<3;x++) ag[x]=pre*r1[x]/r4;
            for (int e=0;e<6;e++){ double s=0; for (int x=0;x<3;x++) s+=Di[i][a][e][x]*ag[x]; f1[e][a]=s; } }
        double complex sl0=0;
        for (int e=0;e<6;e++){ fftd(f1[e],N1,1); double complex s=0;
            for (int a=0;a<N1;a++){ const int K=kfd(a,N1); if (!K) continue; const double kk=K*n[i];
                const double complex ph=cexp(I*K*E[i][1])*f1[e][a]/N1; s+=ph/(I*kk); if (e==0) sl0+=ph/(kk*kk); }
            disp[i][e]+=creal(s); }
        disp[i][1]+=1.5*n[i]*creal(sl0);
    }
    pint_disp6k_t[0]=t1-t0; pint_disp6k_t[1]=tf; pint_disp6k_t[2]=tk;
}
