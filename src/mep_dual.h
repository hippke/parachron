// Forward-mode dual numbers for the derivatives of the osculating Jacobi elements (ln a, lambda, k, h, q, p) with
// respect to the three velocity components at fixed position (exact, no finite-difference round-off).
// Shared by mep_coarse.c (cfg.dual) and mep_disp6.c.
#ifndef MEP_DUAL_H
#define MEP_DUAL_H
#include <math.h>
typedef struct { double x, d[3]; } dn;
static inline dn dc(double x){ dn r={x,{0,0,0}}; return r; }
static inline dn dadd(dn a, dn b){ dn r={a.x+b.x,{a.d[0]+b.d[0],a.d[1]+b.d[1],a.d[2]+b.d[2]}}; return r; }
static inline dn dsub(dn a, dn b){ dn r={a.x-b.x,{a.d[0]-b.d[0],a.d[1]-b.d[1],a.d[2]-b.d[2]}}; return r; }
static inline dn dmul(dn a, dn b){ dn r={a.x*b.x,{a.d[0]*b.x+a.x*b.d[0],a.d[1]*b.x+a.x*b.d[1],a.d[2]*b.x+a.x*b.d[2]}}; return r; }
static inline dn dscale(dn a, double s){ dn r={a.x*s,{a.d[0]*s,a.d[1]*s,a.d[2]*s}}; return r; }
static inline dn dfun(dn a, double f, double df){ dn r={f,{df*a.d[0],df*a.d[1],df*a.d[2]}}; return r; }
static inline dn ddiv(dn a, dn b){ const double ib=1./b.x, q=a.x*ib; dn r={q,{(a.d[0]-q*b.d[0])*ib,(a.d[1]-q*b.d[1])*ib,(a.d[2]-q*b.d[2])*ib}}; return r; }
static inline dn dsqrt(dn a){ const double s=sqrt(a.x); return dfun(a,s,0.5/s); }
static inline dn dsin(dn a){ return dfun(a,sin(a.x),cos(a.x)); }
static inline dn dcos(dn a){ return dfun(a,cos(a.x),-sin(a.x)); }
static inline dn dlog(dn a){ return dfun(a,log(a.x),1./a.x); }
static inline dn dacos_clamped(dn a){ if (a.x>=1. || a.x<=-1.) return dc(a.x>=1.?0.:M_PI); return dfun(a,acos(a.x),-1./sqrt(1.-a.x*a.x)); }
static inline dn datan2(dn y, dn x){ const double r2=x.x*x.x+y.x*y.x; dn r={atan2(y.x,x.x),{0,0,0}};
    for (int c=0;c<3;c++) r.d[c]=(x.x*y.d[c]-y.x*x.d[c])/r2;
    return r; }
static void mep_dedv_dual(double mu, const double* r, const double* v, double D[6][3]){
    dn V[3]; for (int x=0;x<3;x++){ V[x]=dc(v[x]); V[x].d[x]=1.; }
    const double rn=sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
    const dn v2=dadd(dadd(dmul(V[0],V[0]),dmul(V[1],V[1])),dmul(V[2],V[2]));
    const dn a=ddiv(dc(1.),dsub(dc(2./rn),dscale(v2,1./mu)));
    const dn hv[3]={dsub(dscale(V[2],r[1]),dscale(V[1],r[2])), dsub(dscale(V[0],r[2]),dscale(V[2],r[0])), dsub(dscale(V[1],r[0]),dscale(V[0],r[1]))};
    const dn hn=dsqrt(dadd(dadd(dmul(hv[0],hv[0]),dmul(hv[1],hv[1])),dmul(hv[2],hv[2])));
    const dn inc=dacos_clamped(ddiv(hv[2],hn)), Om=datan2(hv[0],dscale(hv[1],-1.));
    const dn vxh[3]={dsub(dmul(V[1],hv[2]),dmul(V[2],hv[1])), dsub(dmul(V[2],hv[0]),dmul(V[0],hv[2])), dsub(dmul(V[0],hv[1]),dmul(V[1],hv[0]))};
    dn ev[3]; for (int x=0;x<3;x++) ev[x]=dsub(dscale(vxh[x],1./mu),dc(r[x]/rn));
    const dn e=dsqrt(dadd(dadd(dmul(ev[0],ev[0]),dmul(ev[1],ev[1])),dmul(ev[2],ev[2])));
    const dn nv[2]={dcos(Om),dsin(Om)}, hu[3]={ddiv(hv[0],hn),ddiv(hv[1],hn),ddiv(hv[2],hn)};
    const dn cx[3]={dmul(nv[1],ev[2]), dscale(dmul(nv[0],ev[2]),-1.), dsub(dmul(nv[0],ev[1]),dmul(nv[1],ev[0]))};
    const dn om=datan2(dadd(dadd(dmul(cx[0],hu[0]),dmul(cx[1],hu[1])),dmul(cx[2],hu[2])), dadd(dmul(nv[0],ev[0]),dmul(nv[1],ev[1])));
    const dn pom=dadd(Om,om);
    const dn rv=dadd(dadd(dscale(V[0],r[0]),dscale(V[1],r[1])),dscale(V[2],r[2]));
    const dn sE=ddiv(rv,dsqrt(dscale(a,mu))), cE=dsub(dc(1.),ddiv(dc(rn),a));
    const dn Ea=datan2(sE,cE), M=dsub(Ea,dmul(e,dsin(Ea)));
    const dn si=dsin(dscale(inc,0.5));
    const dn El[6]={dlog(a), dadd(M,pom), dmul(e,dcos(pom)), dmul(e,dsin(pom)), dmul(si,dcos(Om)), dmul(si,dsin(Om))};
    for (int k=0;k<6;k++) for (int c=0;c<3;c++) D[k][c]=El[k].d[c];
}
#endif
