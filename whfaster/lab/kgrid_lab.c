// Kepler solver grid test (Dagli & Rein 2026, Fig. 2/3/7 style).
// usage: kepler_grid kernel(asm|c) mode(std|highe) out.txt
//  std:   a in [0.1,1] AU (ne x na grid), e in [0,0.9], dt = 5 d, M = 1. Random true anomaly per point.
//  highe: dt/P in 10^[-3.1,-0.1], 1-e in 10^[-8,0], a = 1.
// Per point: relative a-error after 1 Kepler step and after 1000 orbits; a from the double state
// evaluated in __float128 (so only the solver's error is measured).
#include "rebound.h"
#include "integrator_whfast512_data.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <quadmath.h>
void reb_whfast512c_kepler_step_kb(struct simd_data*); void reb_whfast512c_kepler_step_kbcomp(struct simd_data*); void reb_whfast512c_kepler_step_kbre(struct simd_data*); void reb_whfast512c_kepler_step_kbcompre(struct simd_data*); void reb_whfast512c_kepler_step_ko(struct simd_data*); void reb_whfast512c_kepler_step_kbcpcc(struct simd_data*); void reb_whfast512c_kepler_step_kbrpcc(struct simd_data*); void reb_whfast512c_kepler_step_kbrrecc(struct simd_data*); void reb_whfast512c_kepler_step_kprod(struct simd_data*); void reb_whfast512c_kepler_step_kbrre(struct simd_data*); void reb_whfast512c_kepler_step_kbrp0(struct simd_data*); void reb_whfast512c_kepler_step_kbrcc(struct simd_data*); void reb_whfast512c_kepler_step_kbru(struct simd_data*); void reb_whfast512c_kepler_step_kbr(struct simd_data*); void reb_whfast512c_kepler_step_kbx(struct simd_data*); void reb_whfast512c_kepler_step_koub(struct simd_data*); void reb_whfast512c_kepler_step_kbcc(struct simd_data*); void reb_whfast512c_kepler_step_kbrp(struct simd_data*); void reb_whfast512_kepler_step(struct simd_data*);
void reb_whfast512c_kepler_step_avx512(struct simd_data*);
static int use_asm;
static void (*kf)(struct simd_data*); static void kstep(struct simd_data* d){ kf(d); }
static __float128 aq(struct simd_data* d, int l){
    __float128 x=d->x.lanes[l],y=d->y.lanes[l],z=d->z.lanes[l],vx=d->vx.lanes[l],vy=d->vy.lanes[l],vz=d->vz.lanes[l],M=d->M.lanes[l];
    __float128 E = (vx*vx+vy*vy+vz*vz)/2 - M/sqrtq(x*x+y*y+z*z);
    return -M/(2*E);
}
static uint64_t rs=12345; __attribute__((constructor)) static void seed_init(void){ const char* s=getenv("SEED"); if(s) rs=strtoull(s,0,10); } static double ur(){ rs=rs*6364136223846793005ULL+1442695040888963407ULL; return (rs>>11)*(1.0/9007199254740992.0); }
int main(int argc, char** argv){
    if(!strcmp(argv[1],"asm"))kf=reb_whfast512_kepler_step; else if(!strcmp(argv[1],"kb"))kf=reb_whfast512c_kepler_step_kb; else if(!strcmp(argv[1],"kbcomp"))kf=reb_whfast512c_kepler_step_kbcomp; else if(!strcmp(argv[1],"kbre"))kf=reb_whfast512c_kepler_step_kbre; else if(!strcmp(argv[1],"kbcompre"))kf=reb_whfast512c_kepler_step_kbcompre; else if(!strcmp(argv[1],"kbrp"))kf=reb_whfast512c_kepler_step_kbrp; else if(!strcmp(argv[1],"kbcc"))kf=reb_whfast512c_kepler_step_kbcc; else if(!strcmp(argv[1],"koub"))kf=reb_whfast512c_kepler_step_koub; else if(!strcmp(argv[1],"kbx"))kf=reb_whfast512c_kepler_step_kbx; else if(!strcmp(argv[1],"kbr"))kf=reb_whfast512c_kepler_step_kbr; else if(!strcmp(argv[1],"kbru"))kf=reb_whfast512c_kepler_step_kbru; else if(!strcmp(argv[1],"kbrcc"))kf=reb_whfast512c_kepler_step_kbrcc; else if(!strcmp(argv[1],"kbrp0"))kf=reb_whfast512c_kepler_step_kbrp0; else if(!strcmp(argv[1],"kbrre"))kf=reb_whfast512c_kepler_step_kbrre; else if(!strcmp(argv[1],"kprod"))kf=reb_whfast512c_kepler_step_kprod; else if(!strcmp(argv[1],"kbrrecc"))kf=reb_whfast512c_kepler_step_kbrrecc; else if(!strcmp(argv[1],"kbrpcc"))kf=reb_whfast512c_kepler_step_kbrpcc; else if(!strcmp(argv[1],"kbcpcc"))kf=reb_whfast512c_kepler_step_kbcpcc; else kf=reb_whfast512c_kepler_step_ko; int highe = !strcmp(argv[2],"highe"); int phase = !strcmp(argv[2],"phase"); int ssm = !strcmp(argv[2],"ss"); FILE* f=fopen(argv[3],"w");
    const int NA=40, NE=40; const double DAY=2*M_PI/365.25;
    struct simd_data* d = aligned_alloc(64, sizeof(struct simd_data));
    int tot=0;
    static double pts[25600][2];
    if (phase){ for (int i=0;i<16;i++) for (int j=0;j<16;j++) for (int k=0;k<16;k++){ pts[tot][0]=0.1+0.9*(i+0.5)/16; pts[tot][1]=0.9*(j+0.5)/16; tot++; } }
    else for (int i=0;i<NA;i++) for (int j=0;j<NE;j++){
        if (!highe){ pts[i*NE+j][0]=0.1+0.9*(i+0.5)/NA; pts[i*NE+j][1]=0.9*(j+0.5)/NE; }
        else { pts[i*NE+j][0]=pow(10,-3.1+3.0*(i+0.5)/NA); pts[i*NE+j][1]=1-pow(10,-8+8.0*(j+0.5)/NE); }  // [0]=dt/P (highe)
        tot++; }
    if (ssm){
        const double sa[8]={0.387,0.723,1.0,1.524,5.203,9.537,19.19,30.07}, se[8]={0.2056,0.0068,0.0167,0.0934,0.0484,0.0539,0.0473,0.0086};
        for (int rep=0; rep<64; rep++){
            memset(d,0,sizeof(*d)); d->mask=0xFF; __float128 a0[8];
            for (int l=0;l<8;l++){ struct reb_particle prim={.m=1}; struct reb_particle p=reb_particle_from_orbit(1.,prim,0.,sa[l],se[l],0.03,2*M_PI*ur(),2*M_PI*ur(),2*M_PI*ur());
                d->x.lanes[l]=p.x; d->y.lanes[l]=p.y; d->z.lanes[l]=p.z; d->vx.lanes[l]=p.vx; d->vy.lanes[l]=p.vy; d->vz.lanes[l]=p.vz; d->M.lanes[l]=1; d->dt.lanes[l]=6*DAY; }
            for (int l=0;l<8;l++) a0[l]=aq(d,l);
            for (long n=0;n<1000000;n++) kstep(d);
            for (int l=0;l<8;l++) fprintf(f,"%d %d %.6e\n", rep, l, (double)((aq(d,l)-a0[l])/a0[l]));
        }
        fclose(f); return 0;
    }
    for (int b=0;b<tot;b+=8){
        memset(d,0,sizeof(*d)); d->mask=0xFF;
        double a[8], e[8], dt[8], P[8];
        for (int l=0;l<8;l++){
            int k=b+l; double fan=2*M_PI*ur();
            if (!highe){ a[l]=pts[k][0]; e[l]=pts[k][1]; dt[l]=5*DAY; }
            else { a[l]=1; e[l]=pts[k][1]; dt[l]=pts[k][0]*2*M_PI; }
            P[l]=2*M_PI*pow(a[l],1.5);
            struct reb_particle prim={.m=1};
            struct reb_particle p=reb_particle_from_orbit(1., prim, 0., a[l], e[l], 0.3, 0.2, 0.1, fan);
            d->x.lanes[l]=p.x; d->y.lanes[l]=p.y; d->z.lanes[l]=p.z; d->vx.lanes[l]=p.vx; d->vy.lanes[l]=p.vy; d->vz.lanes[l]=p.vz;
            d->M.lanes[l]=1; d->dt.lanes[l]=dt[l];
        }
        // all lanes must share dt for the kernel: in std mode they do; in highe mode run lanes separately
        __float128 a0[8]; for (int l=0;l<8;l++) a0[l]=aq(d,l);
        double err1[8], errN[8];
        if (!highe){
            struct simd_data s=*d; kstep(&s); for (int l=0;l<8;l++) err1[l]=(double)((aq(&s,l)-a0[l])/a0[l]);
            // 1000 orbits: per lane different number of steps -> run max and record each lane at its count
            long nmax=0, nl[8]; for (int l=0;l<8;l++){ nl[l]=lround(1000*P[l]/dt[l]); if (nl[l]>nmax) nmax=nl[l]; }
            s=*d; for (long n=1;n<=nmax;n++){ kstep(&s); for (int l=0;l<8;l++) if (n==nl[l]) errN[l]=(double)((aq(&s,l)-a0[l])/a0[l]); }
        } else {
            for (int l=0;l<8;l++){
                struct simd_data s=*d; for (int m=0;m<8;m++) s.dt.lanes[m]=dt[l];
                kstep(&s); err1[l]=(double)((aq(&s,l)-a0[l])/a0[l]);
                long nl=lround(fmin(1000*P[l]/dt[l], 200000)); s=*d; for (int m=0;m<8;m++) s.dt.lanes[m]=dt[l];
                for (long n=0;n<nl;n++) kstep(&s);
                errN[l]=(double)((aq(&s,l)-a0[l])/a0[l]);
            }
        }
        for (int l=0;l<8;l++) fprintf(f,"%.6e %.10e %.4e %.4e\n", highe?pts[b+l][0]:a[l], highe?1-e[l]:e[l], err1[l], errN[l]);
    }
    fclose(f);
}
