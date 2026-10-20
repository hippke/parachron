// Lab driver: time and compare kernel variants on the Solar System (dt, GR on, nsys=1).
// usage: lab N_steps repeats dt_days
#include "rebound.h"
#include "integrator_whfast512.h"
#include "integrator_whfast512_data.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../tools/ss_ics.h"
typedef int (*kfun)(struct simd_data*, uint64_t*, int, volatile sig_atomic_t*, int, int, int, int);
#define X(name) int reb_whfast512c_full_steps_##name(struct simd_data*, uint64_t*, int, volatile sig_atomic_t*, int, int, int, int);
#include "variants.def"
#undef X
#define X(name) {#name, reb_whfast512c_full_steps_##name},
static struct { const char* name; kfun f; } V[] = {
#include "variants.def"
};
#undef X
uint64_t* lab_counts_cnt(void) __attribute__((weak));
uint64_t* lab_hist_cnt(void) __attribute__((weak));
static double now(){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec; }
int main(int argc, char** argv){
    uint64_t N = argc>1 ? atoll(argv[1]) : 10000000;
    int rep = argc>2 ? atoi(argv[2]) : 5;
    double dtd = argc>3 ? atof(argv[3]) : 6.0;
    int nv = sizeof(V)/sizeof(V[0]);
    struct reb_simulation* r = ss_create(dtd, 0);
    reb_simulation_set_integrator(r, "whfast512");
    struct reb_integrator_whfast512_state* w = r->integrator.state;
    w->gr_potential=1; w->corrector=17; w->kernel=2; w->concatenate_steps=1;
    r->integrator.callbacks.step(r, r->integrator.state); // one step, stays unsynchronized (data allocated)
    struct simd_data* d0 = aligned_alloc(64, sizeof(struct simd_data));
    memcpy(d0, w->data, sizeof(struct simd_data));
    struct simd_data* d = aligned_alloc(64, sizeof(struct simd_data));
    static double times[64][4096]; struct simd_data fin[64];
    volatile sig_atomic_t sig=0;
    for (int k=0;k<rep;k++) for (int v=0; v<nv; v++){
        memcpy(d, d0, sizeof(*d)); uint64_t n=N;
        double t0=now(); V[v].f(d,&n,0,&sig,1,1,0,0); double t1=now();
        times[v][k]=(t1-t0)/N*1e9; fin[v]=*d;
    }
    for (int v=0; v<nv; v++){
        double s[4096]; memcpy(s,times[v],sizeof(double)*rep);
        for(int i=0;i<rep;i++)for(int j=i+1;j<rep;j++) if(s[j]<s[i]){double t=s[i];s[i]=s[j];s[j]=t;}
        // max relative position difference vs variant 0
        double md=0; for(int i=0;i<8;i++){ double dx=fin[v].x.lanes[i]-fin[0].x.lanes[i]; double rr=fabs(fin[0].x.lanes[i])+1e-300; if(fabs(dx)/rr>md) md=fabs(dx)/rr; }
        int bit = memcmp(&fin[v].x,&fin[0].x, 6*sizeof(reb_simd8d))==0;
        double rr[4096]; for(int k=0;k<rep;k++) rr[k]=times[v][k]/times[0][k];
        for(int i=0;i<rep;i++)for(int j=i+1;j<rep;j++) if(rr[j]<rr[i]){double t=rr[i];rr[i]=rr[j];rr[j]=t;}
        printf("%-14s paired ratio vs v0: median %.4f  [q25 %.4f, q75 %.4f]   ", V[v].name, rr[rep/2], rr[rep/4], rr[3*rep/4]);
        printf("%-14s min %6.1f med %6.1f ns/step  %.1f  %s maxreldx=%.1e\n", V[v].name, s[0], s[rep/2], s[0], bit?"BITWISE==v0":"differs", md);
    }
    if (lab_counts_cnt){ uint64_t* c=lab_counts_cnt(); printf("counts: skipped-first-Halley %llu, newton iters 1:%llu 2:%llu 3:%llu 4:%llu 5:%llu fallback:%llu\n", (unsigned long long)c[0],(unsigned long long)c[1],(unsigned long long)c[2],(unsigned long long)c[3],(unsigned long long)c[4],(unsigned long long)c[5],(unsigned long long)c[9]); }
    if (lab_hist_cnt){ uint64_t* h=lab_hist_cnt(); printf("first-Newton |dX| histogram per lane, bins: <1e-16, 1e-16..1e-15, ..., >1e-5\n"); for(int l=0;l<8;l++){ printf("lane %d:",l); for(int b=0;b<13;b++) printf(" %8llu",(unsigned long long)h[l*13+b]); printf("\n"); } }
    return 0;
}
