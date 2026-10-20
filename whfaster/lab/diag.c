#include "rebound.h"
#include "integrator_whfast512_data.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
extern double kb_diag[8][7];
void reb_whfast512c_kepler_step_kbd(struct simd_data*);
void reb_whfast512c_kepler_step_kbdc(struct simd_data*);
int main(int argc,char**argv){
    int comp=atoi(argv[1]);
    struct simd_data* d=aligned_alloc(64,sizeof(*d)); memset(d,0,sizeof(*d)); d->mask=0xFF;
    double es[8]={0.05,0.15,0.25,0.35,0.45,0.55,0.2,0.4}; double as[8]={0.4,0.4,0.4,0.4,0.4,0.4,0.8,0.8};
    for(int l=0;l<8;l++){ struct reb_particle prim={.m=1}; struct reb_particle p=reb_particle_from_orbit(1.,prim,0.,as[l],es[l],0.3,0.2,0.1,1.0+l);
      d->x.lanes[l]=p.x; d->y.lanes[l]=p.y; d->z.lanes[l]=p.z; d->vx.lanes[l]=p.vx; d->vy.lanes[l]=p.vy; d->vz.lanes[l]=p.vz; d->M.lanes[l]=1; d->dt.lanes[l]=5*2*M_PI/365.25; }
    for (long n=0;n<2000000;n++){ if(comp) reb_whfast512c_kepler_step_kbdc(d); else reb_whfast512c_kepler_step_kbd(d); }
    for(int l=0;l<8;l++){ double N=kb_diag[l][6]; printf("a=%.1f e=%.2f: mean rel diff G1 %+.2e (rms %.1e)  G2 %+.2e (rms %.1e)  G3 %+.2e (rms %.1e)\n", as[l],es[l],
        kb_diag[l][0]/N, sqrt(kb_diag[l][3]/N), kb_diag[l][1]/N, sqrt(kb_diag[l][4]/N), kb_diag[l][2]/N, sqrt(kb_diag[l][5]/N)); }
}
