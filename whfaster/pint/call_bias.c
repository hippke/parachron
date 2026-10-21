// per-call energy offset of WHFast512 (iteration 36): same steps in C calls vs 1 call, averaged over M start states
//   ./call_bias kernel corrector gr C M [jacobi_recursive]
#include <math.h>
#include <stdlib.h>
#include "pint_coarse.h"
#include "integrator_whfast512.h"
#include <stdio.h>
#include "../tools/ss_ics.h"
static double en(struct reb_simulation* s){ double E=reb_simulation_energy(s); const double c=10065.32, m0=s->particles[0].m;
  for (unsigned i=1;i<s->N;i++){ const double dx=s->particles[i].x-s->particles[0].x, dy=s->particles[i].y-s->particles[0].y, dz=s->particles[i].z-s->particles[0].z; E+=-3.*m0*m0*s->particles[i].m/(c*c*(dx*dx+dy*dy+dz*dz)); } return E; }
static int JR;
static struct reb_simulation* mk(struct reb_simulation* src, int kern, int corr, int gr){
  struct reb_simulation* s=reb_simulation_create(); for (unsigned i=0;i<src->N;i++) reb_simulation_add(s,src->particles[i]);
  s->dt=6.*DAY; s->exact_finish_time=0; reb_simulation_set_integrator(s,"whfast512");
  struct reb_integrator_whfast512_state* w=s->integrator.state; w->gr_potential=gr; w->corrector=corr; w->kernel=kern; w->jacobi_recursive=JR; return s; }

int main(int argc,char**argv){ JR=argc>6?atoi(argv[6]):0;
  const int kern=atoi(argv[1]), corr=atoi(argv[2]), gr=atoi(argv[3]), C=atoi(argv[4]), M=atoi(argv[5]);
  const int nst=600;   // steps per sample
  struct reb_simulation* r=mk(ss_create(6.,0),kern,corr,gr); ((struct reb_integrator_whfast512_state*)r->integrator.state)->concatenate_steps=997;
  double s1=0,s2=0;
  for (int q=0;q<M;q++){
    reb_simulation_steps(r,1);      // new start state
    struct reb_simulation* a=mk(r,kern,corr,gr); struct reb_simulation* b=mk(r,kern,corr,gr);
    ((struct reb_integrator_whfast512_state*)a->integrator.state)->concatenate_steps=nst; reb_simulation_steps(a,1);
    ((struct reb_integrator_whfast512_state*)b->integrator.state)->concatenate_steps=nst/C; for (int c=0;c<C;c++) reb_simulation_steps(b,1);
    const double d=(en(b)-en(a))/fabs(en(a))/(C-1); s1+=d; s2+=d*d;
    reb_simulation_free(a); reb_simulation_free(b);
  }
  const double m=s1/M, sd=sqrt(s2/M-m*m);
  printf("kernel %d corrector %2d gr %d jacobi_recursive %d: extra energy per call %+.2e +- %.1e\n",kern,corr,gr,JR,m,sd/sqrt(M));
}
