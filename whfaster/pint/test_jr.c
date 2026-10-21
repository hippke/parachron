// jacobi_recursive = 1 vs 0 for N_systems = 1, 2, 4 and the C kernels 3, 4: one call of 10 steps must agree to round-off
#include <math.h>
#include <stdio.h>
#include "rebound.h"
#include "integrator_whfast512.h"
static struct reb_simulation* mk(int nsys, int jr, int kern){
  struct reb_simulation* r=reb_simulation_create();
  const int npl = nsys==1 ? 8 : (nsys==2 ? 4 : 2);
  for (int s=0;s<nsys;s++){
    reb_simulation_add_fmt(r,"m",1.0);
    for (int i=0;i<npl;i++) reb_simulation_add_fmt(r,"m a e inc Omega omega l",1e-5*(1+i+s),1.0+0.4*i+0.1*s,0.05+0.01*i,0.02*i,0.3*i+s,0.7*i,1.1*i+0.5*s);
  }
  r->dt=0.01; r->exact_finish_time=0; reb_simulation_set_integrator(r,"whfast512");
  struct reb_integrator_whfast512_state* w=r->integrator.state; w->N_systems=nsys; w->gr_potential=0; w->corrector=17; w->kernel=kern; w->jacobi_recursive=jr; w->concatenate_steps=10;
  return r;
}
int main(){
  for (int kern=3; kern<=4; kern++) for (int nsys=1;nsys<=4;nsys*=2){
    struct reb_simulation* a=mk(nsys,0,kern); struct reb_simulation* b=mk(nsys,1,kern);
    const double e0=reb_simulation_energy(a);
    for (int q=0;q<1;q++){ reb_simulation_steps(a,1); reb_simulation_steps(b,1); }
    double dx=0; for (unsigned i=0;i<a->N;i++){ double d=fabs(a->particles[i].x-b->particles[i].x)+fabs(a->particles[i].vy-b->particles[i].vy); if (d>dx) dx=d; }
    printf("kernel %d N_systems %d: max |diff| jr=1 vs jr=0 %.2e; dE/E jr0 %+.2e jr1 %+.2e; status %d %d\n",kern,nsys,dx,(reb_simulation_energy(a)-e0)/fabs(e0),(reb_simulation_energy(b)-e0)/fabs(e0),a->status,b->status);
    reb_simulation_free(a); reb_simulation_free(b);
  }
}
