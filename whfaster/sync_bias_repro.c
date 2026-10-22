// Standalone reproduction of the WHFast512 synchronization energy offset (for an upstream report).
// Solar System (Sun + 8 planets), dt = 6 d. For M start states along a trajectory, the same N steps are done in
// one call and in C calls; the extra relative energy change per call is averaged. Uses only the public REBOUND API.
//   cc -O3 -march=native sync_bias_repro.c -I<rebound>/src -L<rebound>/src -lrebound -lm -o sync_bias_repro
#include <stdio.h>
#include <math.h>
#include "rebound.h"
#include "integrator_whfast512.h"
static const double pos[9][3]={{-0.008816286905115728,-0.0010954664916791675,0.0002143249385447027},{-0.05942272929227954,-0.46308699693348293,-0.032897989948949075},{-0.7276101005375593,0.006575003332463933,0.041795901908847084},{-0.5789530452882667,-0.8361530119313055,0.0002611520181901174},{-1.45202492400084,0.827519404194876,0.052981833432457694},{4.492983939852296,2.0661626247490354,-0.10909246996001629},{8.4974210980544,-4.8620394993693585,-0.2537835862373596},{12.959111916929283,14.760785302864473,-0.1130656917933948},{29.787987348666505,-2.51460654509393,-0.6347108842010732}};
static const double vel[9][3]={{0.00014315746073017681,-0.0004912441820893999,8.127678560998346e-07},{1.2978664284760637,-0.09524541469911743,-0.12677574364801253},{-0.019239782390457125,-1.1813975672919448,-0.01509392594251431},{0.8098712561282222,-0.5682496529341624,2.6169897281383047e-05},{-0.37436417754222295,-0.6365841544564991,-0.004143932260467942},{-0.18818907783656452,0.41919544951404614,0.0024710497024424977},{0.14292308496870448,0.2808676923735748,-0.010574288572728728},{-0.1734971049470612,0.14019515029516152,0.0027683484887051457},{0.014142947617173336,0.18292110872737416,-0.004092845767710294}};
static const double mass[9]={0.9999999999950272,1.6601208254808336e-07,2.447838287784771e-06,3.0404326489511185e-06,3.2271560828978514e-07,0.0009547919099366768,0.0002858856700231729,4.366249613200406e-05,5.151383772628957e-05};
static struct reb_simulation* mk(const struct reb_particle* ps){
    struct reb_simulation* s=reb_simulation_create(); s->G=1.;
    for (int i=0;i<9;i++) reb_simulation_add(s,ps[i]);
    s->dt=6.*2.*M_PI/365.25; s->exact_finish_time=0; reb_simulation_set_integrator(s,"whfast512");
    struct reb_integrator_whfast512_state* w=s->integrator.state; w->gr_potential=0; w->corrector=17;
    return s;
}
int main(){
    struct reb_particle ps[9];
    for (int i=0;i<9;i++){ struct reb_particle p={0}; p.m=mass[i]; p.x=pos[i][0]; p.y=pos[i][1]; p.z=pos[i][2]; p.vx=vel[i][0]; p.vy=vel[i][1]; p.vz=vel[i][2]; ps[i]=p; }
    struct reb_simulation* r=mk(ps); reb_simulation_move_to_com(r);
    ((struct reb_integrator_whfast512_state*)r->integrator.state)->concatenate_steps=997;
    const int M=30000, C=10, N=600; double s1=0, s2=0;
    for (int q=0;q<M;q++){
        reb_simulation_steps(r,1);    // next start state
        struct reb_simulation* a=mk(r->particles); struct reb_simulation* b=mk(r->particles);
        ((struct reb_integrator_whfast512_state*)a->integrator.state)->concatenate_steps=N; reb_simulation_steps(a,1);
        ((struct reb_integrator_whfast512_state*)b->integrator.state)->concatenate_steps=N/C; for (int c=0;c<C;c++) reb_simulation_steps(b,1);
        const double ea=reb_simulation_energy(a), d=(reb_simulation_energy(b)-ea)/fabs(ea)/(C-1);
        s1+=d; s2+=d*d; reb_simulation_free(a); reb_simulation_free(b);
    }
    const double m=s1/M, sem=sqrt((s2/M-m*m)/M);
    printf("extra dE/E per WHFast512 call (synchronization): %+.2e +- %.1e (%d start states)\n", m, sem, M);
    return 0;
}
