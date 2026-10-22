// Single-core speed of unmodified upstream REBOUND WHFast512 for the paper's benchmark (arXiv 2609.25339):
// Sun + 8 planets (ICs of examples/whfast512_solar_system), dt = 6 d, GR potential, 17th-order corrector,
// concatenate_steps = 1e6. Runs for a wall-clock budget (argv[1] seconds, default 600) in calls of 1e6 steps and
// projects the time for 5 Gyr. Public API only.
//   cc -O3 -march=native bench_upstream.c -I<rebound>/src -L<rebound>/src -lrebound -lm -o bench_upstream
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "rebound.h"
#include "integrator_whfast512.h"
static const double pos[9][3]={{-0.008816286905115728,-0.0010954664916791675,0.0002143249385447027},{-0.05942272929227954,-0.46308699693348293,-0.032897989948949075},{-0.7276101005375593,0.006575003332463933,0.041795901908847084},{-0.5789530452882667,-0.8361530119313055,0.0002611520181901174},{-1.45202492400084,0.827519404194876,0.052981833432457694},{4.492983939852296,2.0661626247490354,-0.10909246996001629},{8.4974210980544,-4.8620394993693585,-0.2537835862373596},{12.959111916929283,14.760785302864473,-0.1130656917933948},{29.787987348666505,-2.51460654509393,-0.6347108842010732}};
static const double vel[9][3]={{0.00014315746073017681,-0.0004912441820893999,8.127678560998346e-07},{1.2978664284760637,-0.09524541469911743,-0.12677574364801253},{-0.019239782390457125,-1.1813975672919448,-0.01509392594251431},{0.8098712561282222,-0.5682496529341624,2.6169897281383047e-05},{-0.37436417754222295,-0.6365841544564991,-0.004143932260467942},{-0.18818907783656452,0.41919544951404614,0.0024710497024424977},{0.14292308496870448,0.2808676923735748,-0.010574288572728728},{-0.1734971049470612,0.14019515029516152,0.0027683484887051457},{0.014142947617173336,0.18292110872737416,-0.004092845767710294}};
static const double mass[9]={0.9999999999950272,1.6601208254808336e-07,2.447838287784771e-06,3.0404326489511185e-06,3.2271560828978514e-07,0.0009547919099366768,0.0002858856700231729,4.366249613200406e-05,5.151383772628957e-05};
static double now(){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec; }
static int cmp(const void* a, const void* b){ const double x=*(const double*)a, y=*(const double*)b; return (x>y)-(x<y); }
int main(int argc, char** argv){
    const double budget = argc>1 ? atof(argv[1]) : 600.;
    const double dt_days = 6.;
    struct reb_simulation* r=reb_simulation_create(); r->G=1.; r->dt=dt_days/365.25*2.*M_PI; r->exact_finish_time=0;
    reb_simulation_set_integrator(r,"whfast512");
    struct reb_integrator_whfast512_state* w=r->integrator.state; w->gr_potential=1; w->corrector=17; w->concatenate_steps=1000000;
    for (int i=0;i<9;i++){ struct reb_particle p={0}; p.m=mass[i]; p.x=pos[i][0]; p.y=pos[i][1]; p.z=pos[i][2]; p.vx=vel[i][0]; p.vy=vel[i][1]; p.vz=vel[i][2]; reb_simulation_add(r,p); }
    reb_simulation_move_to_com(r);
    const double E0=reb_simulation_energy(r);
    static double ns[100000]; int nc=0; const double t0=now();
    while (now()-t0<budget && nc<100000){ const double a=now(); reb_simulation_steps(r,1); ns[nc++]=(now()-a)/1e6*1e9; }
    const double wall=now()-t0, steps=1e6*nc, nsavg=wall/steps*1e9;
    qsort(ns,nc,sizeof(double),cmp);
    const double steps5=5e9*365.25/dt_days;
    printf("upstream WHFast512, dt %.0f d, GR, corrector 17: %d calls x 1e6 steps = %.3g steps (%.1f Myr) in %.1f s\n", dt_days, nc, steps, steps*dt_days/365.25/1e6, wall);
    printf("ns/step: mean %.2f, median %.2f, min %.2f, max %.2f\n", nsavg, ns[nc/2], ns[0], ns[nc-1]);
    printf("projected 5 Gyr (%.3g steps): %.2f h at the mean rate (%.2f h at the median)\n", steps5, steps5*nsavg*1e-9/3600, steps5*ns[nc/2]*1e-9/3600);
    printf("relative energy error at the end: %.2e (status %d)\n", fabs((reb_simulation_energy(r)-E0)/E0), r->status);
    return 0;
}
