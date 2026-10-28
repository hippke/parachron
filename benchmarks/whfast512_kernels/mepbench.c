// MEP speed-up on this machine: the Solar System (paper setup) over 1600 slices of 60000 steps (1.58 Myr), either
// serially (one WHFast512 call per slice, jacobi_recursive = 1, as the paper's serial reference) or with MEP on
// `cores` slices per window (1600/cores windows).
//   ./mepbench serial <kernel>
//   ./mepbench mep <kernel> <cores> [tol=3e-6]
// Output (one line): wall-clock seconds, mean iterations per window (MEP), kernel used.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "rebound.h"
#include "integrator_whfast512.h"

static const double pos[9][3]={{-0.008816286905115728,-0.0010954664916791675,0.0002143249385447027},{-0.05942272929227954,-0.46308699693348293,-0.032897989948949075},{-0.7276101005375593,0.006575003332463933,0.041795901908847084},{-0.5789530452882667,-0.8361530119313055,0.0002611520181901174},{-1.45202492400084,0.827519404194876,0.052981833432457694},{4.492983939852296,2.0661626247490354,-0.10909246996001629},{8.4974210980544,-4.8620394993693585,-0.2537835862373596},{12.959111916929283,14.760785302864473,-0.1130656917933948},{29.787987348666505,-2.51460654509393,-0.6347108842010732}};
static const double vel[9][3]={{0.00014315746073017681,-0.0004912441820893999,8.127678560998346e-07},{1.2978664284760637,-0.09524541469911743,-0.12677574364801253},{-0.019239782390457125,-1.1813975672919448,-0.01509392594251431},{0.8098712561282222,-0.5682496529341624,2.6169897281383047e-05},{-0.37436417754222295,-0.6365841544564991,-0.004143932260467942},{-0.18818907783656452,0.41919544951404614,0.0024710497024424977},{0.14292308496870448,0.2808676923735748,-0.010574288572728728},{-0.1734971049470612,0.14019515029516152,0.0027683484887051457},{0.014142947617173336,0.18292110872737416,-0.004092845767710294}};
static const double mass[9]={0.9999999999950272,1.6601208254808336e-07,2.447838287784771e-06,3.0404326489511185e-06,3.2271560828978514e-07,0.0009547919099366768,0.0002858856700231729,4.366249613200406e-05,5.151383772628957e-05};
#define SLICES 1600
#define SLICE_STEPS 60000

static double now(void){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec; }

int main(int argc, char** argv){
    if (argc<3){ fprintf(stderr,"usage: mepbench serial <kernel> | mepbench mep <kernel> <cores> [tol]\n"); return 2; }
    const int serial = strcmp(argv[1],"serial")==0;
    const int kernel = atoi(argv[2]);
    struct reb_simulation* r=reb_simulation_create(); r->G=1.; r->dt=6./365.25*2.*3.14159265358979323846; r->exact_finish_time=0;
    reb_simulation_set_integrator(r,"whfast512");
    struct reb_integrator_whfast512_state* w=r->integrator.state; w->gr_potential=1; w->corrector=17; w->kernel=kernel;
    for (int i=0;i<9;i++){ struct reb_particle p={0}; p.m=mass[i]; p.x=pos[i][0]; p.y=pos[i][1]; p.z=pos[i][2]; p.vx=vel[i][0]; p.vy=vel[i][1]; p.vz=vel[i][2]; reb_simulation_add(r,p); }
    reb_simulation_move_to_com(r);
    if (serial){
        w->jacobi_recursive=1; w->concatenate_steps=SLICE_STEPS;
        const double a=now();
        for (int s=0;s<SLICES;s++) reb_simulation_steps(r,1);
        printf("%.3f 0 %u\n", now()-a, w->kernel_used);
        reb_simulation_free(r);
        return 0;
    }
    struct reb_mep_options o; reb_mep_options_default(&o);
    o.cores = argc>3 ? atoi(argv[3]) : 16;
    if (argc>4) o.tol = atof(argv[4]);
    o.slice_steps = SLICE_STEPS; o.checkpoint_every = 0;
    struct reb_mep* m=reb_mep_create(r,&o,NULL);
    if (!m){ fprintf(stderr,"mepbench: MEP setup failed\n"); return 1; }
    const double tmax = (double)SLICES*SLICE_STEPS*r->dt*(1.+1e-12);
    reb_simulation_free(r);
    const double a=now();
    const enum REB_STATUS st=reb_mep_integrate(m,tmax);
    const double wall=now()-a;
    if (st!=REB_STATUS_SUCCESS && st!=REB_STATUS_RUNNING){ fprintf(stderr,"mepbench: status %d %s\n",st,m->message); return 1; }
    double kmean=0.; for (int k=0;k<REB_MEP_KHIST;k++) kmean+=k*m->khist[k];
    printf("%.3f %.3f %u\n", wall, kmean/(m->windows>0?m->windows:1), m->kernel);
    reb_mep_free(m);
    return 0;
}
