// Throughput of the paired WHFast512 kernel: two independent simulations in one instruction stream
// (reb_integrator_whfast512_steps_pair) against the same two stepped one after the other. Paper setup (Sun + 8 planets,
// dt = 6 d, GR potential, corrector 17); the second simulation starts 1e-12 AU away.
// Output (one line): ns per step per simulation for single and paired (medians of ABBA rounds), their ratio, the
// kernel used, and whether the paired results equal the single ones bitwise.
//   cc -O2 pbench.c -I<src> -L<src> -lrebound -lm      ./pbench <kernel> <rounds>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include "rebound.h"
#include "integrator_whfast512.h"

static const double pos[9][3]={{-0.008816286905115728,-0.0010954664916791675,0.0002143249385447027},{-0.05942272929227954,-0.46308699693348293,-0.032897989948949075},{-0.7276101005375593,0.006575003332463933,0.041795901908847084},{-0.5789530452882667,-0.8361530119313055,0.0002611520181901174},{-1.45202492400084,0.827519404194876,0.052981833432457694},{4.492983939852296,2.0661626247490354,-0.10909246996001629},{8.4974210980544,-4.8620394993693585,-0.2537835862373596},{12.959111916929283,14.760785302864473,-0.1130656917933948},{29.787987348666505,-2.51460654509393,-0.6347108842010732}};
static const double vel[9][3]={{0.00014315746073017681,-0.0004912441820893999,8.127678560998346e-07},{1.2978664284760637,-0.09524541469911743,-0.12677574364801253},{-0.019239782390457125,-1.1813975672919448,-0.01509392594251431},{0.8098712561282222,-0.5682496529341624,2.6169897281383047e-05},{-0.37436417754222295,-0.6365841544564991,-0.004143932260467942},{-0.18818907783656452,0.41919544951404614,0.0024710497024424977},{0.14292308496870448,0.2808676923735748,-0.010574288572728728},{-0.1734971049470612,0.14019515029516152,0.0027683484887051457},{0.014142947617173336,0.18292110872737416,-0.004092845767710294}};
static const double mass[9]={0.9999999999950272,1.6601208254808336e-07,2.447838287784771e-06,3.0404326489511185e-06,3.2271560828978514e-07,0.0009547919099366768,0.0002858856700231729,4.366249613200406e-05,5.151383772628957e-05};

static double now(void){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec; }
static int cmp(const void* a, const void* b){ const double x=*(const double*)a, y=*(const double*)b; return (x>y)-(x<y); }
static struct reb_simulation* sim(int kernel, double dx){
    struct reb_simulation* r=reb_simulation_create(); r->G=1.; r->dt=6./365.25*2.*3.14159265358979323846; r->exact_finish_time=0;
    for (int i=0;i<9;i++){ struct reb_particle p={0}; p.m=mass[i]; p.x=pos[i][0]+(i==3?dx:0); p.y=pos[i][1]; p.z=pos[i][2]; p.vx=vel[i][0]; p.vy=vel[i][1]; p.vz=vel[i][2]; reb_simulation_add(r,p); }
    reb_simulation_move_to_com(r);
    reb_simulation_set_integrator(r,"whfast512");
    struct reb_integrator_whfast512_state* w=r->integrator.state; w->gr_potential=1; w->corrector=17; w->kernel=kernel;
    return r;
}
static void single(struct reb_simulation* r, uint64_t n){ ((struct reb_integrator_whfast512_state*)r->integrator.state)->concatenate_steps=n; reb_simulation_steps(r,1); }
int main(int argc, char** argv){
    const int kernel = argc>1 ? atoi(argv[1]) : 2;
    const int rounds = argc>2 ? atoi(argv[2]) : 5;
    const uint64_t N = 1000000;
    struct reb_simulation *a1=sim(kernel,0), *a2=sim(kernel,1e-12), *b1=sim(kernel,0), *b2=sim(kernel,1e-12);
    single(a1,1000); single(a2,1000); reb_integrator_whfast512_steps_pair(b1,b2,1000);      // warm-up, same steps
    double* ts=malloc(sizeof(double)*rounds); double* tp=malloc(sizeof(double)*rounds); double* q=malloc(sizeof(double)*rounds);
    for (int i=0;i<rounds;i++){                       // ABBA: single, paired, paired, single
        double t0=now(); single(a1,N); single(a2,N); double s1=now()-t0;
        t0=now(); reb_integrator_whfast512_steps_pair(b1,b2,N); double p1=now()-t0;
        t0=now(); reb_integrator_whfast512_steps_pair(b1,b2,N); double p2=now()-t0;
        t0=now(); single(a1,N); single(a2,N); double s2=now()-t0;
        ts[i]=(s1+s2)/(4.*N)*1e9; tp[i]=(p1+p2)/(4.*N)*1e9; q[i]=tp[i]/ts[i];
    }
    int same=1;
    for (int i=0;i<9;i++){ if (memcmp(&a1->particles[i],&b1->particles[i],sizeof(double)*6) || memcmp(&a2->particles[i],&b2->particles[i],sizeof(double)*6)) same=0; }
    qsort(ts,rounds,sizeof(double),cmp); qsort(tp,rounds,sizeof(double),cmp); qsort(q,rounds,sizeof(double),cmp);
    const int used=((struct reb_integrator_whfast512_state*)a1->integrator.state)->kernel_used;
    printf("%.2f %.2f %.4f %d %s\n", ts[rounds/2], tp[rounds/2], q[rounds/2], used, same ? "identical" : "DIFFERENT");
    return 0;
}
