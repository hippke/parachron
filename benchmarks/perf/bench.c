/**
 * REBOUND use-case performance benchmarks.
 *
 * Each benchmark reproduces the hot loop of a compute-heavy research use case
 * (see reports/02_use_cases.md) at a size that reaches a steady rate within
 * a fraction of a second. The harness repeatedly calls a benchmark's `work`
 * function (one "chunk") and reports throughput in benchmark-specific units.
 *
 * Usage:
 *   bench --list
 *   bench --run NAME [--warmup S] [--time S] [--window S]   (JSON on stdout)
 *   bench --check NAME                                      (fixed-work state signature)
 *
 * Only the public REBOUND API and integrator state headers are used, so the
 * same file builds against older and newer source trees.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include "rebound.h"
#include <dlfcn.h>
#include "integrator_whfast.h"
#include "integrator_whfast512.h"
#include "integrator_ias15.h"
#include "integrator_eos.h"
#include "obertas_ics.h"

/////////////////////////////////////////////////////////////////////////
// Helpers

static double now(void){
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + 1e-9*(double)t.tv_nsec;
}

// Small deterministic RNG (xorshift64*) so setups don't depend on libc rand().
static uint64_t rng_state = 0x9E3779B97F4A7C15ULL;
static void rng_seed(uint64_t s){ rng_state = s*0x9E3779B97F4A7C15ULL + 1; }
static double rng_uniform(double a, double b){
    rng_state ^= rng_state >> 12; rng_state ^= rng_state << 25; rng_state ^= rng_state >> 27;
    uint64_t x = rng_state * 2685821657736338717ULL;
    return a + (b-a)*((double)(x >> 11) * (1.0/9007199254740992.0));
}
static double rng_loguniform(double a, double b){ return pow(10., rng_uniform(log10(a), log10(b))); }

static struct reb_simulation* new_sim(void){
    struct reb_simulation* r = reb_simulation_create();
    r->rand_seed = 42;  // REBOUND's own RNG (MEGNO, reb_random_*) must be deterministic
    return r;
}

static int sim_stopped(struct reb_simulation* r){
    // Statuses > 0 mean the integration ended (collision, escape, ...).
    return r->status > 0;
}

// Order-independent-enough signature of the full particle state.
static double state_signature(struct reb_simulation* r){
    double s = 0;
    for (size_t i=0; i<r->N; i++){
        const struct reb_particle p = r->particles[i];
        s += (double)(i+1)*(p.x + 2.*p.y + 3.*p.z + 5.*p.vx + 7.*p.vy + 11.*p.vz + 13.*p.m);
    }
    return s;
}

// Solar System initial conditions (Sun + 8 planets), units: G=1, AU, yr/2pi.
// Same values as examples/whfast512_solar_system.
static const double ss_pos[9][3] = {
    {-0.008816286905115728, -0.0010954664916791675, 0.0002143249385447027},
    {-0.05942272929227954, -0.46308699693348293, -0.032897989948949075},
    {-0.7276101005375593, 0.006575003332463933, 0.041795901908847084},
    {-0.5789530452882667, -0.8361530119313055, 0.0002611520181901174},
    {-1.45202492400084, 0.827519404194876, 0.052981833432457694},
    {4.492983939852296, 2.0661626247490354, -0.10909246996001629},
    {8.4974210980544, -4.8620394993693585, -0.2537835862373596},
    {12.959111916929283, 14.760785302864473, -0.1130656917933948},
    {29.787987348666505, -2.51460654509393, -0.6347108842010732}};
static const double ss_vel[9][3] = {
    {0.00014315746073017681, -0.0004912441820893999, 8.127678560998346e-07},
    {1.2978664284760637, -0.09524541469911743, -0.12677574364801253},
    {-0.019239782390457125, -1.1813975672919448, -0.01509392594251431},
    {0.8098712561282222, -0.5682496529341624, 2.6169897281383047e-05},
    {-0.37436417754222295, -0.6365841544564991, -0.004143932260467942},
    {-0.18818907783656452, 0.41919544951404614, 0.0024710497024424977},
    {0.14292308496870448, 0.2808676923735748, -0.010574288572728728},
    {-0.1734971049470612, 0.14019515029516152, 0.0027683484887051457},
    {0.014142947617173336, 0.18292110872737416, -0.004092845767710294}};
static const double ss_mass[9] = {
    0.9999999999950272, 1.6601208254808336e-07, 2.447838287784771e-06, 3.0404326489511185e-06,
    3.2271560828978514e-07, 0.0009547919099366768, 0.0002858856700231729,
    4.366249613200406e-05, 5.151383772628957e-05};

static void add_solar_system(struct reb_simulation* r){
    for (int i=0; i<9; i++){
        struct reb_particle p = {0};
        p.m = ss_mass[i];
        p.x = ss_pos[i][0]; p.y = ss_pos[i][1]; p.z = ss_pos[i][2];
        p.vx = ss_vel[i][0]; p.vy = ss_vel[i][1]; p.vz = ss_vel[i][2];
        reb_simulation_add(r, p);
    }
    reb_simulation_move_to_com(r);
}

// GR potential (Nobili & Roxburgh 1986), as REBOUNDx gr_potential / Brown & Rein 2020.
static void gr_potential_force(struct reb_simulation* r){
    const double C2 = 10065.32*10065.32;  // speed of light^2 in these units
    struct reb_particle* const ps = r->particles;
    const struct reb_particle src = ps[0];
    const double prefac1 = 6.*(r->G*src.m)*(r->G*src.m)/C2;
    for (size_t i=1; i<r->N; i++){
        const double dx = ps[i].x - src.x, dy = ps[i].y - src.y, dz = ps[i].z - src.z;
        const double r2 = dx*dx + dy*dy + dz*dz;
        const double prefac = prefac1/(r2*r2);
        ps[i].ax -= prefac*dx; ps[i].ay -= prefac*dy; ps[i].az -= prefac*dz;
        ps[0].ax += ps[i].m/src.m*prefac*dx; ps[0].ay += ps[i].m/src.m*prefac*dy; ps[0].az += ps[i].m/src.m*prefac*dz;
    }
}

/////////////////////////////////////////////////////////////////////////
// Generic "one long simulation" benchmark context

struct simctx {
    struct reb_simulation* r;
    struct reb_simulation* (*make)(int variant);
    int variant;           // which initial condition (cycled on restart)
    size_t chunk_steps;    // steps per work() call
    double work_per_time;  // if >0: work = simulated time * this (adaptive integrators)
    double signature_acc;  // accumulated over restarts
    double step_mult;      // real time steps per REBOUND step (WHFast512 concatenates steps)
};

static void* simctx_create(struct reb_simulation* (*make)(int), size_t chunk_steps, double work_per_time){
    struct simctx* c = calloc(1, sizeof(*c));
    c->make = make; c->chunk_steps = chunk_steps; c->work_per_time = work_per_time;
    c->step_mult = 1.;
    c->r = make(0);
    return c;
}

static double simctx_work(void* p){
    struct simctx* c = p;
    double t0 = c->r->t;
    reb_simulation_steps(c->r, c->chunk_steps);
    double w = c->work_per_time>0. ? (c->r->t - t0)*c->work_per_time : (double)c->chunk_steps*c->step_mult;
    if (sim_stopped(c->r)){
        // System went unstable / collided: start the next initial condition.
        c->signature_acc += state_signature(c->r);
        reb_simulation_free(c->r);
        c->variant++;
        c->r = c->make(c->variant);
    }
    return w;
}
static double simctx_signature(void* p){ struct simctx* c = p; return c->signature_acc + state_signature(c->r); }
static void simctx_destroy(void* p){ struct simctx* c = p; reb_simulation_free(c->r); free(c); }

/////////////////////////////////////////////////////////////////////////
// 1) Solar System, plain WH + GR (Abbot et al. 2021/2023 ensembles)

static struct reb_simulation* make_ss_wh_gr(int variant){
    struct reb_simulation* r = new_sim();
    reb_simulation_set_integrator(r, "whfast");
    struct reb_integrator_whfast_state* wh = r->integrator.state;
    wh->safe_mode = 0;
    r->dt = 3.16/365.25*2.*M_PI;  // 3.16 days
    r->additional_forces = gr_potential_force;
    r->force_is_velocity_dependent = 0;
    r->exact_finish_time = 0;
    add_solar_system(r);
    r->particles[1].x += 1e-12*variant;
    return r;
}
static void* create_ss_wh_gr(void){ return simctx_create(make_ss_wh_gr, 20000, 0); }

// 2) Solar System, WHCKL (17th-order corrector + lazy kernel) + GR (Brown & Rein 2020, 2022)
static struct reb_simulation* make_ss_whckl_gr(int variant){
    struct reb_simulation* r = new_sim();
    reb_simulation_set_integrator(r, "whfast");
    struct reb_integrator_whfast_state* wh = r->integrator.state;
    wh->safe_mode = 0;
    wh->corrector = 17;
    wh->kernel = REB_INTEGRATOR_WHFAST_KERNEL_LAZY;
    r->dt = 8.062/365.25*2.*M_PI;
    r->additional_forces = gr_potential_force;
    r->force_is_velocity_dependent = 0;
    r->exact_finish_time = 0;
    add_solar_system(r);
    r->particles[1].x += 1e-12*variant;
    return r;
}
// Large chunks: every reb_simulation_steps() call synchronizes and applies the corrector.
static void* create_ss_whckl_gr(void){ return simctx_create(make_ss_whckl_gr, 20000, 0); }

// 3) Solar System, WHFast512 with built-in GR (Javaheri, Rein & Tamayo 2023)
static int avx512_ok = -1;
static int have_avx512(void){
    if (avx512_ok<0){
#if defined(__x86_64__) || defined(__amd64__)
        __builtin_cpu_init();
        avx512_ok = __builtin_cpu_supports("avx512f") ? 1 : 0;
#else
        avx512_ok = 0;
#endif
    }
    return avx512_ok;
}
static struct reb_simulation* make_ss_whfast512(int variant){
    struct reb_simulation* r = new_sim();
    reb_simulation_set_integrator(r, "whfast512");
    struct reb_integrator_whfast512_state* w = r->integrator.state;
    w->gr_potential = 1;
    w->concatenate_steps = 1000;  // one REBOUND step = 1000 time steps
    r->dt = 5.0/365.25*2.*M_PI;
    r->exact_finish_time = 0;
    add_solar_system(r);
    r->particles[1].x += 1e-12*variant;
    return r;
}
static void* create_ss_whfast512(void){
    struct simctx* c = simctx_create(make_ss_whfast512, 100, 0);
    c->step_mult = 1000.;
    return c;
}

// 3b) WHFast512 portable C kernels (integrator.kernel, added Oct 2026). The kernel is set through
// the integrator's field descriptors (by name), so bench.c still builds against older trees; there
// the benchmarks count as unavailable.
static int set_integrator_uint(struct reb_simulation* r, const char* name, unsigned int value){
    const struct reb_binarydata_field_descriptor* fd = r->integrator.callbacks.field_descriptor_list;
    if (!fd) return 0;
    for (; fd->name[0]; fd++){
        if (strcmp(fd->name, name)==0){
            *(unsigned int*)((char*)r->integrator.state + fd->offset) = value;
            return 1;
        }
    }
    return 0;
}
static unsigned int w512c_kernel = 3; // set by the create functions below
// type 0: Solar System (8 planets); 1: three compact planets (one system: four-lane path); 2: 4 systems x 2 planets
static struct reb_simulation* make_w512c_type(int type, int variant){
    struct reb_simulation* r = new_sim();
    reb_simulation_set_integrator(r, "whfast512");
    struct reb_integrator_whfast512_state* w = r->integrator.state;
    w->concatenate_steps = 1000;
    set_integrator_uint(r, "kernel", w512c_kernel);
    r->exact_finish_time = 0;
    if (type==0){
        w->gr_potential = 1;
        r->dt = 5.0/365.25*2.*M_PI;
        add_solar_system(r);
    }else{
        const int nsys = type==1 ? 1 : 4, np = type==1 ? 3 : 2;
        w->N_systems = nsys;
        r->dt = 0.035*2.*M_PI;
        for (int s=0; s<nsys; s++){
            reb_simulation_add_fmt(r, "m", 1.0);
            struct reb_particle star = r->particles[r->N-1];
            for (int k=0; k<np; k++){
                reb_simulation_add_fmt(r, "m a e inc l omega primary", 3e-6, pow(1.12, k)*(1.+0.01*s), 0.01, 0.01*k, 1.3*k+0.2*s, 0.5*k, star);
            }
        }
        if (nsys==1) reb_simulation_move_to_com(r);
    }
    r->particles[1].x += 1e-12*variant;
    return r;
}
static struct reb_simulation* make_w512c_ss(int variant){ return make_w512c_type(0, variant); }
static struct reb_simulation* make_w512c_small(int variant){ return make_w512c_type(1, variant); }
static struct reb_simulation* make_w512c_4x2(int variant){ return make_w512c_type(2, variant); }
static void* w512c_create(struct reb_simulation* (*make)(int), unsigned int kernel){
    w512c_kernel = kernel;
    struct simctx* c = simctx_create(make, 100, 0);
    c->step_mult = 1000.;
    return c;
}
static void* create_w512c_ss_avx2(void){ return w512c_create(make_w512c_ss, 3); }
static void* create_w512c_ss_avx512(void){ return w512c_create(make_w512c_ss, 2); }
static void* create_w512c_small_avx2(void){ return w512c_create(make_w512c_small, 3); }
static void* create_w512c_4x2_avx2(void){ return w512c_create(make_w512c_4x2, 3); }
// Available if the tree has the kernel setting and the CPU runs the kernel (one test step).
static int w512c_kernel_ok(unsigned int kernel){
    struct reb_simulation* r = reb_simulation_create();
    reb_simulation_set_integrator(r, "whfast512");
    int ok = set_integrator_uint(r, "kernel", kernel);
    if (ok){
        struct reb_integrator_whfast512_state* w = r->integrator.state;
        w->concatenate_steps = 1;
        r->exact_finish_time = 0;
        r->dt = 0.01;
        reb_simulation_add_fmt(r, "m", 1.0);
        reb_simulation_add_fmt(r, "m a", 1e-3, 1.0);
        reb_simulation_steps(r, 1);
        ok = (r->status <= 0);
    }
    reb_simulation_free(r);
    return ok;
}
static int have_w512c_avx2(void){ static int ok = -1; if (ok<0) ok = w512c_kernel_ok(3); return ok; }
static int have_w512c_avx512(void){ static int ok = -1; if (ok<0) ok = w512c_kernel_ok(2); return ok; }

// 3c) Two simulations per call (reb_integrator_whfast512_steps_pair, Oct 2026), looked up at run
// time. Trees without it step the two simulations one after the other (same work and results),
// so an A/B comparison measures paired vs separate stepping.
typedef void (*pair_fn_t)(struct reb_simulation*, struct reb_simulation*, uint64_t);
struct pairctx { struct reb_simulation* a; struct reb_simulation* b; struct reb_simulation* (*make)(int); };
static void* pair_create(struct reb_simulation* (*make)(int), unsigned int kernel){
    struct pairctx* c = calloc(1, sizeof(*c));
    w512c_kernel = kernel;
    c->make = make; c->a = make(0); c->b = make(1);
    return c;
}
static double pair_work(void* p){
    struct pairctx* c = p;
    static pair_fn_t fn = NULL; static int looked = 0;
    if (!looked){ fn = (pair_fn_t)dlsym(RTLD_DEFAULT, "reb_integrator_whfast512_steps_pair"); looked = 1; }
    const uint64_t N = 100000;
    if (fn){
        fn(c->a, c->b, N);
    }else{
        struct reb_simulation* rs[2] = {c->a, c->b};
        for (int i=0; i<2; i++){
            struct reb_integrator_whfast512_state* w = rs[i]->integrator.state;
            w->concatenate_steps = N;
            reb_simulation_steps(rs[i], 1);
        }
    }
    return 2.*N;   // simulation-steps
}
static double pair_signature(void* p){ struct pairctx* c = p; return state_signature(c->a) + state_signature(c->b); }
static void pair_destroy(void* p){ struct pairctx* c = p; reb_simulation_free(c->a); reb_simulation_free(c->b); free(c); }
static void* create_pair_ss_avx512(void){ return pair_create(make_w512c_ss, 2); }
static void* create_pair_3p_avx2(void){ return pair_create(make_w512c_small, 3); }

// 3d) Output after every step: reb_simulation_steps(r, 1) in a loop (each call synchronizes),
// WHFast512 with the default kernel. Measures the per-call cost (item 53).
static struct reb_simulation* make_w512_out1(int variant){
    struct reb_simulation* r = new_sim();
    reb_simulation_set_integrator(r, "whfast512");
    struct reb_integrator_whfast512_state* w = r->integrator.state;
    w->gr_potential = 1;
    w->concatenate_steps = 1;
    r->dt = 5.0/365.25*2.*M_PI;
    r->exact_finish_time = 0;
    add_solar_system(r);
    r->particles[1].x += 1e-12*variant;
    return r;
}
static void* create_w512_out1(void){ return simctx_create(make_w512_out1, 1, 0); }
static double work_w512_out1(void* p){
    struct simctx* c = p;
    for (int i=0; i<2000; i++) reb_simulation_steps(c->r, 1);
    return 2000.;
}

// 4) Five Earth-mass planets, WHFast dt=0.05 yr (Obertas et al. 2017)
static struct reb_simulation* make_obertas(int variant){
    const int k = variant % OBERTAS_N_SYSTEMS;
    struct reb_simulation* r = new_sim();
    r->G = 4.*M_PI*M_PI;
    reb_simulation_set_integrator(r, "whfast");
    r->dt = 0.05;
    r->exact_finish_time = 0;
    struct reb_particle star = {0};
    star.m = 1.;
    reb_simulation_add(r, star);
    for (int i=0; i<5; i++){
        reb_simulation_add_fmt(r, "primary m a f", star, obertas_m_planet, obertas_a[k][i], obertas_f[k][i]);
    }
    reb_simulation_move_to_com(r);
    return r;
}
static void* create_obertas(void){ return simctx_create(make_obertas, 20000, 0); }

// 5) Random compact three-planet systems with line collision search (SPOCK / Tamayo 2016, 2020)
//    Distributions follow spock/run_integrations/runfunctions.py:run_random, but spacing is
//    restricted to 15-25 mutual Hill radii so most systems survive the measurement window.
static struct reb_simulation* make_spock(int variant){
    rng_seed(1000 + (uint64_t)variant);
    struct reb_simulation* r = new_sim();
    r->G = 4.*M_PI*M_PI;
    reb_simulation_set_integrator(r, "whfast");
    struct reb_integrator_whfast_state* wh = r->integrator.state;
    wh->safe_mode = 0;
    r->collision = REB_COLLISION_LINE;  // default resolve: halt on collision
    r->exact_finish_time = 0;
    double M1 = rng_loguniform(1e-7, 1e-4), M2 = rng_loguniform(1e-7, 1e-4), M3 = rng_loguniform(1e-7, 1e-4);
    double a1 = 1.;
    double h12 = a1*pow((M1+M2)/3., 1./3.);
    double a2 = a1 + rng_uniform(15., 25.)*h12;
    double h23 = a2*pow((M2+M3)/3., 1./3.);
    double a3 = a2 + rng_uniform(15., 25.)*h23;
    double minhill = fmin(h12, h23);
    double e1 = rng_loguniform(M2/pow((a2-a1)/a1, 2), 0.3*(a2-a1)/a1);
    double e2 = rng_loguniform(fmax(M1/pow((a2-a1)/a2, 2), M3/pow((a3-a2)/a2, 2)), 0.3*fmin((a2-a1)/a2, (a3-a2)/a2));
    double e3 = rng_loguniform(M2/pow((a3-a2)/a3, 2), 0.3*(a3-a2)/a3);
    struct reb_particle star = {0};
    star.m = 1.;
    reb_simulation_add(r, star);
    double Ms[3] = {M1, M2, M3}, as[3] = {a1, a2, a3}, es[3] = {e1, e2, e3};
    for (int i=0; i<3; i++){
        reb_simulation_add_fmt(r, "primary m a e inc pomega Omega f r", star, Ms[i], as[i], es[i],
                rng_loguniform(1e-3, 1e-1), rng_uniform(0, 2*M_PI), rng_uniform(0, 2*M_PI), rng_uniform(0, 2*M_PI), minhill);
    }
    reb_simulation_move_to_com(r);
    r->dt = 2.*sqrt(3.)/100.;  // ~3.5% of the inner period
    return r;
}
static void* create_spock(void){ return simctx_create(make_spock, 20000, 0); }

// 6) IAS15, five Saturn-mass planets in the HL Tau gaps (Rein & Tamayo 2017 test case)
//    Work unit: inner-planet orbits (adaptive step size, so steps/s would hide step-size changes).
static struct reb_simulation* make_hltau(int variant){
    struct reb_simulation* r = new_sim();
    r->G = 4.*M_PI*M_PI;  // AU, Msun, yr
    reb_simulation_set_integrator(r, "ias15");
    rng_seed(7 + (uint64_t)variant);
    struct reb_particle star = {0};
    star.m = 1.;
    reb_simulation_add(r, star);
    const double gaps[5] = {13.2, 32.3, 42.2, 68.8, 74.7};
    for (int i=0; i<5; i++){
        reb_simulation_add_fmt(r, "primary m a inc l", star, 3e-4, gaps[i], rng_uniform(0, 1e-3), rng_uniform(0, 2*M_PI));
    }
    reb_simulation_move_to_com(r);
    r->dt = 0.1;
    return r;
}
static void* create_hltau(void){ return simctx_create(make_hltau, 2000, 1./pow(13.2, 1.5)); }

// 7) IAS15 with 60 planetary embryos and merging collisions (giant-impact phase studies)
static struct reb_simulation* make_embryos(int variant){
    struct reb_simulation* r = new_sim();
    r->G = 4.*M_PI*M_PI;
    reb_simulation_set_integrator(r, "ias15");
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    rng_seed(11 + (uint64_t)variant);
    struct reb_particle star = {0};
    star.m = 1.; star.r = 0.005;
    reb_simulation_add(r, star);
    for (int i=0; i<60; i++){
        double a = 0.5 + 1.0*i/60.;
        reb_simulation_add_fmt(r, "primary m r a e inc l pomega", star, 3.7e-8, 2e-5, a,
                rng_uniform(0, 0.02), rng_uniform(0, 0.01), rng_uniform(0, 2*M_PI), rng_uniform(0, 2*M_PI));
    }
    reb_simulation_move_to_com(r);
    r->dt = 1e-3;
    return r;
}
static void* create_embryos(void){ return simctx_create(make_embryos, 200, 1.); }

/////////////////////////////////////////////////////////////////////////
// 8) MEGNO stability map, Jupiter + Saturn-like grid (ipython_examples/Megno.ipynb;
//    Hadden & Lithwick 2018 and exoplanet discovery papers). Work unit: grid points.

struct gridctx { int i; double acc; };
#define MEGNO_NGRID 40
static void* create_megno(void){ return calloc(1, sizeof(struct gridctx)); }
static double work_megno(void* p){
    struct gridctx* c = p;
    int k = c->i % (MEGNO_NGRID*MEGNO_NGRID);
    double a = 7. + 3.*(k % MEGNO_NGRID)/(MEGNO_NGRID-1.);
    double e = 0.5*(k / MEGNO_NGRID)/(MEGNO_NGRID-1.);
    struct reb_simulation* r = new_sim();
    reb_simulation_set_integrator(r, "whfast");
    r->dt = 5.;
    reb_simulation_add_fmt(r, "m", 1.);
    reb_simulation_add_fmt(r, "m a M omega e", 0.000954, 5.204, 0.600, 0.257, 0.048);
    reb_simulation_add_fmt(r, "m a M omega e", 0.000285, a, 0.871, 1.616, e);
    reb_simulation_move_to_com(r);
    reb_simulation_init_megno(r);
    r->exit_max_distance = 20.;
    r->exact_finish_time = 0;
    reb_simulation_integrate(r, 5e2*2.*M_PI);
    double megno = sim_stopped(r) ? 10. : reb_simulation_megno(r);
    c->acc += megno;
    reb_simulation_free(r);
    c->i++;
    return 1.;
}
static double sig_grid(void* p){ return ((struct gridctx*)p)->acc; }
static void destroy_free(void* p){ free(p); }

// 9) Millions of short posterior integrations (Wang et al. 2018, HR 8799). Work unit: simulations.
static void* create_posterior(void){ return calloc(1, sizeof(struct gridctx)); }
static double work_posterior(void* p){
    struct gridctx* c = p;
    rng_seed(5000 + (uint64_t)c->i);
    struct reb_simulation* r = new_sim();
    r->G = 4.*M_PI*M_PI;
    reb_simulation_set_integrator(r, "whfast");
    struct reb_particle star = {0};
    star.m = 1.52;
    reb_simulation_add(r, star);
    const double a[4] = {16.4, 26.7, 41.4, 71.6}, m[4] = {9.6e-3, 9.6e-3, 9.6e-3, 6.9e-3};
    for (int i=0; i<4; i++){
        reb_simulation_add_fmt(r, "primary m a e inc l pomega", star, m[i]*rng_uniform(0.8, 1.2), a[i]*rng_uniform(0.97, 1.03),
                rng_uniform(0, 0.15), rng_uniform(0, 0.05), rng_uniform(0, 2*M_PI), rng_uniform(0, 2*M_PI));
    }
    reb_simulation_move_to_com(r);
    r->dt = 0.01*sqrt(pow(a[0], 3)/star.m);  // 1% of the inner period (yr, G=4pi^2)
    r->exit_min_distance = 0.5;
    reb_simulation_steps(r, 300);
    c->acc += state_signature(r);
    reb_simulation_free(r);
    c->i++;
    return 1.;
}

/////////////////////////////////////////////////////////////////////////
// 10) Saturn's rings viscous overstability: SEI + tree collisions (Rein & Latter 2013;
//     examples/overstability, same box: 200x5x20 root boxes, ~520 particles).

static double cor_half(const struct reb_simulation* r, double v){ (void)r; (void)v; return 0.5; }
static struct reb_simulation* make_rings(int variant){
    struct reb_simulation* r = new_sim();
    r->OMEGA = 1.;
    r->OMEGAZ = 3.6;
    r->dt = 2e-3*2.*M_PI;
    const double particle_r = 1., tau = 1.64;
    r->coefficient_of_restitution = cor_half;
    reb_simulation_set_integrator(r, "sei");
    r->collision = REB_COLLISION_TREE;
    r->collision_resolve = reb_collision_resolve_hardsphere;
    r->gravity = REB_GRAVITY_NONE;
    r->boundary = REB_BOUNDARY_SHEAR;
    r->root_size = 1;
    r->N_root_x = 200; r->N_root_y = 5; r->N_root_z = 20;
    r->N_ghost_x = 1; r->N_ghost_y = 1; r->N_ghost_z = 0;
    const double bx = 200., by = 5.;
    const double N = tau*bx*by/(M_PI*particle_r*particle_r);
    rng_seed(3 + (uint64_t)variant);
    while (r->N < N){
        struct reb_particle p = {0};
        p.x = rng_uniform(-0.5, 0.5)*bx;
        p.y = rng_uniform(-0.5, 0.5)*by;
        p.z = 10.0*rng_uniform(-0.5, 0.5)*particle_r;
        p.vy = -1.5*p.x;
        p.m = 1.;
        p.r = particle_r;
        reb_simulation_add(r, p);
    }
    return r;
}
static void* create_rings(void){ return simctx_create(make_rings, 1, 0); }

// 10b) Shearing sheet as in examples/shearing_sheet: power-law size distribution (1-4 m, slope -3),
//      tree gravity and tree collisions, Bridges et al. coefficient of restitution (~5,000 particles).
static double cor_bridges(const struct reb_simulation* r, double v){
    (void)r;
    double eps = 0.32*pow(fabs(v)*100.,-0.234);
    if (eps>1) eps=1;
    if (eps<0) eps=0;
    return eps;
}
static double rng_powerlaw(double min, double max, double slope){
    double y = rng_uniform(0., 1.);
    if (slope == -1) return exp(y*log(max/min) + log(min));
    return pow((pow(max,slope+1.)-pow(min,slope+1.))*y+pow(min,slope+1.), 1./(slope+1.));
}
static struct reb_simulation* make_rings_poly(int variant){
    struct reb_simulation* r = new_sim();
    r->opening_angle2 = .5;
    reb_simulation_set_integrator(r, "sei");
    r->boundary = REB_BOUNDARY_SHEAR;
    r->gravity = REB_GRAVITY_TREE;
    r->collision = REB_COLLISION_TREE;
    r->collision_resolve = reb_collision_resolve_hardsphere;
    r->OMEGA = 0.00013143527;
    r->G = 6.67428e-11;
    r->softening = 0.1;
    r->dt = 1e-3*2.*M_PI/r->OMEGA;
    r->root_size = 100.;
    r->N_root_x = 2; r->N_root_y = 2; r->N_root_z = 1;
    r->N_ghost_x = 2; r->N_ghost_y = 2; r->N_ghost_z = 0;
    r->coefficient_of_restitution = cor_bridges;
    r->minimum_collision_velocity = 1.*r->OMEGA*0.001;
    const double bx = 200., by = 200.;
    const double total_mass = 400.*bx*by;
    rng_seed(11 + (uint64_t)variant);
    double mass = 0;
    while (mass < total_mass){
        struct reb_particle p = {0};
        p.x = rng_uniform(-bx/2., bx/2.);
        p.y = rng_uniform(-by/2., by/2.);
        // Box-Muller normal with sigma 1 m
        double u1 = rng_uniform(1e-300, 1.), u2 = rng_uniform(0., 1.);
        p.z = sqrt(-2.*log(u1))*cos(2.*M_PI*u2);
        p.vy = -1.5*p.x*r->OMEGA;
        p.r = rng_powerlaw(1., 4., -3.);
        p.m = 400.*4./3.*M_PI*p.r*p.r*p.r;
        reb_simulation_add(r, p);
        mass += p.m;
    }
    return r;
}
static void* create_rings_poly(void){ return simctx_create(make_rings_poly, 1, 0); }

// 11) Self-gravitating disc with tree gravity (examples/selfgravity_disc), N=10,000
static struct reb_simulation* make_selfgrav(int variant){
    struct reb_simulation* r = new_sim();
    r->rand_seed = 42 + variant;
    reb_simulation_set_integrator(r, "leapfrog");
    r->gravity = REB_GRAVITY_TREE;
    r->boundary = REB_BOUNDARY_OPEN;
    r->opening_angle2 = 1.5;
    r->softening = 0.02;
    r->dt = 3e-2;
    const double boxsize = 10.2;
    r->root_size = boxsize;
    const double disc_mass = 2e-1;
    const int N = 10000;
    struct reb_particle star = {0};
    star.m = 1;
    reb_simulation_add(r, star);
    for (int i=0; i<N; i++){
        struct reb_particle pt = {0};
        double a = reb_random_powerlaw(r, boxsize/10., boxsize/2./1.2, -1.5);
        double phi = reb_random_uniform(r, 0, 2.*M_PI);
        pt.x = a*cos(phi); pt.y = a*sin(phi); pt.z = a*reb_random_normal(r, 0.001);
        double mu = star.m + disc_mass*(pow(a,-1.5)-pow(boxsize/10.,-1.5))/(pow(boxsize/2./1.2,-1.5)-pow(boxsize/10.,-1.5));
        double vkep = sqrt(r->G*mu/a);
        pt.vx = vkep*sin(phi); pt.vy = -vkep*cos(phi);
        pt.m = disc_mass/(double)N;
        reb_simulation_add(r, pt);
    }
    return r;
}
static void* create_selfgrav(void){ return simctx_create(make_selfgrav, 2, 0); }
static struct reb_simulation* make_selfgrav_group(int variant){
    struct reb_simulation* r = make_selfgrav(variant);
    r->tree_group_size = 32;  // opt-in group walk (different, more accurate results at the same opening angle)
    return r;
}
static void* create_selfgrav_group(void){ return simctx_create(make_selfgrav_group, 1, 0); }

// 12) MERCURIUS planetesimal-driven migration (examples/planetesimal_disk_migration, 2,500 planetesimals)
static struct reb_simulation* make_mercurius(int variant){
    struct reb_simulation* r = new_sim();
    r->rand_seed = 12 + variant;
    reb_simulation_set_integrator(r, "mercurius");
    r->testparticle_type = 1;
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    r->boundary = REB_BOUNDARY_OPEN;
    r->root_size = 6; r->N_root_x = 2; r->N_root_y = 2;
    const double m_earth = 3.003e-6, m_neptune = 5.1e-4, a_scat = 1, a_mig = 1.67;
    r->dt = 6.283*pow(a_scat,1.5)/50;
    struct reb_particle star = {0};
    star.m = 1; star.r = 0.005;
    reb_simulation_add(r, star);
    struct reb_particle p = reb_particle_from_orbit(r->G, star, m_neptune, a_scat, 0, reb_random_normal(r, 1e-5), 0, 0, 0);
    p.r = 0.000467;
    reb_simulation_add(r, p);
    p = reb_particle_from_orbit(r->G, star, 2.3*m_earth, a_mig, 0, reb_random_normal(r, 1e-5), 0, 0, 0);
    p.r = 0.0000788215;
    reb_simulation_add(r, p);
    r->N_active = r->N;
    const int Np = 2500;
    const double mp = 2.3*10*m_earth/Np;
    while (r->N < (size_t)Np + r->N_active){
        double a = reb_random_powerlaw(r, a_mig-0.02, a_mig+1, 1);
        double e = reb_random_rayleigh(r, 0.005), inc = reb_random_rayleigh(r, 0.005);
        double Om = reb_random_uniform(r, 0, 2*M_PI), om = reb_random_uniform(r, 0, 2*M_PI), f = reb_random_uniform(r, 0, 2*M_PI);
        struct reb_particle pt = reb_particle_from_orbit(r->G, star, mp, a, e, inc, Om, om, f);
        pt.r = 0.00000934532;
        reb_simulation_add(r, pt);
    }
    reb_simulation_move_to_com(r);
    return r;
}
static void* create_mercurius(void){ return simctx_create(make_mercurius, 20, 0); }

// 13b) Fully active planetesimal disc (all bodies interact, e.g. embryo/planetesimal growth):
//      2 planets + 500 active planetesimals with merging, MERCURIUS (variant 0) or TRACE (variant 1).
static struct reb_simulation* make_disc_active(int variant, int trace){
    struct reb_simulation* r = new_sim();
    r->rand_seed = 12 + variant;
    reb_simulation_set_integrator(r, trace ? "trace" : "mercurius");
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    const double m_earth = 3.003e-6;
    r->dt = 6.283/50;
    struct reb_particle star = {0};
    star.m = 1; star.r = 0.005;
    reb_simulation_add(r, star);
    struct reb_particle p = reb_particle_from_orbit(r->G, star, 5.1e-4, 1.0, 0, reb_random_normal(r, 1e-5), 0, 0, 0);
    p.r = 0.000467;
    reb_simulation_add(r, p);
    p = reb_particle_from_orbit(r->G, star, 2.3*m_earth, 1.67, 0, reb_random_normal(r, 1e-5), 0, 0, 0);
    p.r = 0.0000788215;
    reb_simulation_add(r, p);
    const int Np = 500;
    const double mp = 2.3*10*m_earth/Np;
    for (int i=0; i<Np; i++){
        double a = reb_random_powerlaw(r, 1.65, 2.67, 1);
        double e = reb_random_rayleigh(r, 0.01), inc = reb_random_rayleigh(r, 0.005);
        struct reb_particle pt = reb_particle_from_orbit(r->G, star, mp, a, e, inc, reb_random_uniform(r, 0, 2*M_PI), reb_random_uniform(r, 0, 2*M_PI), reb_random_uniform(r, 0, 2*M_PI));
        pt.r = 0.00000934532*pow(1000./Np, 1./3.);
        reb_simulation_add(r, pt);
    }
    reb_simulation_move_to_com(r);
    return r;
}
static struct reb_simulation* make_mercurius_disc(int variant){ return make_disc_active(variant, 0); }
static struct reb_simulation* make_trace_disc(int variant){ return make_disc_active(variant, 1); }
static void* create_mercurius_disc(void){ return simctx_create(make_mercurius_disc, 1, 0); }
static void* create_trace_disc(void){ return simctx_create(make_trace_disc, 1, 0); }

// 13) TRACE with two planets and 2,000 test particles (debris/Kuiper-belt style)
static struct reb_simulation* make_trace(int variant){
    struct reb_simulation* r = new_sim();
    rng_seed(21 + (uint64_t)variant);
    reb_simulation_add_fmt(r, "m", 1.);
    reb_simulation_add_fmt(r, "m a", 1e-3, 5.);
    reb_simulation_add_fmt(r, "m a", 3e-4, 9.5);
    for (int i=0; i<2000; i++){
        reb_simulation_add_fmt(r, "a e l", rng_uniform(20., 30.), rng_uniform(0., 0.1), rng_uniform(0, 2*M_PI));
    }
    r->N_active = 3;
    r->dt = 0.5;
    reb_simulation_set_integrator(r, "trace");
    return r;
}
static void* create_trace(void){ return simctx_create(make_trace, 5, 0); }

// 14) EOS (PMLF6 kernel with jerk) with two planets and 300 test particles
static struct reb_simulation* make_eos(int variant){
    struct reb_simulation* r = new_sim();
    rng_seed(31 + (uint64_t)variant);
    reb_simulation_add_fmt(r, "m", 1.);
    reb_simulation_add_fmt(r, "m a", 1e-3, 1.);
    reb_simulation_add_fmt(r, "m a", 3e-4, 1.9);
    for (int i=0; i<300; i++){
        reb_simulation_add_fmt(r, "a e l", rng_uniform(2.5, 4.), rng_uniform(0., 0.1), rng_uniform(0, 2*M_PI));
    }
    r->N_active = 3;
    reb_simulation_set_integrator(r, "eos");
    struct reb_integrator_eos_state* eos = r->integrator.state;
    eos->phi0 = REB_INTEGRATOR_EOS_TYPE_PMLF6;
    eos->phi1 = REB_INTEGRATOR_EOS_TYPE_PMLF6;
    r->dt = 0.05;
    return r;
}
static void* create_eos(void){ return simctx_create(make_eos, 10, 0); }

/////////////////////////////////////////////////////////////////////////
// Registry

struct bench {
    const char* name;
    const char* unit;
    const char* usecase;
    void* (*create)(void);
    double (*work)(void*);
    double (*signature)(void*);
    void (*destroy)(void*);
    int (*available)(void);
    int warm_chunks;   // fixed-work mode: chunks run before timing starts
    int meas_chunks;   // fixed-work mode: timed chunks (~0.5 s on the reference machine)
};

static int always(void){ return 1; }

static const struct bench benches[] = {
    // name             unit        use case                                                          create/work/signature/destroy                     available    warm meas
    {"ss_wh_gr",        "steps/s",  "Solar System 5 Gyr ensembles, plain WH + GR (Abbot 2021/2023)",   create_ss_wh_gr,     simctx_work, simctx_signature, simctx_destroy, always,    5,   30},
    {"ss_whckl_gr",     "steps/s",  "Solar System ensembles, WHCKL + GR (Brown & Rein 2020/2022)",      create_ss_whckl_gr,  simctx_work, simctx_signature, simctx_destroy, always,    4,   22},
    {"ss_whfast512",    "steps/s",  "Solar System with WHFast512 + GR (Javaheri et al. 2023)",           create_ss_whfast512, simctx_work, simctx_signature, simctx_destroy, have_avx512,    5,   32},
    {"w512c_ss_avx2",   "steps/s",  "Solar System, WHFast512 C kernel AVX2 + GR (runs without AVX-512)",  create_w512c_ss_avx2,   simctx_work, simctx_signature, simctx_destroy, have_w512c_avx2,   5, 30},
    {"w512c_ss_avx512", "steps/s",  "Solar System, WHFast512 C kernel AVX-512 + GR",                      create_w512c_ss_avx512, simctx_work, simctx_signature, simctx_destroy, have_w512c_avx512, 5, 30},
    {"w512c_3p_avx2",   "steps/s",  "3 compact planets, WHFast512 C AVX2 (four-lane path; SPOCK-like)",   create_w512c_small_avx2, simctx_work, simctx_signature, simctx_destroy, have_w512c_avx2,  5, 40},
    {"w512c_4x2_avx2",  "steps/s",  "4 systems x 2 planets per vector, WHFast512 C AVX2 (ensembles)",     create_w512c_4x2_avx2,  simctx_work, simctx_signature, simctx_destroy, have_w512c_avx2,   5, 40},
    {"w512_out1",       "calls/s",  "Solar System, WHFast512, output after every step (steps(r,1) per call)", create_w512_out1, work_w512_out1, simctx_signature, simctx_destroy, always, 3, 30},
    {"w512c_pair_ss512", "steps/s", "2 Solar Systems per call, WHFast512 C AVX-512 (steps_pair)",          create_pair_ss_avx512, pair_work, pair_signature, pair_destroy, have_w512c_avx512, 2, 16},
    {"w512c_pair_3p_avx2", "steps/s", "2 x 3 compact planets per call, WHFast512 C AVX2 four-lane (steps_pair)", create_pair_3p_avx2, pair_work, pair_signature, pair_destroy, have_w512c_avx2, 2, 24},
    {"obertas_5p",      "steps/s",  "5 Earth-mass planets, WHFast (Obertas et al. 2017)",                create_obertas,      simctx_work, simctx_signature, simctx_destroy, always,    5,   28},
    {"spock_3p_line",   "steps/s",  "Compact 3-planet systems + line collisions (SPOCK, Tamayo 2016/20)", create_spock,       simctx_work, simctx_signature, simctx_destroy, always,   10,   55},
    {"ias15_hltau",     "orbits/s", "IAS15, 5 Saturn-mass planets (Rein & Tamayo 2017)",                create_hltau,        simctx_work, simctx_signature, simctx_destroy, always,    5,   30},
    {"ias15_embryos",   "yr/s",     "IAS15, 60 embryos with merging collisions (giant-impact phase)",    create_embryos,      simctx_work, simctx_signature, simctx_destroy, always,    2,   11},
    {"megno_grid",      "points/s", "MEGNO stability map, 2 planets x 500 yr (Hadden & Lithwick 2018)",  create_megno,        work_megno,  sig_grid,         destroy_free,   always,  100, 1200},
    {"short_posterior", "sims/s",   "Many short WHFast runs of posterior draws (Wang et al. 2018)",      create_posterior,    work_posterior, sig_grid,      destroy_free,   always,  200, 2500},
    {"rings_sei_tree",  "steps/s",  "Saturn's rings, SEI + tree collisions (Rein & Latter 2013)",        create_rings,        simctx_work, simctx_signature, simctx_destroy, always,    3,  400},
    {"rings_poly_tree", "steps/s",  "Shearing sheet with size distribution, tree gravity + collisions (examples/shearing_sheet)", create_rings_poly, simctx_work, simctx_signature, simctx_destroy, always,    2,   20},
    {"selfgrav_tree",   "steps/s",  "Self-gravitating disc, tree gravity, N=10k",                        create_selfgrav,     simctx_work, simctx_signature, simctx_destroy, always,    3,   18},
    {"selfgrav_group", "steps/s",  "Self-gravitating disc, tree gravity with group walk (tree_group_size = 32, opt-in)", create_selfgrav_group, simctx_work, simctx_signature, simctx_destroy, always, 3, 18},
    {"mercurius_plts",  "steps/s",  "MERCURIUS, 2 planets + 2,500 planetesimals (planetesimal migration)", create_mercurius,  simctx_work, simctx_signature, simctx_destroy, always,   10,   50},
    {"trace_tp",        "steps/s",  "TRACE, 2 planets + 2,000 test particles",                           create_trace,        simctx_work, simctx_signature, simctx_destroy, always,   20,  135},
    {"mercurius_disc",  "steps/s",  "MERCURIUS, 2 planets + 500 active planetesimals with merging (embryo/planetesimal growth)", create_mercurius_disc, simctx_work, simctx_signature, simctx_destroy, always, 3, 30},
    {"trace_disc",      "steps/s",  "TRACE, 2 planets + 500 active planetesimals with merging", create_trace_disc, simctx_work, simctx_signature, simctx_destroy, always, 3, 30},
    {"eos_pmlf6_tp",    "steps/s",  "EOS PMLF6, 2 planets + 300 test particles",                         create_eos,          simctx_work, simctx_signature, simctx_destroy, always,    2,   12},
};
static const int N_benches = sizeof(benches)/sizeof(benches[0]);

static const struct bench* find(const char* name){
    for (int i=0; i<N_benches; i++) if (strcmp(benches[i].name, name)==0) return &benches[i];
    return NULL;
}

static void json_str(const char* s){
    putchar('"');
    for (; *s; s++){ if (*s=='"' || *s=='\\') putchar('\\'); putchar(*s); }
    putchar('"');
}

int main(int argc, char* argv[]){
    const char* mode = NULL; const char* name = NULL;
    double warmup = 0.3, total = 1.0, window = 0.05;
    int check_chunks = 3;
    int fixed_work = 0, warm_override = -1, meas_override = -1;
    for (int i=1; i<argc; i++){
        if (!strcmp(argv[i], "--list")) mode = "list";
        else if (!strcmp(argv[i], "--run") && i+1<argc){ mode = "run"; name = argv[++i]; }
        else if (!strcmp(argv[i], "--check") && i+1<argc){ mode = "check"; name = argv[++i]; }
        else if (!strcmp(argv[i], "--warmup") && i+1<argc) warmup = atof(argv[++i]);
        else if (!strcmp(argv[i], "--time") && i+1<argc) total = atof(argv[++i]);
        else if (!strcmp(argv[i], "--window") && i+1<argc) window = atof(argv[++i]);
        else if (!strcmp(argv[i], "--chunks") && i+1<argc) check_chunks = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--fixed") && i+1<argc){ mode = "run"; name = argv[++i]; fixed_work = 1; }
        else if (!strcmp(argv[i], "--warm-chunks") && i+1<argc) warm_override = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--meas-chunks") && i+1<argc) meas_override = atoi(argv[++i]);
        else { fprintf(stderr, "Unknown argument %s\n", argv[i]); return 2; }
    }
    if (!mode){ fprintf(stderr, "Usage: bench --list | --fixed NAME [--warm-chunks N --meas-chunks N] | --run NAME [--time S --window S] | --check NAME [--chunks N]\n"); return 2; }

    if (!strcmp(mode, "list")){
        printf("[");
        for (int i=0; i<N_benches; i++){
            printf(i?",\n {":"{");
            printf("\"name\":"); json_str(benches[i].name);
            printf(",\"unit\":"); json_str(benches[i].unit);
            printf(",\"usecase\":"); json_str(benches[i].usecase);
            printf(",\"available\":%d,\"warm_chunks\":%d,\"meas_chunks\":%d}", benches[i].available(), benches[i].warm_chunks, benches[i].meas_chunks);
        }
        printf("]\n");
        return 0;
    }

    const struct bench* b = find(name);
    if (!b){ fprintf(stderr, "Unknown benchmark %s\n", name); return 2; }
    if (!b->available()){ printf("{\"name\":\"%s\",\"skipped\":\"not available on this CPU\"}\n", b->name); return 0; }

    if (!strcmp(mode, "check")){
        // Fixed amount of work from a fresh start; the signature detects changes in results.
        void* ctx = b->create();
        double w = 0;
        for (int k=0; k<check_chunks; k++) w += b->work(ctx);
        printf("{\"name\":\"%s\",\"work\":%.17g,\"signature\":%.17g}\n", b->name, w, b->signature(ctx));
        b->destroy(ctx);
        return 0;
    }

    double t_setup = now();
    void* ctx = b->create();
    t_setup = now() - t_setup;

    if (fixed_work){
        // Fixed-work run: identical simulated interval every time, independent of speed.
        int warm = warm_override>=0 ? warm_override : b->warm_chunks;
        int meas = meas_override>=0 ? meas_override : b->meas_chunks;
        double t_warm = now();
        for (int k=0; k<warm; k++) b->work(ctx);
        t_warm = now() - t_warm;
        printf("{\"name\":\"%s\",\"unit\":\"%s\",\"mode\":\"fixed\",\"setup_s\":%.6f,\"warmup_s\":%.6f,\"chunks\":[", b->name, b->unit, t_setup, t_warm);
        for (int k=0; k<meas; k++){
            double t0 = now();
            double w = b->work(ctx);
            printf("%s[%.10g,%.9f]", k?",":"", w, now()-t0);
        }
        printf("],\"signature\":%.17g}\n", b->signature(ctx));
        b->destroy(ctx);
        return 0;
    }

    // Time-boxed run (used by `run.py probe`): windowed rates over wall time.
    size_t cap = 4096, n = 0;
    double* rates = malloc(cap*sizeof(double));
    double* times = malloc(cap*sizeof(double));
    const double t_start = now();
    double w_win = 0, t_win = t_start, t_now = t_start;
    size_t chunks = 0;
    while (t_now - t_start < total){
        w_win += b->work(ctx);
        chunks++;
        t_now = now();
        if (t_now - t_win >= window){
            if (n==cap){ cap*=2; rates = realloc(rates, cap*sizeof(double)); times = realloc(times, cap*sizeof(double)); }
            rates[n] = w_win/(t_now - t_win);
            times[n] = t_now - t_start;
            n++;
            w_win = 0; t_win = t_now;
        }
    }
    printf("{\"name\":\"%s\",\"unit\":\"%s\",\"mode\":\"timed\",\"setup_s\":%.6f,\"warmup_s\":%g,\"chunks\":%zu,\"chunk_s\":%.6g,\"windows\":[", b->name, b->unit, t_setup, warmup, chunks, (t_now-t_start)/(double)chunks);
    for (size_t i=0; i<n; i++) printf("%s[%.4f,%.6g]", i?",":"", times[i], rates[i]);
    printf("]}\n");
    b->destroy(ctx);
    free(rates); free(times);
    return 0;
}
