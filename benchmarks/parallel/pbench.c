// Multicore benchmarks for REBOUND (public API only, so it builds against any source tree).
//
//   ./pbench CASE [scale]
//
// Runs a fixed amount of work (fixed simulated time or number of steps) for one case and prints one
// JSON line: {"case":..., "seconds":..., "steps":..., "signature":...}. The thread count is set
// from outside with OMP_NUM_THREADS (OpenMP build) and has no effect on the serial build.
// The signature is a weighted sum over the final particle state: equal signatures for different
// thread counts mean the result does not depend on the number of threads.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include "rebound.h"
#include "integrator_whfast.h"

static double now(void){
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9*ts.tv_nsec;
}

static uint64_t rng_state = 0x9E3779B97F4A7C15ULL;
static void rng_seed(uint64_t s){ rng_state = s*0x9E3779B97F4A7C15ULL + 1; }
static double rng_uniform(double a, double b){
    rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7; rng_state ^= rng_state << 17;
    return a + (b-a)*((rng_state >> 11)*(1.0/9007199254740992.0));
}

static double signature(struct reb_simulation* r){
    double s = 0;
    for (size_t i=0; i<r->N; i++){
        const struct reb_particle p = r->particles[i];
        const double w = 1. + 0.618*(double)i;
        s += w*(p.x + 2.*p.y + 3.*p.z + 5.*p.vx + 7.*p.vy + 11.*p.vz) + 13.*p.m;
    }
    return s;
}

static struct reb_simulation* new_sim(void){
    struct reb_simulation* r = reb_simulation_create();
    r->exit_max_distance = 0; // no escape checks
    return r;
}

// Star plus N-1 equal-mass bodies in a thin, dynamically cold disc.
static void add_disc(struct reb_simulation* r, int N, double m, double a0, double a1, uint64_t seed){
    rng_seed(seed);
    struct reb_particle star = {0};
    star.m = 1.;
    reb_simulation_add(r, star);
    for (int i=1; i<N; i++){
        const double a = rng_uniform(a0, a1);
        reb_simulation_add_fmt(r, "primary m a e inc l pomega Omega", star, m, a,
                rng_uniform(0, 0.02), rng_uniform(0, 0.01), rng_uniform(0, 2*M_PI), rng_uniform(0, 2*M_PI), rng_uniform(0, 2*M_PI));
    }
    reb_simulation_move_to_com(r);
}

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
        p.x = ss_pos[i][0]; p.y = ss_pos[i][1]; p.z = ss_pos[i][2];
        p.vx = ss_vel[i][0]; p.vy = ss_vel[i][1]; p.vz = ss_vel[i][2];
        p.m = ss_mass[i];
        reb_simulation_add(r, p);
    }
    reb_simulation_move_to_com(r);
}

// Returns a simulation and the amount of work: integrate to t_end (adaptive) or n_steps (fixed).
struct work { struct reb_simulation* r; double t_end; long n_steps; };

static struct work make_case(const char* name, double scale){
    struct work w = {0};
    struct reb_simulation* r = new_sim();
    w.r = r;
    r->G = 4.*M_PI*M_PI; // AU, Msun, yr
    if (!strcmp(name, "ias15_hltau")){          // 5 planets: must not get slower
        reb_simulation_set_integrator(r, "ias15");
        rng_seed(7);
        struct reb_particle star = {0}; star.m = 1.;
        reb_simulation_add(r, star);
        const double gaps[5] = {13.2, 32.3, 42.2, 68.8, 74.7};
        for (int i=0; i<5; i++) reb_simulation_add_fmt(r, "primary m a inc l", star, 3e-4, gaps[i], rng_uniform(0, 1e-3), rng_uniform(0, 2*M_PI));
        reb_simulation_move_to_com(r);
        r->dt = 0.1; w.t_end = 80000.*scale;
    }else if (!strcmp(name, "ias15_embryos")){ // 60 embryos + merging (as in benchmarks/perf)
        reb_simulation_set_integrator(r, "ias15");
        r->collision = REB_COLLISION_DIRECT;
        r->collision_resolve = reb_collision_resolve_merge;
        rng_seed(11);
        struct reb_particle star = {0}; star.m = 1.; star.r = 0.005;
        reb_simulation_add(r, star);
        for (int i=0; i<60; i++){
            const double a = 0.5 + 1.0*i/60.;
            reb_simulation_add_fmt(r, "primary m r a e inc l pomega", star, 3.7e-8, 2e-5, a,
                    rng_uniform(0, 0.02), rng_uniform(0, 0.01), rng_uniform(0, 2*M_PI), rng_uniform(0, 2*M_PI));
        }
        reb_simulation_move_to_com(r);
        r->dt = 1e-3; w.t_end = 20.*scale;
    }else if (!strcmp(name, "ias15_n128") || !strcmp(name, "ias15_n512") || !strcmp(name, "ias15_n2048")){
        const int N = atoi(name+7);
        reb_simulation_set_integrator(r, "ias15");
        add_disc(r, N, 1e-7, 1., 3., 31);
        r->dt = 1e-3;
        w.t_end = (N==128 ? 10. : (N==512 ? 0.5 : 0.04))*scale;
    }else if (!strcmp(name, "whfast_ss")){      // 8 planets: must not get slower
        reb_simulation_set_integrator(r, "whfast");
        r->G = 1.;                        // units as in benchmarks/perf (G=1, AU, Msun)
        add_solar_system(r);
        r->dt = 5./365.25*2.*M_PI;        // 5 days
        w.n_steps = (long)(200000*scale);
    }else if (!strcmp(name, "whfast_tp10k")){   // 2 planets + 10,000 test particles
        reb_simulation_set_integrator(r, "whfast");
        rng_seed(21);
        reb_simulation_add_fmt(r, "m", 1.);
        reb_simulation_add_fmt(r, "m a", 1e-3, 5.);
        reb_simulation_add_fmt(r, "m a", 3e-4, 9.5);
        for (int i=0; i<10000; i++) reb_simulation_add_fmt(r, "a e l", rng_uniform(20., 30.), rng_uniform(0., 0.1), rng_uniform(0, 2*M_PI));
        r->N_active = 3;
        r->dt = 0.5; w.n_steps = (long)(400*scale);
    }else if (!strcmp(name, "whfast_n256")){    // 256 massive bodies (planetesimal disc)
        reb_simulation_set_integrator(r, "whfast");
        add_disc(r, 256, 1e-8, 1., 3., 41);
        r->dt = 0.01; w.n_steps = (long)(4000*scale);
    }else{
        fprintf(stderr, "unknown case %s\n", name);
        exit(1);
    }
    return w;
}

int main(int argc, char** argv){
    if (argc<2){
        fprintf(stderr, "usage: %s CASE [scale]\ncases: ias15_hltau ias15_embryos ias15_n128 ias15_n512 ias15_n2048 whfast_ss whfast_tp10k whfast_n256\n", argv[0]);
        return 1;
    }
    const double scale = argc>2 ? atof(argv[2]) : 1.;
    struct work w = make_case(argv[1], scale);
    struct reb_simulation* r = w.r;
    if (r->integrator.state && !strncmp(argv[1], "whfast", 6)){
        ((struct reb_integrator_whfast_state*)r->integrator.state)->safe_mode = 0; // as in long research runs
    }
    // warm-up: one step (allocations, first force evaluation)
    reb_simulation_steps(r, 1);
    const long steps0 = r->steps_done;
    const double t0 = now();
    if (w.n_steps){
        reb_simulation_steps(r, w.n_steps);
    }else{
        reb_simulation_integrate(r, r->t + w.t_end);
    }
    const double t1 = now();
    printf("{\"case\": \"%s\", \"seconds\": %.6f, \"steps\": %ld, \"N\": %zu, \"signature\": %.17g}\n",
            argv[1], t1-t0, (long)(r->steps_done-steps0), r->N, signature(r));
    reb_simulation_free(r);
    return 0;
}
