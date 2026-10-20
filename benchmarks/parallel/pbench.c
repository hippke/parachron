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
#include "integrator_ias15.h"

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

// --- Setups ported from benchmarks/perf/bench.c (same initial conditions) ---
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
// Self-gravitating disc, tree gravity, N = 10,000 (examples/selfgravity_disc)
static void setup_selfgrav(struct reb_simulation* r){
    r->rand_seed = 42;
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
}
// Shearing sheet with size distribution, tree gravity + tree collisions (~5,000 particles)
static void setup_rings_poly(struct reb_simulation* r){
    r->rand_seed = 11;   // collision resolution order is shuffled with rand_seed (random by default)
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
    rng_seed(11);
    double mass = 0;
    while (mass < total_mass){
        struct reb_particle p = {0};
        p.x = rng_uniform(-bx/2., bx/2.);
        p.y = rng_uniform(-by/2., by/2.);
        double u1 = rng_uniform(1e-300, 1.), u2 = rng_uniform(0., 1.);
        p.z = sqrt(-2.*log(u1))*cos(2.*M_PI*u2);
        p.vy = -1.5*p.x*r->OMEGA;
        p.r = rng_powerlaw(1., 4., -3.);
        p.m = 400.*4./3.*M_PI*p.r*p.r*p.r;
        reb_simulation_add(r, p);
        mass += p.m;
    }
}
// MERCURIUS planetesimal-driven migration: 2 planets + 2,500 planetesimals (test particles type 1)
static void setup_mercurius_plts(struct reb_simulation* r){
    r->rand_seed = 12;
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
}
// Fully active planetesimal disc: 2 planets + 500 active planetesimals with merging
static void setup_disc_active(struct reb_simulation* r, int trace){
    r->rand_seed = 12;
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
}
// TRACE with two planets and 2,000 test particles
static void setup_trace_tp(struct reb_simulation* r){
    rng_seed(21);
    reb_simulation_add_fmt(r, "m", 1.);
    reb_simulation_add_fmt(r, "m a", 1e-3, 5.);
    reb_simulation_add_fmt(r, "m a", 3e-4, 9.5);
    for (int i=0; i<2000; i++){
        reb_simulation_add_fmt(r, "a e l", rng_uniform(20., 30.), rng_uniform(0., 0.1), rng_uniform(0, 2*M_PI));
    }
    r->N_active = 3;
    r->dt = 0.5;
    reb_simulation_set_integrator(r, "trace");
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
    }else if (!strcmp(name, "ias15_tp2000")){  // 20 massive bodies + 2,000 test particles (type 0)
        reb_simulation_set_integrator(r, "ias15");
        add_disc(r, 20, 1e-6, 1., 3., 51);
        rng_seed(52);
        const struct reb_particle star = r->particles[0];
        for (int i=0; i<2000; i++) reb_simulation_add_fmt(r, "primary a e inc l", star, rng_uniform(1., 3.), rng_uniform(0, 0.05), rng_uniform(0, 0.02), rng_uniform(0, 2*M_PI));
        r->N_active = 20;
        r->dt = 1e-3; w.t_end = 2.*scale;
    }else if (!strcmp(name, "tree_selfgrav")){  // leapfrog, tree gravity, N = 10,000
        r->G = 1.;
        setup_selfgrav(r);
        w.n_steps = (long)(40*scale);
    }else if (!strcmp(name, "tree_selfgrav_group")){ // same, opt-in group walk (tree_group_size = 32)
        r->G = 1.;
        setup_selfgrav(r);
#ifdef REB_HAS_TREE_GROUP_SIZE
        r->tree_group_size = 32;
#endif
        w.n_steps = (long)(40*scale);
    }else if (!strcmp(name, "tree_rings")){     // SEI, shearing sheet, tree gravity + collisions
        setup_rings_poly(r);
        w.n_steps = (long)(60*scale);
    }else if (!strcmp(name, "mercurius_plts")){ // 2 planets + 2,500 planetesimals (type 1)
        r->G = 1.;
        setup_mercurius_plts(r);
        w.n_steps = (long)(1000*scale);
    }else if (!strcmp(name, "mercurius_disc") || !strcmp(name, "trace_disc")){ // 500 active bodies + merging
        r->G = 1.;
        setup_disc_active(r, name[0]=='t');
        w.n_steps = (long)(200*scale);
    }else if (!strcmp(name, "trace_tp")){       // 2 planets + 2,000 test particles
        r->G = 1.;
        setup_trace_tp(r);
        w.n_steps = (long)(2000*scale);
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

// Ensembles: K independent Solar System clones (positions perturbed at the 1e-12 level), advanced by
// reb_simulations_steps with as many workers as OMP_NUM_THREADS (also in the serial build: the
// ensemble pool uses plain threads). ensemble_ss512: WHFast512 + GR (1000 time steps per REBOUND
// step); ensemble_ss: WHFast. The signature is the sum over all clones.
static int run_ensemble(const char* name, double scale){
    // ensemble_3p512: 3 inner planets (Mercury, Venus, Earth) per system, WHFast512 (the four-lane
    // AVX2 kernel and the AVX-512 kernels can pair two simulations per kernel call). The kernel can
    // be set with PBENCH_W512_KERNEL (0 auto, 1 asm, 2 C AVX-512, 3 C AVX2, 4 C generic).
    const int is3p = !strcmp(name, "ensemble_3p512");
    const int is512 = !strcmp(name, "ensemble_ss512") || is3p;
    if (!is512 && strcmp(name, "ensemble_ss")){ fprintf(stderr, "unknown case %s\n", name); return 1; }
    const char* ek = getenv("PBENCH_W512_KERNEL");
    const unsigned int kernel = ek ? (unsigned int)atoi(ek) : 0;
    const char* e = getenv("OMP_NUM_THREADS");
    const int nthreads = e ? atoi(e) : 1;
    const int K = 32;
    struct reb_simulation* sims[32];
    for (int k=0; k<K; k++){
        struct reb_simulation* r = new_sim();
        r->G = 1.;
        if (is512){
            reb_simulation_set_integrator(r, "whfast512");
            struct reb_integrator_whfast512_state* w512 = r->integrator.state;
            w512->gr_potential = 1;
            w512->concatenate_steps = 1000;
            w512->kernel = kernel;
            r->dt = 5.0/365.25*2.*M_PI;
            r->exact_finish_time = 0;
        }else{
            reb_simulation_set_integrator(r, "whfast");
            r->dt = 5./365.25*2.*M_PI;
        }
        add_solar_system(r);
        if (is3p){ while (r->N > 4) reb_simulation_remove_particle(r, r->N-1); }
        r->particles[1].x += 1e-12*k;
        if (!is512) ((struct reb_integrator_whfast_state*)r->integrator.state)->safe_mode = 0;
        reb_simulation_steps(r, 1); // warm-up
        sims[k] = r;
    }
    const long n_steps = (long)((is512 ? 200 : 25000)*scale);
    const double t0 = now();
    reb_simulations_steps(sims, (size_t)K, (size_t)n_steps, nthreads);
    const double t1 = now();
    double sig = 0;
    long steps = 0;
    for (int k=0; k<K; k++){ sig += signature(sims[k]); steps += (long)sims[k]->steps_done-1; reb_simulation_free(sims[k]); }
    printf("{\"case\": \"%s\", \"seconds\": %.6f, \"steps\": %ld, \"N\": %d, \"signature\": %.17g}\n", name, t1-t0, steps, K, sig);
    return 0;
}

int main(int argc, char** argv){
    if (argc<2){
        fprintf(stderr, "usage: %s CASE [scale]\ncases (IAS15 cases also as CASE+stages): ias15_hltau ias15_embryos ias15_n128 ias15_n512 ias15_n2048 ias15_tp2000 tree_selfgrav tree_selfgrav_group tree_rings mercurius_plts mercurius_disc trace_disc trace_tp ensemble_ss512 ensemble_3p512 ensemble_ss whfast_ss whfast_tp10k whfast_n256\n", argv[0]);
        return 1;
    }
    const double scale = argc>2 ? atof(argv[2]) : 1.;
    if (!strncmp(argv[1], "ensemble_", 9)) return run_ensemble(argv[1], scale);
    // A case name ending in "+stages" runs the IAS15 case with ias15.parallel_stages = 1.
    char name[64];
    snprintf(name, sizeof(name), "%s", argv[1]);
    char* plus = strstr(name, "+stages");
    if (plus) *plus = 0;
    struct work w = make_case(name, scale);
    struct reb_simulation* r = w.r;
    if (plus){
        if (strncmp(name, "ias15", 5)){ fprintf(stderr, "+stages needs an ias15 case\n"); return 1; }
        ((struct reb_integrator_ias15_state*)r->integrator.state)->parallel_stages = 1;
    }
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
