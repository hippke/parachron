// IAS15 stage mode (ias15.parallel_stages = 1) in all supported configurations: prints, per
// configuration, a bitwise hash of the final state in stage mode and the largest position
// difference to the default mode. The hashes must be identical for the serial library and the
// OpenMP library at any OMP_NUM_THREADS.
//   cc -O2 -I src stage_check.c -L src -lrebound -Wl,-rpath,$PWD/src -lm -o stage_check
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "rebound.h"
#include "integrator_ias15.h"

static uint64_t rng;
static double urand(void){ rng ^= rng<<13; rng ^= rng>>7; rng ^= rng<<17; return (rng>>11)*(1.0/9007199254740992.0); }

struct cfg { const char* name; size_t N; size_t N_active; int tptype; int compensated; int ghost; unsigned int ignore; double tmax; };

static struct reb_simulation* make(const struct cfg* c, int stages){
    rng = 88172645463325252ULL + c->N*31 + c->N_active*7 + c->tptype*3 + c->compensated*5 + c->ghost*11 + c->ignore*13;
    struct reb_simulation* r = reb_simulation_create();
    reb_simulation_set_integrator(r, "ias15");
    ((struct reb_integrator_ias15_state*)r->integrator.state)->parallel_stages = stages;
    r->gravity = c->compensated ? REB_GRAVITY_COMPENSATED : REB_GRAVITY_BASIC;
    if (c->ghost){
        // shearing-free periodic box with softening (a planetesimal patch)
        r->softening = 0.02; r->G = 1.;
        r->root_size = 1.; r->N_root_x = 1; r->N_root_y = 1; r->N_root_z = 1;
        r->boundary = REB_BOUNDARY_PERIODIC; r->N_ghost_x = 1; r->N_ghost_y = 1; r->N_ghost_z = 0;
        for (size_t i=0; i<c->N; i++){
            struct reb_particle p = {0};
            p.x = urand()-0.5; p.y = urand()-0.5; p.z = 0.05*(urand()-0.5);
            p.vx = 0.01*(urand()-0.5); p.vy = 0.01*(urand()-0.5);
            p.m = 1e-4*(0.5+urand());
            reb_simulation_add(r, p);
        }
    }else{
        struct reb_particle star = {0}; star.m = 1.;
        reb_simulation_add(r, star);
        for (size_t i=1; i<c->N; i++){
            const double a = 1. + 2.*urand();
            const double m = i<c->N_active ? 1e-6*(0.5+urand()) : 1e-10*urand();
            reb_simulation_add_fmt(r, "primary m a e inc l pomega Omega", star, m, a, 0.05*urand(), 0.02*urand(), 6.283*urand(), 6.283*urand(), 6.283*urand());
        }
        reb_simulation_move_to_com(r);
    }
    if (c->N_active && c->N_active<c->N){ r->N_active = c->N_active; r->testparticle_type = c->tptype; }
    r->gravity_ignore_terms = c->ignore;
    return r;
}

static uint64_t hash(struct reb_simulation* r){
    uint64_t h = 1469598103934665603ULL;
    for (size_t i=0; i<r->N; i++){
        const double v[6] = {r->particles[i].x, r->particles[i].y, r->particles[i].z, r->particles[i].vx, r->particles[i].vy, r->particles[i].vz};
        uint64_t u[6]; memcpy(u, v, sizeof(u));
        for (int k=0; k<6; k++){ h ^= u[k]; h *= 1099511628211ULL; }
    }
    return h;
}

int main(void){
    const struct cfg cfgs[] = {
        {"basic, all active, N=40",           40,   0, 0, 0, 0, 0, 3.},
        {"basic, all active, N=200",         200,   0, 0, 0, 0, 0, 0.5},
        {"basic, 20 active + 100 tp type 0", 120,  20, 0, 0, 0, 0, 2.},
        {"basic, 20 active + 100 tp type 1", 120,  20, 1, 0, 0, 0, 2.},
        {"basic, 150 active + 300 tp type 1",450, 150, 1, 0, 0, 0, 0.3},
        {"basic, 3 active + 2000 tp type 0",2003,   3, 0, 0, 0, 0, 0.3},
        {"basic, ignore terms 0-1, N=80",     80,   0, 0, 0, 0, 1, 1.},
        {"basic, ignore terms of 0, N=80",    80,   0, 0, 0, 0, 2, 1.},
        {"basic, ghost boxes, N=60",          60,   0, 0, 0, 1, 0, 0.5},
        {"basic, ghost boxes, N=300",        300,   0, 0, 0, 1, 0, 0.2},
        {"compensated, N=40",                 40,   0, 0, 1, 0, 0, 3.},
        {"compensated, N=200",               200,   0, 0, 1, 0, 0, 0.5},
        {"compensated, 30 active + 90 tp t1",120,  30, 1, 1, 0, 0, 2.},
        {"compensated, ignore terms 0-1",     80,   0, 0, 1, 0, 1, 1.},
    };
    for (size_t k=0; k<sizeof(cfgs)/sizeof(cfgs[0]); k++){
        const struct cfg* c = &cfgs[k];
        struct reb_simulation* a = make(c, 0);
        struct reb_simulation* b = make(c, 1);
        reb_simulation_integrate(a, c->tmax);
        reb_simulation_integrate(b, c->tmax);
        double dmax = 0, scale = 0;
        for (size_t i=0; i<a->N; i++){
            const double d = fabs(a->particles[i].x-b->particles[i].x)+fabs(a->particles[i].y-b->particles[i].y)+fabs(a->particles[i].z-b->particles[i].z);
            if (d>dmax) dmax = d;
            const double s = fabs(a->particles[i].x)+fabs(a->particles[i].y);
            if (s>scale) scale = s;
        }
        printf("%-36s steps %5llu/%5llu  max|dx|/scale %.1e  stage hash %016llx\n", c->name,
               (unsigned long long)a->steps_done, (unsigned long long)b->steps_done, dmax/scale, (unsigned long long)hash(b));
        reb_simulation_free(a); reb_simulation_free(b);
    }
    return 0;
}
