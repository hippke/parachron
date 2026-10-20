// Bitwise check of the gravity routines: prints one hash of the exact acceleration bits per
// configuration. Run it against the serial library and against the OpenMP library with several
// OMP_NUM_THREADS values; all outputs must be identical. "gravity_check blocks" does the same for
// basic gravity with gravity_block_order = 1.
//
//   cc -O2 -I src gravity_check.c -L src -lrebound -Wl,-rpath,$PWD/src -lm -o gravity_check
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "rebound.h"

static uint64_t rng = 88172645463325252ULL;
static double urand(void){ rng ^= rng<<13; rng ^= rng>>7; rng ^= rng<<17; return (rng>>11)*(1.0/9007199254740992.0); }

static uint64_t hash_acc(struct reb_simulation* r){
    uint64_t h = 1469598103934665603ULL;
    for (size_t i=0; i<r->N; i++){
        const double a[3] = {r->particles[i].ax, r->particles[i].ay, r->particles[i].az};
        uint64_t u[3];
        memcpy(u, a, sizeof(u));
        for (int c=0; c<3; c++){ h ^= u[c]; h *= 1099511628211ULL; }
    }
    return h;
}

int main(int argc, char** argv){
    // With an argument "blocks": basic gravity with gravity_block_order = 1 (compensated skipped).
    const int blocks = argc>1 && !strcmp(argv[1], "blocks");
    const size_t Ns[] = {3, 17, 40, 200, 700, 1100, 1600};
    const char* gravities[] = {"basic", "compensated"};
    for (int g=0; g<(blocks?1:2); g++){
    for (int in=0; in<7; in++){
    for (int tp=0; tp<3; tp++){          // 0: no test particles, 1: type 0, 2: type 1
    for (int ign=0; ign<3; ign++){       // gravity_ignore_terms
    for (int ghost=0; ghost<2; ghost++){
        if (g==1 && ghost) continue;     // compensated gravity has no ghost boxes
        rng = 88172645463325252ULL + 7919*in + 104729*tp + 13*ign + 31*ghost;
        struct reb_simulation* r = reb_simulation_create();
        reb_simulation_set_integrator(r, "leapfrog");
        r->gravity = g==0 ? REB_GRAVITY_BASIC : REB_GRAVITY_COMPENSATED;
        r->gravity_block_order = blocks;
        r->softening = ghost ? 0.01 : 0.;
        const size_t N = Ns[in];
        for (size_t i=0; i<N; i++){
            struct reb_particle p = {0};
            p.x = urand()-0.5; p.y = urand()-0.5; p.z = 0.1*(urand()-0.5);
            p.m = (i==0) ? 1. : 1e-3*urand();
            reb_simulation_add(r, p);
        }
        if (tp){
            r->N_active = N/2>2 ? N/2 : 2;
            r->testparticle_type = tp==2;
        }
        r->gravity_ignore_terms = ign;
        if (ghost){
            r->root_size = 1.; r->N_root_x = 1; r->N_root_y = 1; r->N_root_z = 1;
            r->boundary = REB_BOUNDARY_PERIODIC;
            r->N_ghost_x = 1; r->N_ghost_y = 1; r->N_ghost_z = 0;
        }
        reb_simulation_update_acceleration(r);
        printf("%s N=%4zu tp=%d ignore=%d ghost=%d  %016llx\n", gravities[g], N, tp, ign, ghost, (unsigned long long)hash_acc(r));
        reb_simulation_free(r);
    }}}}}
    return 0;
}
