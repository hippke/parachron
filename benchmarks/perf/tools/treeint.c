// Self-gravitating disc, leapfrog, tree gravity: energy drift and position deviation vs a
// direct-summation run (same initial conditions), after n steps.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "rebound.h"
static uint32_t s=11; static double u(void){ s=s*1664525u+1013904223u; return (s+0.5)/4294967296.0; }
static struct reb_simulation* make(void){
    s=11;
    struct reb_simulation* r = reb_simulation_create();
    reb_simulation_set_integrator(r, "leapfrog"); r->dt = 1e-2;
    r->boundary = REB_BOUNDARY_OPEN; r->softening = 0.02; r->root_size = 10.2*4; r->G = 1;
    struct reb_particle star = {0}; star.m = 1; reb_simulation_add(r, star);
    for (int i=0;i<3000;i++){
        double a = pow(pow(1.,-1.5) + u()*(pow(5.,-1.5)-pow(1.,-1.5)), -1./1.5);
        double phi = 2*M_PI*u(), vk = sqrt(1./a);
        struct reb_particle p = {0};
        p.x = a*cos(phi); p.y = a*sin(phi); p.z = a*0.05*(u()-0.5);
        p.vx = -vk*sin(phi); p.vy = vk*cos(phi); p.m = 0.2/3000.;
        reb_simulation_add(r, p);
    }
    reb_simulation_move_to_com(r);
    return r;
}
int main(int argc, char** argv){
    int nsteps = atoi(argv[1]); double theta = atof(argv[2]);
    struct reb_simulation* d = make(); d->gravity = REB_GRAVITY_BASIC;
    struct reb_simulation* t = make(); t->gravity = REB_GRAVITY_TREE; t->opening_angle2 = theta*theta;
    double E0 = reb_simulation_energy(d);
    clock_t c0=clock(); reb_simulation_steps(t, nsteps); double tt=(double)(clock()-c0)/CLOCKS_PER_SEC;
    reb_simulation_steps(d, nsteps);
    double Ed = reb_simulation_energy(d), Et = reb_simulation_energy(t);
    double dev=0; for (size_t i=1;i<d->N;i++){ double dx=d->particles[i].x-t->particles[i].x, dy=d->particles[i].y-t->particles[i].y, dz=d->particles[i].z-t->particles[i].z; dev += sqrt(dx*dx+dy*dy+dz*dz); }
    printf("theta=%.2f steps=%d  dE/E direct=%.2e tree=%.2e  mean |x_tree-x_direct|=%.2e  tree time %.2fs\n", theta, nsteps, fabs((Ed-E0)/E0), fabs((Et-E0)/E0), dev/(d->N-1), tt);
    return 0;
}
