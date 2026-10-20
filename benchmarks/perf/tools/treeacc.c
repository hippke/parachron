// Tree gravity accuracy vs cost: relative force error against direct summation, time per evaluation.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "rebound.h"
extern int reb_tree_mac_bmax;
static double now(){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec; }
static uint32_t s=11; static double u(void){ s=s*1664525u+1013904223u; return (s+0.5)/4294967296.0; }
static struct reb_simulation* make(int kind){
    struct reb_simulation* r = reb_simulation_create();
    if (kind==0){ // self-gravitating disc, as benchmarks/perf selfgrav_tree (examples/selfgravity_disc)
        r->boundary = REB_BOUNDARY_OPEN; r->softening = 0.02; r->root_size = 10.2; r->G = 1;
        struct reb_particle star = {0}; star.m = 1; reb_simulation_add(r, star);
        for (int i=0;i<10000;i++){
            double a = pow(pow(0.1,-1.5) + u()*(pow(5.,-1.5)-pow(0.1,-1.5)), -1./1.5);
            double phi = 2*M_PI*u(), vk = sqrt(1./a);
            struct reb_particle p = {0};
            p.x = a*cos(phi); p.y = a*sin(phi); p.z = a*0.05*(u()-0.5);
            p.vx = -vk*sin(phi); p.vy = vk*cos(phi); p.m = 0.2/10000.;
            reb_simulation_add(r, p);
        }
    }else{ // shearing sheet with size distribution (examples/shearing_sheet), ghost boxes
        r->boundary = REB_BOUNDARY_SHEAR; r->OMEGA = 0.00013143527; r->G = 6.67428e-11; r->softening = 0.1;
        r->root_size = 100.; r->N_root_x = 2; r->N_root_y = 2; r->N_root_z = 1;
        r->N_ghost_x = 2; r->N_ghost_y = 2; r->N_ghost_z = 0;
        double mass = 0;
        while (mass < 400.*200.*200.){
            struct reb_particle p = {0};
            p.x = (u()-0.5)*200.; p.y = (u()-0.5)*200.;
            double u1=u(), u2=u(); p.z = sqrt(-2.*log(u1))*cos(2.*M_PI*u2);
            double y=u(); p.r = pow((pow(4.,-2.)-pow(1.,-2.))*y+pow(1.,-2.), -0.5);
            p.m = 400.*4./3.*M_PI*p.r*p.r*p.r; reb_simulation_add(r, p); mass += p.m;
        }
    }
    return r;
}
int main(int argc, char** argv){
    int kind = atoi(argv[1]);
    struct reb_simulation* r = make(kind);
    size_t N = r->N;
    double* ad = malloc(sizeof(double)*3*N);
    r->gravity = REB_GRAVITY_BASIC;
    reb_simulation_update_acceleration(r);
    for (size_t i=0;i<N;i++){ ad[3*i]=r->particles[i].ax; ad[3*i+1]=r->particles[i].ay; ad[3*i+2]=r->particles[i].az; }
    r->gravity = REB_GRAVITY_TREE;
    double thetas[] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2};
    double* err = malloc(sizeof(double)*N);
    for (int mac=0; mac<2; mac++){
        reb_tree_mac_bmax = mac;
        for (int t=0;t<10;t++){
            r->opening_angle2 = thetas[t]*thetas[t];
            reb_simulation_update_acceleration(r);  // warm up
            int rep = 3; double t0 = now();
            for (int k=0;k<rep;k++) reb_simulation_update_acceleration(r);
            double dt = (now()-t0)/rep;
            for (size_t i=0;i<N;i++){
                double ex=r->particles[i].ax-ad[3*i], ey=r->particles[i].ay-ad[3*i+1], ez=r->particles[i].az-ad[3*i+2];
                double an = sqrt(ad[3*i]*ad[3*i]+ad[3*i+1]*ad[3*i+1]+ad[3*i+2]*ad[3*i+2]);
                err[i] = sqrt(ex*ex+ey*ey+ez*ez)/an;
            }
            // skip the star in the disc case (its force is dominated by cancellation)
            size_t i0 = (kind==0) ? 1 : 0;
            int cmp(const void* a, const void* b){ double x=*(double*)a, y=*(double*)b; return x<y?-1:x>y; }
            qsort(err+i0, N-i0, sizeof(double), cmp);
            printf("%s %s theta=%.1f  time=%7.2f ms  err_median=%.2e  err_99=%.2e\n", kind?"sheet":"disc", mac?"bmax ":"width", thetas[t], dt*1e3, err[i0+(N-i0)/2], err[i0+(size_t)(0.99*(N-i0))]);
        }
    }
    return 0;
}
