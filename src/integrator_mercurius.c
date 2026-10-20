/**
 * integrator_mercurius.c:  The hybrid symplectic MERCURIUS integrator.
 * 
 * Copyright (c) 2019 Hanno Rein, Dan Tamayo
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "rebound.h"
#include "rebound_internal.h"
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "gravity.h"
#include "integrator_mercurius.h"
#include "integrator_whfast.h"
#include "collision.h"
#include "binarydata.h"
#define MIN(a, b) ((a) > (b) ? (b) : (a))    ///< Returns the minimum of a and b
#define MAX(a, b) ((a) > (b) ? (a) : (b))    ///< Returns the maximum of a and b


void reb_integrator_mercurius_step(struct reb_simulation* r, void* state);
void* reb_integrator_mercurius_create();
void reb_integrator_mercurius_free(void* state);
void reb_integrator_mercurius_synchronize(struct reb_simulation* r, void* state);
void reb_integrator_mercurius_did_add_particle(struct reb_simulation* r);
void reb_integrator_mercurius_will_remove_particle(struct reb_simulation* r, size_t index);
const struct reb_binarydata_field_descriptor reb_integrator_mercurius_field_descriptor_list[];


const struct reb_integrator reb_integrator_mercurius = {
    .documentation = 
    "MERCURIUS is a hybrid symplectic integrator very similar to MERCURY "
    "by [Chambers (1999)]. It uses WHFast for long term integrations but "
    "switches over smoothly to IAS15 for close encounters. The MERCURIUS "
    "implementation is described in [Rein et al (2019)]."
    "\n\n"
    "[Rein et al (2019)]: https://ui.adsabs.harvard.edu/abs/2019MNRAS.485.5490R/abstract\n"
    "[Chambers (1999)]: https://ui.adsabs.harvard.edu/abs/1999MNRAS.304..793C/abstract\n" 
    ,
    .step = reb_integrator_mercurius_step,
    .create = reb_integrator_mercurius_create,
    .free = reb_integrator_mercurius_free,
    .synchronize = reb_integrator_mercurius_synchronize,
    .did_add_particle = reb_integrator_mercurius_did_add_particle,
    .will_remove_particle = reb_integrator_mercurius_will_remove_particle,
    .field_descriptor_list = reb_integrator_mercurius_field_descriptor_list,
};

const struct reb_binarydata_field_descriptor reb_integrator_mercurius_field_descriptor_list[] = {
    { "The critical switchover radii of particles are calculated automatically "
        "based on multiple criteria. One criterion calculates the Hill radius of "
        "particles and then multiplies it with the `r_crit_hill` parameter. The "
        "parameter is in units of the Hill radius. The default value is 3.",
        REB_DOUBLE,      "r_crit_hill",     offsetof(struct reb_integrator_mercurius_state, r_crit_hill), 0, 0, 0},
    { "If this flag is set to 1 (the default), the integrator will recalculate "
        "heliocentric coordinates and synchronize after every timestep to avoid "
        "problems with outputs or particle modifications between timesteps. Setting "
        "this flag to 0 will result in a speedup, but care must be taken to synchronize "
        "and recalculate coordinates manually if needed.",
        REB_UINT,        "safe_mode",       offsetof(struct reb_integrator_mercurius_state, safe_mode), 0, 0, 0},
    { "This is a function pointer to the force switching function. " 
        "By default the MERCURY switching function will be used. "
        "The argument `d` is the distance between two particles. "
        "The argument `dcrit` is the maximum of the critical distances of the two particles. "
        "The return value is a scalar between 0 and 1. "
        "If this function always returns 1, then the integrator effectively becomes the standard Wisdom-Holman integrator.",
        REB_FUNCTIONPOINTER,"L",            offsetof(struct reb_integrator_mercurius_state, L), 0, 0, 0},
    { "", REB_POINTER,     "dcrit",           offsetof(struct reb_integrator_mercurius_state, dcrit), offsetof(struct reb_integrator_mercurius_state, N_allocated_dcrit), sizeof(double), 0},
    { "", REB_VEC3D,       "com_pos",         offsetof(struct reb_integrator_mercurius_state, com_pos), 0, 0, 0},
    { "", REB_VEC3D,       "com_vel",         offsetof(struct reb_integrator_mercurius_state, com_vel), 0, 0, 0},
    { 0 }, // Null terminated list
};


void* reb_integrator_mercurius_create(){
    struct reb_integrator_mercurius_state* mercurius = calloc(sizeof(struct reb_integrator_mercurius_state),1);
    mercurius->mode = REB_INTEGRATOR_MERCURIUS_MODE_WH;
    mercurius->r_crit_hill = 3;
    mercurius->tponly_encounter = 0;
    mercurius->safe_mode = 1;
    mercurius->L = reb_integrator_mercurius_L_mercury;
    return mercurius;
}

void reb_integrator_mercurius_free(void* p){
    struct reb_integrator_mercurius_state* mercurius = p;
    free(mercurius->particles_backup);
    free(mercurius->particles_backup_xv);
    free(mercurius->sap_box);
    free(mercurius->sap_order);
    free(mercurius->sap_pairs);
    free(mercurius->enc_gather);
    free(mercurius->particles_backup_additional_forces);
    free(mercurius->encounter_map);
    free(mercurius->dcrit);
    free(mercurius);
}


void reb_integrator_mercurius_calculate_acceleration_mode_encounter(struct reb_simulation* r);

double reb_integrator_mercurius_L_mercury(const struct reb_simulation* const r, double d, double dcrit){
    // This is the changeover function used by the Mercury integrator.
    double y = (d-0.1*dcrit)/(0.9*dcrit);
    if (y<0.){
        return 0.;
    }else if (y>1.){
        return 1.;
    }else{
        return 10.*(y*y*y) - 15.*(y*y*y*y) + 6.*(y*y*y*y*y);
    }
    (void)r; // The simulation object is not used here.
}

double reb_integrator_mercurius_L_C4(const struct reb_simulation* const r, double d, double dcrit){
    // This is the changeover function C4 proposed by Hernandez (2019)
    double y = (d-0.1*dcrit)/(0.9*dcrit);
    if (y<0.){
        return 0.;
    }else if (y>1.){
        return 1.;
    }else{
        return (70.*y*y*y*y -315.*y*y*y +540.*y*y -420.*y +126.)*y*y*y*y*y;
    }
    (void)r; // The simulation object is not used here.
}

double reb_integrator_mercurius_L_C5(const struct reb_simulation* const r, double d, double dcrit){
    // This is the changeover function C5 proposed by Hernandez (2019)
    double y = (d-0.1*dcrit)/(0.9*dcrit);
    if (y<0.){
        return 0.;
    }else if (y>1.){
        return 1.;
    }else{
        return (-252.*y*y*y*y*y +1386.*y*y*y*y -3080.*y*y*y +3465.*y*y -1980.*y +462.)*y*y*y*y*y*y;
    }
    (void)r; // The simulation object is not used here.
}

static double f(double x){
    if (x<0) return 0;
    return exp(-1./x);
}

double reb_integrator_mercurius_L_infinity(const struct reb_simulation* const r, double d, double dcrit){
    // Infinitely differentiable function.
    double y = (d-0.1*dcrit)/(0.9*dcrit);
    if (y<0.){
        return 0.;
    }else if (y>1.){
        return 1.;
    }else{
        return f(y) /(f(y) + f(1.-y));
    }
    (void)r; // The simulation object is not used here.
}


void reb_integrator_mercurius_inertial_to_dh(struct reb_simulation* r, struct reb_integrator_mercurius_state* mercurius){
    struct reb_particle* restrict const particles = r->particles;
    struct reb_vec3d com_pos = {0};
    struct reb_vec3d com_vel = {0};
    double mtot = 0.;
    const size_t N_active = (r->N_active==SIZE_MAX || r->testparticle_type==1)?r->N:r->N_active;
    const size_t N = r->N;
    for (size_t i=0;i<N_active;i++){
        double m = particles[i].m;
        com_pos.x += m * particles[i].x;
        com_pos.y += m * particles[i].y;
        com_pos.z += m * particles[i].z;
        com_vel.x += m * particles[i].vx;
        com_vel.y += m * particles[i].vy;
        com_vel.z += m * particles[i].vz;
        mtot += m; 
    }
    com_pos.x /= mtot; com_pos.y /= mtot; com_pos.z /= mtot;
    com_vel.x /= mtot; com_vel.y /= mtot; com_vel.z /= mtot;
    // Particle 0 is also changed to allow for easy collision detection
    struct reb_particle p0 = particles[0];
    for (size_t i=0;i<N;i++){ 
        particles[i].x -= p0.x;
        particles[i].y -= p0.y;
        particles[i].z -= p0.z;
        particles[i].vx -= com_vel.x;
        particles[i].vy -= com_vel.y;
        particles[i].vz -= com_vel.z;
    }
    mercurius->com_pos = com_pos;
    mercurius->com_vel = com_vel;
}

void reb_integrator_mercurius_dh_to_inertial(struct reb_simulation* r, struct reb_integrator_mercurius_state* mercurius){
    struct reb_particle* restrict const particles = r->particles;
    struct reb_particle temp = {0};
    const size_t N = r->N;
    const size_t N_active = (r->N_active==SIZE_MAX || r->testparticle_type==1)?r->N:r->N_active;
    for (size_t i=1;i<N_active;i++){
        double m = particles[i].m;
        temp.x += m * particles[i].x;
        temp.y += m * particles[i].y;
        temp.z += m * particles[i].z;
        temp.vx += m * particles[i].vx;
        temp.vy += m * particles[i].vy;
        temp.vz += m * particles[i].vz;
        temp.m += m;
    }
    temp.m += r->particles[0].m;
    temp.x /= temp.m; 
    temp.y /= temp.m;
    temp.z /= temp.m;
    temp.vx /= particles[0].m; 
    temp.vy /= particles[0].m;
    temp.vz /= particles[0].m;
    // Use com to calculate central object's position.
    // This ignores previous values stored in particles[0].
    // Should not matter unless collisions occurred.
    particles[0].x = mercurius->com_pos.x - temp.x; 
    particles[0].y = mercurius->com_pos.y - temp.y; 
    particles[0].z = mercurius->com_pos.z - temp.z; 

    for (size_t i=1;i<N;i++){
        particles[i].x += particles[0].x;
        particles[i].y += particles[0].y;
        particles[i].z += particles[0].z;
        particles[i].vx += mercurius->com_vel.x;
        particles[i].vy += mercurius->com_vel.y;
        particles[i].vz += mercurius->com_vel.z;
    }
    particles[0].vx = mercurius->com_vel.x - temp.vx; 
    particles[0].vy = mercurius->com_vel.y - temp.vy; 
    particles[0].vz = mercurius->com_vel.z - temp.vz; 
}


// Close-encounter test for one pair (i < j), unchanged from the original double loop.
static inline void reb_mercurius_encounter_predict_pair(struct reb_integrator_mercurius_state* const mercurius, const struct reb_particle* const particles, const struct reb_vec6d* const particles_backup, const double* const dcrit, const double dt, const size_t N_active, const size_t i, const size_t j){
    const double dxn = particles[i].x - particles[j].x;
    const double dyn = particles[i].y - particles[j].y;
    const double dzn = particles[i].z - particles[j].z;
    const double dvxn = particles[i].vx - particles[j].vx;
    const double dvyn = particles[i].vy - particles[j].vy;
    const double dvzn = particles[i].vz - particles[j].vz;
    const double rn = (dxn*dxn + dyn*dyn + dzn*dzn);
    const double dxo = particles_backup[i].x - particles_backup[j].x;
    const double dyo = particles_backup[i].y - particles_backup[j].y;
    const double dzo = particles_backup[i].z - particles_backup[j].z;
    const double dvxo = particles_backup[i].vx - particles_backup[j].vx;
    const double dvyo = particles_backup[i].vy - particles_backup[j].vy;
    const double dvzo = particles_backup[i].vz - particles_backup[j].vz;
    const double ro = (dxo*dxo + dyo*dyo + dzo*dzo);

    const double drndt = (dxn*dvxn+dyn*dvyn+dzn*dvzn)*2.;
    const double drodt = (dxo*dvxo+dyo*dvyo+dzo*dvzo)*2.;

    double dcritmax2 = MAX(dcrit[i],dcrit[j]);
    dcritmax2 *= 1.21*dcritmax2;
    // Conservative early exit: the cubic Hermite interpolant of r^2 used below is
    // h00*ro + h01*rn + h10*dt*drodt + h11*dt*drndt with h00+h01=1 (both in [0,1]) and
    // |h10|,|h11| <= 4/27 on [0,1], so it never drops below this bound. If the bound is
    // above twice the threshold, the test at the end cannot succeed.
    const double r2_lower = MIN(rn,ro) - 0.15*(fabs(dt*drodt) + fabs(dt*drndt));
    if (r2_lower > 2.*dcritmax2) return;

    const double a = 6.*(ro-rn)+3.*dt*(drodt+drndt); 
    const double b = 6.*(rn-ro)-2.*dt*(2.*drodt+drndt); 
    const double c = dt*drodt; 

    double rmin = MIN(rn,ro);

    const double s = b*b-4.*a*c;
    const double sr = sqrt(MAX(0.,s));
    const double tmin1 = (-b + sr)/(2.*a); 
    const double tmin2 = (-b - sr)/(2.*a); 
    if (tmin1>0. && tmin1<1.){
        const double rmin1 = (1.-tmin1)*(1.-tmin1)*(1.+2.*tmin1)*ro
        + tmin1*tmin1*(3.-2.*tmin1)*rn
        + tmin1*(1.-tmin1)*(1.-tmin1)*dt*drodt
        - tmin1*tmin1*(1.-tmin1)*dt*drndt;
        rmin = MIN(MAX(rmin1,0.),rmin);
    }
    if (tmin2>0. && tmin2<1.){
        const double rmin2 = (1.-tmin2)*(1.-tmin2)*(1.+2.*tmin2)*ro
        + tmin2*tmin2*(3.-2.*tmin2)*rn
        + tmin2*(1.-tmin2)*(1.-tmin2)*dt*drodt
        - tmin2*tmin2*(1.-tmin2)*dt*drndt;
        rmin = MIN(MAX(rmin2,0.),rmin);
    }

    if (rmin < dcritmax2){
        if (mercurius->encounter_map[i]==0){
            mercurius->encounter_map[i] = i;
            mercurius->encounter_N++;
        }
        if (mercurius->encounter_map[j]==0){
            mercurius->encounter_map[j] = j;
            mercurius->encounter_N++;
        }
        if (j<N_active){ // Two massive particles have a close encounter
            mercurius->tponly_encounter = 0;
        }
    }
}

// Broadphase (default for many active bodies): the Hermite cubic f(t) of the squared separation
// used below is bounded from below by its Bernstein control values
//   b0 = ro, b1 = ro + dt*drodt/3, b2 = rn - dt*drndt/3, b3 = rn,
// and b1 = |P1_i - P1_j|^2 - (dt/3)^2 |v0_i - v0_j|^2 with P1 = x0 + v0 dt/3 (b2 likewise with
// P2 = x1 - v1 dt/3). Each particle gets the bounding box of x0, P1, P2, x1, padded by
// p = 1.1 (1+1e-6) dcrit + (1+1e-6) (dt/3) max(|v0|,|v1|) + 1e-7 (|x0|+|x1|). If two boxes are
// disjoint along one axis, every control value exceeds (1.1 max(dcrit_i,dcrit_j))^2 = dcritmax2,
// so the pair test cannot trigger. The outcome (flags and counts) does not depend on the order in
// which pairs are tested, so the results are bitwise identical.
#define REB_MERCURIUS_SAP_MIN_ACTIVE 32
static size_t reb_mercurius_encounter_broadphase(struct reb_simulation* const r, struct reb_integrator_mercurius_state* const mercurius, const size_t N_active){
    const size_t N = r->N;
    const struct reb_particle* const ps = r->particles;
    const struct reb_vec6d* const b0 = mercurius->particles_backup_xv;
    const double* const dcrit = mercurius->dcrit;
    const double dt3 = r->dt/3.;
    if (mercurius->N_allocated_sap < N){
        mercurius->sap_box = realloc(mercurius->sap_box, sizeof(double)*6*N);
        mercurius->sap_order = realloc(mercurius->sap_order, sizeof(size_t)*N);
        mercurius->N_allocated_sap = N;
        mercurius->N_sap_order = 0;
    }
    double* const box = mercurius->sap_box; // xmin, xmax, ymin, ymax, zmin, zmax
    size_t* const order = mercurius->sap_order;
    for (size_t k=0; k<N; k++){
        const double x0[3] = {b0[k].x, b0[k].y, b0[k].z};
        const double v0[3] = {b0[k].vx, b0[k].vy, b0[k].vz};
        const double x1[3] = {ps[k].x, ps[k].y, ps[k].z};
        const double v1[3] = {ps[k].vx, ps[k].vy, ps[k].vz};
        const double vm0 = sqrt(v0[0]*v0[0]+v0[1]*v0[1]+v0[2]*v0[2]);
        const double vm1 = sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);
        const double R = sqrt(x0[0]*x0[0]+x0[1]*x0[1]+x0[2]*x0[2]) + sqrt(x1[0]*x1[0]+x1[1]*x1[1]+x1[2]*x1[2]);
        const double pad = 1.1*(1.+1e-6)*dcrit[k] + (1.+1e-6)*dt3*(vm0 > vm1 ? vm0 : vm1) + 1e-7*R;
        for (int a=0; a<3; a++){
            const double p1 = x0[a] + v0[a]*dt3;
            const double p2 = x1[a] - v1[a]*dt3;
            double lo = x0[a], hi = x0[a];
            if (p1 < lo) lo = p1;
            if (p1 > hi) hi = p1;
            if (p2 < lo) lo = p2;
            if (p2 > hi) hi = p2;
            if (x1[a] < lo) lo = x1[a];
            if (x1[a] > hi) hi = x1[a];
            box[6*k+2*a] = lo - pad;
            box[6*k+2*a+1] = hi + pad;
        }
    }
    if (mercurius->N_sap_order != N){
        for (size_t k=0; k<N; k++) order[k] = k;
        mercurius->N_sap_order = N;
    }
    for (size_t a=1; a<N; a++){ // insertion sort by xmin (nearly sorted)
        const size_t ka = order[a];
        const double la = box[6*ka];
        size_t b = a;
        while (b>0 && box[6*order[b-1]] > la){
            order[b] = order[b-1];
            b--;
        }
        order[b] = ka;
    }
    size_t n = 0;
    for (size_t a=0; a<N; a++){
        const size_t i = order[a];
        const double* const bi = &box[6*i];
        for (size_t b=a+1; b<N; b++){
            const size_t j = order[b];
            const double* const bj = &box[6*j];
            if (bj[0] > bi[1]) break;
            if (i >= N_active && j >= N_active) continue;
            if (bj[2] > bi[3] || bi[2] > bj[3]) continue;
            if (bj[4] > bi[5] || bi[4] > bj[5]) continue;
            if (2*(n+1) > mercurius->N_allocated_sap_pairs){
                mercurius->N_allocated_sap_pairs = mercurius->N_allocated_sap_pairs ? 2*mercurius->N_allocated_sap_pairs : 1024;
                mercurius->sap_pairs = realloc(mercurius->sap_pairs, sizeof(size_t)*mercurius->N_allocated_sap_pairs);
            }
            mercurius->sap_pairs[2*n]   = i < j ? i : j;
            mercurius->sap_pairs[2*n+1] = i < j ? j : i;
            n++;
        }
    }
    return n;
}

static void reb_mercurius_encounter_predict(struct reb_simulation* const r, struct reb_integrator_mercurius_state* mercurius){
    // This function predicts close encounters during the timestep
    // It makes use of the old and new position and velocities obtained
    // after the Kepler step.
    struct reb_particle* const particles = r->particles;
    const struct reb_vec6d* const particles_backup = mercurius->particles_backup_xv;
    const double* const dcrit = mercurius->dcrit;
    const size_t N = r->N;
    const size_t N_active = r->N_active==SIZE_MAX?r->N:r->N_active;
    const double dt = r->dt;
    mercurius->encounter_N = 1;
    mercurius->encounter_map[0] = 1;
    if (r->testparticle_type==1){
        mercurius->tponly_encounter = 0; // testparticles affect massive particles
    }else{
        mercurius->tponly_encounter = 1;
    }
    for (size_t i=1; i<N; i++){
        mercurius->encounter_map[i] = 0;
    }
    if (N_active >= REB_MERCURIUS_SAP_MIN_ACTIVE){
        const size_t n = reb_mercurius_encounter_broadphase(r, mercurius, N_active);
        for (size_t c=0; c<n; c++){
            reb_mercurius_encounter_predict_pair(mercurius, particles, particles_backup, dcrit, dt, N_active, mercurius->sap_pairs[2*c], mercurius->sap_pairs[2*c+1]);
        }
        return;
    }
    for (size_t i=0; i<N_active; i++){
        for (size_t j=i+1; j<N; j++){
            reb_mercurius_encounter_predict_pair(mercurius, particles, particles_backup, dcrit, dt, N_active, i, j);
        }
    }
}

void reb_integrator_mercurius_interaction_step(struct reb_simulation* const r, double dt){
    struct reb_particle* restrict const particles = r->particles;
    const size_t N = r->N;
    for (size_t i=1;i<N;i++){
        particles[i].vx += dt*particles[i].ax;
        particles[i].vy += dt*particles[i].ay;
        particles[i].vz += dt*particles[i].az;
    }
}

void reb_integrator_mercurius_jump_step(struct reb_simulation* const r, double dt){
    struct reb_particle* restrict const particles = r->particles;
    const size_t N_active = r->N_active==SIZE_MAX?r->N:r->N_active;
    const size_t N = r->testparticle_type==0 ? N_active: r->N;
    double px=0., py=0., pz=0.;
    for (size_t i=1;i<N;i++){
        px += r->particles[i].vx*r->particles[i].m; // in dh
        py += r->particles[i].vy*r->particles[i].m; 
        pz += r->particles[i].vz*r->particles[i].m;
    }
    px /= r->particles[0].m;
    py /= r->particles[0].m;
    pz /= r->particles[0].m;
    const size_t N_all = r->N;
    for (size_t i=1;i<N_all;i++){
        particles[i].x += dt*px;
        particles[i].y += dt*py;
        particles[i].z += dt*pz;
    }
}

static void reb_integrator_mercurius_kepler_step(struct reb_simulation* const r, double dt){
    struct reb_particle* restrict const particles = r->particles;
    const size_t N = r->N;
    reb_integrator_whfast_kepler_solver_range(particles, 1, N, r->G*particles[0].m, dt, r); // in dh
}

static void reb_mercurius_encounter_step(struct reb_simulation* const r, struct reb_integrator_mercurius_state* mercurius, const double _dt){
    // Only particles having a close encounter are integrated by IAS15.
    if (mercurius->encounter_N<2){
        return; // If there are no particles (other than the star) having a close encounter, then there is nothing to do.
    }
    size_t N_active = (r->N_active==SIZE_MAX)?r->N:r->N_active;
    size_t i_enc = 0;
    mercurius->encounter_N_active = 0;
    for (size_t i=0; i<r->N; i++){
        if(mercurius->encounter_map[i]){  
            struct reb_particle tmp = r->particles[i];      // Copy for potential use for tponly_encounter
            const struct reb_vec6d xv = mercurius->particles_backup_xv[i]; // Use coordinates before whfast step
            r->particles[i].x = xv.x; r->particles[i].y = xv.y; r->particles[i].z = xv.z;
            r->particles[i].vx = xv.vx; r->particles[i].vy = xv.vy; r->particles[i].vz = xv.vz;
            mercurius->encounter_map[i_enc] = i;
            i_enc++;
            if (i<N_active){
                mercurius->encounter_N_active++;
                if (mercurius->tponly_encounter){
                    mercurius->particles_backup[i] = tmp;         // Make copy of particles after the kepler step.
                                                                  // used to restore the massive objects' states in the case
                                                                  // of only massless test-particle encounters
                }
            }
        }
    }

    mercurius->mode = REB_INTEGRATOR_MERCURIUS_MODE_ENCOUNTER;

    // run
    const double old_dt = r->dt;
    const double dtsign = old_dt>=0.?1.:-1.;
    const double old_t = r->t;
    double t_needed = r->t + _dt; 

    struct reb_integrator_ias15_state* ias15 = reb_integrator_ias15.create();
    r->map = mercurius->encounter_map;
    r->N_map = mercurius->encounter_N;
    r->N_targets = SIZE_MAX; // Search for any possible collisions between N_map particles.

    r->dt = 0.0001*_dt; // start with a small timestep.

    // No additional forces during encounter
    r->gravity = REB_GRAVITY_CUSTOM;
    r->gravity_custom = reb_integrator_mercurius_calculate_acceleration_mode_encounter;

    while(dtsign*r->t < dtsign*t_needed && fabs(r->dt/old_dt)>1e-14 && r->status<=0){
        struct reb_particle star = r->particles[0]; // backup velocity
        r->particles[0].vx = 0; // star does not move in dh 
        r->particles[0].vy = 0;
        r->particles[0].vz = 0;
        reb_integrator_ias15.step(r, ias15);
        r->particles[0].vx = star.vx; // restore every timestep for collisions
        r->particles[0].vy = star.vy;
        r->particles[0].vz = star.vz;

        if (dtsign*(r->t+r->dt) > dtsign*t_needed){
            r->dt = t_needed-r->t;
        }

        // Search and resolve collisions
        reb_collision_search(r);

        struct reb_particle p0 = r->particles[0];
        star.vx = p0.vx; // keep track of changed star velocity for later collisions
        star.vy = p0.vy;
        star.vz = p0.vz;
        if (p0.x !=0 || p0.y !=0 || p0.z !=0){
            // Collision with star occurred
            // Shift all particles back to heliocentric coordinates
            // Ignore stars velocity:
            //   - will not be used after this
            //   - com velocity is unchanged. this velocity will be used
            //     to reconstruct star's velocity later.
            for (size_t i=0; i<r->N; i++){
                r->particles[i].x -= p0.x;
                r->particles[i].y -= p0.y;
                r->particles[i].z -= p0.z;
            }
        }
    }

    // if only test particles encountered massive bodies, reset the
    // massive body coordinates to their post Kepler step state
    if(mercurius->tponly_encounter){
        for (size_t i=1;i<mercurius->encounter_N_active;i++){
            size_t mi = mercurius->encounter_map[i];
            r->particles[mi] = mercurius->particles_backup[mi];
        }
    }

    reb_integrator_ias15.free(ias15);

    r->t = old_t;
    r->dt = old_dt;
    mercurius->mode = REB_INTEGRATOR_MERCURIUS_MODE_WH;
    r->map = NULL;
    r->N_map = 0;
}

double reb_integrator_mercurius_calculate_dcrit_for_particle(struct reb_simulation* r, struct reb_integrator_mercurius_state* mercurius, size_t i){
    const double m0 = r->particles[0].m;
    const double dx  = r->particles[i].x;  // in dh
    const double dy  = r->particles[i].y;
    const double dz  = r->particles[i].z;
    const double dvx = r->particles[i].vx - r->particles[0].vx; 
    const double dvy = r->particles[i].vy - r->particles[0].vy; 
    const double dvz = r->particles[i].vz - r->particles[0].vz; 
    const double _r = sqrt(dx*dx + dy*dy + dz*dz);
    const double v2 = dvx*dvx + dvy*dvy + dvz*dvz;

    const double GM = r->G*(m0+r->particles[i].m);
    const double a = GM*_r / (2.*GM - _r*v2);
    const double vc = sqrt(GM/fabs(a));
    double dcrit = 0;
    // Criteria 1: average velocity
    dcrit = MAX(dcrit, vc*0.4*r->dt);
    // Criteria 2: current velocity
    dcrit = MAX(dcrit, sqrt(v2)*0.4*r->dt);
    // Criteria 3: Hill radius
    dcrit = MAX(dcrit, mercurius->r_crit_hill*a*cbrt(r->particles[i].m/(3.*r->particles[0].m)));
    // Criteria 4: physical radius
    dcrit = MAX(dcrit, 2.*r->particles[i].r);
    return dcrit;
}

#define REB_MERCURIUS_BLOCK 32
void reb_integrator_mercurius_calculate_acceleration_mode_encounter(struct reb_simulation* r){
    struct reb_integrator_mercurius_state* mercurius = r->integrator.state;
    struct reb_particle* const particles = r->particles;
    const double G = r->G;
    const double softening2 = r->softening*r->softening;
    const int _testparticle_type   = r->testparticle_type;
    double (*_L) (const struct reb_simulation* const r, double d, double dcrit) = mercurius->L;
    const double m0 = r->particles[0].m;
    const double* const dcrit = mercurius->dcrit;
    const size_t encounter_N = mercurius->encounter_N;
    const size_t encounter_N_active = mercurius->encounter_N_active;
    size_t* map = mercurius->encounter_map;
    // Pairs farther apart than the changeover region have L = 1 exactly for the default
    // changeover function (see the WH-mode acceleration), so their force (1-L) is exactly zero.
    const int L_default = (_L == reb_integrator_mercurius_L_mercury);
#ifdef OPENMP
    // Serial code unless the parallel loop below pays off (same results, see rebound_internal.h).
    if (!reb_omp_parallel_onesided((double)encounter_N*(double)encounter_N_active))
#endif // OPENMP
    {
    particles[0].ax = 0; // map[0] is always 0 
    particles[0].ay = 0; 
    particles[0].az = 0; 
    // Acceleration due to star
    for (size_t i=1; i<encounter_N; i++){
        size_t mi = map[i];
        const double x = particles[mi].x;
        const double y = particles[mi].y;
        const double z = particles[mi].z;
        const double _r = sqrt(x*x + y*y + z*z + softening2);
        double prefact = -G/(_r*_r*_r)*m0;
        particles[mi].ax    = prefact*x;
        particles[mi].ay    = prefact*y;
        particles[mi].az    = prefact*z;
    }
    // We're in a heliocentric coordinate system.
    // The star feels no acceleration
    // Contiguous copies of the encounter particles' positions and critical radii (same values).
    const int gathered = L_default && encounter_N >= 16;
    double *gx = NULL, *gy = NULL, *gz = NULL, *gdc = NULL;
    if (gathered){
        if (mercurius->N_allocated_enc_gather < encounter_N){
            mercurius->enc_gather = realloc(mercurius->enc_gather, sizeof(double)*4*encounter_N);
            mercurius->N_allocated_enc_gather = encounter_N;
        }
        gx = mercurius->enc_gather;
        gy = gx + encounter_N;
        gz = gy + encounter_N;
        gdc = gz + encounter_N;
        for (size_t k=0; k<encounter_N; k++){
            const size_t mk = map[k];
            gx[k] = particles[mk].x;
            gy[k] = particles[mk].y;
            gz[k] = particles[mk].z;
            gdc[k] = dcrit[mk];
        }
    }
    // Interactions between active-active
    for (size_t i=2; i<encounter_N_active; i++){
        size_t mi = map[i];
        if (gathered){
            // Far-pair test on contiguous copies, a block at a time (vectorisable); the remaining
            // pairs are handled as below, in the original order.
            const double xi = gx[i], yi = gy[i], zi = gz[i], dci = gdc[i];
            for (size_t j0=1; j0<i; j0+=REB_MERCURIUS_BLOCK){
                const size_t nb = (i-j0) < REB_MERCURIUS_BLOCK ? (i-j0) : REB_MERCURIUS_BLOCK;
                double bdx[REB_MERCURIUS_BLOCK], bdy[REB_MERCURIUS_BLOCK], bdz[REB_MERCURIUS_BLOCK], br2[REB_MERCURIUS_BLOCK];
                int near[REB_MERCURIUS_BLOCK];
                int any = 0;
                for (size_t k=0; k<nb; k++){
                    const double dx = xi - gx[j0+k];
                    const double dy = yi - gy[j0+k];
                    const double dz = zi - gz[j0+k];
                    const double dcritmax = MAX(dci,gdc[j0+k]);
                    const double r2 = dx*dx + dy*dy + dz*dz + softening2;
                    bdx[k] = dx; bdy[k] = dy; bdz[k] = dz; br2[k] = r2;
                    near[k] = !(r2 > 1.0002000100*dcritmax*dcritmax);
                    any |= near[k];
                }
                if (!any) continue;
                for (size_t k=0; k<nb; k++){
                    if (!near[k]) continue;
                    const size_t mj = map[j0+k];
                    const double dx = bdx[k], dy = bdy[k], dz = bdz[k];
                    const double dcritmax = MAX(dci,gdc[j0+k]);
                    const double _r = sqrt(br2[k]);
                    const double L = _L(r,_r,dcritmax);
                    double prefact = G*(1.-L)/(_r*_r*_r);
                    double prefactj = -prefact*particles[mj].m;
                    double prefacti = prefact*particles[mi].m;
                    particles[mi].ax    += prefactj*dx;
                    particles[mi].ay    += prefactj*dy;
                    particles[mi].az    += prefactj*dz;
                    particles[mj].ax    += prefacti*dx;
                    particles[mj].ay    += prefacti*dy;
                    particles[mj].az    += prefacti*dz;
                }
            }
            continue;
        }
        for (size_t j=1; j<i; j++){
            size_t mj = map[j];
            const double dx = particles[mi].x - particles[mj].x;
            const double dy = particles[mi].y - particles[mj].y;
            const double dz = particles[mi].z - particles[mj].z;
            const double dcritmax = MAX(dcrit[mi],dcrit[mj]);
            const double r2 = dx*dx + dy*dy + dz*dz + softening2;
            if (L_default && r2 > 1.0002000100*dcritmax*dcritmax) continue; // (1.0001 dcritmax)^2: no force
            const double _r = sqrt(r2);
            const double L = _L(r,_r,dcritmax);
            double prefact = G*(1.-L)/(_r*_r*_r);
            double prefactj = -prefact*particles[mj].m;
            double prefacti = prefact*particles[mi].m;
            particles[mi].ax    += prefactj*dx;
            particles[mi].ay    += prefactj*dy;
            particles[mi].az    += prefactj*dz;
            particles[mj].ax    += prefacti*dx;
            particles[mj].ay    += prefacti*dy;
            particles[mj].az    += prefacti*dz;
        }
    }
    // Interactions between active-testparticle
    const size_t startitestp = MAX(encounter_N_active,2);
    for (size_t i=startitestp; i<encounter_N; i++){
        size_t mi = map[i];
        for (size_t j=1; j<encounter_N_active; j++){
            size_t mj = map[j];
            const double dx = particles[mi].x - particles[mj].x;
            const double dy = particles[mi].y - particles[mj].y;
            const double dz = particles[mi].z - particles[mj].z;
            const double dcritmax = MAX(dcrit[mi],dcrit[mj]);
            const double r2 = dx*dx + dy*dy + dz*dz + softening2;
            if (L_default && r2 > 1.0002000100*dcritmax*dcritmax) continue; // (1.0001 dcritmax)^2: no force
            const double _r = sqrt(r2);
            const double L = _L(r,_r,dcritmax);
            double prefact = G*(1.-L)/(_r*_r*_r);
            double prefactj = -prefact*particles[mj].m;
            particles[mi].ax    += prefactj*dx;
            particles[mi].ay    += prefactj*dy;
            particles[mi].az    += prefactj*dz;
            if (_testparticle_type){
                double prefacti = prefact*particles[mi].m;
                particles[mj].ax    += prefacti*dx;
                particles[mj].ay    += prefacti*dy;
                particles[mj].az    += prefacti*dz;
            }
        }
    }
    }
#ifdef OPENMP
    else {
    // One-sided parallel loop, bitwise identical to the serial loop above: the serial half-pair
    // loop adds to each encounter particle the star term first, then its partners in ascending
    // order (own row j < k, then later rows i > k), then (type 1) the test-particle back-reaction;
    // a partner-row term (prefact*m_i)*(x_i - x_k) equals (-prefact*m_i)*(x_k - x_i) exactly.
    // The same far-pair test skips the same pairs.
    particles[0].ax = 0; // map[0] is always 0
    particles[0].ay = 0;
    particles[0].az = 0;
    const size_t startitestp = MAX(encounter_N_active,2);
#pragma omp parallel for schedule(guided) num_threads(reb_omp_threads())
    for (size_t k=1; k<encounter_N; k++){
        const size_t mk = map[k];
        const double xk = particles[mk].x, yk = particles[mk].y, zk = particles[mk].z;
        double ax, ay, az;
        {
            const double _r = sqrt(xk*xk + yk*yk + zk*zk + softening2);
            double prefact = -G/(_r*_r*_r)*m0;
            ax = prefact*xk; ay = prefact*yk; az = prefact*zk;
        }
        if (k<encounter_N_active || k<startitestp){
            if (k<encounter_N_active){
                for (size_t f=1; f<encounter_N_active; f++){
                    if (f==k) continue;
                    const size_t mf = map[f];
                    const double dx = xk - particles[mf].x;
                    const double dy = yk - particles[mf].y;
                    const double dz = zk - particles[mf].z;
                    const double dcritmax = MAX(dcrit[mk],dcrit[mf]);
                    const double r2 = dx*dx + dy*dy + dz*dz + softening2;
                    if (L_default && r2 > 1.0002000100*dcritmax*dcritmax) continue;
                    const double _r = sqrt(r2);
                    const double L = _L(r,_r,dcritmax);
                    double prefact = G*(1.-L)/(_r*_r*_r);
                    double prefactj = -prefact*particles[mf].m;
                    ax += prefactj*dx;
                    ay += prefactj*dy;
                    az += prefactj*dz;
                }
                if (_testparticle_type){
                    for (size_t f=startitestp; f<encounter_N; f++){
                        const size_t mf = map[f];
                        const double dx = particles[mf].x - xk;
                        const double dy = particles[mf].y - yk;
                        const double dz = particles[mf].z - zk;
                        const double dcritmax = MAX(dcrit[mf],dcrit[mk]);
                        const double r2 = dx*dx + dy*dy + dz*dz + softening2;
                        if (L_default && r2 > 1.0002000100*dcritmax*dcritmax) continue;
                        const double _r = sqrt(r2);
                        const double L = _L(r,_r,dcritmax);
                        double prefact = G*(1.-L)/(_r*_r*_r);
                        double prefacti = prefact*particles[mf].m;
                        ax += prefacti*dx;
                        ay += prefacti*dy;
                        az += prefacti*dz;
                    }
                }
            }
        }else{
            for (size_t f=1; f<encounter_N_active; f++){
                const size_t mf = map[f];
                const double dx = xk - particles[mf].x;
                const double dy = yk - particles[mf].y;
                const double dz = zk - particles[mf].z;
                const double dcritmax = MAX(dcrit[mk],dcrit[mf]);
                const double r2 = dx*dx + dy*dy + dz*dz + softening2;
                if (L_default && r2 > 1.0002000100*dcritmax*dcritmax) continue;
                const double _r = sqrt(r2);
                const double L = _L(r,_r,dcritmax);
                double prefact = G*(1.-L)/(_r*_r*_r);
                double prefactj = -prefact*particles[mf].m;
                ax += prefactj*dx;
                ay += prefactj*dy;
                az += prefactj*dz;
            }
        }
        particles[mk].ax = ax; particles[mk].ay = ay; particles[mk].az = az;
    }
    }
#endif // OPENMP
}

static void reb_integrator_mercurius_calculate_acceleration_mode_wh(struct reb_simulation* r, struct reb_integrator_mercurius_state* mercurius){
    struct reb_particle* const particles = r->particles;
    const size_t N = r->N;
    const size_t N_active = ((r->N_active==SIZE_MAX)?N:r->N_active);
    const double G = r->G;
    const double softening2 = r->softening*r->softening;
    const int _testparticle_type   = r->testparticle_type;
    double (*_L) (const struct reb_simulation* const r, double d, double dcrit) = mercurius->L;
    const double* const dcrit = mercurius->dcrit;
    // The default changeover function returns exactly 1 for d > dcrit; with a small margin
    // (y >= 1.0001 in reb_integrator_mercurius_L_mercury) the call can be skipped.
    const int L_default = (_L == reb_integrator_mercurius_L_mercury);
#ifdef OPENMP
    // Serial code unless the parallel loop below pays off (same results, see rebound_internal.h).
    if (!reb_omp_parallel_onesided((double)N*(double)N_active))
#endif // OPENMP
    {
    for (size_t i=0; i<N; i++){
        particles[i].ax = 0; 
        particles[i].ay = 0; 
        particles[i].az = 0; 
    }
    for (size_t i=2; i<N_active; i++){
        if (reb_sigint > 1) return;
        if (i-1 >= 16){
            // Blocks of partners: separations and distances first (vectorisable square roots),
            // then the changeover function and the accumulation in the original order.
            // Same operations as the plain loop below, so the results are identical.
            const double xi = particles[i].x, yi = particles[i].y, zi = particles[i].z, mi = particles[i].m, dci = dcrit[i];
            for (size_t j0=1; j0<i; j0+=REB_MERCURIUS_BLOCK){
                const size_t nb = (i-j0) < REB_MERCURIUS_BLOCK ? (i-j0) : REB_MERCURIUS_BLOCK;
                double bdx[REB_MERCURIUS_BLOCK], bdy[REB_MERCURIUS_BLOCK], bdz[REB_MERCURIUS_BLOCK], br[REB_MERCURIUS_BLOCK], bpf[REB_MERCURIUS_BLOCK];
                for (size_t k=0; k<nb; k++){
                    bdx[k] = xi - particles[j0+k].x;
                    bdy[k] = yi - particles[j0+k].y;
                    bdz[k] = zi - particles[j0+k].z;
                }
                for (size_t k=0; k<nb; k++){
                    const double _r = sqrt(bdx[k]*bdx[k] + bdy[k]*bdy[k] + bdz[k]*bdz[k] + softening2);
                    br[k] = _r;
                    bpf[k] = G/(_r*_r*_r); // = G*L/r^3 for L = 1 (G*1 is exactly G)
                }
                for (size_t k=0; k<nb; k++){
                    const size_t j = j0+k;
                    const double _r = br[k];
                    const double dcritmax = MAX(dci,dcrit[j]);
                    double prefact;
                    if (L_default && _r > 1.0001*dcritmax){
                        prefact = bpf[k];
                    }else{
                        const double L = _L(r,_r,dcritmax);
                        prefact = G*L/(_r*_r*_r);
                    }
                    const double prefactj = -prefact*particles[j].m;
                    const double prefacti = prefact*mi;
                    particles[i].ax    += prefactj*bdx[k];
                    particles[i].ay    += prefactj*bdy[k];
                    particles[i].az    += prefactj*bdz[k];
                    particles[j].ax    += prefacti*bdx[k];
                    particles[j].ay    += prefacti*bdy[k];
                    particles[j].az    += prefacti*bdz[k];
                }
            }
            continue;
        }
        for (size_t j=1; j<i; j++){
            const double dx = particles[i].x - particles[j].x;
            const double dy = particles[i].y - particles[j].y;
            const double dz = particles[i].z - particles[j].z;
            const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
            const double dcritmax = MAX(dcrit[i],dcrit[j]);
            const double L = (L_default && _r > 1.0001*dcritmax) ? 1. : _L(r,_r,dcritmax);
            const double prefact = G*L/(_r*_r*_r);
            const double prefactj = -prefact*particles[j].m;
            const double prefacti = prefact*particles[i].m;
            particles[i].ax    += prefactj*dx;
            particles[i].ay    += prefactj*dy;
            particles[i].az    += prefactj*dz;
            particles[j].ax    += prefacti*dx;
            particles[j].ay    += prefacti*dy;
            particles[j].az    += prefacti*dz;
        }
    }
    const size_t startitestp = MAX(N_active,2);
    for (size_t i=startitestp; i<N; i++){
        if (reb_sigint > 1) return;
        for (size_t j=1; j<N_active; j++){
            const double dx = particles[i].x - particles[j].x;
            const double dy = particles[i].y - particles[j].y;
            const double dz = particles[i].z - particles[j].z;
            const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
            const double dcritmax = MAX(dcrit[i],dcrit[j]);
            const double L = (L_default && _r > 1.0001*dcritmax) ? 1. : _L(r,_r,dcritmax);
            const double prefact = G*L/(_r*_r*_r);
            const double prefactj = -prefact*particles[j].m;
            particles[i].ax    += prefactj*dx;
            particles[i].ay    += prefactj*dy;
            particles[i].az    += prefactj*dz;
            if (_testparticle_type){
                const double prefacti = prefact*particles[i].m;
                particles[j].ax    += prefacti*dx;
                particles[j].ay    += prefacti*dy;
                particles[j].az    += prefacti*dz;
            }
        }
    }
    }
#ifdef OPENMP
    else {
    // One-sided parallel loop, bitwise identical to the serial loop above (same terms, same
    // order per particle: partners in ascending order, then the type-1 back-reaction of the test
    // particles; a partner-row term (prefact*m_i)*(x_i - x_k) equals (-prefact*m_i)*(x_k - x_i)).
    particles[0].ax = 0;
    particles[0].ay = 0;
    particles[0].az = 0;
    const size_t startitestp = MAX(N_active,2);
#pragma omp parallel for schedule(guided) num_threads(reb_omp_threads())
    for (size_t k=1; k<N; k++){
        const double xk = particles[k].x, yk = particles[k].y, zk = particles[k].z;
        double ax = 0., ay = 0., az = 0.;
        const int active = k<N_active;
        if (!active && k<startitestp){
            particles[k].ax = 0.; particles[k].ay = 0.; particles[k].az = 0.;
            continue;
        }
        for (size_t j=1; j<N_active; j++){
            if (j==k) continue;
            const double dx = xk - particles[j].x;
            const double dy = yk - particles[j].y;
            const double dz = zk - particles[j].z;
            const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
            const double dcritmax = MAX(dcrit[k],dcrit[j]);
            const double L = (L_default && _r > 1.0001*dcritmax) ? 1. : _L(r,_r,dcritmax);
            const double prefact = G*L/(_r*_r*_r);
            const double prefactj = -prefact*particles[j].m;
            ax += prefactj*dx;
            ay += prefactj*dy;
            az += prefactj*dz;
        }
        if (active && _testparticle_type){
            for (size_t i=startitestp; i<N; i++){
                const double dx = particles[i].x - xk;
                const double dy = particles[i].y - yk;
                const double dz = particles[i].z - zk;
                const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
                const double dcritmax = MAX(dcrit[i],dcrit[k]);
                const double L = (L_default && _r > 1.0001*dcritmax) ? 1. : _L(r,_r,dcritmax);
                const double prefact = G*L/(_r*_r*_r);
                const double prefacti = prefact*particles[i].m;
                ax += prefacti*dx;
                ay += prefacti*dy;
                az += prefacti*dz;
            }
        }
        particles[k].ax = ax; particles[k].ay = ay; particles[k].az = az;
    }
    }
#endif // OPENMP

    if (r->additional_forces){
        // Additional forces are only calculated in the kick step, not during close encounter
        // shift pos and velocity so that external forces are calculated in inertial frame
        // Note: Copying avoids degrading floating point performance
        if(r->N>mercurius->N_allocated_additional_forces){
            mercurius->particles_backup_additional_forces = realloc(mercurius->particles_backup_additional_forces, r->N*sizeof(struct reb_particle));
            mercurius->N_allocated_additional_forces = r->N;
        }
        memcpy(mercurius->particles_backup_additional_forces,r->particles,r->N*sizeof(struct reb_particle)); 
        reb_integrator_mercurius_dh_to_inertial(r, mercurius);

        r->additional_forces(r);

        struct reb_particle* restrict const particles = r->particles;
        struct reb_particle* restrict const backup = mercurius->particles_backup_additional_forces;
        for (size_t i=0;i<r->N;i++){
            particles[i].x = backup[i].x;
            particles[i].y = backup[i].y;
            particles[i].z = backup[i].z;
            particles[i].vx = backup[i].vx;
            particles[i].vy = backup[i].vy;
            particles[i].vz = backup[i].vz;
        }
    }
}

void reb_integrator_mercurius_did_add_particle(struct reb_simulation* r){
    struct reb_integrator_mercurius_state* mercurius = r->integrator.state;
    switch (mercurius->mode){
        case REB_INTEGRATOR_MERCURIUS_MODE_WH:
            // Nothing to do here. r->did_modify_particles will get set automatically
            break;
        case REB_INTEGRATOR_MERCURIUS_MODE_ENCOUNTER:
            if (mercurius->N_allocated_dcrit<r->N){
                mercurius->dcrit              = realloc(mercurius->dcrit, sizeof(double)*r->N);
                mercurius->N_allocated_dcrit = r->N;
            }
            mercurius->dcrit[r->N-1] = reb_integrator_mercurius_calculate_dcrit_for_particle(r, mercurius, r->N-1);
            if (mercurius->N_allocated<r->N){
                mercurius->particles_backup   = realloc(mercurius->particles_backup,sizeof(struct reb_particle)*r->N);
                mercurius->particles_backup_xv = realloc(mercurius->particles_backup_xv,sizeof(struct reb_vec6d)*r->N);
                mercurius->encounter_map      = realloc(mercurius->encounter_map,sizeof(size_t)*r->N);
                r->map = mercurius->encounter_map;
                mercurius->N_allocated = r->N;
            }
            mercurius->encounter_map[mercurius->encounter_N] = r->N-1;
            mercurius->encounter_N++;
            r->N_map++;
            if (r->N_active==SIZE_MAX){ 
                // If global N_active is not set, then all particles are active, so the new one as well.
                // Otherwise, assume we're adding non active particle. 
                mercurius->encounter_N_active++;
            }
            break;
    }
}

void reb_integrator_mercurius_will_remove_particle(struct reb_simulation* r, size_t index){
    struct reb_integrator_mercurius_state* mercurius = r->integrator.state;
    if (mercurius->N_allocated_dcrit>0 && index<mercurius->N_allocated_dcrit){
        for (size_t i=0;i<r->N-1;i++){
            if (i>=index){
                mercurius->dcrit[i] = mercurius->dcrit[i+1];
            }
        }
    }
    if (mercurius->mode==REB_INTEGRATOR_MERCURIUS_MODE_ENCOUNTER){
        int after_to_be_removed_particle = 0;
        size_t encounter_index = SIZE_MAX;
        for (size_t i=0;i<mercurius->encounter_N;i++){
            if (after_to_be_removed_particle == 1){
                mercurius->encounter_map[i-1] = mercurius->encounter_map[i] - 1; 
            }
            if (mercurius->encounter_map[i]==index){
                encounter_index = i;
                after_to_be_removed_particle = 1;
            }
        }
        if (encounter_index == SIZE_MAX){
            reb_simulation_error(r,"Cannot find particle in encounter map.");
            return;
        }
        if (encounter_index<mercurius->encounter_N_active){
            mercurius->encounter_N_active--;
        }
        mercurius->encounter_N--;
        r->N_map--;
    }
}

void reb_integrator_mercurius_step(struct reb_simulation* r, void* state){
    struct reb_integrator_mercurius_state* mercurius = state;
    if (r->N_var_config){
        reb_simulation_warning(r,"Mercurius does not work with variational equations.");
    }

    const size_t N = r->N;

    if (mercurius->N_allocated_dcrit<N){
        // Need to safe these arrays in Simulationarchive
        mercurius->dcrit              = realloc(mercurius->dcrit, sizeof(double)*N);
        mercurius->N_allocated_dcrit = N;
        // Heliocentric coordinates were never calculated.
        // This will get triggered on first step only (not when loaded from archive)
        r->did_modify_particles = 1;
    }
    if (mercurius->N_allocated<N){
        // These arrays are only used within one timestep. 
        // Can be recreated without loosing bit-wise reproducibility
        mercurius->particles_backup   = realloc(mercurius->particles_backup,sizeof(struct reb_particle)*N);
        mercurius->particles_backup_xv = realloc(mercurius->particles_backup_xv,sizeof(struct reb_vec6d)*N);
        mercurius->encounter_map      = realloc(mercurius->encounter_map,sizeof(size_t)*N);
        mercurius->N_allocated = N;
    }
    if (mercurius->safe_mode || r->did_modify_particles){
        if (r->is_synchronized==0){
            reb_integrator_mercurius_synchronize(r, state);
            reb_simulation_warning(r,"Particles were modified while simulation was not synchronized.");
        }
        reb_integrator_mercurius_inertial_to_dh(r, mercurius);
    }

    if (r->did_modify_particles){
        if (r->is_synchronized==0){
            reb_integrator_mercurius_synchronize(r, mercurius);
            reb_integrator_mercurius_inertial_to_dh(r, mercurius);
            reb_simulation_warning(r,"MERCURIUS: Recalculating dcrit but pos/vel were not synchronized before.");
        }
        mercurius->dcrit[0] = 2.*r->particles[0].r; // central object only uses physical radius
        for (size_t i=1;i<N;i++){
            mercurius->dcrit[i] = reb_integrator_mercurius_calculate_dcrit_for_particle(r, mercurius, i);
        }
    }

    // Calculate collisions only with DIRECT method
    if (r->collision != REB_COLLISION_NONE && r->collision != REB_COLLISION_DIRECT){
        reb_simulation_warning(r,"Mercurius only works with a direct collision search.");
    }

    // Calculate gravity with special function
    if (r->gravity != REB_GRAVITY_BASIC && r->gravity != REB_GRAVITY_CUSTOM){
        reb_simulation_warning(r,"Mercurius has its own gravity routine. Gravity routine set by the user will be ignored.");
    }
    mercurius->mode = REB_INTEGRATOR_MERCURIUS_MODE_WH;

    reb_integrator_mercurius_calculate_acceleration_mode_wh(r, mercurius);

    if (r->is_synchronized){
        reb_integrator_mercurius_interaction_step(r,r->dt/2.);
    }else{
        reb_integrator_mercurius_interaction_step(r,r->dt);
    }
    reb_integrator_mercurius_jump_step(r,r->dt/2.);

    // COM step
    mercurius->com_pos.x += r->dt*mercurius->com_vel.x;
    mercurius->com_pos.y += r->dt*mercurius->com_vel.y;
    mercurius->com_pos.z += r->dt*mercurius->com_vel.z;

    // Make copy of particles before the kepler step.
    // Then evolve all particles in kepler step.
    // Result will be used in encounter prediction.
    // Particles having a close encounter will be overwritten 
    // later by encounter step.
    // Only positions and velocities are copied: the Kepler step and the encounter
    // prediction do not change any other field.
    {
        const struct reb_particle* const p = r->particles;
        struct reb_vec6d* const b = mercurius->particles_backup_xv;
        for (size_t i=0; i<N; i++){
            b[i].x = p[i].x; b[i].y = p[i].y; b[i].z = p[i].z;
            b[i].vx = p[i].vx; b[i].vy = p[i].vy; b[i].vz = p[i].vz;
        }
    }
    reb_integrator_mercurius_kepler_step(r, r->dt);

    reb_mercurius_encounter_predict(r, mercurius);

    reb_mercurius_encounter_step(r, mercurius, r->dt);

    reb_integrator_mercurius_jump_step(r, r->dt/2.);

    r->is_synchronized = 0;
    if (mercurius->safe_mode){
        reb_integrator_mercurius_synchronize(r, state);
    }

    r->t+=r->dt;
    r->dt_last_done = r->dt;
    r->N_targets = 1; // Only search for collisions with star in-between timesteps.
}

void reb_integrator_mercurius_synchronize(struct reb_simulation* r, void* state){
    struct reb_integrator_mercurius_state* mercurius = state;
    if (r->is_synchronized == 0){
        mercurius->mode = REB_INTEGRATOR_MERCURIUS_MODE_WH;
        reb_integrator_mercurius_calculate_acceleration_mode_wh(r, mercurius);
        reb_integrator_mercurius_interaction_step(r, r->dt/2.);

        reb_integrator_mercurius_dh_to_inertial(r, mercurius);

        r->is_synchronized = 1;
    }
}

