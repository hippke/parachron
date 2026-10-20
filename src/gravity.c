/**
 * @file     gravity.c
 * @brief     Direct gravity calculation, O(N^2).
 * @author     Hanno Rein <hanno@hanno-rein.de>
 *
 * @details     This is the crudest implementation of an N-body code
 * which sums up every pair of particles. It is only useful very small 
 * particle numbers (N<~100) as it scales as O(N^2). Note that the MPI
 * implementation is not well tested and only works for very specific
 * problems. This should be resolved in the future. 
 *
 * 
 * @section LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
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
#include <math.h>
#include "particle.h"
#include "output.h"
#include "tree.h"
#include "boundary.h"
#define MAX(a, b) ((a) > (b) ? (a) : (b))    ///< Returns the maximum of a and b

#ifdef MPI
#include "communication_mpi.h"
#endif

// Main Gravity Routine


// Group walk (opt-in, r->tree_group_size = G > 0; Barnes 1990, as in PKDGRAV and Bonsai).
// The tree is cut into groups: the largest cells with at most G particles. For each group one
// walk collects an interaction list: a cell is accepted if w < theta*(d - R), with d the distance
// from the group centre to the cell's centre of mass and R the group radius (so the usual
// criterion w < theta*d holds for every particle in the group), otherwise it is opened; leaves
// are taken as particles. Every particle of the group then sums over the same list. One walk
// per group instead of per particle, and the inner loop over the group's particles vectorises.
// The results differ from the per-particle walk (deterministic, but not bitwise identical).
static void reb_tree_group_collect(const struct reb_treecell* node, size_t* idx, size_t* n){
    if (node->pt >= 0){ idx[(*n)++] = (size_t)node->pt; return; }
    for (int o=0; o<8; o++) if (node->oct[o]) reb_tree_group_collect(node->oct[o], idx, n);
}
static void reb_tree_group_find(const struct reb_treecell* node, const size_t G, const struct reb_treecell** gl, size_t* ng){
    const size_t cnt = node->pt >= 0 ? 1 : (size_t)(-node->pt);
    if (cnt <= G){ gl[(*ng)++] = node; return; }
    for (int o=0; o<8; o++) if (node->oct[o]) reb_tree_group_find(node->oct[o], G, gl, ng);
}
#define REB_TREE_GROUP_MAX 256
#define REB_TREE_GROUP_STACK 1024
static void reb_gravity_tree_group_walk(struct reb_simulation* const r){
    struct reb_particle* const particles = r->particles;
    const size_t N = r->N;
    const double G = r->G;
    const double softening2 = r->softening*r->softening;
    const double theta = sqrt(r->opening_angle2);
    const size_t N_root = r->N_root_x*r->N_root_y*r->N_root_z;
    size_t Gs = r->tree_group_size > REB_TREE_GROUP_MAX ? REB_TREE_GROUP_MAX : r->tree_group_size;
    // Groups and their particles (in depth-first order), once per force evaluation.
    const struct reb_treecell** gl = malloc(sizeof(*gl)*(N+1));
    size_t ng = 0;
    for (size_t i=0; i<N_root; i++) if (r->tree_root[i]) reb_tree_group_find(r->tree_root[i], Gs, gl, &ng);
    size_t* gidx = malloc(sizeof(size_t)*(N+1));
    size_t* goff = malloc(sizeof(size_t)*(ng+1));
    size_t n = 0;
    for (size_t g=0; g<ng; g++){ goff[g] = n; reb_tree_group_collect(gl[g], gidx, &n); }
    goff[ng] = n;
    unsigned char* ingroup = calloc(N>0?N:1, 1);
    size_t capc = 1024, capp = 1024;
    double* cl = malloc(sizeof(double)*4*capc);   // accepted cells: mx, my, mz, m
#ifdef QUADRUPOLE
    double* cq = malloc(sizeof(double)*6*capc);   // and their quadrupole moments
#endif
    size_t* pl = malloc(sizeof(size_t)*capp);     // particles (leaves) outside the group
    size_t* plg = malloc(sizeof(size_t)*(REB_TREE_GROUP_MAX+1)); // leaves inside the group
    double gx[REB_TREE_GROUP_MAX], gy[REB_TREE_GROUP_MAX], gz[REB_TREE_GROUP_MAX];
    double gax[REB_TREE_GROUP_MAX], gay[REB_TREE_GROUP_MAX], gaz[REB_TREE_GROUP_MAX];
    const struct reb_treecell* stack[REB_TREE_GROUP_STACK];
    for (int gbx=-r->N_ghost_x; gbx<=r->N_ghost_x; gbx++){
    for (int gby=-r->N_ghost_y; gby<=r->N_ghost_y; gby++){
    for (int gbz=-r->N_ghost_z; gbz<=r->N_ghost_z; gbz++){
        const struct reb_vec6d gb = reb_boundary_get_ghostbox(r, gbx,gby,gbz);
        for (size_t g=0; g<ng; g++){
            if (reb_sigint > 1) goto done;
            const size_t* const idx = &gidx[goff[g]];
            const size_t k = goff[g+1]-goff[g];
            double c0=0, c1=0, c2=0;
            for (size_t a=0; a<k; a++){
                const struct reb_particle* const p = &particles[idx[a]];
                gx[a] = p->x + gb.x; gy[a] = p->y + gb.y; gz[a] = p->z + gb.z;
                gax[a] = 0; gay[a] = 0; gaz[a] = 0;
                c0 += gx[a]; c1 += gy[a]; c2 += gz[a];
                ingroup[idx[a]] = 1;
            }
            c0 /= (double)k; c1 /= (double)k; c2 /= (double)k;
            double R = 0;
            for (size_t a=0; a<k; a++){
                const double dx = gx[a]-c0, dy = gy[a]-c1, dz = gz[a]-c2;
                const double d = sqrt(dx*dx+dy*dy+dz*dz);
                if (d > R) R = d;
            }
            size_t nc = 0, np = 0, npg = 0;
            for (size_t i=0; i<N_root; i++){
                if (!r->tree_root[i]) continue;
                int sp = 0;
                stack[sp++] = r->tree_root[i];
                while (sp > 0){
                    const struct reb_treecell* const node = stack[--sp];
                    if (node->pt >= 0){
                        const size_t j = (size_t)node->pt;
                        if (ingroup[j]){ // (also in ghost boxes: like the per-particle walk, no interaction with the own image)
                            plg[npg++] = j;
                        }else{
                            if (np == capp){ capp *= 2; pl = realloc(pl, sizeof(size_t)*capp); }
                            pl[np++] = j;
                        }
                        continue;
                    }
                    const double dx = c0-node->mx, dy = c1-node->my, dz = c2-node->mz;
                    const double d = sqrt(dx*dx+dy*dy+dz*dz);
                    if (d > R && node->w < theta*(d-R)){
                        if (nc == capc){
                            capc *= 2;
                            cl = realloc(cl, sizeof(double)*4*capc);
#ifdef QUADRUPOLE
                            cq = realloc(cq, sizeof(double)*6*capc);
#endif
                        }
                        cl[4*nc] = node->mx; cl[4*nc+1] = node->my; cl[4*nc+2] = node->mz; cl[4*nc+3] = node->m;
#ifdef QUADRUPOLE
                        cq[6*nc] = node->mxx; cq[6*nc+1] = node->mxy; cq[6*nc+2] = node->mxz;
                        cq[6*nc+3] = node->myy; cq[6*nc+4] = node->myz; cq[6*nc+5] = node->mzz;
#endif
                        nc++;
                    }else{
                        if (sp + 8 > REB_TREE_GROUP_STACK){
                            reb_simulation_error(r, "Tree too deep for the group walk.");
                            goto done;
                        }
                        for (int o=7; o>=0; o--) if (node->oct[o]) stack[sp++] = node->oct[o];
                    }
                }
            }
            // Cells (vectorised over the group's particles).
            for (size_t b=0; b<nc; b++){
                const double mx = cl[4*b], my = cl[4*b+1], mz = cl[4*b+2], m = cl[4*b+3];
#ifdef QUADRUPOLE
                const double qxx = cq[6*b], qxy = cq[6*b+1], qxz = cq[6*b+2], qyy = cq[6*b+3], qyz = cq[6*b+4], qzz = cq[6*b+5];
#endif
                for (size_t a=0; a<k; a++){
                    const double dx = gx[a]-mx, dy = gy[a]-my, dz = gz[a]-mz;
                    const double r2 = dx*dx + dy*dy + dz*dz;
                    const double _r = sqrt(r2 + softening2);
                    const double prefact = -G/(_r*_r*_r)*m;
#ifdef QUADRUPOLE
                    double qprefact = G/(_r*_r*_r*_r*_r);
                    gax[a] += qprefact*(dx*qxx + dy*qxy + dz*qxz);
                    gay[a] += qprefact*(dx*qxy + dy*qyy + dz*qyz);
                    gaz[a] += qprefact*(dx*qxz + dy*qyz + dz*qzz);
                    const double mrr = dx*dx*qxx + dy*dy*qyy + dz*dz*qzz + 2.*dx*dy*qxy + 2.*dx*dz*qxz + 2.*dy*dz*qyz;
                    qprefact *= -5.0/(2.0*_r*_r)*mrr;
                    gax[a] += (qprefact + prefact)*dx;
                    gay[a] += (qprefact + prefact)*dy;
                    gaz[a] += (qprefact + prefact)*dz;
#else
                    gax[a] += prefact*dx;
                    gay[a] += prefact*dy;
                    gaz[a] += prefact*dz;
#endif
                }
            }
            // Particles outside the group (vectorised over the group's particles).
            for (size_t b=0; b<np; b++){
                const struct reb_particle* const pj = &particles[pl[b]];
                const double px = pj->x, py = pj->y, pz = pj->z, m = pj->m;
                for (size_t a=0; a<k; a++){
                    const double dx = gx[a]-px, dy = gy[a]-py, dz = gz[a]-pz;
                    const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
                    const double prefact = -G/(_r*_r*_r)*m;
                    gax[a] += prefact*dx;
                    gay[a] += prefact*dy;
                    gaz[a] += prefact*dz;
                }
            }
            // Particles of the group itself (no self-interaction).
            for (size_t b=0; b<npg; b++){
                const struct reb_particle* const pj = &particles[plg[b]];
                for (size_t a=0; a<k; a++){
                    if (idx[a] == plg[b]) continue;
                    const double dx = gx[a]-pj->x, dy = gy[a]-pj->y, dz = gz[a]-pj->z;
                    const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
                    const double prefact = -G/(_r*_r*_r)*pj->m;
                    gax[a] += prefact*dx;
                    gay[a] += prefact*dy;
                    gaz[a] += prefact*dz;
                }
            }
            for (size_t a=0; a<k; a++){
                struct reb_particle* const p = &particles[idx[a]];
                p->ax += gax[a];
                p->ay += gay[a];
                p->az += gaz[a];
                ingroup[idx[a]] = 0;
            }
        }
    }}}
done:
    free(gl); free(gidx); free(goff); free(ingroup); free(cl); free(pl); free(plg);
#ifdef QUADRUPOLE
    free(cq);
#endif
}

void reb_gravity_tree_calculate_acceleration(struct reb_simulation* r){
    PROFILING_START();
    // Construct tree first.
    // Note the boundary_check and distribute_particles can change the number of particles on the current node. 
    // This is not compatible with some integrators. Also the reason this is up here, rather than in 
    // the switch statement below.

    // Check if particles are in box 
    PROFILING_START();
    reb_boundary_check(r);     
    PROFILING_STOP(PROFILING_CAT_BOUNDARY);
#ifdef MPI
    // Check if particles are in local rootbox, if not distribute 
    reb_communication_mpi_distribute_particles(r);
#endif // MPI

    reb_tree_construct(r);
    if (!r->tree_root){
        reb_simulation_error(r,"Tree does not exist. Cannot calculate accelerations.");
        return;
    }

    // Update center of mass and quadrupole moments in tree in preparation of force calculation.
    // Also distributed essential tree if MPI is used.
    reb_tree_calculate_gravity_data(r); 


    struct reb_particle* const particles = r->particles;
    const size_t N = r->N;

#pragma omp parallel for schedule(guided)
    for (size_t i=0; i<N; i++){
        particles[i].ax = 0; 
        particles[i].ay = 0; 
        particles[i].az = 0; 
    }
#ifndef MPI
    if (r->tree_group_size > 0){
        reb_gravity_tree_group_walk(r);
        reb_tree_delete(r);
        PROFILING_STOP(PROFILING_CAT_GRAVITY);
        return;
    }
#endif // MPI
    // Visit particles in the depth-first leaf order of the tree: consecutive particles are
    // spatially close and walk mostly the same nodes, which keeps those nodes in cache.
    // Each particle's own walk is unchanged, so the results are identical.
    size_t* order = malloc(sizeof(size_t)*(N>0?N:1));
    if (reb_tree_particle_order(r, order) != N){
        for (size_t i=0; i<N; i++) order[i] = i; // e.g. particles missing from the tree (MPI): original order
    }
    // Summing over all Ghost Boxes
    for (int gbx=-r->N_ghost_x; gbx<=r->N_ghost_x; gbx++){
        for (int gby=-r->N_ghost_y; gby<=r->N_ghost_y; gby++){
            for (int gbz=-r->N_ghost_z; gbz<=r->N_ghost_z; gbz++){
                // The ghost box offset is the same for all particles.
                const struct reb_vec6d gb0 = reb_boundary_get_ghostbox(r, gbx,gby,gbz);
                // Summing over all particle pairs
#pragma omp parallel for schedule(guided)
                for (size_t k=0; k<N; k++){
                    const size_t i = order[k];
#ifndef OPENMP
                    if (reb_sigint > 1){ free(order); return; }
#endif // OPENMP
                    struct reb_vec6d gb = gb0;
                    // Precalculated shifted position
                    gb.x += particles[i].x;
                    gb.y += particles[i].y;
                    gb.z += particles[i].z;
                    reb_tree_calculate_acceleration_for_particle(r, i, gb);
                }
            }
        }
    }
    free(order);
    // Delete tree (if it exists)    
    reb_tree_delete(r);
    PROFILING_STOP(PROFILING_CAT_GRAVITY);
}

void reb_gravity_jacobi_calculate_acceleration(struct reb_simulation* r){
    PROFILING_START();
    struct reb_particle* const particles = r->particles;
    const size_t N = r->N;
    const double G = r->G;
    double Rjx = 0.;
    double Rjy = 0.;
    double Rjz = 0.;
    double Mj = 0.;
    for (size_t j=0; j<N; j++){
        particles[j].ax = 0; 
        particles[j].ay = 0; 
        particles[j].az = 0; 
        for (size_t i=0; i<j+1; i++){
            if (j>1){
                ////////////////
                // Jacobi Term
                // Note: ignoring j==1 term here and below as they cancel
                const double Qjx = particles[j].x - Rjx/Mj; 
                const double Qjy = particles[j].y - Rjy/Mj;
                const double Qjz = particles[j].z - Rjz/Mj;
                const double dr = sqrt(Qjx*Qjx + Qjy*Qjy + Qjz*Qjz);
                double dQjdri = Mj; 
                if (i<j){
                    dQjdri = -particles[j].m; //rearranged such that m==0 does not diverge
                }
                const double prefact = G*dQjdri/(dr*dr*dr);
                particles[i].ax    += prefact*Qjx;
                particles[i].ay    += prefact*Qjy;
                particles[i].az    += prefact*Qjz;
            }
            if (i!=j && (i!=0 || j!=1)){
                ////////////////
                // Direct Term
                // Note: ignoring i==0 && j==1 term here and above as they cancel 
                const double dx = particles[i].x - particles[j].x;
                const double dy = particles[i].y - particles[j].y;
                const double dz = particles[i].z - particles[j].z;
                const double dr = sqrt(dx*dx + dy*dy + dz*dz);
                const double prefact = G /(dr*dr*dr);
                const double prefacti = prefact*particles[i].m;
                const double prefactj = prefact*particles[j].m;

                particles[i].ax    -= prefactj*dx;
                particles[i].ay    -= prefactj*dy;
                particles[i].az    -= prefactj*dz;
                particles[j].ax    += prefacti*dx;
                particles[j].ay    += prefacti*dy;
                particles[j].az    += prefacti*dz;
            }
        }
        Rjx += particles[j].m*particles[j].x;
        Rjy += particles[j].m*particles[j].y;
        Rjz += particles[j].m*particles[j].z;
        Mj += particles[j].m;
    }
    PROFILING_STOP(PROFILING_CAT_GRAVITY);
}

#define REB_GRAVITY_BLOCK 32      // partners per block in the two-pass inner loop
#define REB_GRAVITY_BLOCK_MIN 16  // rows with fewer partners use the plain loop
// AVX2 clone selected at load time (REB_TARGET_CLONES, rebound_internal.h); identical results.
REB_TARGET_CLONES void reb_gravity_basic_calculate_acceleration(struct reb_simulation* r){
    PROFILING_START();
    struct reb_particle* const particles = r->particles;
    const size_t N = r->N;
    const double G = r->G;
    const double softening2 = r->softening*r->softening;
    const unsigned int gravity_ignore_terms = r->gravity_ignore_terms;
    const size_t N_active = ((r->N_active==SIZE_MAX)?N:r->N_active);
    const int _testparticle_type   = r->testparticle_type;
    const int N_ghost_x = r->N_ghost_x;
    const int N_ghost_y = r->N_ghost_y;
    const int N_ghost_z = r->N_ghost_z;
#ifndef OPENMP // OPENMP off
    const size_t starti = (gravity_ignore_terms==REB_GRAVITY_IGNORE_TERMS_NONE)?1:2;
    const size_t startj = (gravity_ignore_terms==REB_GRAVITY_IGNORE_TERMS_INVOLVING_0)?1:0;
#endif // OPENMP
#pragma omp parallel for 
    for (size_t i=0; i<N; i++){
        particles[i].ax = 0; 
        particles[i].ay = 0; 
        particles[i].az = 0; 
    }
    // Summing over all Ghost Boxes
    for (int gbx=-N_ghost_x; gbx<=N_ghost_x; gbx++){
        for (int gby=-N_ghost_y; gby<=N_ghost_y; gby++){
            for (int gbz=-N_ghost_z; gbz<=N_ghost_z; gbz++){
                struct reb_vec6d gb = reb_boundary_get_ghostbox(r, gbx,gby,gbz);
                // All active particle pairs
#ifndef OPENMP // OPENMP off, do O(1/2*N^2)
                for (size_t i=starti; i<N_active; i++){
                    if (reb_sigint > 1) return;
                    // Particle i's position and acceleration live in registers during the inner loop
                    // (the compiler cannot prove that the writes to particles[j] do not alias them).
                    // The additions happen in the same order as before, so results are unchanged.
                    const double xi = gb.x+particles[i].x;
                    const double yi = gb.y+particles[i].y;
                    const double zi = gb.z+particles[i].z;
                    const double mi = particles[i].m;
                    double axi = particles[i].ax;
                    double ayi = particles[i].ay;
                    double azi = particles[i].az;
                    if (i-startj < REB_GRAVITY_BLOCK_MIN){
                        for (size_t j=startj; j<i; j++){
                            const double dx = xi - particles[j].x;
                            const double dy = yi - particles[j].y;
                            const double dz = zi - particles[j].z;
                            const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
                            const double prefact = G/(_r*_r*_r);
                            const double prefactj = -prefact*particles[j].m;
                            const double prefacti = prefact*mi;

                            axi    += prefactj*dx;
                            ayi    += prefactj*dy;
                            azi    += prefactj*dz;
                            particles[j].ax    += prefacti*dx;
                            particles[j].ay    += prefacti*dy;
                            particles[j].az    += prefacti*dz;
                        }
                    }else{
                        // Blocks of up to REB_GRAVITY_BLOCK partners: first the separations and 1/r^3
                        // factors (independent per pair, vectorisable: square roots and divisions run
                        // two or more at a time), then the accumulation in the original order.
                        for (size_t j0=startj; j0<i; j0+=REB_GRAVITY_BLOCK){
                            const size_t nb = (i-j0) < REB_GRAVITY_BLOCK ? (i-j0) : REB_GRAVITY_BLOCK;
                            double bdx[REB_GRAVITY_BLOCK], bdy[REB_GRAVITY_BLOCK], bdz[REB_GRAVITY_BLOCK], bpf[REB_GRAVITY_BLOCK];
                            for (size_t k=0; k<nb; k++){
                                const struct reb_particle* const pj = &particles[j0+k];
                                bdx[k] = xi - pj->x;
                                bdy[k] = yi - pj->y;
                                bdz[k] = zi - pj->z;
                            }
                            for (size_t k=0; k<nb; k++){
                                const double _r = sqrt(bdx[k]*bdx[k] + bdy[k]*bdy[k] + bdz[k]*bdz[k] + softening2);
                                bpf[k] = G/(_r*_r*_r);
                            }
                            for (size_t k=0; k<nb; k++){
                                struct reb_particle* const pj = &particles[j0+k];
                                const double dx = bdx[k], dy = bdy[k], dz = bdz[k];
                                const double prefact = bpf[k];
                                const double prefactj = -prefact*pj->m;
                                const double prefacti = prefact*mi;

                                axi    += prefactj*dx;
                                ayi    += prefactj*dy;
                                azi    += prefactj*dz;
                                pj->ax    += prefacti*dx;
                                pj->ay    += prefacti*dy;
                                pj->az    += prefacti*dz;
                            }
                        }
                    }
                    particles[i].ax = axi;
                    particles[i].ay = ayi;
                    particles[i].az = azi;
                }
#else // OPENMP on, do O(N^2)
#pragma omp parallel for
                for (size_t i=0; i<N; i++){
                    for (size_t j=0; j<N_active; j++){
                        if (gravity_ignore_terms==REB_GRAVITY_IGNORE_TERMS_BETWEEN_0_AND_1 && ((j==1 && i==0) || (i==1 && j==0) )) continue;
                        if (gravity_ignore_terms==REB_GRAVITY_IGNORE_TERMS_INVOLVING_0 && ((j==0 || i==0) )) continue;
                        if (i==j) continue;
                        const double dx = (gb.x+particles[i].x) - particles[j].x;
                        const double dy = (gb.y+particles[i].y) - particles[j].y;
                        const double dz = (gb.z+particles[i].z) - particles[j].z;
                        const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
                        const double prefact = -G/(_r*_r*_r)*particles[j].m;

                        particles[i].ax    += prefact*dx;
                        particles[i].ay    += prefact*dy;
                        particles[i].az    += prefact*dz;
                    }
                }
#endif // OPENMP
       // Interactions of test particles with active particles
#ifndef OPENMP // OPENMP off
                const size_t startitestp = MAX(N_active, starti);
                for (size_t i=startitestp; i<N; i++){
                    if (reb_sigint > 1) return;
                    const double xi = gb.x+particles[i].x;
                    const double yi = gb.y+particles[i].y;
                    const double zi = gb.z+particles[i].z;
                    const double mi = particles[i].m;
                    double axi = particles[i].ax;
                    double ayi = particles[i].ay;
                    double azi = particles[i].az;
                    for (size_t j=startj; j<N_active; j++){
                        const double dx = xi - particles[j].x;
                        const double dy = yi - particles[j].y;
                        const double dz = zi - particles[j].z;
                        const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
                        const double prefact = G/(_r*_r*_r);
                        const double prefactj = -prefact*particles[j].m;

                        axi    += prefactj*dx;
                        ayi    += prefactj*dy;
                        azi    += prefactj*dz;
                        if (_testparticle_type){
                            const double prefacti = prefact*mi;
                            particles[j].ax    += prefacti*dx;
                            particles[j].ay    += prefacti*dy;
                            particles[j].az    += prefacti*dz;
                        }
                    }
                    particles[i].ax = axi;
                    particles[i].ay = ayi;
                    particles[i].az = azi;
                }
#else // OPENMP on
                if (_testparticle_type){
#pragma omp parallel for
                    for (size_t i=0; i<N_active; i++){
                        for (size_t j=N_active; j<N; j++){
                            if (gravity_ignore_terms==REB_GRAVITY_IGNORE_TERMS_BETWEEN_0_AND_1 && ((j==1 && i==0) )) continue;
                            if (gravity_ignore_terms==REB_GRAVITY_IGNORE_TERMS_INVOLVING_0 && ((j==0 || i==0) )) continue;
                            const double dx = (gb.x+particles[i].x) - particles[j].x;
                            const double dy = (gb.y+particles[i].y) - particles[j].y;
                            const double dz = (gb.z+particles[i].z) - particles[j].z;
                            const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
                            const double prefact = -G/(_r*_r*_r)*particles[j].m;

                            particles[i].ax    += prefact*dx;
                            particles[i].ay    += prefact*dy;
                            particles[i].az    += prefact*dz;
                        }
                    }
                }
#endif // OPENMP
            }
        }
    }
    PROFILING_STOP(PROFILING_CAT_GRAVITY);
}

void reb_gravity_compensated_calculate_acceleration(struct reb_simulation* r){
    PROFILING_START();
    struct reb_particle* const particles = r->particles;
    const size_t N = r->N;
    const double G = r->G;
    const double softening2 = r->softening*r->softening;
    const unsigned int gravity_ignore_terms = r->gravity_ignore_terms;
    const size_t N_active = ((r->N_active==SIZE_MAX)?N:r->N_active);
    const int _testparticle_type   = r->testparticle_type;
    if (r->N_allocated_gravity_cs<N){
        r->gravity_cs = realloc(r->gravity_cs,N*sizeof(struct reb_vec3d));
        r->N_allocated_gravity_cs = N;
    }
    struct reb_vec3d* restrict const cs = r->gravity_cs;
#pragma omp parallel for schedule(guided)
    for (size_t i=0; i<N; i++){
        particles[i].ax = 0.; 
        particles[i].ay = 0.; 
        particles[i].az = 0.; 
        cs[i].x = 0.;
        cs[i].y = 0.;
        cs[i].z = 0.;
    }
    // Summing over all massive particle pairs
#ifdef OPENMP
#pragma omp parallel for schedule(guided)
    for (size_t i=0; i<N_active; i++){
        for (size_t j=0; j<N_active; j++){
            if (gravity_ignore_terms==REB_GRAVITY_IGNORE_TERMS_BETWEEN_0_AND_1 && ((j==1 && i==0) || (i==1 && j==0))) continue;
            if (gravity_ignore_terms==REB_GRAVITY_IGNORE_TERMS_INVOLVING_0 && ((j==0 || i==0))) continue;
            if (i==j) continue;
            const double dx = particles[i].x - particles[j].x;
            const double dy = particles[i].y - particles[j].y;
            const double dz = particles[i].z - particles[j].z;
            const double r2 = dx*dx + dy*dy + dz*dz + softening2;
            const double r = sqrt(r2);
            const double prefact  = G/(r2*r);
            const double prefactj = -prefact*particles[j].m;

            {
                double ix = prefactj*dx;
                double yx = ix - cs[i].x;
                double tx = particles[i].ax + yx;
                cs[i].x = (tx - particles[i].ax) - yx;
                particles[i].ax = tx;

                double iy = prefactj*dy;
                double yy = iy- cs[i].y;
                double ty = particles[i].ay + yy;
                cs[i].y = (ty - particles[i].ay) - yy;
                particles[i].ay = ty;

                double iz = prefactj*dz;
                double yz = iz - cs[i].z;
                double tz = particles[i].az + yz;
                cs[i].z = (tz - particles[i].az) - yz;
                particles[i].az = tz;
            }
        }
    }

    // Testparticles
#pragma omp parallel for schedule(guided)
    for (size_t i=N_active; i<N; i++){
        for (size_t j=0; j<N_active; j++){
            if (gravity_ignore_terms==REB_GRAVITY_IGNORE_TERMS_BETWEEN_0_AND_1 && ((j==1 && i==0) || (i==1 && j==0))) continue;
            if (gravity_ignore_terms==REB_GRAVITY_IGNORE_TERMS_INVOLVING_0 && ((j==0 || i==0))) continue;
            const double dx = particles[i].x - particles[j].x;
            const double dy = particles[i].y - particles[j].y;
            const double dz = particles[i].z - particles[j].z;
            const double r2 = dx*dx + dy*dy + dz*dz + softening2;
            const double r = sqrt(r2);
            const double prefact  = G/(r2*r);
            const double prefactj = -prefact*particles[j].m;

            {
                double ix = prefactj*dx;
                double yx = ix - cs[i].x;
                double tx = particles[i].ax + yx;
                cs[i].x = (tx - particles[i].ax) - yx;
                particles[i].ax = tx;

                double iy = prefactj*dy;
                double yy = iy- cs[i].y;
                double ty = particles[i].ay + yy;
                cs[i].y = (ty - particles[i].ay) - yy;
                particles[i].ay = ty;

                double iz = prefactj*dz;
                double yz = iz - cs[i].z;
                double tz = particles[i].az + yz;
                cs[i].z = (tz - particles[i].az) - yz;
                particles[i].az = tz;
            }
        }
    }
    if (_testparticle_type){
#pragma omp parallel for schedule(guided)
        for (size_t j=0; j<N_active; j++){
            for (size_t i=N_active; i<N; i++){
                if (gravity_ignore_terms==REB_GRAVITY_IGNORE_TERMS_BETWEEN_0_AND_1 && ((j==1 && i==0) || (i==1 && j==0))) continue;
                if (gravity_ignore_terms==REB_GRAVITY_IGNORE_TERMS_INVOLVING_0 && ((j==0 || i==0))) continue;
                const double dx = particles[i].x - particles[j].x;
                const double dy = particles[i].y - particles[j].y;
                const double dz = particles[i].z - particles[j].z;
                const double r2 = dx*dx + dy*dy + dz*dz + softening2;
                const double r = sqrt(r2);
                const double prefact  = G/(r2*r);
                const double prefacti = prefact*particles[i].m;
                {
                    double ix = prefacti*dx;
                    double yx = ix - cs[j].x;
                    double tx = particles[j].ax + yx;
                    cs[j].x = (tx - particles[j].ax) - yx;
                    particles[j].ax = tx;

                    double iy = prefacti*dy;
                    double yy = iy - cs[j].y;
                    double ty = particles[j].ay + yy;
                    cs[j].y = (ty - particles[j].ay) - yy;
                    particles[j].ay = ty;

                    double iz = prefacti*dz;
                    double yz = iz - cs[j].z;
                    double tz = particles[j].az + yz;
                    cs[j].z = (tz - particles[j].az) - yz;
                    particles[j].az = tz;
                }
            }
        }
    }
#else // OPENMP
    for (size_t i=0; i<N_active; i++){
        if (reb_sigint > 1) return;
        for (size_t j=i+1; j<N_active; j++){
            if (gravity_ignore_terms==REB_GRAVITY_IGNORE_TERMS_BETWEEN_0_AND_1 && ((j==1 && i==0) || (i==1 && j==0))) continue;
            if (gravity_ignore_terms==REB_GRAVITY_IGNORE_TERMS_INVOLVING_0 && ((j==0 || i==0))) continue;
            const double dx = particles[i].x - particles[j].x;
            const double dy = particles[i].y - particles[j].y;
            const double dz = particles[i].z - particles[j].z;
            const double r2 = dx*dx + dy*dy + dz*dz + softening2;
            const double r = sqrt(r2);
            const double prefact  = G/(r2*r);
            const double prefacti = prefact*particles[i].m;
            const double prefactj = -prefact*particles[j].m;

            {
                double ix = prefactj*dx;
                double yx = ix - cs[i].x;
                double tx = particles[i].ax + yx;
                cs[i].x = (tx - particles[i].ax) - yx;
                particles[i].ax = tx;

                double iy = prefactj*dy;
                double yy = iy- cs[i].y;
                double ty = particles[i].ay + yy;
                cs[i].y = (ty - particles[i].ay) - yy;
                particles[i].ay = ty;

                double iz = prefactj*dz;
                double yz = iz - cs[i].z;
                double tz = particles[i].az + yz;
                cs[i].z = (tz - particles[i].az) - yz;
                particles[i].az = tz;
            }

            {
                double ix = prefacti*dx;
                double yx = ix - cs[j].x;
                double tx = particles[j].ax + yx;
                cs[j].x = (tx - particles[j].ax) - yx;
                particles[j].ax = tx;

                double iy = prefacti*dy;
                double yy = iy - cs[j].y;
                double ty = particles[j].ay + yy;
                cs[j].y = (ty - particles[j].ay) - yy;
                particles[j].ay = ty;

                double iz = prefacti*dz;
                double yz = iz - cs[j].z;
                double tz = particles[j].az + yz;
                cs[j].z = (tz - particles[j].az) - yz;
                particles[j].az = tz;
            }
        }
    }

    // Testparticles
    for (size_t i=N_active; i<N; i++){
        if (reb_sigint > 1) return;
        for (size_t j=0; j<N_active; j++){
            if (gravity_ignore_terms==REB_GRAVITY_IGNORE_TERMS_BETWEEN_0_AND_1 && ((j==1 && i==0) || (i==1 && j==0))) continue;
            if (gravity_ignore_terms==REB_GRAVITY_IGNORE_TERMS_INVOLVING_0 && ((j==0 || i==0))) continue;
            const double dx = particles[i].x - particles[j].x;
            const double dy = particles[i].y - particles[j].y;
            const double dz = particles[i].z - particles[j].z;
            const double r2 = dx*dx + dy*dy + dz*dz + softening2;
            const double r = sqrt(r2);
            const double prefact  = G/(r2*r);
            const double prefactj = -prefact*particles[j].m;

            {
                double ix = prefactj*dx;
                double yx = ix - cs[i].x;
                double tx = particles[i].ax + yx;
                cs[i].x = (tx - particles[i].ax) - yx;
                particles[i].ax = tx;

                double iy = prefactj*dy;
                double yy = iy- cs[i].y;
                double ty = particles[i].ay + yy;
                cs[i].y = (ty - particles[i].ay) - yy;
                particles[i].ay = ty;

                double iz = prefactj*dz;
                double yz = iz - cs[i].z;
                double tz = particles[i].az + yz;
                cs[i].z = (tz - particles[i].az) - yz;
                particles[i].az = tz;
            }
            if (_testparticle_type){
                const double prefacti = prefact*particles[i].m;
                {
                    double ix = prefacti*dx;
                    double yx = ix - cs[j].x;
                    double tx = particles[j].ax + yx;
                    cs[j].x = (tx - particles[j].ax) - yx;
                    particles[j].ax = tx;

                    double iy = prefacti*dy;
                    double yy = iy - cs[j].y;
                    double ty = particles[j].ay + yy;
                    cs[j].y = (ty - particles[j].ay) - yy;
                    particles[j].ay = ty;

                    double iz = prefacti*dz;
                    double yz = iz - cs[j].z;
                    double tz = particles[j].az + yz;
                    cs[j].z = (tz - particles[j].az) - yz;
                    particles[j].az = tz;
                }
            }
        }
    }
#endif // OPENMP
    PROFILING_STOP(PROFILING_CAT_GRAVITY);
}

void reb_gravity_basic_calculate_acceleration_var(struct reb_simulation* r){
    PROFILING_START();
    struct reb_particle* const particles = r->particles;
    struct reb_particle* const particles_var = r->particles_var;
    const double G = r->G;
    const unsigned int gravity_ignore_terms = r->gravity_ignore_terms;
    const int _testparticle_type   = r->testparticle_type;
    const size_t N = r->N;
    const size_t N_active = ((r->N_active==SIZE_MAX)?N:r->N_active);
    const size_t starti = gravity_ignore_terms==REB_GRAVITY_IGNORE_TERMS_NONE?1:2;
    const size_t startj = gravity_ignore_terms==REB_GRAVITY_IGNORE_TERMS_INVOLVING_0?1:0;
    for (size_t v=0;v<r->N_var_config;v++){
        struct reb_variational_configuration const vc = r->var_config[v];
        if (vc.order==1){
            //////////////////
            /// 1st order  ///
            //////////////////
            struct reb_particle* const particles_var1 = particles_var + vc.index;
            if (vc.testparticle<0){
                for (size_t i=0; i<N; i++){
                    particles_var1[i].ax = 0.; 
                    particles_var1[i].ay = 0.; 
                    particles_var1[i].az = 0.; 
                }
                for (size_t i=starti; i<N_active; i++){
                    for (size_t j=startj; j<i; j++){
                        const double dx = particles[i].x - particles[j].x;
                        const double dy = particles[i].y - particles[j].y;
                        const double dz = particles[i].z - particles[j].z;
                        const double r2 = dx*dx + dy*dy + dz*dz;
                        const double _r  = sqrt(r2);
                        const double r3inv = 1./(r2*_r);
                        const double r5inv = 3.*r3inv/r2;
                        const double ddx = particles_var1[i].x - particles_var1[j].x;
                        const double ddy = particles_var1[i].y - particles_var1[j].y;
                        const double ddz = particles_var1[i].z - particles_var1[j].z;
                        const double Gmi = G * particles[i].m;
                        const double Gmj = G * particles[j].m;

                        // Variational equations
                        const double dxdx = dx*dx*r5inv - r3inv;
                        const double dydy = dy*dy*r5inv - r3inv;
                        const double dzdz = dz*dz*r5inv - r3inv;
                        const double dxdy = dx*dy*r5inv;
                        const double dxdz = dx*dz*r5inv;
                        const double dydz = dy*dz*r5inv;
                        const double dax =   ddx * dxdx + ddy * dxdy + ddz * dxdz;
                        const double day =   ddx * dxdy + ddy * dydy + ddz * dydz;
                        const double daz =   ddx * dxdz + ddy * dydz + ddz * dzdz;

                        // Variational mass contributions
                        const double dGmi = G*particles_var1[i].m;
                        const double dGmj = G*particles_var1[j].m;

                        particles_var1[i].ax += Gmj * dax - dGmj*r3inv*dx;
                        particles_var1[i].ay += Gmj * day - dGmj*r3inv*dy;
                        particles_var1[i].az += Gmj * daz - dGmj*r3inv*dz;

                        particles_var1[j].ax -= Gmi * dax - dGmi*r3inv*dx;
                        particles_var1[j].ay -= Gmi * day - dGmi*r3inv*dy;
                        particles_var1[j].az -= Gmi * daz - dGmi*r3inv*dz; 
                    }
                }
                for (size_t i=N_active; i<N; i++){
                    for (size_t j=startj; j<N_active; j++){
                        const double dx = particles[i].x - particles[j].x;
                        const double dy = particles[i].y - particles[j].y;
                        const double dz = particles[i].z - particles[j].z;
                        const double r2 = dx*dx + dy*dy + dz*dz;
                        const double _r  = sqrt(r2);
                        const double r3inv = 1./(r2*_r);
                        const double r5inv = 3.*r3inv/r2;
                        const double ddx = particles_var1[i].x - particles_var1[j].x;
                        const double ddy = particles_var1[i].y - particles_var1[j].y;
                        const double ddz = particles_var1[i].z - particles_var1[j].z;
                        const double Gmi = G * particles[i].m;
                        const double Gmj = G * particles[j].m;

                        // Variational equations
                        const double dxdx = dx*dx*r5inv - r3inv;
                        const double dydy = dy*dy*r5inv - r3inv;
                        const double dzdz = dz*dz*r5inv - r3inv;
                        const double dxdy = dx*dy*r5inv;
                        const double dxdz = dx*dz*r5inv;
                        const double dydz = dy*dz*r5inv;
                        const double dax =   ddx * dxdx + ddy * dxdy + ddz * dxdz;
                        const double day =   ddx * dxdy + ddy * dydy + ddz * dydz;
                        const double daz =   ddx * dxdz + ddy * dydz + ddz * dzdz;

                        // Variational mass contributions
                        const double dGmi = G*particles_var1[i].m;
                        const double dGmj = G*particles_var1[j].m;

                        particles_var1[i].ax += Gmj * dax - dGmj*r3inv*dx;
                        particles_var1[i].ay += Gmj * day - dGmj*r3inv*dy;
                        particles_var1[i].az += Gmj * daz - dGmj*r3inv*dz;
                        if (_testparticle_type){
                            // Warning! This does not make sense when the mass is varied!
                            particles_var1[j].ax -= Gmi * dax - dGmi*r3inv*dx;
                            particles_var1[j].ay -= Gmi * day - dGmi*r3inv*dy;
                            particles_var1[j].az -= Gmi * daz - dGmi*r3inv*dz; 
                        }
                    }
                }
            }else{ //testparticle
                size_t i = vc.testparticle;
                particles_var1[0].ax = 0.; 
                particles_var1[0].ay = 0.; 
                particles_var1[0].az = 0.; 
                for (size_t j=0; j<N; j++){
                    if (i==j) continue;
                    if (gravity_ignore_terms==REB_GRAVITY_IGNORE_TERMS_BETWEEN_0_AND_1 && ((j==1 && i==0) || (i==1 && j==0))) continue;
                    if (gravity_ignore_terms==REB_GRAVITY_IGNORE_TERMS_INVOLVING_0 && ((j==0 || i==0))) continue;
                    const double dx = particles[i].x - particles[j].x;
                    const double dy = particles[i].y - particles[j].y;
                    const double dz = particles[i].z - particles[j].z;
                    const double r2 = dx*dx + dy*dy + dz*dz;
                    const double _r  = sqrt(r2);
                    const double r3inv = 1./(r2*_r);
                    const double r5inv = 3.*r3inv/r2;
                    const double ddx = particles_var1[0].x;
                    const double ddy = particles_var1[0].y;
                    const double ddz = particles_var1[0].z;
                    const double Gmj = G * particles[j].m;

                    // Variational equations
                    const double dxdx = dx*dx*r5inv - r3inv;
                    const double dydy = dy*dy*r5inv - r3inv;
                    const double dzdz = dz*dz*r5inv - r3inv;
                    const double dxdy = dx*dy*r5inv;
                    const double dxdz = dx*dz*r5inv;
                    const double dydz = dy*dz*r5inv;
                    const double dax =   ddx * dxdx + ddy * dxdy + ddz * dxdz;
                    const double day =   ddx * dxdy + ddy * dydy + ddz * dydz;
                    const double daz =   ddx * dxdz + ddy * dydz + ddz * dzdz;

                    // No variational mass contributions for test particles!

                    particles_var1[0].ax += Gmj * dax;
                    particles_var1[0].ay += Gmj * day;
                    particles_var1[0].az += Gmj * daz;

                }
            }
        }else if (vc.order==2){
            if (_testparticle_type){
                reb_simulation_error(r,"testparticletype=1 not implemented for second order variational equations.");
            }
            if (gravity_ignore_terms!=REB_GRAVITY_IGNORE_TERMS_NONE){
                reb_simulation_error(r,"Second order variational equations do not support gravity_ignore_terms");
            }
            //////////////////
            /// 2nd order  ///
            //////////////////
            struct reb_particle* const particles_var2 = particles_var + vc.index;
            struct reb_particle* const particles_var1a = particles_var + vc.index_1st_order_a;
            struct reb_particle* const particles_var1b = particles_var + vc.index_1st_order_b;
            if (vc.testparticle<0){
                for (size_t i=0; i<N; i++){
                    particles_var2[i].ax = 0.; 
                    particles_var2[i].ay = 0.; 
                    particles_var2[i].az = 0.; 
                }
                for (size_t i=0; i<N; i++){
                    for (size_t j=i+1; j<N; j++){
                        const double dx = particles[i].x - particles[j].x;
                        const double dy = particles[i].y - particles[j].y;
                        const double dz = particles[i].z - particles[j].z;
                        const double r2 = dx*dx + dy*dy + dz*dz;
                        const double r  = sqrt(r2);
                        const double r3inv = 1./(r2*r);
                        const double r5inv = r3inv/r2;
                        const double r7inv = r5inv/r2;
                        const double ddx = particles_var2[i].x - particles_var2[j].x;
                        const double ddy = particles_var2[i].y - particles_var2[j].y;
                        const double ddz = particles_var2[i].z - particles_var2[j].z;
                        const double Gmi = G * particles[i].m;
                        const double Gmj = G * particles[j].m;
                        const double ddGmi = G*particles_var2[i].m;
                        const double ddGmj = G*particles_var2[j].m;

                        // Variational equations
                        // delta^(2) terms
                        double dax =         ddx * ( 3.*dx*dx*r5inv - r3inv )
                        + ddy * ( 3.*dx*dy*r5inv )
                        + ddz * ( 3.*dx*dz*r5inv );
                        double day =         ddx * ( 3.*dy*dx*r5inv )
                        + ddy * ( 3.*dy*dy*r5inv - r3inv )
                        + ddz * ( 3.*dy*dz*r5inv );
                        double daz =         ddx * ( 3.*dz*dx*r5inv )
                        + ddy * ( 3.*dz*dy*r5inv )
                        + ddz * ( 3.*dz*dz*r5inv - r3inv );

                        // delta^(1) delta^(1) terms
                        const double dk1dx = particles_var1a[i].x - particles_var1a[j].x;
                        const double dk1dy = particles_var1a[i].y - particles_var1a[j].y;
                        const double dk1dz = particles_var1a[i].z - particles_var1a[j].z;
                        const double dk2dx = particles_var1b[i].x - particles_var1b[j].x;
                        const double dk2dy = particles_var1b[i].y - particles_var1b[j].y;
                        const double dk2dz = particles_var1b[i].z - particles_var1b[j].z;

                        const double rdk1 =  dx*dk1dx + dy*dk1dy + dz*dk1dz;
                        const double rdk2 =  dx*dk2dx + dy*dk2dy + dz*dk2dz;
                        const double dk1dk2 =  dk1dx*dk2dx + dk1dy*dk2dy + dk1dz*dk2dz;
                        dax     +=        3.* r5inv * dk2dx * rdk1
                        + 3.* r5inv * dk1dx * rdk2
                        + 3.* r5inv    * dx * dk1dk2  
                        - 15.      * dx * r7inv * rdk1 * rdk2;
                        day     +=        3.* r5inv * dk2dy * rdk1
                        + 3.* r5inv * dk1dy * rdk2
                        + 3.* r5inv    * dy * dk1dk2  
                        - 15.      * dy * r7inv * rdk1 * rdk2;
                        daz     +=        3.* r5inv * dk2dz * rdk1
                        + 3.* r5inv * dk1dz * rdk2
                        + 3.* r5inv    * dz * dk1dk2  
                        - 15.      * dz * r7inv * rdk1 * rdk2;

                        const double dk1Gmi = G * particles_var1a[i].m;
                        const double dk1Gmj = G * particles_var1a[j].m;
                        const double dk2Gmi = G * particles_var1b[i].m;
                        const double dk2Gmj = G * particles_var1b[j].m;

                        particles_var2[i].ax += Gmj * dax 
                        - ddGmj*r3inv*dx 
                        - dk2Gmj*r3inv*dk1dx + 3.*dk2Gmj*r5inv*dx*rdk1
                        - dk1Gmj*r3inv*dk2dx + 3.*dk1Gmj*r5inv*dx*rdk2;
                        particles_var2[i].ay += Gmj * day 
                        - ddGmj*r3inv*dy
                        - dk2Gmj*r3inv*dk1dy + 3.*dk2Gmj*r5inv*dy*rdk1
                        - dk1Gmj*r3inv*dk2dy + 3.*dk1Gmj*r5inv*dy*rdk2;
                        particles_var2[i].az += Gmj * daz 
                        - ddGmj*r3inv*dz
                        - dk2Gmj*r3inv*dk1dz + 3.*dk2Gmj*r5inv*dz*rdk1
                        - dk1Gmj*r3inv*dk2dz + 3.*dk1Gmj*r5inv*dz*rdk2;

                        particles_var2[j].ax -= Gmi * dax 
                        - ddGmi*r3inv*dx
                        - dk2Gmi*r3inv*dk1dx + 3.*dk2Gmi*r5inv*dx*rdk1
                        - dk1Gmi*r3inv*dk2dx + 3.*dk1Gmi*r5inv*dx*rdk2;
                        particles_var2[j].ay -= Gmi * day 
                        - ddGmi*r3inv*dy
                        - dk2Gmi*r3inv*dk1dy + 3.*dk2Gmi*r5inv*dy*rdk1
                        - dk1Gmi*r3inv*dk2dy + 3.*dk1Gmi*r5inv*dy*rdk2;
                        particles_var2[j].az -= Gmi * daz 
                        - ddGmi*r3inv*dz
                        - dk2Gmi*r3inv*dk1dz + 3.*dk2Gmi*r5inv*dz*rdk1
                        - dk1Gmi*r3inv*dk2dz + 3.*dk1Gmi*r5inv*dz*rdk2;
                    }
                }
            }else{ //testparticle
                size_t i = vc.testparticle;
                particles_var2[0].ax = 0.; 
                particles_var2[0].ay = 0.; 
                particles_var2[0].az = 0.; 
                for (size_t j=0; j<N; j++){
                    if (i==j) continue;
                    const double dx = particles[i].x - particles[j].x;
                    const double dy = particles[i].y - particles[j].y;
                    const double dz = particles[i].z - particles[j].z;
                    const double r2 = dx*dx + dy*dy + dz*dz;
                    const double r  = sqrt(r2);
                    const double r3inv = 1./(r2*r);
                    const double r5inv = r3inv/r2;
                    const double r7inv = r5inv/r2;
                    const double ddx = particles_var2[0].x;
                    const double ddy = particles_var2[0].y;
                    const double ddz = particles_var2[0].z;
                    const double Gmj = G * particles[j].m;

                    // Variational equations
                    // delta^(2) terms
                    double dax =         ddx * ( 3.*dx*dx*r5inv - r3inv )
                    + ddy * ( 3.*dx*dy*r5inv )
                    + ddz * ( 3.*dx*dz*r5inv );
                    double day =         ddx * ( 3.*dy*dx*r5inv )
                    + ddy * ( 3.*dy*dy*r5inv - r3inv )
                    + ddz * ( 3.*dy*dz*r5inv );
                    double daz =         ddx * ( 3.*dz*dx*r5inv )
                    + ddy * ( 3.*dz*dy*r5inv )
                    + ddz * ( 3.*dz*dz*r5inv - r3inv );

                    // delta^(1) delta^(1) terms
                    const double dk1dx = particles_var1a[0].x;
                    const double dk1dy = particles_var1a[0].y;
                    const double dk1dz = particles_var1a[0].z;
                    const double dk2dx = particles_var1b[0].x;
                    const double dk2dy = particles_var1b[0].y;
                    const double dk2dz = particles_var1b[0].z;

                    const double rdk1 =  dx*dk1dx + dy*dk1dy + dz*dk1dz;
                    const double rdk2 =  dx*dk2dx + dy*dk2dy + dz*dk2dz;
                    const double dk1dk2 =  dk1dx*dk2dx + dk1dy*dk2dy + dk1dz*dk2dz;
                    dax     +=        3.* r5inv * dk2dx * rdk1
                    + 3.* r5inv * dk1dx * rdk2
                    + 3.* r5inv    * dx * dk1dk2  
                    - 15.      * dx * r7inv * rdk1 * rdk2;
                    day     +=        3.* r5inv * dk2dy * rdk1
                    + 3.* r5inv * dk1dy * rdk2
                    + 3.* r5inv    * dy * dk1dk2  
                    - 15.      * dy * r7inv * rdk1 * rdk2;
                    daz     +=        3.* r5inv * dk2dz * rdk1
                    + 3.* r5inv * dk1dz * rdk2
                    + 3.* r5inv    * dz * dk1dk2  
                    - 15.      * dz * r7inv * rdk1 * rdk2;

                    // No variational mass contributions for test particles!

                    particles_var2[0].ax += Gmj * dax; 
                    particles_var2[0].ay += Gmj * day;
                    particles_var2[0].az += Gmj * daz;
                }
            }
        }
    }
    PROFILING_STOP(PROFILING_CAT_GRAVITY);
}

void reb_gravity_basic_calculate_and_apply_jerk(struct reb_simulation* r, const double v){
    PROFILING_START();
    struct reb_particle* const particles = r->particles;
    const size_t N = r->N;
    const double G = r->G;
    const size_t N_active = ((r->N_active==SIZE_MAX)?N:r->N_active);
    const int _testparticle_type   = r->testparticle_type;
    const size_t starti = (r->gravity_ignore_terms==REB_GRAVITY_IGNORE_TERMS_NONE)?1:2;
    const size_t startj = (r->gravity_ignore_terms==REB_GRAVITY_IGNORE_TERMS_INVOLVING_0)?1:0;
    // All interactions between active particles
#pragma omp parallel for
    for (size_t i=starti; i<N_active; i++){
#ifndef OPENMP
        if (reb_sigint > 1) return;
#endif // OPENMP
        for (size_t j=startj; j<i; j++){
            const double dx = particles[i].x - particles[j].x; 
            const double dy = particles[i].y - particles[j].y; 
            const double dz = particles[i].z - particles[j].z; 

            const double dax = particles[i].ax - particles[j].ax; 
            const double day = particles[i].ay - particles[j].ay; 
            const double daz = particles[i].az - particles[j].az; 

            const double dr = sqrt(dx*dx + dy*dy + dz*dz);
            const double alphasum = dax*dx+day*dy+daz*dz;
            const double prefact2 = 2.*v*G /(dr*dr*dr);
            const double prefact2i = prefact2*particles[j].m;
            const double prefact2j = prefact2*particles[i].m;
            const double prefact1 = alphasum*prefact2/dr *3./dr;
            const double prefact1i = prefact1*particles[j].m;
            const double prefact1j = prefact1*particles[i].m;
            particles[i].vx    += dx*prefact1i - dax*prefact2i;
            particles[i].vy    += dy*prefact1i - day*prefact2i;
            particles[i].vz    += dz*prefact1i - daz*prefact2i;
            particles[j].vx    += dax*prefact2j - dx*prefact1j;
            particles[j].vy    += day*prefact2j - dy*prefact1j;
            particles[j].vz    += daz*prefact2j - dz*prefact1j;
        }
    }
    // Interactions between active particles and test particles
#pragma omp parallel for
    for (size_t i=N_active; i<N; i++){
#ifndef OPENMP
        if (reb_sigint > 1) return;
#endif // OPENMP
        for (size_t j=startj; j<N_active; j++){ // Only active particles exert forces on test particles.
            const double dx = particles[i].x - particles[j].x; 
            const double dy = particles[i].y - particles[j].y; 
            const double dz = particles[i].z - particles[j].z; 

            const double dax = particles[i].ax - particles[j].ax; 
            const double day = particles[i].ay - particles[j].ay; 
            const double daz = particles[i].az - particles[j].az; 

            const double dr = sqrt(dx*dx + dy*dy + dz*dz);
            const double alphasum = dax*dx+day*dy+daz*dz;
            const double prefact2 = 2.*v*G /(dr*dr*dr);
            const double prefact1 = alphasum*prefact2/dr *3./dr;
            const double prefact1i = prefact1*particles[j].m;
            const double prefact2i = prefact2*particles[j].m;
            particles[i].vx    += dx*prefact1i - dax*prefact2i;
            particles[i].vy    += dy*prefact1i - day*prefact2i;
            particles[i].vz    += dz*prefact1i - daz*prefact2i;
            if (_testparticle_type){
                const double prefact1j = prefact1*particles[i].m;
                const double prefact2j = prefact2*particles[i].m;
                particles[j].vx    += dax*prefact2j - dx*prefact1j;
                particles[j].vy    += day*prefact2j - dy*prefact1j;
                particles[j].vz    += daz*prefact2j - dz*prefact1j;
            }
        }
    }
    PROFILING_STOP(PROFILING_CAT_GRAVITY);
}


