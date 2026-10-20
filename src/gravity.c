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

#pragma omp parallel for schedule(static) if(reb_omp_parallel(0.5*(double)N)) num_threads(reb_omp_threads())
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
#pragma omp parallel for schedule(guided) if(reb_omp_parallel(64.*(double)N)) num_threads(reb_omp_threads())
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

// Structure-of-arrays kernels shared by the serial-arrays entry point (all builds) and the
// OpenMP kernels below.
struct reb_gravity_soa {
    const double* x; const double* y; const double* z; const double* m;
    double* ax; double* ay; double* az;
};
// f[l] = G/r^3 with r = sqrt(dx^2+dy^2+dz^2+softening^2), the serial loop's expression; a helper
// with non-aliasing arrays so that the loop is vectorised (square roots, divisions).
static inline __attribute__((always_inline)) void reb_gravity_rinv3(const size_t n, const double* const REB_RESTRICT dx, const double* const REB_RESTRICT dy, const double* const REB_RESTRICT dz, const double G, const double softening, double* const REB_RESTRICT f){
    const double softening2 = softening*softening; // squared here: provably >= 0, so sqrt vectorises
    for (size_t l=0; l<n; l++){
        const double _r = sqrt(dx[l]*dx[l] + dy[l]*dy[l] + dz[l]*dz[l] + softening2);
        f[l] = G/(_r*_r*_r);
    }
}
// The serial half-pair loop for rows i in [ifirst, i1), partners j in [jlo, min(i, jhi)); the
// partner terms (to a_j) are added only if partner is set (test-particle rows of type 0 do not
// act on their partners).
static inline __attribute__((always_inline)) void reb_gravity_rows_ex(const struct reb_gravity_soa p, const size_t ifirst, const size_t i1, const size_t jlo, const size_t jhi, const int partner, const struct reb_vec6d gb, const double G, const double softening){
    double* const REB_RESTRICT pax = p.ax;
    double* const REB_RESTRICT pay = p.ay;
    double* const REB_RESTRICT paz = p.az;
    double bdx[REB_GRAVITY_BLOCK], bdy[REB_GRAVITY_BLOCK], bdz[REB_GRAVITY_BLOCK], bf[REB_GRAVITY_BLOCK];
    for (size_t i=ifirst; i<i1; i++){
        const double xi = gb.x+p.x[i];
        const double yi = gb.y+p.y[i];
        const double zi = gb.z+p.z[i];
        const double mi = p.m[i];
        double axi = pax[i], ayi = pay[i], azi = paz[i];
        const size_t jend = i<jhi ? i : jhi;
        for (size_t b0=jlo; b0<jend; b0+=REB_GRAVITY_BLOCK){
            const size_t nb = (jend-b0) < REB_GRAVITY_BLOCK ? (jend-b0) : REB_GRAVITY_BLOCK;
            for (size_t l=0; l<nb; l++){
                bdx[l] = xi - p.x[b0+l];
                bdy[l] = yi - p.y[b0+l];
                bdz[l] = zi - p.z[b0+l];
            }
            reb_gravity_rinv3(nb, bdx, bdy, bdz, G, softening, bf);
            for (size_t l=0; l<nb; l++){
                const double prefact = bf[l];
                const double prefactj = -prefact*p.m[b0+l];
                const double prefacti = prefact*mi;
                axi += prefactj*bdx[l];
                ayi += prefactj*bdy[l];
                azi += prefactj*bdz[l];
                if (partner){
                    pax[b0+l] += prefacti*bdx[l];
                    pay[b0+l] += prefacti*bdy[l];
                    paz[b0+l] += prefacti*bdz[l];
                }
            }
        }
        pax[i] = axi; pay[i] = ayi; paz[i] = azi;
    }
}
static inline __attribute__((always_inline)) void reb_gravity_rows(const struct reb_gravity_soa p, const size_t ifirst, const size_t i1, const size_t jlo, const struct reb_vec6d gb, const double G, const double softening){
    reb_gravity_rows_ex(p, ifirst, i1, jlo, SIZE_MAX, 1, gb, G, softening);
}
// Diagonal tile: the serial half-pair loop for rows i in [i0, i1), partners j in [max(i0,startj), i).
REB_TARGET_CLONES static void reb_gravity_tile_diagonal(const struct reb_gravity_soa p, const size_t i0, const size_t i1, const size_t starti, const size_t startj, const struct reb_vec6d gb, const double G, const double softening){
    reb_gravity_rows(p, i0>starti?i0:starti, i1, i0>startj?i0:startj, gb, G, softening);
}
// Rows i in [i0, i1) of the serial half-pair loop with all their partners j < i (all particles
// active, no ghost boxes): adds -f*m_j*dx to particle i and f*m_i*dx to particle j. Used by IAS15's
// stage-parallel mode, which splits a force evaluation into such row blocks, each accumulated
// into its own zeroed arrays and then summed in a fixed order.
REB_TARGET_CLONES void reb_gravity_basic_rows_arrays(const double* x, const double* y, const double* z, const double* m, double* ax, double* ay, double* az, const size_t i0, const size_t i1, const double G, const double softening){
    const struct reb_gravity_soa p = {x, y, z, m, ax, ay, az};
    const struct reb_vec6d gb = {0};
    reb_gravity_rows(p, i0, i1, 0, gb, G, softening);
}
// General row block for IAS15's stage mode: rows [i0, i1) (from at least starti), partners
// [startj, min(i, jhi)), partner terms only if partner, ghost-box offset gb.
REB_TARGET_CLONES void reb_gravity_basic_rows_arrays_ex(const double* x, const double* y, const double* z, const double* m, double* ax, double* ay, double* az, const size_t i0, const size_t i1, const size_t starti, const size_t startj, const size_t jhi, const int partner, const struct reb_vec6d gb, const double G, const double softening){
    const struct reb_gravity_soa p = {x, y, z, m, ax, ay, az};
    const size_t ifirst = i0>starti ? i0 : starti;
    if (partner){
        reb_gravity_rows_ex(p, ifirst, i1, startj, jhi, 1, gb, G, softening);
    }else{
        reb_gravity_rows_ex(p, ifirst, i1, startj, jhi, 0, gb, G, softening);
    }
}
// Compensated gravity for particles [k0, k1) on plain arrays: the same additions in the same order
// (Kahan summation over partners in ascending index order, then the test-particle back-reaction)
// as reb_gravity_compensated_calculate_acceleration, so the result is bitwise identical to it at
// the same positions. Each particle only receives its own terms: chunks are independent.
REB_TARGET_CLONES void reb_gravity_compensated_range_arrays(const double* x, const double* y, const double* z, const double* m, double* ax, double* ay, double* az, const size_t k0, const size_t k1, const size_t N, const size_t N_active, const int testparticle_type, const unsigned int gravity_ignore_terms, const double G, const double softening){
    const double softening2 = softening*softening;
    for (size_t i=k0; i<k1; i++){
        double axi = 0., ayi = 0., azi = 0., csx = 0., csy = 0., csz = 0.;
        const size_t jn = N_active<N ? N_active : N;
        for (size_t j=0; j<jn; j++){
            if (gravity_ignore_terms==REB_GRAVITY_IGNORE_TERMS_BETWEEN_0_AND_1 && ((j==1 && i==0) || (i==1 && j==0))) continue;
            if (gravity_ignore_terms==REB_GRAVITY_IGNORE_TERMS_INVOLVING_0 && ((j==0 || i==0))) continue;
            if (i==j) continue;
            const double dx = x[i] - x[j];
            const double dy = y[i] - y[j];
            const double dz = z[i] - z[j];
            const double r2 = dx*dx + dy*dy + dz*dz + softening2;
            const double r = sqrt(r2);
            const double prefact  = G/(r2*r);
            const double prefactj = -prefact*m[j];
            double t, yv;
            yv = prefactj*dx - csx; t = axi + yv; csx = (t - axi) - yv; axi = t;
            yv = prefactj*dy - csy; t = ayi + yv; csy = (t - ayi) - yv; ayi = t;
            yv = prefactj*dz - csz; t = azi + yv; csz = (t - azi) - yv; azi = t;
        }
        if (testparticle_type && i<N_active){
            for (size_t k=N_active; k<N; k++){
                // back-reaction of test particle k on active particle i (serial: j = i, i = k)
                if (gravity_ignore_terms==REB_GRAVITY_IGNORE_TERMS_BETWEEN_0_AND_1 && ((i==1 && k==0) || (k==1 && i==0))) continue;
                if (gravity_ignore_terms==REB_GRAVITY_IGNORE_TERMS_INVOLVING_0 && ((i==0 || k==0))) continue;
                const double dx = x[k] - x[i];
                const double dy = y[k] - y[i];
                const double dz = z[k] - z[i];
                const double r2 = dx*dx + dy*dy + dz*dz + softening2;
                const double r = sqrt(r2);
                const double prefact  = G/(r2*r);
                const double prefacti = prefact*m[k];
                double t, yv;
                yv = prefacti*dx - csx; t = axi + yv; csx = (t - axi) - yv; axi = t;
                yv = prefacti*dy - csy; t = ayi + yv; csy = (t - ayi) - yv; ayi = t;
                yv = prefacti*dz - csz; t = azi + yv; csz = (t - azi) - yv; azi = t;
            }
        }
        ax[i] = axi; ay[i] = ayi; az[i] = azi;
    }
}
// Serial half-pair gravity on plain arrays (all particles active, no ghost boxes, no ignored
// terms): the same additions in the same order as the serial loop in
// reb_gravity_basic_calculate_acceleration. Accelerations must be zero on entry. Used by
// IAS15's stage-parallel mode, which evaluates the forces of several substeps at once.
void reb_gravity_basic_arrays(const size_t N, const double* x, const double* y, const double* z, const double* m, double* ax, double* ay, double* az, const double G, const double softening){
    const struct reb_gravity_soa p = {x, y, z, m, ax, ay, az};
    const struct reb_vec6d gb = {0};
    reb_gravity_tile_diagonal(p, 0, N, 1, 0, gb, G, softening);
}
#ifdef OPENMP
// One-sided basic gravity for particle k (OpenMP builds): adds to particle k's acceleration the
// contributions of all its partners, in exactly the order and with exactly the floating-point
// expressions of the serial half-pair loop in reb_gravity_basic_calculate_acceleration:
//   1) k's own row (serial: i = k): partners j in [startj, min(k, N_active)), dx = (gb+x_k) - x_j,
//      a_k += (-f*m_j)*dx;
//   2) later active rows (serial: i > k adds to particle j = k): dx = (gb+x_i) - x_k, a_k += (f*m_i)*dx;
//   3) back-reaction of test particles (testparticle_type = 1), same form as 2), i >= N_active.
// with f = G/r^3, r^2 = dx^2+dy^2+dz^2+softening^2. Each particle only writes its own acceleration,
// so particles can be distributed over threads freely and the result is bitwise identical to the
// serial loop for any number of threads (also with ghost boxes). Costs N^2 instead of N^2/2 pairs.
// Separations and 1/r^3 are computed in blocks (vectorisable), the accumulation stays in order.
static inline __attribute__((always_inline)) void reb_gravity_basic_onesided_segment(const struct reb_particle* const particles, const size_t k, const size_t j0, const size_t j1, const int own_row, const struct reb_vec6d gb, const double G, const double softening, double* const ax, double* const ay, double* const az){
    const double softening2 = softening*softening; // squared here: provably >= 0, so sqrt vectorises
    const double xk = particles[k].x;
    const double yk = particles[k].y;
    const double zk = particles[k].z;
    const double gxk = gb.x+xk, gyk = gb.y+yk, gzk = gb.z+zk;
    double axk = *ax, ayk = *ay, azk = *az;
    double bdx[REB_GRAVITY_BLOCK], bdy[REB_GRAVITY_BLOCK], bdz[REB_GRAVITY_BLOCK], bpf[REB_GRAVITY_BLOCK];
    for (size_t b0=j0; b0<j1; b0+=REB_GRAVITY_BLOCK){
        const size_t nb = (j1-b0) < REB_GRAVITY_BLOCK ? (j1-b0) : REB_GRAVITY_BLOCK;
        if (own_row){
            for (size_t l=0; l<nb; l++){
                const struct reb_particle* const pj = &particles[b0+l];
                bdx[l] = gxk - pj->x;
                bdy[l] = gyk - pj->y;
                bdz[l] = gzk - pj->z;
            }
        }else{
            for (size_t l=0; l<nb; l++){
                const struct reb_particle* const pi = &particles[b0+l];
                bdx[l] = (gb.x+pi->x) - xk;
                bdy[l] = (gb.y+pi->y) - yk;
                bdz[l] = (gb.z+pi->z) - zk;
            }
        }
        for (size_t l=0; l<nb; l++){
            const double _r = sqrt(bdx[l]*bdx[l] + bdy[l]*bdy[l] + bdz[l]*bdz[l] + softening2);
            bpf[l] = G/(_r*_r*_r);
        }
        if (own_row){
            for (size_t l=0; l<nb; l++){
                const double prefactj = -bpf[l]*particles[b0+l].m;
                axk += prefactj*bdx[l];
                ayk += prefactj*bdy[l];
                azk += prefactj*bdz[l];
            }
        }else{
            for (size_t l=0; l<nb; l++){
                const double prefacti = bpf[l]*particles[b0+l].m;
                axk += prefacti*bdx[l];
                ayk += prefacti*bdy[l];
                azk += prefacti*bdz[l];
            }
        }
    }
    *ax = axk; *ay = ayk; *az = azk;
}
static inline __attribute__((always_inline)) void reb_gravity_basic_onesided(struct reb_particle* const particles, const size_t k, const size_t N, const size_t N_active, const size_t starti, const size_t startj, const int testparticle_type, const struct reb_vec6d gb, const double G, const double softening){
    double ax = particles[k].ax, ay = particles[k].ay, az = particles[k].az;
    if (k>=starti){ // own row (active particles and test particles)
        const size_t jend = k<N_active ? k : N_active;
        if (jend>startj) reb_gravity_basic_onesided_segment(particles, k, startj, jend, 1, gb, G, softening, &ax, &ay, &az);
    }
    if (k<N_active && k>=startj){ // later rows: active partners, then test particles (type 1)
        const size_t i0 = (k+1>starti) ? k+1 : starti;
        if (N_active>i0) reb_gravity_basic_onesided_segment(particles, k, i0, N_active, 0, gb, G, softening, &ax, &ay, &az);
        if (testparticle_type){
            const size_t t0 = (N_active>starti) ? N_active : starti;
            if (N>t0) reb_gravity_basic_onesided_segment(particles, k, t0, N, 0, gb, G, softening, &ax, &ay, &az);
        }
    }
    particles[k].ax = ax; particles[k].ay = ay; particles[k].az = az;
}
// Particles [k0, k1); own function (not inlined into the outlined OpenMP region) so that it
// gets the AVX2 clone, called once per chunk rather than once per particle.
REB_TARGET_CLONES static void reb_gravity_basic_onesided_range(struct reb_particle* const particles, const size_t k0, const size_t k1, const size_t N, const size_t N_active, const size_t starti, const size_t startj, const int testparticle_type, const struct reb_vec6d gb, const double G, const double softening){
    for (size_t k=k0; k<k1; k++){
        reb_gravity_basic_onesided(particles, k, N, N_active, starti, startj, testparticle_type, gb, G, softening);
    }
}
#define REB_GRAVITY_ONESIDED_CHUNK 64

// Two-pass parallel basic gravity (OpenMP builds), bitwise identical to the serial loop.
//
// The serial loop adds to particle k, in this order: its own row (partners j < k, ascending), then
// one term from every later row i > k (ascending), then the test-particle rows (testparticle_type 1).
// Here the active particles are split into row blocks [s_b, s_{b+1}). Each block is handled by one
// thread at a time, and every particle still receives exactly the serial terms in the serial order:
//   pass 1, block b: (a) for each row i of the block, partners j < s_b: compute f = G/r^3, add the
//            own-row term a_i += (-f*m_j)*((gb+x_i) - x_j) and store f in the table (row-major,
//            sequential writes); (b) the diagonal tile, partners s_b <= j < i: the serial half-pair
//            loop restricted to the block (own term to a_i, partner term to a_j), no table needed;
//   pass 2, block b: the later rows i >= s_{b+1} (and test-particle rows), ascending:
//            a_k += (f*m_i)*((gb+x_i) - x_k) for all k of the block, from the table (a contiguous
//            segment per row; vectorised across the block, independent per-particle sums).
// Test particles (type 0) only have their own row, computed in pass 1 without storing it.
// The square roots and divisions are done once per pair, as in the serial code. The block layout
// does not change the result, so it is also independent of the number of threads.
// Positions, masses and accelerations are kept in contiguous arrays during the computation.
#ifndef REB_GRAVITY_TABLE_MAX
#define REB_GRAVITY_TABLE_MAX (1u<<23)  // table entries (64 MB); larger systems use the one-sided loop
#endif
#define REB_GRAVITY_MAX_BLOCKS 256
#define REB_GRAVITY_PASS2_CHUNK 64
// Off-diagonal part of row i: partners j in [j0, j1); stores f into frow[j-j0] if frow != NULL.
REB_TARGET_CLONES static void reb_gravity_tile_row(const struct reb_gravity_soa p, const size_t i, const size_t j0, const size_t j1, const struct reb_vec6d gb, const double G, const double softening, double* const REB_RESTRICT frow){
    const double xi = gb.x+p.x[i];
    const double yi = gb.y+p.y[i];
    const double zi = gb.z+p.z[i];
    const double* const REB_RESTRICT px = p.x;
    const double* const REB_RESTRICT py = p.y;
    const double* const REB_RESTRICT pz = p.z;
    const double* const REB_RESTRICT pm = p.m;
    double axi = p.ax[i], ayi = p.ay[i], azi = p.az[i];
    double bdx[REB_GRAVITY_BLOCK], bdy[REB_GRAVITY_BLOCK], bdz[REB_GRAVITY_BLOCK], bf[REB_GRAVITY_BLOCK];
    for (size_t b0=j0; b0<j1; b0+=REB_GRAVITY_BLOCK){
        const size_t nb = (j1-b0) < REB_GRAVITY_BLOCK ? (j1-b0) : REB_GRAVITY_BLOCK;
        for (size_t l=0; l<nb; l++){
            bdx[l] = xi - px[b0+l];
            bdy[l] = yi - py[b0+l];
            bdz[l] = zi - pz[b0+l];
        }
        reb_gravity_rinv3(nb, bdx, bdy, bdz, G, softening, bf);
        if (frow){
            double* const REB_RESTRICT fb = frow + (b0-j0);
            for (size_t l=0; l<nb; l++) fb[l] = bf[l];
        }
        for (size_t l=0; l<nb; l++){
            const double prefactj = -bf[l]*pm[b0+l];
            axi += prefactj*bdx[l];
            ayi += prefactj*bdy[l];
            azi += prefactj*bdz[l];
        }
    }
    p.ax[i] = axi; p.ay[i] = ayi; p.az[i] = azi;
}
// One later row for n partners: a_k += (f_k*m_i)*(gx_i - x_k) (separate function so that the
// compiler sees independent, non-aliasing arrays and vectorises the loop).
static inline __attribute__((always_inline)) void reb_gravity_later_row(const size_t n, const double* const REB_RESTRICT f, const double mi, const double gxi, const double gyi, const double gzi, const double* const REB_RESTRICT xk, const double* const REB_RESTRICT yk, const double* const REB_RESTRICT zk, double* const REB_RESTRICT ax, double* const REB_RESTRICT ay, double* const REB_RESTRICT az){
    for (size_t l=0; l<n; l++){
        const double prefacti = f[l]*mi;
        ax[l] += prefacti*(gxi - xk[l]);
        ay[l] += prefacti*(gyi - yk[l]);
        az[l] += prefacti*(gzi - zk[l]);
    }
}
// Pass 2 for block [k0, k1): later rows i in [i0, i1), row i's factors for partners
// [kbase, ...) at table + rowoff(i); partners k in [max(k0,startj), k1).
REB_TARGET_CLONES static void reb_gravity_tile_later(const struct reb_gravity_soa p, const size_t k0, const size_t k1, const size_t i0, const size_t i1, const double* const table, const size_t* const rowoff, const size_t startj, const struct reb_vec6d gb){
    const size_t ks = k0>startj ? k0 : startj;
    if (ks>=k1) return;
    double* const REB_RESTRICT pax = p.ax + ks;
    double* const REB_RESTRICT pay = p.ay + ks;
    double* const REB_RESTRICT paz = p.az + ks;
    const double* const REB_RESTRICT xk = p.x + ks;
    const double* const REB_RESTRICT yk = p.y + ks;
    const double* const REB_RESTRICT zk = p.z + ks;
    const size_t n = k1-ks;
    for (size_t i=i0; i<i1; i++){
        reb_gravity_later_row(n, table + rowoff[i] + (ks-startj), p.m[i], gb.x+p.x[i], gb.y+p.y[i], gb.z+p.z[i], xk, yk, zk, pax, pay, paz);
    }
}
// Threads for the parallel basic-gravity kernels: about N_active/64 when the active pairs dominate
// (more threads than that cost more in synchronisation and cross-core traffic than they save, as
// measured on a 20-thread desktop CPU); all threads when test-particle rows dominate. Results do not
// depend on the thread count.
static int reb_gravity_omp_threads(const size_t N, const size_t N_active){
    const int T = reb_omp_threads();
    if (N-N_active >= N_active) return T;
    size_t c = N_active/REB_OMP_GRAVITY_ROWS_PER_THREAD;
    if (c<2) c = 2;
    return c<(size_t)T ? (int)c : T;
}
// Returns 0 if the table would be too large (the caller then uses the one-sided loop).
static int reb_gravity_basic_two_pass(struct reb_simulation* const r, const size_t N, const size_t N_active, const size_t starti, const size_t startj, const int testparticle_type, const struct reb_vec6d gb, const double G, const double softening){
    if (starti>N_active || startj>=starti) return 0;
    // Row blocks of the active particles, boundaries at equal pass-1 work (~ s^2), at least 8 rows.
    const int nthreads = reb_gravity_omp_threads(N, N_active);
    size_t nb = 4*(size_t)nthreads;
    if (nb>REB_GRAVITY_MAX_BLOCKS) nb = REB_GRAVITY_MAX_BLOCKS;
    if (nb>N_active/8) nb = N_active/8;
    if (nb<2) return 0;
    size_t s[REB_GRAVITY_MAX_BLOCKS+1];
    for (size_t b=0; b<=nb; b++){
        size_t v = (size_t)((double)N_active*sqrt((double)b/(double)nb));
        if (b==nb) v = N_active;
        s[b] = v;
    }
    for (size_t b=1; b<=nb; b++){ if (s[b]<s[b-1]+1) s[b] = s[b-1]+1; } // non-empty blocks
    if (s[nb]>N_active) return 0;
    // Table: row i (active, in block b) stores partners [startj, s_b); test-particle rows (type 1)
    // store [startj, N_active). rowoff[i] = offset of row i.
    const size_t Nrows = testparticle_type ? N : N_active;
    const size_t need_rowoff = Nrows;
    size_t size = 0;
    // first pass over rows to get the size (O(N))
    for (size_t b=0; b<nb; b++){
        const size_t w = s[b]>startj ? s[b]-startj : 0;
        size += (s[b+1]-s[b])*w;
    }
    if (testparticle_type) size += (N-N_active)*(N_active-startj);
    if (size>REB_GRAVITY_TABLE_MAX) return 0;
    const size_t need = size + 7*N + (need_rowoff*sizeof(size_t)+sizeof(double)-1)/sizeof(double);
    if (r->N_allocated_gravity_pair_table<need){
        free(r->gravity_pair_table);
        r->gravity_pair_table = malloc(sizeof(double)*need);
        r->N_allocated_gravity_pair_table = r->gravity_pair_table ? need : 0;
        if (!r->gravity_pair_table) return 0;
    }
    double* const table = r->gravity_pair_table;
    double* const soa = table + size;
    size_t* const rowoff = (size_t*)(soa + 7*N);
    const struct reb_gravity_soa p = {soa, soa+N, soa+2*N, soa+3*N, soa+4*N, soa+5*N, soa+6*N};
    {
        size_t off = 0;
        for (size_t b=0; b<nb; b++){
            const size_t w = s[b]>startj ? s[b]-startj : 0;
            for (size_t i=s[b]; i<s[b+1]; i++){ rowoff[i] = off; off += w; }
        }
        for (size_t i=N_active; i<Nrows; i++){ rowoff[i] = off; off += N_active-startj; }
    }
    // Pass-2 chunks, ordered by block (block 0 = most later rows first).
    size_t chunk_first[REB_GRAVITY_MAX_BLOCKS+1];
    size_t nchunks = 0;
    for (size_t b=0; b<nb; b++){
        chunk_first[b] = nchunks;
        nchunks += (s[b+1]-s[b]+REB_GRAVITY_PASS2_CHUNK-1)/REB_GRAVITY_PASS2_CHUNK;
    }
    chunk_first[nb] = nchunks;
    unsigned short chunk_block_buf[1024];
    unsigned short* const chunk_block = nchunks<=1024 ? chunk_block_buf : malloc(sizeof(unsigned short)*nchunks);
    for (size_t b=0; b<nb; b++){
        for (size_t c=chunk_first[b]; c<chunk_first[b+1]; c++) chunk_block[c] = (unsigned short)b;
    }
    struct reb_particle* const particles = r->particles;
    // Full team (same size as every other region); threads 0..nthreads-1 do the work.
    size_t next_block = 0, next_chunk = 0;
    const size_t ntp = N>N_active ? N-N_active : 0;
#pragma omp parallel num_threads(reb_omp_threads())
    {
        const int tid = omp_get_thread_num();
        const int nw = nthreads<omp_get_num_threads() ? nthreads : omp_get_num_threads();
        const int worker = tid<nw;
        size_t lo = 0, hi = 0;
        if (worker){
            reb_omp_static_range(N, nw, tid, &lo, &hi);
            for (size_t i=lo; i<hi; i++){
                soa[i] = particles[i].x; soa[N+i] = particles[i].y; soa[2*N+i] = particles[i].z; soa[3*N+i] = particles[i].m;
                soa[4*N+i] = particles[i].ax; soa[5*N+i] = particles[i].ay; soa[6*N+i] = particles[i].az;
            }
        }
#pragma omp barrier
        if (worker){
            // Pass 1, active blocks (largest work last in index, so hand them out first).
            for (size_t bb; (bb = reb_omp_next(&next_block)) < nb; ){
                const size_t b = nb-1-bb;
                const size_t i0 = s[b], i1 = s[b+1];
                if (s[b]>startj){
                    for (size_t i=(i0>starti?i0:starti); i<i1; i++){
                        reb_gravity_tile_row(p, i, startj, s[b], gb, G, softening, table + rowoff[i]);
                    }
                }
                reb_gravity_tile_diagonal(p, i0, i1, starti, startj, gb, G, softening);
            }
            // Pass 1, test-particle rows (own row only; stored for type 1).
            reb_omp_static_range(ntp, nw, tid, &lo, &hi);
            for (size_t i=N_active+lo; i<N_active+hi; i++){
                if (i<starti || N_active<=startj) continue;
                reb_gravity_tile_row(p, i, startj, N_active, gb, G, softening, testparticle_type ? table + rowoff[i] : NULL);
            }
        }
#pragma omp barrier
        if (worker){
            // Pass 2: later rows for each active block, in chunks of REB_GRAVITY_PASS2_CHUNK
            // particles (independent: each particle only receives its own terms). Early blocks
            // have the most later rows, so chunks are handed out dynamically, largest first.
            for (size_t c; (c = reb_omp_next(&next_chunk)) < nchunks; ){
                const size_t b = chunk_block[c];
                const size_t k0 = s[b] + (c-chunk_first[b])*REB_GRAVITY_PASS2_CHUNK;
                const size_t k1 = (k0+REB_GRAVITY_PASS2_CHUNK<s[b+1]) ? k0+REB_GRAVITY_PASS2_CHUNK : s[b+1];
                const size_t ifirst = s[b+1]>starti ? s[b+1] : starti;
                reb_gravity_tile_later(p, k0, k1, ifirst, Nrows, table, rowoff, startj, gb);
            }
        }
#pragma omp barrier
        if (worker){
            reb_omp_static_range(N, nw, tid, &lo, &hi);
            for (size_t i=lo; i<hi; i++){
                particles[i].ax = soa[4*N+i]; particles[i].ay = soa[5*N+i]; particles[i].az = soa[6*N+i];
            }
        }
    }
    if (chunk_block!=chunk_block_buf) free(chunk_block);
    return 1;
}
#endif // OPENMP

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
    const size_t starti = (gravity_ignore_terms==REB_GRAVITY_IGNORE_TERMS_NONE)?1:2;
    const size_t startj = (gravity_ignore_terms==REB_GRAVITY_IGNORE_TERMS_INVOLVING_0)?1:0;
#ifdef OPENMP
    // Work = number of pair interactions. Pairs between active particles cost ~1.5x more in the
    // parallel kernel than in the serial loop; test-particle rows cost the same (each pair is
    // computed once either way), so they pay off with less work.
    const int parallel = reb_omp_parallel_gravity((double)N_active*(double)N_active)
                      || reb_omp_parallel((double)(N-N_active)*(double)N_active);
#endif // OPENMP
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
#ifdef OPENMP
                if (!parallel)
#endif // OPENMP
                // Serial: O(1/2*N^2), each pair once (Newton's third law)
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
#ifdef OPENMP
                // Parallel: one-sided, each particle sums all its partners in the serial order
                // (bitwise identical to the serial loop, see reb_gravity_basic_onesided).
                else if (!reb_gravity_basic_two_pass(r, N, N_active, starti, startj, _testparticle_type, gb, G, r->softening)){
                    const size_t nchunk = (N+REB_GRAVITY_ONESIDED_CHUNK-1)/REB_GRAVITY_ONESIDED_CHUNK;
                    const int nthreads = reb_gravity_omp_threads(N, N_active);
                    // Full team, threads 0..nthreads-1 work (see reb_omp_static_range).
#pragma omp parallel num_threads(reb_omp_threads())
                    {
                        const int tid = omp_get_thread_num();
                        const int nw = nthreads<omp_get_num_threads() ? nthreads : omp_get_num_threads();
                        if (tid<nw){
                            size_t c0, c1;
                            reb_omp_static_range(nchunk, nw, tid, &c0, &c1);
                            for (size_t c=c0; c<c1; c++){
                                const size_t k0 = c*REB_GRAVITY_ONESIDED_CHUNK;
                                const size_t k1 = (k0+REB_GRAVITY_ONESIDED_CHUNK<N) ? k0+REB_GRAVITY_ONESIDED_CHUNK : N;
                                reb_gravity_basic_onesided_range(particles, k0, k1, N, N_active, starti, startj, _testparticle_type, gb, G, r->softening);
                            }
                        }
                    }
                }
#endif // OPENMP
       // Interactions of test particles with active particles
#ifdef OPENMP
                if (!parallel)   // the parallel branch above handled test particles too
#endif // OPENMP
                {
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
                }
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
#pragma omp parallel for schedule(static) if(reb_omp_parallel(0.5*(double)N)) num_threads(reb_omp_threads())
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
    // Work = pair interactions in the one-sided loops (each particle sums all partners in index
    // order, the same order as the serial loop below, so results are identical).
    if (reb_omp_parallel_gravity((double)N*(double)N_active)){
#pragma omp parallel for schedule(guided) num_threads(reb_omp_threads())
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
#pragma omp parallel for schedule(guided) num_threads(reb_omp_threads())
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
#pragma omp parallel for schedule(guided) num_threads(reb_omp_threads())
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
    }else
#endif // OPENMP
    {
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
    }
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
    // Not parallel: each pair also updates particle j, which another thread may own (data race).
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
#pragma omp parallel for schedule(static) if(!_testparticle_type && reb_omp_parallel((double)(N-N_active)*(double)N_active)) num_threads(reb_omp_threads())
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


