/**
 * @file 	simulation.c
 * @brief 	Functions which create, free, and operate on simulations.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section 	LICENSE
 * Copyright (c) 2026 Hanno Rein
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
#include <stdio.h>
#include <math.h>
#include <string.h>
#ifdef MPI
#include <mpi.h>
#include "communication_mpi.h"
#endif // MPI
#ifdef OPENMP
#include <omp.h>
#endif // OPENMP
#include "tree.h"
#include "particle.h"
#include "simulation.h"
#include "simulationarchive.h"
#include "output.h"
#include "server.h"
#include "display.h"
#include "boundary.h"
#include "binarydata.h"
#include "fmemopen.h" // own implementation of fmemopen
#include "collision.h"
#include "tools.h"
#include "gravity.h"
#include "integrator_bs.h"
#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#include <emscripten/html5.h>
#endif // __EMSCRIPTEN__

#define MAX(a, b) ((a) < (b) ? (b) : (a))       ///< Returns the maximum of a and b

// Used for visualization
struct reb_thread_info {
    struct reb_simulation* r;
    double tmax;
};

// For unit tests to check if python struct has same size as c struct.
size_t reb_simulation_struct_size(){
    return sizeof(struct reb_simulation);
}

// User triggered stop.
void reb_simulation_stop(struct reb_simulation* const r){
    r->status = REB_STATUS_USER;
}

// Print or save warning message. r can be NULL.
void reb_simulation_warning(struct reb_simulation* const r, const char* const msg){
    int save_messages = 0;
    if (r != NULL && r->save_messages) save_messages = 1;
    reb_message(&r->messages, save_messages, REB_MESSAGE_TYPE_WARNING, msg);
}

// Print or save error message. r can be NULL.
void reb_simulation_error(struct reb_simulation* const r, const char* const msg){
    int save_messages = 0;
    if (r != NULL && r->save_messages) save_messages = 1;
    reb_message(&r->messages, save_messages, REB_MESSAGE_TYPE_ERROR, msg);
}

// Print or save info message. r can be NULL.
void reb_simulation_info(struct reb_simulation* const r, const char* const msg){
    int save_messages = 0;
    if (r != NULL && r->save_messages) save_messages = 1;
    reb_message(&r->messages, save_messages, REB_MESSAGE_TYPE_INFO, msg);
}

// One bare timestep. Without heartbeat.
static void reb_simulation_step(struct reb_simulation* const r);

struct reb_simulation* reb_simulation_create(){
    struct reb_simulation* r = calloc(1,sizeof(struct reb_simulation));
    // Init/reset values to default
    // Note: Only need to set non-zero values because
    //       we memset everything to zero above.
    r->rand_seed = reb_tools_get_rand_seed();
    r->G                    = 1;
    r->dt                   = 0.001;
    r->gravity              = REB_GRAVITY_BASIC;
    r->root_size            = -1;
    r->opening_angle2       = 0.25;
    r->OMEGAZ               = -1.0;
    r->N_root_x             = 1;
    r->N_root_y             = 1;
    r->N_root_z             = 1;
    r->N_active             = SIZE_MAX; 
    r->N_targets            = SIZE_MAX; 
    r->exact_finish_time    = 1;
    r->output_timing_last   = -1;
    r->simulationarchive_version = 5;
    r->is_synchronized = 1;
    reb_simulation_set_integrator(r, "ias15");

    return r;
}

void reb_simulation_free(struct reb_simulation* const r){
    // Free everything in the simulation first
    if (r->integrator.callbacks.free){
        r->integrator.callbacks.free(r->integrator.state);
    }
    r->gravity = REB_GRAVITY_BASIC; // Some integrators set their own gravity routine. Resetting.
    r->gravity_ignore_terms = REB_GRAVITY_IGNORE_TERMS_NONE;
    if (r->simulationarchive_filename){
        free(r->simulationarchive_filename);
    }
    if(r->display_settings){
        free(r->display_settings);
    }
#ifdef OPENGL
    if(r->display_data){
        // Waiting for visualization to shut down.
        if (r->display_data->window){ // Not needed under normal circumstances
            usleep(100);
        }
        if (r->display_data->window){ // still running?
            reb_simulation_info(NULL, "Waiting for OpenGL visualization to shut down...");
            while(r->display_data->window){
                usleep(100);
            }
        }
        pthread_mutex_destroy(&(r->display_data->mutex));
        if (r->display_data->r_copy){
            reb_simulation_free(r->display_data->r_copy);
        }
        if (r->display_data->particle_data){
            free(r->display_data->particle_data);
        }
        if (r->display_data->orbit_data){
            free(r->display_data->orbit_data);
        }
        free(r->display_data);
    }
#endif //OPENGL
#ifdef SERVER
    reb_simulation_stop_server(r);
#endif // SERVER
    free(r->gravity_cs);
    free(r->gravity_pair_table);
    free(r->collisions);
    reb_tree_delete(r);
    reb_tree_free_cells(r);
    free(r->tree_root);
    if(r->free_particle_ap){
        for(size_t i=0; i<r->N; i++){
            r->free_particle_ap(&r->particles[i]);
        }
    }
    free(r->particles);
    free(r->particles_var);
    for (size_t i=0;i<r->N_name_list;i++){
        free(r->name_list[i]);
    }
    free(r->name_list);
    free(r->name_hash_table);
    if (r->messages){
        for (size_t i=0;i<reb_messages_max_N;i++){
            free(r->messages[i]);
        }
    }
    if (r->messages){
        free(r->messages);
    }
    if (r->extras_cleanup){
        r->extras_cleanup(r);
    }
    if (r->var_config){
        free(r->var_config);
    }
    free(r->odes);

    // Free simulation itself
    free(r);
}


static void* set_integrator(struct reb_simulation* r, const char* name, struct reb_integrator integrator){
    // Memory for name is constant. Does not need to be released.
    if (r->integrator.callbacks.free && r->integrator.state){
        r->integrator.callbacks.free(r->integrator.state);
    }
    r->integrator.callbacks = integrator;
    r->integrator.name = name;
    if (r->integrator.callbacks.create){
        r->integrator.state = r->integrator.callbacks.create();
    }else{
        r->integrator.state = NULL;
    }
    return r->integrator.state; // for convenience, return pointer to data
}


void* reb_simulation_set_integrator(struct reb_simulation* r, const char* name){
    if (r->integrator.name && !r->is_synchronized){
        reb_simulation_warning(r, "Changing integrators while simulation is not synchronized results in undefined behaviour.");
    }
    // All built-in integrators
#define X(iname) if (strcmp(name, #iname)==0) { return set_integrator(r, #iname, reb_integrator_##iname); }
    REB_BUILTIN_INTEGRATORS
#undef X
    if (reb_integrator_configurations_custom){
        int i=0;
        while(reb_integrator_configurations_custom[i].name){ // {0} terminated array
            const struct reb_integrator_configuration integrator = reb_integrator_configurations_custom[i];
            if (strcmp(integrator.name, name)==0){
                return set_integrator(r, integrator.name, integrator.callbacks);
            }
            i++;
        }
    }

    reb_simulation_error(r, "Integrator not found.");
    return NULL;
}

// Heartbeat wrapper. Runs actual heartbeat and does exit checks.
static void run_heartbeat(struct reb_simulation* const r){
    if (r->heartbeat){ r->heartbeat(r); }               // Heartbeat
    if (r->exit_max_distance){
        // Check for escaping particles
        const double max2 = r->exit_max_distance * r->exit_max_distance;
        const struct reb_particle* const particles = r->particles;
        const size_t N = r->N;
        for (size_t i=0;i<N;i++){
            struct reb_particle p = particles[i];
            double r2 = p.x*p.x + p.y*p.y + p.z*p.z;
            if (r2>max2){
                r->status = REB_STATUS_ESCAPE;
            }
        }
    }
    if (r->exit_min_distance){
        // Check for close encounters
        const double min2 = r->exit_min_distance * r->exit_min_distance;
        const struct reb_particle* const particles = r->particles;
        const size_t N = r->N;
        for (size_t i=0;i<N;i++){
            struct reb_particle pi = particles[i];
            for (size_t j=0;j<i;j++){
                struct reb_particle pj = particles[j];
                const double x = pi.x-pj.x;
                const double y = pi.y-pj.y;
                const double z = pi.z-pj.z;
                const double r2 = x*x + y*y + z*z;
                if (r2<min2){
                    r->status = REB_STATUS_ENCOUNTER;
                }
            }
        }
    }
}

static int error_message_waiting(struct reb_simulation* const r){
    if (r->messages){
        for (size_t i=0;i<reb_messages_max_N;i++){
            if (r->messages[i]!=NULL){
                if (r->messages[i][0]=='e'){
                    return 1;
                }
            }
        }
    }
    return 0;
}


static int reb_check_exit(struct reb_simulation* const r, const double tmax, double* last_full_dt){
    if(r->status <= REB_STATUS_SINGLE_STEP){
        if(r->status == REB_STATUS_SINGLE_STEP){
            r->status = REB_STATUS_PAUSED;
        }else{
            // This allows an arbitrary number of steps before the simulation is paused
            r->status++;
        }
    }
    while(r->status == REB_STATUS_PAUSED || r->status == REB_STATUS_SCREENSHOT){
        // Wait for user to disable paused simulation
#ifdef __EMSCRIPTEN__
        emscripten_sleep(100);
#else
        usleep(1000);
#endif 
        if (reb_sigint){ // cancel while paused
            r->status = REB_STATUS_SIGINT;
        }
    }
    const double dtsign = copysign(1.,r->dt);   // Used to determine integration direction
    if (error_message_waiting(r)){
        r->status = REB_STATUS_GENERIC_ERROR;
    }
    if (r->status>=0){
        // Exit now.
    }else if(tmax!=INFINITY){
        if(r->exact_finish_time==1){
            if ((r->t+r->dt)*dtsign>=tmax*dtsign){  // Next step would overshoot
                if (r->t==tmax){
                    r->status = REB_STATUS_SUCCESS;
                }else if(r->status == REB_STATUS_LAST_STEP){
                    double tscale = 1e-12*fabs(tmax);   // Find order of magnitude for time
                    if (tscale<1e-200){     // Failsafe if tmax==0.
                        tscale = 1e-12;
                    }
                    if (fabs(r->t-tmax)<tscale){
                        r->status = REB_STATUS_SUCCESS;
                    }else{
                        // not there yet, do another step.
                        reb_simulation_synchronize(r);  // Need to synchronize because step size changed.
                        r->dt = tmax-r->t;
                    }
                }else{
                    r->status = REB_STATUS_LAST_STEP;   // Do one small step, then exit.
                    reb_simulation_synchronize(r);      // Need to synchronize because step size changed.
                    if (r->dt_last_done!=0.){           // If first timestep is also last, do not use dt_last_done (which would be 0.)
                        *last_full_dt = r->dt_last_done; // store last full dt before decreasing the timestep to match finish time
                    }
                    r->dt = tmax-r->t;
                }
            }else{
                if (r->status == REB_STATUS_LAST_STEP){
                    // This will get executed if an adaptive integrator reduces
                    // the timestep in what was supposed to be the last timestep.
                    r->status = REB_STATUS_RUNNING;
                }
            }
        }else{
            if (r->t*dtsign>=tmax*dtsign){  // Past the integration time
                r->status = REB_STATUS_SUCCESS; // Exit now.
            }
        }
    }
#ifndef MPI
    if (!r->N && !r->N_odes){
        reb_simulation_warning(r,"No particles in simulation. Will exit.");
        r->status = REB_STATUS_NO_PARTICLES; // Exit now.
    }
#else
    int status_max = 0;
    MPI_Allreduce(&(r->status), &status_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); 
    if (status_max>=0){
        r->status = status_max;
    }

#endif // MPI
    return r->status;
}


// Wall-time bookkeeping for the integration loops. Reading the clock costs tens of
// nanoseconds, comparable to a whole WHFast step for small systems, so the clock is read
// only every `stride` steps; the stride adapts so that reads are >= ~50 us apart (slow
// steps are still timed individually). The elapsed time is spread evenly over the steps
// of a batch. r->walltime is exact at every read and lags by at most ~50 us in between.
struct reb_walltime_tracker {
    struct reb_timeval t_last;
    uint64_t stride;
    uint64_t countdown;
};

static void reb_walltime_tracker_start(struct reb_walltime_tracker* const w){
    gettimeofday(&w->t_last, NULL);
    w->stride = 1;
    w->countdown = 1;
}

static void reb_walltime_tracker_record(struct reb_simulation* const r, struct reb_walltime_tracker* const w, const uint64_t N_steps){
    struct reb_timeval now;
    gettimeofday(&now, NULL);
    const double dt = now.tv_sec-w->t_last.tv_sec+(now.tv_usec-w->t_last.tv_usec)/1e6;
    w->t_last = now;
    if (N_steps==0) return;
    r->walltime_last_step = dt/(double)N_steps;
    r->walltime_last_steps_sum += dt;
    r->walltime_last_steps_N += (int)N_steps;
    if (r->walltime_last_steps_sum > 0.1){
        r->walltime_last_steps = r->walltime_last_steps_sum/r->walltime_last_steps_N;
        r->walltime_last_steps_sum =0;
        r->walltime_last_steps_N = 0;
    }
    r->walltime += dt;
    // Adapt the stride
    if (dt < 50e-6 && w->stride < 65536){
        w->stride *= 2;
    }else if (dt > 1e-3 && w->stride > 1){
        w->stride /= 2;
    }
}

// Call after every step.
static inline void reb_walltime_tracker_step(struct reb_simulation* const r, struct reb_walltime_tracker* const w){
    if (--w->countdown) return;
    reb_walltime_tracker_record(r, w, w->stride);
    w->countdown = w->stride;
}

// Call when the loop ends: accounts for the steps since the last clock read.
static void reb_walltime_tracker_finish(struct reb_simulation* const r, struct reb_walltime_tracker* const w){
    reb_walltime_tracker_record(r, w, w->stride - w->countdown);
}

static void* reb_simulation_integrate_raw(void* args){
    reb_integration_begin();
    struct reb_thread_info* thread_info = (struct reb_thread_info*)args;
    struct reb_simulation* const r = thread_info->r;

    if (thread_info->tmax != r->t){
        int dt_sign = (thread_info->tmax > r->t) ? 1.0 : -1.0; // determine integration direction
        r->dt = copysign(r->dt, dt_sign);
    }

    double last_full_dt = r->dt; // need to store r->dt in case timestep gets artificially shrunk to meet exact_finish_time=1
    r->dt_last_done = 0.; // Reset in case first timestep attempt will fail

    if (r->testparticle_hidewarnings==0 && reb_particle_check_testparticles(r)){
        reb_simulation_warning(r,"At least one test particle (type 0) has finite mass. This might lead to unexpected behaviour. Set testparticle_hidewarnings=1 to hide this warning.");
    }
    if (r->status != REB_STATUS_PAUSED && r->status != REB_STATUS_SCREENSHOT){ // Allow simulation to be paused initially
        r->status = REB_STATUS_RUNNING;
    }
    run_heartbeat(r);
#ifdef __EMSCRIPTEN__
    double t0 = emscripten_performance_now();
#endif
    struct reb_walltime_tracker walltime_tracker;
    reb_walltime_tracker_start(&walltime_tracker);
    while(reb_check_exit(r,thread_info->tmax,&last_full_dt)<0){
#ifdef __EMSCRIPTEN__
        double t1 = emscripten_performance_now();
        if (t1-t0>1000./120.){ // max framerate 120Hz
            t0 = t1;
            emscripten_sleep(0); // allow drawing and event handling
        }

#endif 
#ifdef OPENGL
        if (r->display_data){
            // Note: Mutex is not FIFO.
            // Allow time for mutex to lock in display.c before it is relocked here.
            while (r->display_data->need_copy == 1){
                usleep(10);
            }
            pthread_mutex_lock(&(r->display_data->mutex)); 
        }
#endif //OPENGL
#ifdef SERVER
        if (r->server_data){
            // Note: Mutex is not FIFO.
            // Allow time for mutex to lock in display.c before it is relocked here.
            while (r->server_data->need_copy == 1){
                usleep(10);
            }
#ifdef _WIN32
            WaitForSingleObject(r->server_data->mutex, INFINITE);
#else // _WIN32
            pthread_mutex_lock(&(r->server_data->mutex)); 
#endif // _WIN32
            r->server_data->mutex_locked_by_integrate = 1;
        }
#endif //SERVER
        if (r->simulationarchive_filename){ reb_simulationarchive_heartbeat(r);}
        reb_simulation_step(r); 
        reb_walltime_tracker_step(r, &walltime_tracker);
        run_heartbeat(r);
        if (reb_sigint){
            r->status = REB_STATUS_SIGINT;
        }
#ifdef OPENGL
        if (r->display_data){
            pthread_mutex_unlock(&(r->display_data->mutex));
        }
#endif //OPENGL
#ifdef SERVER
        if (r->server_data){
#ifdef _WIN32
            ReleaseMutex(r->server_data->mutex);
#else // _WIN32
            pthread_mutex_unlock(&(r->server_data->mutex));
#endif // _WIN32
            r->server_data->mutex_locked_by_integrate = 0;
        }
#endif //SERVER
        if (r->usleep > 0){
            usleep(r->usleep);
        }
    }
    reb_walltime_tracker_finish(r, &walltime_tracker);
    reb_simulation_synchronize(r);
    if(r->exact_finish_time==1){ // if finish_time = 1, r->dt could have been shrunk, so set to the last full timestep
        r->dt = last_full_dt; 
    }
    if (r->simulationarchive_filename){ reb_simulationarchive_heartbeat(r);}
    reb_integration_end();
    return NULL;
}


enum REB_STATUS reb_simulation_integrate(struct reb_simulation* const r, double tmax){
    struct reb_thread_info thread_info = {
        .r = r,
        .tmax = tmax, 
    };

#ifdef OPENGL
#ifdef __EMSCRIPTEN__
    if (r->display_data==NULL){
        r->display_data = calloc(sizeof(struct reb_display_data),1);
        r->display_data->r = r;
        // Setup windows, compile shaders, etc.
        reb_display_init(r); // Will return. Display routines running in animation_loop.
    }
    reb_simulation_integrate_raw(&thread_info);
#else // __EMSCRIPTEN__
    if (r->display_data==NULL){
        r->display_data = calloc(sizeof(struct reb_display_data),1);
        r->display_data->r = r;
        if (pthread_mutex_init(&(r->display_data->mutex), NULL)){
            reb_simulation_error(r,"Mutex creation failed.");
        }
    }

    if (pthread_create(&r->display_data->compute_thread,NULL,reb_simulation_integrate_raw,&thread_info)){
        reb_simulation_error(r, "Error creating compute thread.");
    }

    reb_display_init(r); // Display routines need to run on main thread. Will not return until r->status<0.
    if (pthread_join(r->display_data->compute_thread,NULL)){
        reb_simulation_error(r, "Error joining compute thread.");
    }
#endif // __EMSCRIPTEN__
#else // OPENGL
    reb_simulation_integrate_raw(&thread_info);
#endif // OPENGL
    return r->status;
}

// Integrate (or step) several independent simulations with a pool of threads. Each simulation is
// integrated by one thread from start to end, exactly as reb_simulation_integrate would do it, so
// the results are bitwise identical to integrating them one after another.
struct reb_simulations_batch {
    struct reb_simulation** sims;
    size_t N;
    const double* tmax;     // NULL: step mode
    size_t N_steps;
    volatile long next;     // next simulation to start
};

// Two WHFast512 simulations advanced by N_steps REBOUND steps together: the integrator steps use
// the paired kernel (reb_integrator_whfast512_step_pair), everything else is done per simulation
// exactly as in reb_simulation_steps / reb_simulation_step, so each result is bitwise identical
// to reb_simulation_steps(r, N_steps). Only for simulations without per-step extras.
static int reb_simulations_pairable(struct reb_simulation* const r1, struct reb_simulation* const r2){
    struct reb_simulation* rs[2] = {r1, r2};
    for (int i=0; i<2; i++){
        const struct reb_simulation* r = rs[i];
        if (strcmp(r->integrator.name, "whfast512")) return 0;
        if (r->pre_timestep_modifications || r->post_timestep_modifications || r->N_odes || r->N_var
                || r->collision != REB_COLLISION_NONE || r->simulationarchive_filename) return 0;
    }
    return reb_integrator_whfast512_can_pair(r1, r2);
}
static void reb_simulations_steps_pair(struct reb_simulation* const r1, struct reb_simulation* const r2, const size_t N_steps){
    struct reb_simulation* rs[2] = {r1, r2};
    struct reb_walltime_tracker wt[2];
    for (int j=0; j<2; j++){
        rs[j]->status = REB_STATUS_RUNNING;
        run_heartbeat(rs[j]);
        reb_walltime_tracker_start(&wt[j]);
    }
    size_t i = 0;
    for (; i<N_steps && r1->status<0 && r2->status<0; i++){
        reb_integrator_whfast512_step_pair(r1, r2);
        for (int j=0; j<2; j++){
            struct reb_simulation* const r = rs[j];
            // the remainder of reb_simulation_step (no timestep modifications, ODEs, variational
            // particles or collisions: see reb_simulations_pairable)
            r->did_modify_particles = 0;
            reb_boundary_check(r);
            r->steps_done++;
            reb_walltime_tracker_step(r, &wt[j]);
            run_heartbeat(r);
        }
    }
    for (int j=0; j<2; j++){
        struct reb_simulation* const r = rs[j];
        for (size_t k=i; k<N_steps && r->status<0; k++){ // the other one stopped: continue alone
            reb_simulation_step(r);
            reb_walltime_tracker_step(r, &wt[j]);
            run_heartbeat(r);
        }
        reb_walltime_tracker_finish(r, &wt[j]);
        reb_simulation_synchronize(r);
        if (r->status <= 0) r->status = REB_STATUS_SUCCESS;
    }
}

static void* reb_simulations_batch_worker(void* args){
    struct reb_simulations_batch* b = args;
    while (1){
#if defined(_MSC_VER) && !defined(__clang__)
        const long k = _InterlockedExchangeAdd(&b->next, 1);
#else
        const long k = __atomic_fetch_add(&b->next, 1, __ATOMIC_SEQ_CST);
#endif
        if (k < 0 || (size_t)k >= b->N) break;
        struct reb_simulation* const r = b->sims[k];
        if (!r) continue;
        if (b->tmax){
            struct reb_thread_info thread_info = { .r = r, .tmax = b->tmax[k] };
            reb_simulation_integrate_raw(&thread_info);
        }else{
            // Step mode: pair this simulation with the next unclaimed one if both are WHFast512
            // simulations that can share the paired kernel (same results, ~1.2-1.5x faster).
            struct reb_simulation* r2 = NULL;
            if (!strcmp(r->integrator.name, "whfast512")){
#if defined(_MSC_VER) && !defined(__clang__)
                const long k2 = _InterlockedExchangeAdd(&b->next, 1);
#else
                const long k2 = __atomic_fetch_add(&b->next, 1, __ATOMIC_SEQ_CST);
#endif
                if (k2 >= 0 && (size_t)k2 < b->N) r2 = b->sims[k2];
            }
            if (r2 && reb_simulations_pairable(r, r2)){
                reb_simulations_steps_pair(r, r2, b->N_steps);
            }else{
                reb_simulation_steps(r, b->N_steps);
                if (r2) reb_simulation_steps(r2, b->N_steps);
            }
        }
    }
    return NULL;
}

#if defined(_WIN32)
#include <windows.h>
static void* reb_simulations_batch_worker_omp1(void* b);
static DWORD WINAPI reb_simulations_batch_worker_win(LPVOID args){ reb_simulations_batch_worker_omp1(args); return 0; }
#elif !defined(__EMSCRIPTEN__)
#include <pthread.h>
#include <unistd.h>
#endif

// Ensemble workers run their simulations with one OpenMP thread each: the cores are already used
// by the workers, and nested teams would oversubscribe them (OpenMP builds; results do not depend
// on the thread count).
static void* reb_simulations_batch_worker_omp1(void* b){
#ifdef OPENMP
    omp_set_num_threads(1); // affects only this (worker) thread
#endif // OPENMP
    return reb_simulations_batch_worker(b);
}

static int reb_simulations_default_threads(void){
#if defined(_WIN32)
    SYSTEM_INFO si;
    GetSystemInfo(&si);
    return (int)si.dwNumberOfProcessors;
#elif defined(__EMSCRIPTEN__)
    return 1;
#else
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    return n > 0 ? (int)n : 1;
#endif
}

static void reb_simulations_run_batch(struct reb_simulations_batch* b, int N_threads){
    if (N_threads <= 0) N_threads = reb_simulations_default_threads();
    if ((size_t)N_threads > b->N) N_threads = (int)b->N;
    if (N_threads <= 1){
        reb_simulations_batch_worker(b);
        return;
    }
#ifdef OPENMP
    // OpenMP builds: the workers are the threads of an OpenMP team. They then follow the user's
    // thread placement (OMP_PROC_BIND / OMP_PLACES); plain threads would inherit the binding of
    // the calling thread, i.e. all run on one core when binding is enabled. REBOUND's parallel
    // regions inside a worker are nested and run single-threaded.
    if (!omp_in_parallel()){
#pragma omp parallel num_threads(N_threads)
        reb_simulations_batch_worker(b);
        return;
    }
#endif // OPENMP
#if defined(_WIN32)
    HANDLE* th = malloc(sizeof(HANDLE)*(size_t)N_threads);
    int started = 0;
    for (int t=0; t<N_threads; t++){
        th[t] = CreateThread(NULL, 0, reb_simulations_batch_worker_win, b, 0, NULL);
        if (th[t]) started++;
    }
#ifdef OPENMP
    const int omp_threads_caller = omp_get_max_threads();
    omp_set_num_threads(1);
#endif // OPENMP
    reb_simulations_batch_worker(b); // the calling thread works too (also covers failed starts)
#ifdef OPENMP
    omp_set_num_threads(omp_threads_caller);
#endif // OPENMP
    for (int t=0; t<N_threads; t++){
        if (th[t]){ WaitForSingleObject(th[t], INFINITE); CloseHandle(th[t]); }
    }
    free(th);
    (void)started;
#elif defined(__EMSCRIPTEN__)
    reb_simulations_batch_worker(b);
#else
    pthread_t* th = malloc(sizeof(pthread_t)*(size_t)N_threads);
    int* ok = calloc((size_t)N_threads, sizeof(int));
    for (int t=1; t<N_threads; t++){
        ok[t] = (pthread_create(&th[t], NULL, reb_simulations_batch_worker_omp1, b) == 0);
    }
#ifdef OPENMP
    const int omp_threads_caller = omp_get_max_threads();
    omp_set_num_threads(1);
#endif // OPENMP
    reb_simulations_batch_worker(b); // the calling thread is worker 0 (also covers failed starts)
#ifdef OPENMP
    omp_set_num_threads(omp_threads_caller); // restore the caller's setting
#endif // OPENMP
    for (int t=1; t<N_threads; t++){
        if (ok[t]) pthread_join(th[t], NULL);
    }
    free(ok);
    free(th);
#endif
}

void reb_simulations_integrate(struct reb_simulation** sims, size_t N, const double* tmax, int N_threads){
    struct reb_simulations_batch b = { .sims = sims, .N = N, .tmax = tmax, .N_steps = 0, .next = 0 };
    reb_simulations_run_batch(&b, N_threads);
}

void reb_simulations_steps(struct reb_simulation** sims, size_t N, size_t N_steps, int N_threads){
    struct reb_simulations_batch b = { .sims = sims, .N = N, .tmax = NULL, .N_steps = N_steps, .next = 0 };
    reb_simulations_run_batch(&b, N_threads);
}

enum REB_STATUS reb_simulation_steps(struct reb_simulation* const r, size_t N_steps){
    r->status = REB_STATUS_RUNNING; 
    run_heartbeat(r);
    struct reb_walltime_tracker walltime_tracker;
    reb_walltime_tracker_start(&walltime_tracker);
    for (size_t i=0;i<N_steps && r->status<0;i++){
        if (r->simulationarchive_filename){ reb_simulationarchive_heartbeat(r);}
        reb_simulation_step(r);
        reb_walltime_tracker_step(r, &walltime_tracker);
        run_heartbeat(r);
    }
    reb_walltime_tracker_finish(r, &walltime_tracker);
    reb_simulation_synchronize(r);
    if (r->simulationarchive_filename){ reb_simulationarchive_heartbeat(r);}
    if (r->status <= 0) r->status = REB_STATUS_SUCCESS; // No error occurred. Success.
    return r->status;
}
static void reb_simulation_step(struct reb_simulation* const r){
    reb_omp_select_cores(r); // one OpenMP team size for the whole step

    // A 'DKD'-like integrator will do the first 'D' part.
    if (r->pre_timestep_modifications){
        reb_simulation_synchronize(r);
        r->pre_timestep_modifications(r);
    }

    PROFILING_START();

    if (r->integrator.callbacks.step){
        r->integrator.callbacks.step(r, r->integrator.state);
    }

    // Integrate other ODEs
    if (r->N_odes && strcmp(r->integrator.name, "bs")!=0){
        double dt = r->dt_last_done;
        double t = r->t - r->dt_last_done; // Note: floating point inaccuracy
        double forward = (dt>0.) ? 1. : -1.;
        struct reb_integrator_bs_state* bs = reb_integrator_bs.create();
        while(t*forward < r->t*forward && fabs((r->t - t)/(fabs(r->t)+1e-16))>1e-15){
            if (reb_sigint > 1){
                r->status = REB_STATUS_SIGINT;
                return;
            }
            if (bs->dt_proposed !=0.){
                double max_dt = fabs(r->t - t);
                dt = fabs(bs->dt_proposed);
                if (dt > max_dt){ // Don't overshoot N-body timestep
                    dt = max_dt;
                    bs->first_or_last_step = 1;
                }
                dt *= forward;
            }
            int success = reb_integrator_bs_step_odes(r, bs, dt);
            if (success){
                t += dt;
            }
        }
        reb_integrator_bs.free(bs);
    }

    PROFILING_STOP(PROFILING_CAT_INTEGRATOR);

    if (r->post_timestep_modifications){
        reb_simulation_synchronize(r);
        r->post_timestep_modifications(r);
    }

    // Reset tainted particle flag
    r->did_modify_particles = 0;


    if (r->N_var){
        reb_simulation_rescale_var(r);
    }

    PROFILING_START();
    reb_boundary_check(r);     
    PROFILING_STOP(PROFILING_CAT_BOUNDARY);

    if (r->collision != REB_COLLISION_NONE){
#ifdef MPI
        reb_communication_mpi_distribute_particles(r);
#endif // MPI

        PROFILING_START();
        reb_collision_search(r);
        PROFILING_STOP(PROFILING_CAT_COLLISION);
    }


    // Update step counter
    r->steps_done++; // This also counts failed IAS15 steps
}

void reb_simulation_copy_with_messages(struct reb_simulation* r_copy,  struct reb_simulation* r, enum REB_BINARYDATA_ERROR_CODE* warnings){
    // Expect r_copy to be a newly created and initialized simulation.
    char* bufp;
    size_t sizep;
    reb_binarydata_simulation_to_stream(r, &bufp,&sizep);

    FILE* fin = reb_fmemopen(bufp, sizep, "r");
    reb_binarydata_input_fields(r_copy, fin, warnings);
    fclose(fin);

    free(bufp);
}

char* reb_simulation_diff_char(struct reb_simulation* r1, struct reb_simulation* r2){
    char* bufp1;
    char* bufp2;
    char* bufp;
    size_t sizep1, sizep2, size;
    reb_binarydata_simulation_to_stream(r1, &bufp1,&sizep1);
    reb_binarydata_simulation_to_stream(r2, &bufp2,&sizep2);

    reb_binarydata_diff(bufp1, sizep1, bufp2, sizep2, &bufp, &size, REB_BINARYDATA_OUTPUT_BUFFER);

    free(bufp1);
    free(bufp2);
    return bufp;
}

void reb_simulation_synchronize(struct reb_simulation* r){
    if (r->integrator.callbacks.synchronize){
        r->integrator.callbacks.synchronize(r, r->integrator.state);
    }
}


void reb_simulation_update_acceleration(struct reb_simulation* r){
    // The OpenMP core policy is chosen once per step (reb_simulation_step), not here: IAS15 calls
    // this function ~20 times per step, and the policy only affects speed, never results.
    if (r->gravity==REB_GRAVITY_CUSTOM){
        if (!r->gravity_custom){
            reb_simulation_error(r, "REB_GRAVITY_CUSTOM selected, but r->gravity_custom function pointer not provided.");
        }
        // Let the user provided routine do everything, then return.
        r->gravity_custom(r);
        return;
    }
    switch (r->gravity){
        case REB_GRAVITY_NONE: // Do nothing.
            for (size_t j=0; j<r->N; j++){
                r->particles[j].ax = 0;  
                r->particles[j].ay = 0;  
                r->particles[j].az = 0;  
            }  
            break;
        case REB_GRAVITY_TREE:
            reb_gravity_tree_calculate_acceleration(r);
            break;
        case REB_GRAVITY_JACOBI:
            reb_gravity_jacobi_calculate_acceleration(r);
            break;
        case REB_GRAVITY_BASIC:
            reb_gravity_basic_calculate_acceleration(r);
            break;
        case REB_GRAVITY_COMPENSATED:
            reb_gravity_compensated_calculate_acceleration(r);
            break;
        default:
            reb_simulation_error(r, "Gravity module not found.");
            return;
    }
    if (r->N_var){
        switch (r->gravity){
            case REB_GRAVITY_BASIC:
                reb_gravity_basic_calculate_acceleration_var(r);
                break;
            case REB_GRAVITY_NONE:
                break;
            default:
                reb_simulation_error(r, "Variational gravity calculation not implemented in selected gravity module. Please use REB_GRAVITY_BASIC.");
                return;
        }
    }

    if (r->additional_forces){
        r->additional_forces(r);
    }
}

int reb_simulation_diff(struct reb_simulation* r1, struct reb_simulation* r2, int output_option){
    if (output_option!=REB_BINARYDATA_OUTPUT_PRINT && output_option!=REB_BINARYDATA_OUTPUT_NONE){
        reb_simulation_error(r1, "Invalid output_option in reb_simulation_diff().");
        return -1;
    }
    char* bufp1;
    char* bufp2;
    size_t sizep1, sizep2;
    reb_binarydata_simulation_to_stream(r1, &bufp1,&sizep1);
    reb_binarydata_simulation_to_stream(r2, &bufp2,&sizep2);

    int ret = reb_binarydata_diff(bufp1, sizep1, bufp2, sizep2, NULL, NULL, output_option);

    free(bufp1);
    free(bufp2);
    return ret;
}

struct reb_simulation* reb_simulation_copy(struct reb_simulation* r){
    struct reb_simulation* r_copy = reb_simulation_create();
    enum REB_BINARYDATA_ERROR_CODE warnings = REB_BINARYDATA_WARNING_NONE;
    reb_simulation_copy_with_messages(r_copy,r,&warnings);
    r = reb_binarydata_process_warnings(r, warnings);
    return r_copy;
}

// Finds the two largest particles in the simulation. *p1 and *p2 will be set to the indices of the particles.
void reb_simulation_two_largest_particles(struct reb_simulation* r, size_t* p1, size_t* p2) {
    struct reb_particle* particles = r->particles;
    *p1 = SIZE_MAX;
    *p2 = SIZE_MAX;
    double largest1 = -1.0;
    double largest2 = -1.0;
    // Serial in all builds: a linear scan is cheap next to the collision search that calls it,
    // and it keeps the tie-breaking (lowest index first) independent of the thread count.
    for (size_t i=0; i<r->N; i++) {
        if (particles[i].r > largest1) {
            largest2 = largest1;
            *p2 = *p1;
            largest1 = particles[i].r;
            *p1 = i;
        }else{
            if (particles[i].r > largest2) {
                largest2 = particles[i].r;
                *p2 = i;
            }
        }
    }
}

void reb_simulation_get_serialized_particle_data(struct reb_simulation* r, double* m, double* radius, double (*xyz)[3], double (*vxvyvz)[3], double (*xyzvxvyvz)[6]){
    const size_t N = r->N;
    struct reb_particle* restrict const particles = r->particles;
    for (size_t i=0;i<N;i++){
        if (m){
            m[i] = particles[i].m;
        }
        if (radius){
            radius[i] = particles[i].r;
        }
        if (xyz){
            xyz[i][0] = particles[i].x;
            xyz[i][1] = particles[i].y;
            xyz[i][2] = particles[i].z;
        }
        if (vxvyvz){
            vxvyvz[i][0] = particles[i].vx;
            vxvyvz[i][1] = particles[i].vy;
            vxvyvz[i][2] = particles[i].vz;
        }
        if (xyzvxvyvz){
            xyzvxvyvz[i][0] = particles[i].x;
            xyzvxvyvz[i][1] = particles[i].y;
            xyzvxvyvz[i][2] = particles[i].z;
            xyzvxvyvz[i][3] = particles[i].vx;
            xyzvxvyvz[i][4] = particles[i].vy;
            xyzvxvyvz[i][5] = particles[i].vz;
        }
    }
}

void reb_simulation_set_serialized_particle_data(struct reb_simulation* r, double* m, double* radius, double (*xyz)[3], double (*vxvyvz)[3], double (*xyzvxvyvz)[6]){
    const size_t N = r->N;
    struct reb_particle* restrict const particles = r->particles;
    for (size_t i=0;i<N;i++){
        if (m){
            particles[i].m = m[i];
        }
        if (radius){
            particles[i].r = radius[i] ;
        }
        if (xyz){
            particles[i].x = xyz[i][0];
            particles[i].y = xyz[i][1];
            particles[i].z = xyz[i][2];
        }
        if (vxvyvz){
            particles[i].vx = vxvyvz[i][0];
            particles[i].vy = vxvyvz[i][1];
            particles[i].vz = vxvyvz[i][2];
        }
        if (xyzvxvyvz){
            particles[i].x = xyzvxvyvz[i][0];
            particles[i].y = xyzvxvyvz[i][1];
            particles[i].z = xyzvxvyvz[i][2];
            particles[i].vx = xyzvxvyvz[i][3];
            particles[i].vy = xyzvxvyvz[i][4];
            particles[i].vz = xyzvxvyvz[i][5];
        }
    }
}

void reb_simulation_add_display_settings(struct reb_simulation*r){
    if (r->display_settings){
        reb_simulation_error(r,"Simulation already has display settings.");
        return;
    }
    r->display_settings = calloc(1,sizeof(struct reb_display_settings));
    reb_display_settings_init(r, r->display_settings);
}

