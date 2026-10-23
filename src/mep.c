/**
 * @file    mep.c
 * @brief   Mean-element parareal (MEP): parallel-in-time integration of one WHFast512 planetary system on many cores.
 * @details The time axis is cut into windows of N slices; one thread per slice runs the fine WHFast512 integration of
 *          its slice while a cheap coarse propagator in mean elements (mep_coarse.c) links the slices (pipelined
 *          parareal, adaptive number of iterations, energy projection). The converged result is the serial WHFast512
 *          trajectory to within the tolerance. Designed and validated for the Sun and eight planets of the Solar
 *          System (four inner planets treated by perturbation theory, four giants integrated directly).
 *          The run state (current Jacobi elements, time, statistics) can be saved to a small checkpoint file and
 *          resumed later; a resumed run is bitwise identical to an uninterrupted one.
 *          Derived from the research code of the MEP paper. Linux only.
 *
 * @section LICENSE
 * Copyright (c) 2026 Michael Hippke
 *
 * This file is part of Parachron (https://github.com/hippke/parachron), a modified version of rebound.
 *
 * Parachron is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Parachron is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Parachron.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#if defined(__linux__) && !defined(__EMSCRIPTEN__)
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <semaphore.h>
#include <sched.h>
#include <time.h>
#include <unistd.h>
#include <stdatomic.h>
#include "rebound.h"
#include "rebound_internal.h"
#include "integrator_whfast512.h"
#include "mep.h"

#define NMAX 128                // maximum number of slices (threads) per window
#define KMAX REB_MEP_KHIST      // iterations are counted in khist[0..KMAX-1]
#define DAY (2.*M_PI/365.25)
typedef double el[MEP_NP][6];
#define CEL(x) ((const double (*)[6])(x))     // non-const -> const 2-D array (ISO C before C2X needs the cast)

static double now(void){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec; }
static double wrapd(double x){ return x - 2.*M_PI*floor((x+M_PI)/(2.*M_PI)); }
static void set_error(char* buf, const char* msg){ snprintf(buf, REB_MEP_MSG, "%s", msg); }
static int el_finite(const el E){ for (int i=0;i<MEP_NP;i++) for (int e=0;e<6;e++) if (!isfinite(E[i][e])) return 0; return 1; }

// ------------------------------------------------------------------------------------------------ run state (one run at a time)
static pthread_mutex_t run_lock = PTHREAD_MUTEX_INITIALIZER;
static struct reb_mep* M_;                  // the run being integrated
static int N, K, Kmax, nchunk=2;
static double tol;
static int64_t nst;                            // fine steps per slice
static struct mep_cfg cfg;
static int cpu_of_slice[NMAX];              // logical CPU of slice thread n (-1: not pinned)
static el U[2][KMAX+1][NMAX+1];             // boundary states per window parity, iteration, boundary
static sem_t ready[2][KMAX+1][NMAX+1], win_start[2], win_done[2];
static atomic_int stop_at[2];               // iteration at which the window stopped (-1: running)
static double upd_max[2][KMAX+1];
static pthread_mutex_t upd_mx=PTHREAD_MUTEX_INITIALIZER;
static el Ecur;                             // start state of the next window
static atomic_long W_eff;                   // windows of this call (lowered when the run stops early)
static int64_t W_done;                         // windows accepted in this call
static double t_start_wall, walltime_limit;
static enum REB_STATUS stop_status;
static atomic_int fine_interrupted;

// ------------------------------------------------------------------------------------------------ helpers
static double energy_sim(struct reb_simulation* s){
    double E=reb_simulation_energy(s);
    if (M_->gr_potential){
        const double c=10065.32, m0=s->particles[0].m;
        for (size_t i=1;i<s->N;i++){
            const double dx=s->particles[i].x-s->particles[0].x, dy=s->particles[i].y-s->particles[0].y, dz=s->particles[i].z-s->particles[0].z;
            E += -3.*m0*m0*s->particles[i].m/(c*c*(dx*dx+dy*dy+dz*dz));
        }
    }
    return E;
}
static double energy_el(const el E){ struct reb_simulation* s=reb_mep_el_to_sim(M_->masses,E); double e=energy_sim(s); reb_simulation_free(s); return e; }
// Projection of a corrected boundary state onto the target energy Et and angular-momentum vector Lt (targets: the
// start state's values plus the increments of one fine solve, so each slice carries the energy and angular-momentum
// changes of the fine integrator, and the parareal fixed point is unchanged). The barycentric velocities change by
// dv_i = alpha v_i + beta x r_i, the minimal mass-weighted change that reaches both targets; linearised,
//   dE = 2K alpha + L.beta,   dL = alpha L + I beta,   I = sum m (r^2 1 - r r^T);
// two Newton passes, applied as increments (keeps the full relative precision of the tiny corrections). The total
// momentum stays zero. Without the angular-momentum part, |L| drifted by -2.6e-10 per Gyr (16 Solar System runs).
static void L_sim(struct reb_simulation* s, double L[3]){ struct reb_vec3d l=reb_simulation_angular_momentum(s); L[0]=l.x; L[1]=l.y; L[2]=l.z; }
static void L_el(const el E, double L[3]){ struct reb_simulation* s=reb_mep_el_to_sim(M_->masses,E); L_sim(s,L); reb_simulation_free(s); }
static int solve4(double A[4][4], double b[4]){
    for (int c=0;c<4;c++){
        int p=c; for (int r=c+1;r<4;r++) if (fabs(A[r][c])>fabs(A[p][c])) p=r;
        if (A[p][c]==0.) return 1;
        if (p!=c){ for (int q=0;q<4;q++){ const double t=A[c][q]; A[c][q]=A[p][q]; A[p][q]=t; } const double t=b[c]; b[c]=b[p]; b[p]=t; }
        for (int r=c+1;r<4;r++){ const double f=A[r][c]/A[c][c]; for (int q=c;q<4;q++) A[r][q]-=f*A[c][q]; b[r]-=f*b[c]; }
    }
    for (int c=3;c>=0;c--){ double v=b[c]; for (int q=c+1;q<4;q++) v-=A[c][q]*b[q]; b[c]=v/A[c][c]; }
    return 0;
}
static void project_EL(el E, double Et, const double Lt[3]){
    struct reb_simulation* s=reb_mep_el_to_sim(M_->masses,CEL(E));
    for (int it=0; it<2; it++){
        const double e=energy_sim(s); double L[3]; L_sim(s,L);
        double K2=0, I[3][3]={{0}};
        for (size_t i=0;i<s->N;i++){
            const struct reb_particle p=s->particles[i]; const double r[3]={p.x,p.y,p.z};
            K2+=p.m*(p.vx*p.vx+p.vy*p.vy+p.vz*p.vz);
            const double r2=r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
            for (int a=0;a<3;a++) for (int c=0;c<3;c++) I[a][c]+=p.m*((a==c?r2:0.)-r[a]*r[c]);
        }
        double A[4][4]={{K2,L[0],L[1],L[2]},{L[0],I[0][0],I[0][1],I[0][2]},{L[1],I[1][0],I[1][1],I[1][2]},{L[2],I[2][0],I[2][1],I[2][2]}};
        double x[4]={Et-e, Lt[0]-L[0], Lt[1]-L[1], Lt[2]-L[2]};
        if (solve4(A,x)) break;
        for (size_t i=0;i<s->N;i++){
            struct reb_particle* p=&s->particles[i];
            const double bx=x[2]*p->z-x[3]*p->y, by=x[3]*p->x-x[1]*p->z, bz=x[1]*p->y-x[2]*p->x;   // beta x r
            p->vx+=p->vx*x[0]+bx; p->vy+=p->vy*x[0]+by; p->vz+=p->vz*x[0]+bz;
        }
    }
    reb_mep_sim_to_el(s,E); reb_simulation_free(s);
}
static struct reb_simulation* fine_sim(const el Ein){
    struct reb_simulation* s=reb_mep_el_to_sim(M_->masses,Ein);
    s->dt=M_->dt; s->exact_finish_time=0;
    reb_simulation_set_integrator(s,"whfast512");
    struct reb_integrator_whfast512_state* w=s->integrator.state;
    w->gr_potential=M_->gr_potential; w->corrector=M_->corrector; w->jacobi_recursive=1; w->kernel=M_->kernel;
    return s;
}
// fine slice; returns 0 if abandoned because the window stopped at an iteration < kneed
static int fine(const el Ein, el Eout, int kneed, atomic_int* stp){
    struct reb_simulation* s=fine_sim(Ein);
    struct reb_integrator_whfast512_state* w=s->integrator.state;
    const int nch = (nchunk>1 && nst%nchunk==0) ? nchunk : 1;
    w->concatenate_steps=nst/nch;
    for (int q=0;q<nch;q++){
        const int st=atomic_load(stp);
        if (nch>1 && st>=0 && st<kneed){ reb_simulation_free(s); return 0; }
        reb_simulation_steps(s,1);
        if (s->status!=REB_STATUS_RUNNING && s->status!=REB_STATUS_SUCCESS) atomic_store(&fine_interrupted,1);
    }
    reb_mep_sim_to_el(s,Eout); reb_simulation_free(s);
    return 1;
}
static double poserr_el(const el A, const el B){
    struct reb_simulation* a=reb_mep_el_to_sim(M_->masses,A); struct reb_simulation* b=reb_mep_el_to_sim(M_->masses,B);
    double m=0; for (int i=1;i<=MEP_NP;i++){ const double d[3]={a->particles[i].x-b->particles[i].x,a->particles[i].y-b->particles[i].y,a->particles[i].z-b->particles[i].z}; for(int x=0;x<3;x++) if (fabs(d[x])>m) m=fabs(d[x]); }
    reb_simulation_free(a); reb_simulation_free(b); return m;
}
// prepare the state of a window (parity p) that starts from Ecur
static void window_init(int p){
    for (int k=0;k<=KMAX;k++){ memcpy(U[p][k][0],Ecur,sizeof(el)); upd_max[p][k]=0; for (int n=0;n<=N;n++){ sem_destroy(&ready[p][k][n]); sem_init(&ready[p][k][n],0,0); } }
    sem_destroy(&win_done[p]); sem_init(&win_done[p],0,0);
    atomic_store(&stop_at[p],-1);
}

// ------------------------------------------------------------------------------------------------ checkpoints
#define MEP_MAGIC "REBOUND-MEP-CHECKPOINT"
#define MEP_VERSION 1
static int ckpt_write(const struct reb_mep* m){
    if (m->checkpoint[0]==0) return 0;
    char tmp[sizeof(m->checkpoint)+8]; snprintf(tmp,sizeof tmp,"%s.tmp",m->checkpoint);
    FILE* f=fopen(tmp,"wb"); if (!f) return 1;
    const char magic[32]=MEP_MAGIC; const int32_t version=MEP_VERSION;
    int ok = fwrite(magic,1,32,f)==32 && fwrite(&version,4,1,f)==1
    && fwrite(&m->opt,sizeof m->opt,1,f)==1
    && fwrite(&m->t,8,1,f)==1 && fwrite(&m->t_start,8,1,f)==1 && fwrite(&m->dt,8,1,f)==1
    && fwrite(&m->windows,8,1,f)==1 && fwrite(m->khist,8,KMAX,f)==KMAX
    && fwrite(m->masses,8,MEP_NP+1,f)==MEP_NP+1 && fwrite(m->E,8,MEP_NP*6,f)==MEP_NP*6
    && fwrite(&m->gr_potential,4,1,f)==1 && fwrite(&m->corrector,4,1,f)==1 && fwrite(&m->kernel,4,1,f)==1;
    ok = (fclose(f)==0) && ok;
    if (!ok || rename(tmp,m->checkpoint)!=0) return 1;     // rename is atomic: a crash leaves the previous checkpoint
    return 0;
}
static int ckpt_read(struct reb_mep* m, const char* fn){
    FILE* f=fopen(fn,"rb"); if (!f) return 1;
    char magic[32]; int32_t version=0;
    int ok = fread(magic,1,32,f)==32 && strncmp(magic,MEP_MAGIC,32)==0 && fread(&version,4,1,f)==1 && version==MEP_VERSION
    && fread(&m->opt,sizeof m->opt,1,f)==1
    && fread(&m->t,8,1,f)==1 && fread(&m->t_start,8,1,f)==1 && fread(&m->dt,8,1,f)==1
    && fread(&m->windows,8,1,f)==1 && fread(m->khist,8,KMAX,f)==KMAX
    && fread(m->masses,8,MEP_NP+1,f)==MEP_NP+1 && fread(m->E,8,MEP_NP*6,f)==MEP_NP*6
    && fread(&m->gr_potential,4,1,f)==1 && fread(&m->corrector,4,1,f)==1 && fread(&m->kernel,4,1,f)==1;
    fclose(f);
    return ok ? 0 : 2;
}

// ------------------------------------------------------------------------------------------------ slice threads
static void* slice_thread(void* arg){
    const int n=(int)(long)arg;
    if (cpu_of_slice[n]>=0){ cpu_set_t cs; CPU_ZERO(&cs); CPU_SET(cpu_of_slice[n], &cs); pthread_setaffinity_np(pthread_self(), sizeof cs, &cs); }
    el Fprev, Gold;
    const int kend = Kmax;
    for (long w=0; w<atomic_load(&W_eff); w++){
        const int p=w&1; el (*Up)[NMAX+1]=U[p]; sem_t (*rdy)[NMAX+1]=ready[p]; atomic_int* stp=&stop_at[p];
        if (n==0 && w>0) sem_wait(&win_start[p]);
        double EFprev=0, EUprev=0, LFprev[3]={0}, LUprev[3]={0}; int have_F=0;
        for (int k=0;k<=kend;k++){
            { const int st=atomic_load(stp); if (st>=0 && k>st) break; }
            if (n>0) sem_wait(&rdy[k][n]);
            { const int st=atomic_load(stp); if (st>=0 && k>st) break; }
            el g;
            // G is not needed for boundaries that are already exact (n <= k-1)
            if (k==0 || n>k-1) reb_mep_coarse_par(&cfg,(const double (*)[6])Up[k][n],g,M_->opt.helpers>0);
            if (k==0 || n<k-1){
                if (k==0) memcpy(Up[0][n+1],g,sizeof(el)); else memcpy(Up[k][n+1],Up[k-1][n+1],sizeof(el));
            } else if (n==k-1){
                memcpy(Up[k][n+1],Fprev,sizeof(el));    // exact
            } else {
                if (!have_F) break;        // only if the window has already stopped (abandoned speculative fine solve)
                for (int i=0;i<MEP_NP;i++) for (int e=0;e<6;e++){
                    double d=Fprev[i][e]-Gold[i][e]; if (e==1) d=wrapd(d);
                    Up[k][n+1][i][e]=g[i][e]+d;
                }
                for (int i=0;i<MEP_NP;i++) Up[k][n+1][i][1]=wrapd(Up[k][n+1][i][1]);
                double Lt[3]; L_el(CEL(Up[k][n]),Lt); for (int a=0;a<3;a++) Lt[a]+=LFprev[a]-LUprev[a];
                project_EL(Up[k][n+1], energy_el(CEL(Up[k][n])) + (EFprev-EUprev), Lt);
            }
            memcpy(Gold,g,sizeof(el));
            if (k>0){
                const double u = n>=k-1 ? poserr_el(CEL(Up[k][n+1]),CEL(Up[k-1][n+1])) : 0.;
                pthread_mutex_lock(&upd_mx); if (u>upd_max[p][k]) upd_max[p][k]=u; pthread_mutex_unlock(&upd_mx);
            }
            if (n<N-1) sem_post(&rdy[k][n+1]);
            else {
                // last slice: all updates of iteration k are known -> decide
                const int interrupted = reb_sigint>0 || atomic_load(&fine_interrupted);
                const int converged = (k>=K && upd_max[p][k]<=tol) || k==Kmax;
                const int nonfinite = !el_finite(CEL(Up[k][N]));   // never accept a non-finite state
                if (converged || interrupted || nonfinite){
                    long Wnext = atomic_load(&W_eff);
                    if (nonfinite && !interrupted){
                        set_error(M_->message,"MEP stopped: non-finite state in a window (the last completed window is kept).");
                        stop_status = REB_STATUS_GENERIC_ERROR; Wnext = w+1;
                    } else if (interrupted){
                        // some fine or coarse integration of this window may have been cut short: discard the window
                        stop_status = REB_STATUS_SIGINT; Wnext = w+1;
                    } else {
                        struct reb_mep* m=M_;
                        memcpy(Ecur,Up[k][N],sizeof(el));
                        m->khist[k]++; m->windows++; W_done++;
                        m->t = m->t_start + (double)m->windows*N*nst*m->dt;
                        memcpy(m->E,Ecur,sizeof(el));
                        if (walltime_limit>0. && now()-t_start_wall>=walltime_limit && w+1<Wnext){ stop_status = REB_STATUS_USER; Wnext = w+1; }
                        if (m->opt.checkpoint_every>0 && m->windows%m->opt.checkpoint_every==0) ckpt_write(m);
                    }
                    atomic_store(&W_eff,Wnext);
                    // every thread has entered this window, so none uses the other parity's state any more
                    if (w+1<Wnext) window_init(p^1);
                    atomic_store(stp,k);
                    for (int q=0;q<N;q++) sem_post(&win_done[p]);
                    for (int kk=k+1; kk<=KMAX; kk++) for (int q=1;q<=N;q++) sem_post(&rdy[kk][q]);   // release waiting threads
                    if (w+1<Wnext) sem_post(&win_start[p^1]);
                }
            }
            // fine slice for the next iteration (speculative: abandoned if the window stops)
            have_F=0;
            if (k<kend && n>=k){
                if (fine((const double (*)[6])Up[k][n],Fprev,k+1,stp)){ EUprev=energy_el(CEL(Up[k][n])); EFprev=energy_el(CEL(Fprev)); L_el(CEL(Up[k][n]),LUprev); L_el(CEL(Fprev),LFprev); have_F=1; }
            }
        }
        sem_wait(&win_done[p]);    // a thread can finish k = Kmax before the window is decided
    }
    return NULL;
}

// ------------------------------------------------------------------------------------------------ CPU topology
// logical CPUs that are the first hardware thread of a physical core, within this process' affinity mask
static int physical_cores(int* list, int max){
    cpu_set_t cs; CPU_ZERO(&cs);
    if (sched_getaffinity(0,sizeof cs,&cs)!=0) return 0;
    int n=0;
    for (int q=0;q<CPU_SETSIZE && n<max;q++){
        if (!CPU_ISSET(q,&cs)) continue;
        char fn[128]; snprintf(fn,sizeof fn,"/sys/devices/system/cpu/cpu%d/topology/thread_siblings_list",q);
        FILE* f=fopen(fn,"r"); int first=q;
        if (f){ char buf[64]={0}; if (fgets(buf,sizeof buf,f)) sscanf(buf,"%d",&first); fclose(f); }
        if (first==q) list[n++]=q;
    }
    return n;
}

// ------------------------------------------------------------------------------------------------ public API
void reb_mep_options_default(struct reb_mep_options* o){
    memset(o,0,sizeof *o);
    o->cores = 0; o->helpers = -1; o->pin = 1; o->k_min = 3; o->k_max = 6; o->tol = 3e-6; o->slice_steps = 60000;
    o->checkpoint_every = 100; o->walltime = 0.;
}

struct reb_mep* reb_mep_create(struct reb_simulation* r, const struct reb_mep_options* o, const char* checkpoint){
    struct reb_mep* m=calloc(1,sizeof *m);
    if (o) m->opt=*o; else reb_mep_options_default(&m->opt);
    if (checkpoint) snprintf(m->checkpoint,sizeof m->checkpoint,"%s",checkpoint);
    if (r->N!=MEP_NP+1 || r->N_var!=0){ reb_simulation_error(r,"MEP needs a central object and exactly 8 planets (no test particles or variational particles)."); free(m); return NULL; }
    if (r->integrator.state==NULL || strcmp(r->integrator.name,"whfast512")!=0){ reb_simulation_error(r,"MEP needs the WHFast512 integrator (sim.integrator = 'whfast512')."); free(m); return NULL; }
    if (r->G!=1.){ reb_simulation_error(r,"MEP needs G = 1 (units AU, Msun, yr/2pi), as WHFast512."); free(m); return NULL; }
    if (!(r->dt>0.)){ reb_simulation_error(r,"MEP needs a positive time step."); free(m); return NULL; }
    struct reb_integrator_whfast512_state* w=r->integrator.state;
    if (w->N_systems!=1){ reb_simulation_error(r,"MEP needs N_systems = 1."); free(m); return NULL; }
    if (m->opt.cores<=0){ int cores[NMAX]; m->opt.cores=physical_cores(cores,NMAX); }
    if (m->opt.cores<2 || m->opt.cores>NMAX){ reb_simulation_error(r,"MEP needs between 2 and 128 cores (slices per window)."); free(m); return NULL; }
    m->dt=r->dt; m->t=r->t; m->t_start=r->t;
    m->gr_potential=w->gr_potential; m->corrector=w->corrector; m->kernel=w->kernel;
    for (int i=0;i<=MEP_NP;i++) m->masses[i]=r->particles[i].m;
    reb_mep_sim_to_el(r,m->E);
    // resolve the fine kernel now (kernel 0 = fastest available), so that a resumed run uses the same one
    if (m->kernel==0){
        struct reb_simulation* s=reb_simulation_copy(r);
        struct reb_integrator_whfast512_state* ws=s->integrator.state; ws->concatenate_steps=1;
        reb_simulation_steps(s,1); m->kernel=ws->kernel_used; reb_simulation_free(s);
    }
    return m;
}
struct reb_mep* reb_mep_create_from_checkpoint(const char* checkpoint, const struct reb_mep_options* o){
    struct reb_mep* m=calloc(1,sizeof *m);
    const int err=ckpt_read(m,checkpoint);
    if (err){ free(m); return NULL; }
    snprintf(m->checkpoint,sizeof m->checkpoint,"%s",checkpoint);
    if (o){  // options that do not change the result may be overridden
        m->opt.helpers=o->helpers; m->opt.pin=o->pin; m->opt.checkpoint_every=o->checkpoint_every; m->opt.walltime=o->walltime;
    }
    return m;
}
void reb_mep_free(struct reb_mep* m){ free(m); }

enum REB_STATUS reb_mep_integrate(struct reb_mep* m, double tmax){
    if (pthread_mutex_trylock(&run_lock)!=0){ set_error(m->message,"Another MEP integration is running in this process."); return REB_STATUS_GENERIC_ERROR; }
    m->message[0]=0;
    int cores[NMAX]; const int ncore=physical_cores(cores,NMAX);
    N = m->opt.cores;
    if (N<2 || N>NMAX){ set_error(m->message,"MEP needs between 2 and 128 cores (slices per window)."); pthread_mutex_unlock(&run_lock); return REB_STATUS_GENERIC_ERROR; }
    if (m->opt.slice_steps<2 || m->opt.slice_steps%2){ set_error(m->message,"slice_steps must be an even number >= 2."); pthread_mutex_unlock(&run_lock); return REB_STATUS_GENERIC_ERROR; }
    if (m->opt.k_min<1 || m->opt.k_max<m->opt.k_min || m->opt.k_max>KMAX-1){ set_error(m->message,"Need 1 <= k_min <= k_max <= 11."); pthread_mutex_unlock(&run_lock); return REB_STATUS_GENERIC_ERROR; }
    const double Tw = (double)N*m->opt.slice_steps*m->dt;
    const long W = (long)floor((tmax-m->t)/Tw*(1.+1e-12)+1e-9);
    if (W<1){ pthread_mutex_unlock(&run_lock); return REB_STATUS_SUCCESS; }   // less than one window left: nothing to do

    M_=m; K=m->opt.k_min; Kmax=m->opt.k_max; tol=m->opt.tol; nst=m->opt.slice_steps;
    reb_mep_default_cfg(&cfg, m->masses, nst*m->dt/DAY, m->gr_potential);
    const int helpers = m->opt.helpers>=0 ? m->opt.helpers : (N/2<8 ? N/2 : 8);
    m->opt.helpers=helpers;
    const int pin = m->opt.pin && N<=ncore;
    for (int n=0;n<N;n++) cpu_of_slice[n] = pin ? cores[n] : -1;
    if (helpers>0) reb_mep_set_threads(helpers, pin);
    memcpy(Ecur,m->E,sizeof(el));
    atomic_store(&W_eff,W); W_done=0; stop_status=REB_STATUS_SUCCESS; atomic_store(&fine_interrupted,0);
    walltime_limit=m->opt.walltime; t_start_wall=now();
    reb_integration_begin();
    for (int p=0;p<2;p++){ sem_init(&win_start[p],0,0); sem_init(&win_done[p],0,0); for (int k=0;k<=KMAX;k++) for (int n=0;n<=N;n++) sem_init(&ready[p][k][n],0,0); }
    window_init(0);
    pthread_t th[NMAX];
    for (int n=0;n<N;n++) pthread_create(&th[n],NULL,slice_thread,(void*)(long)n);
    for (int n=0;n<N;n++) pthread_join(th[n],NULL);
    for (int p=0;p<2;p++){ sem_destroy(&win_start[p]); sem_destroy(&win_done[p]); for (int k=0;k<=KMAX;k++) for (int n=0;n<=N;n++) sem_destroy(&ready[p][k][n]); }
    reb_integration_end();
    m->walltime_last = now()-t_start_wall;
    if (ckpt_write(m)) set_error(m->message,"Could not write the checkpoint file.");
    M_=NULL;
    pthread_mutex_unlock(&run_lock);
    return stop_status;
}

void reb_mep_update_simulation(const struct reb_mep* m, struct reb_simulation* r){
    struct reb_simulation* s=reb_mep_el_to_sim(m->masses,(const double (*)[6])m->E);
    reb_simulation_remove_all_particles(r);
    for (size_t i=0;i<s->N;i++) reb_simulation_add(r,s->particles[i]);
    reb_simulation_free(s);
    r->t=m->t; r->dt=m->dt;
    if (r->integrator.state==NULL || strcmp(r->integrator.name,"whfast512")!=0) reb_simulation_set_integrator(r,"whfast512");
    struct reb_integrator_whfast512_state* w=r->integrator.state;
    w->gr_potential=m->gr_potential; w->corrector=m->corrector; w->kernel=m->kernel;
}

#else // not Linux (or WebAssembly): MEP is not available
#include <stdlib.h>
#include <string.h>
#include "rebound.h"
void reb_mep_options_default(struct reb_mep_options* o){ memset(o,0,sizeof *o); }
struct reb_mep* reb_mep_create(struct reb_simulation* r, const struct reb_mep_options* o, const char* checkpoint){
    (void)o; (void)checkpoint; reb_simulation_error(r,"MEP (mean-element parareal) is only available on Linux."); return NULL; }
    struct reb_mep* reb_mep_create_from_checkpoint(const char* checkpoint, const struct reb_mep_options* o){ (void)checkpoint; (void)o; return NULL; }
    void reb_mep_free(struct reb_mep* m){ free(m); }
    enum REB_STATUS reb_mep_integrate(struct reb_mep* m, double tmax){ (void)m; (void)tmax; return REB_STATUS_GENERIC_ERROR; }
    void reb_mep_update_simulation(const struct reb_mep* m, struct reb_simulation* r){ (void)m; (void)r; }
#endif
