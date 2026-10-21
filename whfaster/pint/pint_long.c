// Long parareal integration in consecutive windows of N slices (no reference needed), plus an optional serial
// fine run of the same span, both writing the state at every slice boundary for the accuracy contract
// (energy incl. the GR potential, Jacobi elements).
//   usage: pint_long out_prefix windows N slice_days threads k [serial=1] [extra_k_tol_AU=0]
// k: fixed number of iterations per window; if extra_k_tol > 0, iterate further while the last update
// (max position change at any boundary) exceeds extra_k_tol (max k+3).
// Output files: <prefix>_pr.bin, <prefix>_ser.bin: records of (t, E_tot, E[8][6]) doubles, and <prefix>_pr.log.
#include "pint_coarse.h"
#include "integrator_whfast512.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <time.h>
#include "../tools/ss_ics.h"
#define NMAX 256
static double now(){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec; }
typedef double el[PINT_NP][6];
static double slice_days;
static double energy_sim(struct reb_simulation* s){
    // REBOUND energy + WHFast512 GR potential -3 m0^2 m_i / (c^2 r_h^2) (heliocentric r_h)
    double E=reb_simulation_energy(s);
    const double c=10065.32, m0=s->particles[0].m;
    for (unsigned i=1;i<s->N;i++){
        const double dx=s->particles[i].x-s->particles[0].x, dy=s->particles[i].y-s->particles[0].y, dz=s->particles[i].z-s->particles[0].z;
        E += -3.*m0*m0*s->particles[i].m/(c*c*(dx*dx+dy*dy+dz*dz));
    }
    return E;
}
static double energy_el(const el E){ struct reb_simulation* s=pint_el_to_sim(ss_mass,E); double e=energy_sim(s); reb_simulation_free(s); return e; }
static void fine(const el Ein, el Eout){
    struct reb_simulation* s=pint_el_to_sim(ss_mass,Ein);
    s->dt=6.*DAY; s->exact_finish_time=0;
    reb_simulation_set_integrator(s,"whfast512");
    struct reb_integrator_whfast512_state* w=s->integrator.state;
    w->gr_potential=1; w->corrector=17; w->kernel=2; w->concatenate_steps=lround(slice_days/6.);
    reb_simulation_steps(s,1);
    pint_sim_to_el(s,Eout);
    reb_simulation_free(s);
}
static double poserr(const el A, const el B){
    struct reb_simulation* a=pint_el_to_sim(ss_mass,A); struct reb_simulation* b=pint_el_to_sim(ss_mass,B);
    double m=0; for (int i=1;i<=PINT_NP;i++){
        double d[3]={a->particles[i].x-b->particles[i].x,a->particles[i].y-b->particles[i].y,a->particles[i].z-b->particles[i].z};
        for (int x=0;x<3;x++) if (fabs(d[x])>m) m=fabs(d[x]); }
    reb_simulation_free(a); reb_simulation_free(b); return m;
}
static double wrapd(double x){ return x - 2.*M_PI*floor((x+M_PI)/(2.*M_PI)); }
static el *gU, *gF; static int gN, gP, gfirst;
static void* worker(void* arg){ const int t=(int)(long)arg; for (int n=gfirst+t;n<gN;n+=gP) fine(gU[n],gF[n]); return NULL; }
static void fine_sweep(el* U, el* Fv, int N, int first, int P){
    gU=U; gF=Fv; gN=N; gP=P; gfirst=first; pthread_t th[256];
    for (int t=0;t<P;t++) pthread_create(&th[t],NULL,worker,(void*)(long)t);
    for (int t=0;t<P;t++) pthread_join(th[t],NULL);
}
// scale barycentric velocities so that the energy becomes Et
static void project_energy(el E, double Et){
    struct reb_simulation* s=pint_el_to_sim(ss_mass,E);
    for (int it=0; it<2; it++){
        const double e=energy_sim(s);
        double K=0; for (unsigned i=0;i<s->N;i++){ struct reb_particle p=s->particles[i]; K+=0.5*p.m*(p.vx*p.vx+p.vy*p.vy+p.vz*p.vz); }
        const double f=sqrt(1.+(Et-e)/K);
        for (unsigned i=0;i<s->N;i++){ s->particles[i].vx*=f; s->particles[i].vy*=f; s->particles[i].vz*=f; }
    }
    pint_sim_to_el(s,E); reb_simulation_free(s);
}
static void write_rec(FILE* f, double t, const el E){ double e=energy_el(E); fwrite(&t,8,1,f); fwrite(&e,8,1,f); fwrite(E,8,48,f); }
int main(int argc, char** argv){
    if (argc<7){ fprintf(stderr,"usage: pint_long out_prefix windows N slice_days threads k [serial=1] [extra_k_tol=0]\n"); return 1; }
    const char* pre=argv[1]; const int W=atoi(argv[2]), N=atoi(argv[3]); slice_days=atof(argv[4]);
    const int P=atoi(argv[5]), K=atoi(argv[6]); const int ser = argc>7 ? atoi(argv[7]) : 1;
    const double xtol = argc>8 ? atof(argv[8]) : 0.;
    const int diag = argc>9 ? atoi(argv[9]) : 0;
    const int skern = argc>10 ? atoi(argv[10]) : 2;   // kernel of the serial run
    const int proj = argc>11 ? atoi(argv[11]) : 0;    // 1: energy projection of corrected boundary states     // 1: measure boundary jumps (extra fine sweeps, kernels 2 and 1)
    FILE* fj=NULL; if (diag){ char f2[512]; snprintf(f2,sizeof f2,"%s_jumps.txt",argv[1]); fj=fopen(f2,"w"); }
    struct pint_cfg c; pint_default_cfg(&c, ss_mass, slice_days);
    char fn[512];
    struct reb_simulation* r0=ss_create(6.,0); el E0; pint_sim_to_el(r0,E0); reb_simulation_free(r0);
    const double Tsl=slice_days*DAY;
    double tser=0;
    if (ser){
        snprintf(fn,sizeof fn,"%s_ser%d.bin",pre,skern); FILE* fs=fopen(fn,"wb");
        struct reb_simulation* s=pint_el_to_sim(ss_mass,E0);
        s->dt=6.*DAY; s->exact_finish_time=0;
        reb_simulation_set_integrator(s,"whfast512");
        struct reb_integrator_whfast512_state* w=s->integrator.state;
        w->gr_potential=1; w->corrector=17; w->kernel=skern; w->concatenate_steps=lround(slice_days/6.);
        el E; pint_sim_to_el(s,E); write_rec(fs,0.,E);
        double t0=now();
        for (long q=0;q<(long)W*N;q++){ reb_simulation_steps(s,1); pint_sim_to_el(s,E); write_rec(fs,(q+1)*Tsl,E); }
        tser=now()-t0; fclose(fs); reb_simulation_free(s);
        fprintf(stderr,"serial: %.2f s\n",tser);
    }
    snprintf(fn,sizeof fn,"%s_pr.bin",pre); FILE* fp=fopen(fn,"wb");
    snprintf(fn,sizeof fn,"%s_pr.log",pre); FILE* fl=fopen(fn,"w");
    static el U[NMAX+1], Gold[NMAX+1], Fv[NMAX+1], Un[NMAX+1];
    memcpy(U[0],E0,sizeof(el)); write_rec(fp,0.,E0);
    const double t0=now(); double tGtot=0; long nG=0, nF=0;
    for (int w=0; w<W; w++){
        for (int n=0;n<N;n++){ double a=now(); pint_coarse(&c,U[n],U[n+1]); tGtot+=now()-a; nG++; memcpy(Gold[n+1],U[n+1],sizeof(el)); }
        int k; double upd=0;
        for (k=1; k<=K+3; k++){
            if (k>K && !(xtol>0 && upd>xtol)) break;
            static double EU[NMAX+1], EF[NMAX+1];
            if (proj) for (int n=k-1;n<N;n++) EU[n]=energy_el(U[n]);
            fine_sweep(U,Fv,N,k-1,P); nF+=N-(k-1);
            if (proj) for (int n=k-1;n<N;n++) EF[n]=energy_el(Fv[n]);
            memcpy(Un[0],U[0],sizeof(el));
            for (int n=0;n<k;n++) memcpy(Un[n],U[n],sizeof(el));
            memcpy(Un[k],Fv[k-1],sizeof(el));
            for (int n=k;n<N;n++){
                el g; double a=now(); pint_coarse(&c,Un[n],g); tGtot+=now()-a; nG++;
                for (int i=0;i<PINT_NP;i++) for (int e=0;e<6;e++){
                    double d=Fv[n][i][e]-Gold[n+1][i][e]; if (e==1) d=wrapd(d);
                    Un[n+1][i][e]=g[i][e]+d;
                }
                for (int i=0;i<PINT_NP;i++) Un[n+1][i][1]=wrapd(Un[n+1][i][1]);
                if (proj) project_energy(Un[n+1], energy_el(Un[n]) + (EF[n]-EU[n]));
                memcpy(Gold[n+1],g,sizeof(el));
            }
            upd=0; for (int n=1;n<=N;n++){ double e=poserr(Un[n],U[n]); if (e>upd) upd=e; }
            memcpy(U,Un,sizeof(el)*(N+1));
        }
        fprintf(fl,"%d %d %.3e\n",w,k-1,upd);
        if (diag){
            // jump at boundary n: U_n - F(U_{n-1}) (energy, position); reference: F with kernel 1 vs kernel 2
            static el Fd[NMAX+1];
            fine_sweep(U,Fd,N,0,P);
            for (int n=0;n<N;n++){
                el F1; struct reb_simulation* s=pint_el_to_sim(ss_mass,U[n]);
                s->dt=6.*DAY; s->exact_finish_time=0; reb_simulation_set_integrator(s,"whfast512");
                struct reb_integrator_whfast512_state* ww=s->integrator.state;
                ww->gr_potential=1; ww->corrector=17; ww->kernel=1; ww->concatenate_steps=lround(slice_days/6.);
                reb_simulation_steps(s,1); pint_sim_to_el(s,F1); reb_simulation_free(s);
                const double eF=energy_el(Fd[n]), eU=energy_el(U[n+1]), e1=energy_el(F1);
                fprintf(fj,"%d %d %.6e %.6e %.6e %.6e\n",w,n,(eU-eF)/fabs(eF),(e1-eF)/fabs(eF),poserr(U[n+1],Fd[n]),poserr(F1,Fd[n]));
            }
        }
        for (int n=1;n<=N;n++) write_rec(fp,((long)w*N+n)*Tsl,U[n]);
        memcpy(U[0],U[N],sizeof(el));
    }
    const double tpr=now()-t0;
    fclose(fp); fclose(fl); if (fj) fclose(fj);
    fprintf(stderr,"parareal: %.2f s on %d threads (serial %.2f s); fine slices %ld (%.2f per slice of the span), coarse steps %ld, coarse time %.2f s\n",
            tpr,P,tser,nF,(double)nF/(W*N),nG,tGtot);
    return 0;
}
