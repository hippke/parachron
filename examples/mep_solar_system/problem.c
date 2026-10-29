/**
 * Parallel-in-time integration of the Solar System with mean-element parareal (MEP)
 *
 * This example integrates the Sun and the eight planets with WHFast512 on many cores at
 * the same time (one thread per time slice). It writes a small checkpoint file and can be
 * stopped (Ctrl-C or the wall-clock limit) and restarted at any time: if the checkpoint
 * exists, the run is resumed from it, bitwise identically to an uninterrupted run.
 *
 * Usage: ./rebound [tmax in years] [wall-clock limit in seconds]
 * Linux only.
 *
 * Copyright (c) 2026 Michael Hippke. Part of Parachron (https://github.com/hippke/parachron), a modified
 * version of REBOUND; GNU GPL v3 or later (see LICENSE).
 */
#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Initial conditions for the Solar System from NASA Horizons (as in examples/whfast512_solar_system)
double all_ss_pos[9][3] = {
    {-0.008816286905115728, -0.0010954664916791675, 0.0002143249385447027},
    {-0.05942272929227954, -0.46308699693348293, -0.032897989948949075},
    {-0.7276101005375593, 0.006575003332463933, 0.041795901908847084},
    {-0.5789530452882667, -0.8361530119313055, 0.0002611520181901174},
    {-1.45202492400084, 0.827519404194876, 0.052981833432457694},
    {4.492983939852296, 2.0661626247490354, -0.10909246996001629},
    {8.4974210980544, -4.8620394993693585, -0.2537835862373596},
    {12.959111916929283, 14.760785302864473, -0.1130656917933948},
    {29.787987348666505, -2.51460654509393, -0.6347108842010732}
};

double all_ss_vel[9][3] = {
    {0.00014315746073017681, -0.0004912441820893999, 8.127678560998346e-07},
    {1.2978664284760637, -0.09524541469911743, -0.12677574364801253},
    {-0.019239782390457125, -1.1813975672919448, -0.01509392594251431},
    {0.8098712561282222, -0.5682496529341624, 2.6169897281383047e-05},
    {-0.37436417754222295, -0.6365841544564991, -0.004143932260467942},
    {-0.18818907783656452, 0.41919544951404614, 0.0024710497024424977},
    {0.14292308496870448, 0.2808676923735748, -0.010574288572728728},
    {-0.1734971049470612, 0.14019515029516152, 0.0027683484887051457},
    {0.014142947617173336, 0.18292110872737416, -0.004092845767710294}
};

double all_ss_mass[9] = {
    0.9999999999950272,
    1.6601208254808336e-07,
    2.447838287784771e-06,
    3.0404326489511185e-06,
    3.2271560828978514e-07,
    0.0009547919099366768,
    0.0002858856700231729,
    4.366249613200406e-05,
    5.151383772628957e-05
};


int main(int argc, char* argv[]){
    const double tmax = (argc>1 ? atof(argv[1]) : 1e6)*2.*M_PI;     // default: 1 Myr
    const double walltime = argc>2 ? atof(argv[2]) : 0.;             // default: no limit
    const char* checkpoint = "checkpoint.mep";
    struct reb_mep* mep;

    FILE* f = fopen(checkpoint, "rb");
    if (f){
        fclose(f);
        mep = reb_mep_create_from_checkpoint(checkpoint, NULL);
        if (!mep){ fprintf(stderr, "Cannot read %s.\n", checkpoint); return 1; }
        printf("Resuming from %s at t = %.0f yr (%lld windows done).\n", checkpoint, mep->t/(2.*M_PI), (long long)mep->windows);
    }else{
        struct reb_simulation* r = reb_simulation_create();
        for (int i=0; i<9; i++){
            struct reb_particle p = {.m=all_ss_mass[i], .x=all_ss_pos[i][0], .y=all_ss_pos[i][1], .z=all_ss_pos[i][2],
                                     .vx=all_ss_vel[i][0], .vy=all_ss_vel[i][1], .vz=all_ss_vel[i][2]};
            reb_simulation_add(r, p);
        }
        reb_simulation_move_to_com(r);
        r->dt = 6./365.25*2.*M_PI;      // 6 days
        r->exact_finish_time = 0;
        reb_simulation_set_integrator(r, "whfast512");
        struct reb_integrator_whfast512_state* w = r->integrator.state;
        w->gr_potential = 1;
        w->corrector = 17;

        struct reb_mep_options o;
        reb_mep_options_default(&o);    // cores: all physical cores; tol = 2e-6 AU; slices of 60000 steps
        mep = reb_mep_create(r, &o, checkpoint);
        reb_simulation_free(r);         // MEP keeps its own copy of the state
        if (!mep){ fprintf(stderr, "Cannot set up MEP.\n"); return 1; }
        printf("New run on %d cores, %.0f yr per window.\n", mep->opt.cores, mep->opt.cores*mep->opt.slice_steps*mep->dt/(2.*M_PI));
    }

    mep->opt.walltime = walltime;
    enum REB_STATUS status = reb_mep_integrate(mep, tmax);
    double kmean = 0.; for (int k=0; k<REB_MEP_KHIST; k++) kmean += k*mep->khist[k];
    printf("t = %.0f yr after %lld windows (mean %.2f iterations per window), %.1f s wall clock.\n",
           mep->t/(2.*M_PI), (long long)mep->windows, kmean/(mep->windows>0?mep->windows:1), mep->walltime_last);
    if (status==REB_STATUS_USER) printf("Stopped by the wall-clock limit. Run again to continue.\n");
    if (status==REB_STATUS_SIGINT) printf("Interrupted. Run again to continue from the last completed window.\n");
    if (status==REB_STATUS_GENERIC_ERROR) printf("Error: %s\n", mep->message);

    // The current state as a regular simulation, e.g. for output
    struct reb_simulation* r = reb_simulation_create();
    reb_mep_update_simulation(mep, r);
    struct reb_particle sun = r->particles[0];
    struct reb_orbit o = reb_orbit_from_particle(r->G, r->particles[1], sun);
    printf("Mercury: a = %.6f AU, e = %.6f\n", o.a, o.e);
    reb_simulation_free(r);
    reb_mep_free(mep);
    return 0;
}
