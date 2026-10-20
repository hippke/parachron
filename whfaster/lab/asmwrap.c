#include <stdint.h>
#include <signal.h>
#include "integrator_whfast512_data.h"
int reb_whfast512_full_steps_gr1_n1_encounter0_escape0(struct simd_data* data, uint64_t* N_steps, int skip_first_kepler_step, volatile sig_atomic_t* sigint);
int reb_whfast512c_full_steps_asm(struct simd_data* d, uint64_t* N, int skip, volatile sig_atomic_t* sig, int gr, int nsys, int enc, int esc){
    (void)gr;(void)nsys;(void)enc;(void)esc; return reb_whfast512_full_steps_gr1_n1_encounter0_escape0(d,N,skip,sig);
}
