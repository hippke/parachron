#include <stdint.h>
#include <string.h>
#include <signal.h>
#include <math.h>
#include "integrator_whfast512_data.h"
#define REB_W512_TARGET
#define W_INL static inline __attribute__((always_inline))
#define REB_W512_NOINLINE_FALLBACK static __attribute__((noinline))
#define REB_W512_STEP_FUNC static __attribute__((noinline))
#include "ops_generic.h"
#define REB_W512_SUFFIX SUFFIX
#include "kernel_lab.h"
