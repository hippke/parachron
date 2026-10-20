/**
 * @file    rebound_internal.h
 * @brief   Various functions used internally and implemented in rebound.c
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section     LICENSE
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

#ifndef _REBOUND_INTERNAL_H
#define _REBOUND_INTERNAL_H

// Definitions for functions that we need to implement ourselves on Windows.
#ifdef _WIN32
#define _NO_CRT_STDIO_INLINE // WIN32 to use _vsprintf_s
#include <WinSock2.h>
#define _WINSOCKAPI_ //stops windows.h including winsock.h
#include <windows.h>
typedef struct reb_timeval {
    int64_t tv_sec;
    int64_t tv_usec;
} reb_timeval;
REB_API int gettimeofday(struct reb_timeval * tp, void* tzp);
#include <stdarg.h>
int asprintf(char **strp, const char *fmt, ...);
int vasprintf(char **strp, const char *fmt, va_list ap);
int rand_r (unsigned int *seed);
void usleep(__int64 usec);
#include <io.h>
#include <process.h> // for getpid()
#define REB_STR_RED
#define REB_STR_RED_BOLD
#define REB_STR_YELLOW_BOLD
#define REB_STR_BOLD
#define REB_STR_GREEN
#define REB_STR_RESET
#else // Linux and MacOS
#define reb_timeval timeval
#include <sys/time.h>
#include <unistd.h>
#include <pthread.h>
#define REB_STR_RED "\033[31m"
#define REB_STR_RED_BOLD "\033[1;31m"
#define REB_STR_YELLOW_BOLD "\033[1;33m"
#define REB_STR_BOLD "\033[1m"
#define REB_STR_GREEN "\033[32m"
#define REB_STR_RESET "\033[0m"
#endif // _WIN32

// Githash should be provided as a command line argument to the compiler. 
// If not, use this dummy.
#ifndef GITHASH
#define GITHASH notavailable0000000000000000000000000001 
#endif

// Graceful global interrupt handler 
#include <signal.h>
extern volatile sig_atomic_t reb_sigint;  
void reb_sigint_handler(int signum);


// Global constants and variables
REB_API extern const char* reb_build_str;   ///< Date and time build string.
REB_API extern const char* reb_version_str; ///< Version string.
REB_API extern const char* reb_githash_str; ///< Current git hash.
REB_API extern const char* reb_logo[26];    ///< Logo of rebound. 
REB_API extern const unsigned char reb_favicon_png[]; /// < Favicon in PNG format.
REB_API extern const unsigned int reb_favicon_len;
REB_API extern const size_t reb_messages_max_N;
REB_API extern const uint32_t reb_string_size_max; // defined in rebound.c

// Alligned alloc/free. Linux/Max/Windows
REB_API void* reb_aligned_alloc(size_t alignment, size_t size);
REB_API void reb_aligned_free(void *ptr);
// Free any pointer. Alias for free(). Used by python.
REB_API void reb_free(void* p);
// Get the next stored warning message. Used only if save_messages==1. Return value is 0 if no messages are present, 1 otherwise.
REB_API int reb_pop_message(char** messages, char* const buf); 
// Print or save a message. 
enum REB_MESSAGE_TYPE {
    REB_MESSAGE_TYPE_INFO = 'i',
    REB_MESSAGE_TYPE_ERROR = 'e',
    REB_MESSAGE_TYPE_WARNING = 'w',
};
void reb_message(char*** messages, int save_messages, enum REB_MESSAGE_TYPE type, const char* const msg);
// Returns 1 if floating point contraction are enabled. 0 otherwise. For unit testing.
REB_API int reb_check_fp_contract(); 
// Returns 1 if the CPU supports the AVX512.
REB_API int reb_avx512_available();
#ifdef _WIN32
int reb_avx2_fma_available(void); // Windows only; elsewhere the compiler builtins are used
#endif
// Wrapper method to set number of OpenMP threads from python.
REB_API void reb_omp_set_num_threads(int num_threads);
// Helper function for comparing strings
int reb_strcmp_ignore_whitespace(const char *s1, const char *s2);

// Lists to find available integrators
REB_API char** reb_integrators_registered(); // only used by python
extern struct reb_integrator_configuration* reb_integrator_configurations_custom;


// Function multiversioning for a few hot kernels: on x86-64 Linux/glibc with GCC or clang, an
// extra AVX2 version of the function is compiled and chosen at load time (ifunc) on CPUs that
// support it; other platforms and compilers get the plain function. Results are bitwise
// identical: floating-point contraction stays off (-std=c99 / ISO mode; the clone is built
// without "fma"), and +, *, /, sqrt are correctly rounded at any vector width; only the number
// of lanes per instruction and the VEX encoding change. Define REB_NO_TARGET_CLONES to disable.
// Also used in OpenMP builds: the serial code paths there are the same functions, and the
// bodies of parallel regions (outlined by the compiler) simply use the default target.
#if defined(__has_attribute)
#if __has_attribute(target_clones)
#define REB_HAS_TARGET_CLONES_ATTR 1
#endif
#endif
// Shared-memory parallelisation (OPENMP builds). A parallel region costs roughly a microsecond
// (fork, barrier, join), as much as a whole step of a small system. Loops therefore run in
// parallel only when more than one thread is available and the loop has at least
// REB_OMP_MIN_WORK units of work (pair interactions or particles, as documented at each call);
// otherwise the serial code path runs, which is the same code as in non-OpenMP builds.
#ifdef OPENMP
#include <omp.h>
#ifndef REB_OMP_MIN_WORK
#define REB_OMP_MIN_WORK 2048
#endif
static inline int reb_omp_parallel(const double work){
    return work >= REB_OMP_MIN_WORK && omp_get_max_threads() > 1 && !omp_in_parallel();
}
// Pair-force loops (basic and compensated gravity): the parallel kernels do more work than the
// serial half-pair loop and move data between cores (see gravity.c), so they need more work.
// Work = N*N_active pair interactions. Measured on an Intel i5-13500 (6 P + 8 E cores): the
// parallel kernel pays off from N ~ 256 (2^16 pairs) on, and the best thread count is about
// N_active/64 (REB_OMP_GRAVITY_ROWS_PER_THREAD; see reb_gravity_omp_threads in gravity.c).
#ifndef REB_OMP_MIN_WORK_GRAVITY
#define REB_OMP_MIN_WORK_GRAVITY 65536.
#endif
#ifndef REB_OMP_GRAVITY_ROWS_PER_THREAD
#define REB_OMP_GRAVITY_ROWS_PER_THREAD 64
#endif
static inline int reb_omp_parallel_gravity(const double work){
    // Cheap test first: small systems never call into the OpenMP runtime here.
    if (work < REB_OMP_MIN_WORK || work < REB_OMP_MIN_WORK_GRAVITY) return 0;
    return reb_omp_parallel(work);
}
#endif // OPENMP

#if defined(REB_HAS_TARGET_CLONES_ATTR) && defined(__x86_64__) && defined(__linux__) && defined(__GLIBC__) && !defined(REB_NO_TARGET_CLONES)
#define REB_TARGET_CLONES __attribute__((target_clones("avx2","default")))
#else
#define REB_TARGET_CLONES
#endif

#endif // _REBOUND_INTERNAL_H
