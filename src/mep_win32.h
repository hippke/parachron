/**
 * @file    mep_win32.h
 * @brief   Windows (clang-cl) support of MEP: the pthread calls of src/mep.c and src/mep_coarse.c on Win32 threads,
 *          slim reader/writer locks, condition variables and one-time initialisation, and the C99 complex arithmetic
 *          of the coarse step without <complex.h> (the Microsoft header has no C99 complex type).
 * @author  Michael Hippke
 *
 * @section LICENSE
 * Copyright (c) 2026 Michael Hippke
 *
 * This file is part of Parachron, a modified version of REBOUND (https://github.com/hannorein/rebound).
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
#ifndef MEP_WIN32_H
#define MEP_WIN32_H
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES       // M_PI in the Microsoft <math.h>
#endif
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

// ---- complex numbers: clang's built-in _Complex type with the C99 names. Multiplication and division of two complex
// numbers are calls to __muldc3 and __divdc3 (compiler-rt), which the Microsoft linker does not provide: mep_coarse.c
// defines them.
#undef complex                  // the Microsoft <math.h> defines complex as its struct _complex
#define complex _Complex
#define I (__builtin_complex(0.0, 1.0))
#define creal(z) __builtin_creal(z)
#define cimag(z) __builtin_cimag(z)
#define conj(z) __builtin_conj(z)
static inline double mep_cabs(double _Complex z){ return hypot(__builtin_creal(z), __builtin_cimag(z)); }
static inline double _Complex mep_cexp(double _Complex z){
    const double e = exp(__builtin_creal(z));
    return __builtin_complex(e*cos(__builtin_cimag(z)), e*sin(__builtin_cimag(z)));
}
#define cabs mep_cabs
#define cexp mep_cexp
#endif // MEP_WIN32_H

#if defined(MEP_WIN32_THREADS) && !defined(MEP_WIN32_THREADS_H)
#define MEP_WIN32_THREADS_H
#include <WinSock2.h>           // before windows.h, as in rebound_internal.h
#define _WINSOCKAPI_
#include <windows.h>
#include <process.h>
#undef near                     // empty macros of windows.h that are variable names in the coarse step
#undef far
typedef HANDLE pthread_t;
typedef SRWLOCK pthread_mutex_t;
typedef CONDITION_VARIABLE pthread_cond_t;
typedef INIT_ONCE pthread_once_t;
#define PTHREAD_MUTEX_INITIALIZER SRWLOCK_INIT
#define PTHREAD_COND_INITIALIZER CONDITION_VARIABLE_INIT
#define PTHREAD_ONCE_INIT INIT_ONCE_STATIC_INIT
static inline int pthread_mutex_init(pthread_mutex_t* m, const void* a){ (void)a; InitializeSRWLock(m); return 0; }
static inline int pthread_mutex_destroy(pthread_mutex_t* m){ (void)m; return 0; }
static inline int pthread_mutex_lock(pthread_mutex_t* m){ AcquireSRWLockExclusive(m); return 0; }
static inline int pthread_mutex_unlock(pthread_mutex_t* m){ ReleaseSRWLockExclusive(m); return 0; }
static inline int pthread_mutex_trylock(pthread_mutex_t* m){ return TryAcquireSRWLockExclusive(m) ? 0 : 16; }   // EBUSY
static inline int pthread_cond_init(pthread_cond_t* c, const void* a){ (void)a; InitializeConditionVariable(c); return 0; }
static inline int pthread_cond_destroy(pthread_cond_t* c){ (void)c; return 0; }
static inline int pthread_cond_wait(pthread_cond_t* c, pthread_mutex_t* m){ SleepConditionVariableSRW(c, m, INFINITE, 0); return 0; }
static inline int pthread_cond_signal(pthread_cond_t* c){ WakeConditionVariable(c); return 0; }
static inline int pthread_cond_broadcast(pthread_cond_t* c){ WakeAllConditionVariable(c); return 0; }
static BOOL CALLBACK mep_once_cb(PINIT_ONCE o, PVOID fn, PVOID* ctx){ (void)o; (void)ctx; ((void (*)(void))fn)(); return TRUE; }
static inline int pthread_once(pthread_once_t* o, void (*fn)(void)){ InitOnceExecuteOnce(o, mep_once_cb, (PVOID)fn, NULL); return 0; }
struct mep_thread_start { void* (*fn)(void*); void* arg; };
static unsigned __stdcall mep_thread_main(void* p){
    struct mep_thread_start s = *(struct mep_thread_start*)p; free(p); s.fn(s.arg); return 0;
}
static inline int pthread_create(pthread_t* t, const void* attr, void* (*fn)(void*), void* arg){
    (void)attr;
    struct mep_thread_start* s = malloc(sizeof *s); if (!s) return 12;      // ENOMEM
    s->fn = fn; s->arg = arg;
    const uintptr_t h = _beginthreadex(NULL, 0, mep_thread_main, s, 0, NULL);
    if (!h){ free(s); return 11; }                                          // EAGAIN
    *t = (HANDLE)h; return 0;
}
static inline int pthread_join(pthread_t t, void** ret){ (void)ret; WaitForSingleObject(t, INFINITE); CloseHandle(t); return 0; }
static inline int pthread_detach(pthread_t t){ CloseHandle(t); return 0; }
#endif // MEP_WIN32_THREADS
