typedef __m512d V;
typedef __mmask8 K;

#if defined(_MSC_VER) && !defined(__clang__)
#define W_INL static __forceinline
#define REB_W512_NOINLINE_FALLBACK static __declspec(noinline)
#else
#define W_INL static inline __attribute__((always_inline))
#define REB_W512_NOINLINE_FALLBACK static __attribute__((noinline))
#endif

W_INL V w_set1(double a){ return _mm512_set1_pd(a); }
W_INL V w_zero(void){ return _mm512_setzero_pd(); }
W_INL V w_load(const double* p){ return _mm512_load_pd(p); }
W_INL void w_store(double* p, V a){ _mm512_store_pd(p, a); }
W_INL V w_add(V a, V b){ return _mm512_add_pd(a, b); }
W_INL V w_sub(V a, V b){ return _mm512_sub_pd(a, b); }
W_INL V w_mul(V a, V b){ return _mm512_mul_pd(a, b); }
W_INL V w_div(V a, V b){ return _mm512_div_pd(a, b); }
W_INL V w_sqrt(V a){ return _mm512_sqrt_pd(a); }
W_INL V w_abs(V a){ return _mm512_castsi512_pd(_mm512_and_si512(_mm512_castpd_si512(a), _mm512_set1_epi64(0x7FFFFFFFFFFFFFFFLL))); }
W_INL V w_floor(V a){ return _mm512_roundscale_pd(a, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC); }
W_INL V w_fmadd(V a, V b, V c){ return _mm512_fmadd_pd(a, b, c); }
W_INL V w_fnmadd(V a, V b, V c){ return _mm512_fnmadd_pd(a, b, c); }
W_INL V w_fmsub(V a, V b, V c){ return _mm512_fmsub_pd(a, b, c); }
W_INL K w_klt(V a, V b){ return _mm512_cmp_pd_mask(a, b, _CMP_LT_OQ); }
W_INL K w_kgt(V a, V b){ return _mm512_cmp_pd_mask(a, b, _CMP_GT_OQ); }
W_INL K w_kle(V a, V b){ return _mm512_cmp_pd_mask(a, b, _CMP_LE_OQ); }
W_INL K w_kand(K a, K b){ return (K)(a & b); }
W_INL K w_kandnot(K a, K b){ return (K)(~a & b); }
W_INL K w_knot(K a){ return (K)(~a); }
W_INL int w_kany(K a){ return a != 0; }
W_INL K w_kall(void){ return (K)0xFF; }
W_INL K w_kbits(unsigned int bits){ return (K)bits; }
W_INL V w_sel(K k, V a, V b){ return _mm512_mask_blend_pd(k, b, a); }
W_INL V w_maskz(K k, V a){ return _mm512_maskz_mov_pd(k, a); }
// Bit-trick seed for 1/sqrt(a): integer 0x5FE6EB50C7B537A9 - (bits(a) >> 1) (deterministic, same on every machine).
W_INL V w_rsqrt_seed(V a){ return _mm512_castsi512_pd(_mm512_sub_epi64(_mm512_set1_epi64(0x5FE6EB50C7B537A9LL), _mm512_srli_epi64(_mm512_castpd_si512(a), 1))); }
W_INL V w_perm_3201(V a){ return _mm512_permutex_pd(a, 0x4B); }
W_INL V w_perm_2310(V a){ return _mm512_permutex_pd(a, 0x1E); }
W_INL V w_perm_3012(V a){ return _mm512_permutex_pd(a, 0x93); }
W_INL V w_swap_pairs(V a){ return _mm512_permute_pd(a, 0x55); }
W_INL V w_perm_b3(V a){ return _mm512_permutexvar_pd(_mm512_set_epi64(0,3,2,1,7,6,5,4), a); }
W_INL V w_perm_b4(V a){ return _mm512_permutexvar_pd(_mm512_set_epi64(1,0,3,2,4,7,6,5), a); }
W_INL V w_perm_b34(V a){ return _mm512_permutexvar_pd(_mm512_set_epi64(3,2,1,0,6,5,4,7), a); }

// Broadcast lane k (compile-time constant) within the register.
#define w_bcast_lane(a, k) _mm512_permutexvar_pd(_mm512_set1_epi64(k), (a))
