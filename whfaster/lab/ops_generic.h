typedef struct { double v[8]; } V;
typedef struct { unsigned char b; } K;   // bit i set: lane i true
#define W_LOOP for (int i=0; i<8; i++)
W_INL V w_set1(double a){ V r; W_LOOP r.v[i] = a; return r; }
W_INL V w_zero(void){ return w_set1(0.0); }
W_INL V w_load(const double* p){ V r; W_LOOP r.v[i] = p[i]; return r; }
W_INL void w_store(double* p, V a){ W_LOOP p[i] = a.v[i]; }
W_INL V w_add(V a, V b){ V r; W_LOOP r.v[i] = a.v[i] + b.v[i]; return r; }
W_INL V w_sub(V a, V b){ V r; W_LOOP r.v[i] = a.v[i] - b.v[i]; return r; }
W_INL V w_mul(V a, V b){ V r; W_LOOP r.v[i] = a.v[i] * b.v[i]; return r; }
W_INL V w_div(V a, V b){ V r; W_LOOP r.v[i] = a.v[i] / b.v[i]; return r; }
W_INL V w_sqrt(V a){ V r; W_LOOP r.v[i] = sqrt(a.v[i]); return r; }
W_INL V w_abs(V a){ V r; W_LOOP r.v[i] = fabs(a.v[i]); return r; }
W_INL V w_floor(V a){ V r; W_LOOP r.v[i] = floor(a.v[i]); return r; }
W_INL V w_fmadd(V a, V b, V c){ V r; W_LOOP r.v[i] = fma(a.v[i], b.v[i], c.v[i]); return r; }
W_INL V w_fnmadd(V a, V b, V c){ V r; W_LOOP r.v[i] = fma(-a.v[i], b.v[i], c.v[i]); return r; }
W_INL V w_fmsub(V a, V b, V c){ V r; W_LOOP r.v[i] = fma(a.v[i], b.v[i], -c.v[i]); return r; }
// Ordered, quiet comparisons: false if either operand is NaN (as _CMP_LT_OQ etc.).
W_INL K w_klt(V a, V b){ K k = {0}; W_LOOP if (a.v[i] < b.v[i]) k.b |= (unsigned char)(1u<<i); return k; }
W_INL K w_kgt(V a, V b){ K k = {0}; W_LOOP if (a.v[i] > b.v[i]) k.b |= (unsigned char)(1u<<i); return k; }
W_INL K w_kle(V a, V b){ K k = {0}; W_LOOP if (a.v[i] <= b.v[i]) k.b |= (unsigned char)(1u<<i); return k; }
W_INL K w_kand(K a, K b){ K k; k.b = (unsigned char)(a.b & b.b); return k; }
W_INL K w_kandnot(K a, K b){ K k; k.b = (unsigned char)(~a.b & b.b); return k; }
W_INL K w_knot(K a){ K k; k.b = (unsigned char)(~a.b); return k; }
W_INL int w_kany(K a){ return a.b != 0; }
W_INL K w_kall(void){ K k; k.b = 0xFF; return k; }
W_INL K w_kbits(unsigned int bits){ K k; k.b = (unsigned char)bits; return k; }
W_INL V w_sel(K k, V a, V b){ V r; W_LOOP r.v[i] = ((k.b >> i) & 1) ? a.v[i] : b.v[i]; return r; }
W_INL V w_maskz(K k, V a){ V r; W_LOOP r.v[i] = ((k.b >> i) & 1) ? a.v[i] : 0.0; return r; }
// Bit-trick seed for 1/sqrt(a): integer 0x5FE6EB50C7B537A9 - (bits(a) >> 1) (deterministic, same on every machine).
W_INL V w_rsqrt_seed(V a){ V r; W_LOOP { uint64_t b; memcpy(&b, &a.v[i], 8); b = 0x5FE6EB50C7B537A9ULL - (b >> 1); memcpy(&r.v[i], &b, 8); } return r; }
W_INL V w_permute(V a, const int* idx){ V r; W_LOOP r.v[i] = a.v[idx[i]]; return r; }
static const int w_idx_3201[8] = {3,2,0,1, 7,6,4,5};
static const int w_idx_2310[8] = {2,3,1,0, 6,7,5,4};
static const int w_idx_3012[8] = {3,0,1,2, 7,4,5,6};
static const int w_idx_swap[8] = {1,0,3,2, 5,4,7,6};
static const int w_idx_b3[8]   = {4,5,6,7, 1,2,3,0};
static const int w_idx_b4[8]   = {5,6,7,4, 2,3,0,1};
static const int w_idx_b34[8]  = {7,4,5,6, 0,1,2,3};
W_INL V w_perm_3201(V a){ return w_permute(a, w_idx_3201); }
W_INL V w_perm_2310(V a){ return w_permute(a, w_idx_2310); }
W_INL V w_perm_3012(V a){ return w_permute(a, w_idx_3012); }
W_INL V w_swap_pairs(V a){ return w_permute(a, w_idx_swap); }
W_INL V w_perm_b3(V a){ return w_permute(a, w_idx_b3); }
W_INL V w_perm_b4(V a){ return w_permute(a, w_idx_b4); }
W_INL V w_perm_b34(V a){ return w_permute(a, w_idx_b34); }

W_INL V w_bcast_lane(V a, int k){ return w_set1(a.v[k]); }
