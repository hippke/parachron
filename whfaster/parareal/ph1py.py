# Python coarse variants (mgc, live tables, Heun) on 5 epochs
import sys, numpy as np, mgc, ph1eval
EPS=[0,1030,3400,6000,9300]
def pair_rule(hi):
    def f(i,j):
        if j>5: return 0
        if j>=4: return 16
        A=np.exp(ph1eval.D['ref0'][0][:,0]); return hi if (hi and A[i]/A[j]>=0.6) else 32
    return f
for spec in sys.argv[1:]:
    comps=None if 'all' in spec else [0]
    hi=64 if 'hi' in spec else 0
    CG=mgc.Coarse(360000,Ng=32,dt_out=360,live=1,comps=comps,pair_ng=pair_rule(hi),heun=1)
    ph1eval.evaluate(lambda E: CG(E,0), spec, eps=EPS)
