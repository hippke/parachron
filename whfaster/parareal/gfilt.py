# Round 3, iteration 21: coarse step with numerically filtered mean elements (HMM / "stroboscopic averaging" idea).
# The first-order transform (all six components, live tables) leaves a pointwise residual (second-order and missing
# terms, periods 0.3-25 yr) which dominates the per-slice coarse error of Venus/Earth/Mars. A short fine run of L yr
# from the slice start, transformed point by point with the first-order theory and Hann-averaged (analytic drift
# removed), gives the mean elements at t0 without that residual. At the end, the osculating state reconstructed from
# the drifted mean elements is corrected by the same measurement on a short fine run started from it.
#   python3 gfilt.py L_yr [quality|eval] [end=0/1]
import numpy as np, sys, rebound
from ptlib import *
import meanel, mgc
m=np.array(mass); Mint=np.cumsum(m)[:-1]; mu=Mint+m[1:]; NI=4
def pair_rule(i,j):
    if j>5: return 0
    if j>=4: return 16
    return 64 if (i,j) in ((1,2),(2,3)) else 32

def sec_rule(i,j):
    if j>=4: return 16
    return 64 if (i,j) in ((1,2),(2,3)) else 32

class GF:
    def __init__(self, L_yr=10., ds=3, end=1, days=360000, comps=None, Ng=32):
        self.L=L_yr; self.ds=ds; self.end=end; self.days=days; self.T=days*DAY
        self.core=mgc.Coarse(days,Ng=Ng,dt_out=360,live=1,comps=comps,pair_ng=pair_rule,heun=1)
        self.nsamp=int(L_yr*365.25/(6*ds)); self.w=np.hanning(self.nsamp+2)[1:-1]
    def filtered(self, Es, th, rates):
        """Hann-weighted mean of (first-order mean elements - analytic drift) over a fine run of L yr from Es"""
        s=sim_from(from_el(Es)); s.integrator="whfast512"; s.dt=6*DAY; s.exact_finish_time=0
        s.integrator.gr_potential=1; s.integrator.corrector=17; s.integrator.concatenate_steps=self.ds
        acc=np.zeros((NI,6)); m0=th.mean(Es)[:NI]
        for q in range(self.nsamp):
            s.steps(1); E=to_el(state(s)); y=th.mean(E)[:NI]; dt=s.t
            y[:,1]=m0[:,1]+wrap(y[:,1]-m0[:,1]-rates[:,1]*dt)
            y[:,2:]-=rates[:,2:]*dt
            acc+=self.w[q]*y
        return acc/self.w.sum()
    def rates(self, M, th):
        R=th.secular_rates()[:NI].copy(); R[:,1]+=np.sqrt(mu[:NI]/np.exp(3*M[:NI,0])); R[:,0]=0; return R
    def __call__(self, E):
        c=self.core
        th=c.table(0,E); M=c.mean_inner(E,0,th)
        if self.nsamp: M[:NI]=self.filtered(E,th,self.rates(M,th))
        M2=c.drift(M,0,th); M2[NI:]=c.outerWH(E)
        th1=c.table(1,M2)
        R0=th.secular_rates(); R1=th1.secular_rates(); nbar=np.sqrt(mu/np.exp(3*M[:,0]))
        M2[:NI,1]=wrap(M[:NI,1]+(nbar[:NI]+0.5*(R0[:NI,1]+R1[:NI,1]))*self.T)
        M2[:NI,2:]=M[:NI,2:]+0.5*(R0[:NI,2:]+R1[:NI,2:])*self.T
        Ee=M2.copy()
        for _ in range(3): Ee[:NI]=M2[:NI]+th1.disp(Ee)[:NI]
        for _ in range(self.end if self.nsamp else 0):
            th2=c.table(1,Ee)
            Mt=self.filtered(Ee,th2,self.rates(M2,th2))
            d=M2[:NI]-Mt; d[:,1]=wrap(d[:,1]); Ee[:NI]+=d; Ee[:NI,1]=wrap(Ee[:NI,1])
        return Ee

if __name__=='__main__':
    L=float(sys.argv[1]); mode=sys.argv[2]; end=int(sys.argv[3]) if len(sys.argv)>3 else 1
    G=GF(L,end=end)
    if mode=='quality':
        import slice_quality
        slice_quality.quality(G,1000,int(sys.argv[4]) if len(sys.argv)>4 else 20,label=f'GF L={L} end={end}')
    else:
        import ph1eval
        ph1eval.evaluate(G,f'GF L={L} end={end}',eps=[0,1030,3400,6000,9300] if mode=='eval5' else None)

class GF2:
    """filtered mean elements on top of an mgc.Coarse core (nsub path): start mean from a short fine run
    (first-order transform with the core's start tables), end state corrected by the same measurement"""
    def __init__(self, core, L_yr=10., ds=3, end=1):
        self.core=core; self.f=GF.__new__(GF); self.f.L=L_yr; self.f.ds=ds
        self.f.nsamp=int(L_yr*365.25/(6*ds)); self.f.w=np.hanning(self.f.nsamp+2)[1:-1]
        self.end=end; core.mean_hook=self.mean_hook; core.end_hook=self.end_hook if end else None
    def mean_hook(self, E, th, M):
        self.th=th; M=M.copy(); M[:NI]=self.f.filtered(E,th,self.f.rates(M,th)); return M
    def end_hook(self, Ee, M2):
        c=self.core; th2=c.table(1,M2)
        Mt=self.f.filtered(Ee,th2,self.f.rates(M2,th2))
        d=M2[:NI]-Mt; d[:,1]=wrap(d[:,1]); Ee=Ee.copy(); Ee[:NI]+=d; Ee[:NI,1]=wrap(Ee[:NI,1]); return Ee
    def __call__(self, E): return self.core(E,0)
