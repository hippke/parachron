# Round 3: value of a perfect mean semi-major axis. For each Phase-0 reference slice n, the slice average of the
# first-order mean ln a (C transform with ng_hi=64, fine run sampled every 30 d) gives the "true" mean m_bar_n; the
# pointwise transform error at the boundary is m(ref_n) - m_bar_n. G_or uses m_bar_n inside the slice and the true
# end-point offset (oracle; shows what a better transform -- second order -- can buy).
#   python3 oracle_lna.py compute EP...    -> oracle_lna_<ep>.npz
#   python3 oracle_lna.py eval
import numpy as np, sys
from ptlib import *
import ptc, ph1eval
D=ph1eval.D; N=10; days=360000
Gt=ptc.CoarseC(days,ng_hi=64)
def m_of(E): return E[:4,0]-ptc.lna_disp(Gt,E)
if sys.argv[1]=='compute':
    for ep in map(int,sys.argv[2:]):
        ref=D[f'ref{ep}']; mbar=[]
        for n in range(N):
            s=sim_from(from_el(ref[n])); s.integrator="whfast512"; s.dt=6*DAY; s.exact_finish_time=0
            s.integrator.gr_potential=1; s.integrator.corrector=17; s.integrator.concatenate_steps=5
            acc=np.zeros(4); q=0
            for _ in range(days//30):
                s.steps(1); acc+=m_of(to_el(state(s))); q+=1
            mbar.append(acc/q)
        np.savez(f'oracle_lna_{ep}.npz',mbar=np.array(mbar),m=np.array([m_of(r) for r in ref]))
        print('epoch',ep,'done'); sys.stdout.flush()
else:
    G0=ptc.CoarseC(days,heun=1,ng_hi=64)
    OR={ep:np.load(f'oracle_lna_{ep}.npz') for ep in ph1eval.EP}
    cur={}
    def make(ep, mode):
        ref=D[f'ref{ep}']; o=OR[ep]
        def G(E):
            n=int(np.argmin([np.sum(ediff(E,r)[4:,1]**2) for r in ref[:N]]))
            Ei=E.copy(); Ei[:4,0]+=o['mbar'][n]-o['m'][n]
            g=G0(Ei)
            if mode==2: g[:4,0]+=o['m'][n+1]-o['mbar'][n]
            else: g[:4,0]-=o['mbar'][n]-o['m'][n]
            return g
        return G
    for mode,lab in ((1,'oracle start mean only'),(2,'oracle start+end')):
        R=[]
        for ep in ph1eval.EP: R.append(ph1eval.evaluate(make(ep,mode),lab,eps=[ep],verbose=True)[0])
        R=np.array(R); print(lab,"median",np.median(R,0),"worst",R.max(0),"k2<=2:",np.sum(R[:,1]<=2))
