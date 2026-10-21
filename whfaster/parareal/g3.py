# Coarse propagator G3 = config B (C) + second-order resonant terms (so2.py), live.
import numpy as np
from ptlib import *
import ptc, so2
TERMS=[{2:-4,3:8,4:-3}]          # Earth-Mars-Jupiter
class G3:
    def __init__(self, days=360000, terms=TERMS, Ng=32, base=None, **kw):
        self.G0=ptc.CoarseC(days,heun=1,**kw); self.base=base or self.G0; self.T=days*DAY; self.terms=terms; self.Ng=Ng
    def mean_motions(self, E, g):
        d,R=ptc.lna_disp_R(self.G0,E)
        n=np.sqrt(so2.MU/np.exp(3*E[:,0]))
        for p in range(4): n[p]=np.sqrt(so2.MU[p]/np.exp(3*(E[p,0]-d[p])))+R[p,1]
        for p in range(4,8):     # giants: mean motion over the slice from their WH run
            dl=wrap(g[p,1]-E[p,1]); m=np.round((n[p]*self.T-dl)/(2*np.pi)); n[p]=(dl+2*np.pi*m)/self.T
        return n
    def __call__(self, E, nslice=0):
        g=self.base(E); n=self.mean_motions(E,g); T=self.T
        for K in self.terms:
            U=so2.resonant_U2(E,K,Ng=self.Ng)
            nu=sum(k*n[p] for p,k in K.items())
            p0=sum(k*E[p,1] for p,k in K.items()); e0=np.exp(1j*p0); e1=np.exp(1j*(p0+nu*T))
            for p in K:
                if p>=4: continue
                A=2*U[p][0]/(1j*nu)
                g[p,0]+=np.real(A*(e1-e0))
                dl=-1.5*n[p]*T*np.real(A*((e1-e0)/(1j*nu*T)-e0)) + np.real(2*U[p][1]/(1j*nu)*(e1-e0))
                g[p,1]=wrap(g[p,1]+dl)
                for c in range(2,6): g[p,c]+=np.real(2*U[p][c]/(1j*nu)*(e1-e0))
        return g
if __name__=='__main__':
    import ph1eval, sys, time
    G=G3(); E=ph1eval.D['ref0'][0]; t=time.time(); G(E); print(f"G3 call {time.time()-t:.2f}s")
    EPS=[int(x) for x in sys.argv[1].split(',')] if len(sys.argv)>1 else ph1eval.EP
    ph1eval.evaluate(G,'G3 (+E-M-J 2nd order)',eps=EPS)
