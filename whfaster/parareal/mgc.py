# Mean-element coarse propagator (hybrid): inner planets (Mercury..Mars) by a first-order numerical mean-element
# theory (meanel.py) + Kepler/secular drift; outer planets by WHFast512 at a large step in a Sun(+inner mass)+4
# planet subsystem. Fourier tables are cached per slice boundary.
import numpy as np, rebound
from ptlib import *
import meanel
m=np.array(mass); Mint=np.cumsum(m)[:-1]; mu=Mint+m[1:]
class Coarse:
    def __init__(self, slice_days, Ng=64, dt_out=180, sec=1, ninner=4, sec2=0, live=0):
        self.days=slice_days; self.T=slice_days*DAY; self.Ng=Ng; self.dt_out=dt_out; self.sec=sec
        self.ni=ninner; self.tables={}; self.sec2=sec2; self.live=live
    def table(self, n, E):
        if self.live: return meanel.Theory(m,E,Ng=self.Ng,inner=self.ni)
        if n not in self.tables: self.tables[n]=meanel.Theory(m,E,Ng=self.Ng,inner=self.ni)
        return self.tables[n]
    def outerWH(self, E):
        ni=self.ni
        s=rebound.Simulation(); s.add(m=Mint[ni])
        for mm,e in zip(m[ni+1:],E[ni:]):
            ecc=np.hypot(e[2],e[3]); pom=np.arctan2(e[3],e[2]); si=np.hypot(e[4],e[5]); Om=np.arctan2(e[5],e[4])
            s.add(m=mm, a=np.exp(e[0]), l=e[1], e=ecc, pomega=pom, inc=2*np.arcsin(si), Omega=Om)
        s.move_to_com()
        s.integrator="whfast512"; s.dt=self.dt_out*DAY; s.exact_finish_time=0
        s.integrator.gr_potential=1; s.integrator.corrector=17; s.integrator.kernel=2
        nst=int(round(self.days/self.dt_out)); assert abs(nst*self.dt_out-self.days)<1e-6
        s.integrator.concatenate_steps=nst; s.steps(1)
        out=[]
        for o in s.orbits():
            si=np.sin(o.inc/2)
            out.append([np.log(o.a), o.l, o.e*np.cos(o.pomega), o.e*np.sin(o.pomega), si*np.cos(o.Omega), si*np.sin(o.Omega)])
        return np.array(out)
    def mean_inner(self, E, n, th=None):
        th=th or self.table(n,E); M=E.copy(); M[:self.ni]=th.mean(E)[:self.ni]; return M
    def osc_inner(self, M, n):
        th=self.table(n,M); Ee=M.copy()
        for _ in range(3): Ee[:self.ni]=M[:self.ni]+th.disp(Ee)[:self.ni]
        return Ee
    def drift(self, M, n, th0=None):
        """Kepler + first-order secular drift of the inner mean elements over one slice"""
        th0=th0 or self.table(n,M); ni=self.ni
        nbar=np.sqrt(mu/np.exp(3*M[:,0]))
        R=th0.secular_rates() if self.sec else np.zeros((8,6))
        M2=M.copy()
        M2[:ni,1]=wrap(M[:ni,1]+(nbar[:ni]+R[:ni,1])*self.T)
        M2[:ni,2:]+=R[:ni,2:]*self.T
        return M2
    def __call__(self, E, n):
        th=self.table(n,E)
        M=self.mean_inner(E,n,th)
        M2=self.drift(M,n,th)
        M2[self.ni:]=self.outerWH(E)
        return self.osc_inner(M2,n+1)
