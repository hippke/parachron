# Mean-element coarse propagator (hybrid): inner planets (Mercury..Mars) by a first-order numerical mean-element
# theory (meanel.py) + Kepler/secular drift; outer planets by WHFast512 at a large step in a Sun(+inner mass)+4
# planet subsystem. Fourier tables are cached per slice boundary.
import numpy as np, rebound
from ptlib import *
import meanel
m=np.array(mass); Mint=np.cumsum(m)[:-1]; mu=Mint+m[1:]
class Coarse:
    def __init__(self, slice_days, Ng=64, dt_out=180, sec=1, ninner=4, sec2=0, live=0, live_sec=0, Ng_sec=None, nsub=0, comps=None, pair_ng=None, heun=0, meann=0, giant_mean=0, nsec=0, sec_ng=None, gavg=0., w2=None):
        self.days=slice_days; self.T=slice_days*DAY; self.Ng=Ng; self.dt_out=dt_out; self.sec=sec
        self.ni=ninner; self.tables={}; self.sec2=sec2; self.live=live; self.live_sec=live_sec; self.Ng_sec=Ng if Ng_sec is None else Ng_sec; self.nsub=nsub; self.comps=comps; self.pair_ng=pair_ng; self.heun=heun; self.meann=meann; self.giant_mean=giant_mean; self.nsec=nsec; self.sec_ng=sec_ng; self.gavg=gavg; self.w2=w2; self._WS=None; self._U2=None
    def table(self, n, E):
        if self.live: th=meanel.Theory(m,E,Ng=self.Ng,inner=self.ni,pair_ng=self.pair_ng); th.comps=self.comps; return th
        if n not in self.tables:
            self.tables[n]=meanel.Theory(m,E,Ng=self.Ng,inner=self.ni,pair_ng=self.pair_ng); self.tables[n].comps=self.comps
        return self.tables[n]
    def outerWH(self, E, nseg=1):
        """outer planets by WH; returns list of element arrays at nseg+1 equally spaced times (incl. start)"""
        ni=self.ni
        s=rebound.Simulation(); s.add(m=Mint[ni])
        for mm,e in zip(m[ni+1:],E[ni:]):
            ecc=np.hypot(e[2],e[3]); pom=np.arctan2(e[3],e[2]); si=np.hypot(e[4],e[5]); Om=np.arctan2(e[5],e[4])
            s.add(m=mm, a=np.exp(e[0]), l=e[1], e=ecc, pomega=pom, inc=2*np.arcsin(si), Omega=Om)
        s.move_to_com()
        s.integrator="whfast512"; s.dt=self.dt_out*DAY; s.exact_finish_time=0
        s.integrator.gr_potential=1; s.integrator.corrector=17; s.integrator.kernel=2
        nst=int(round(self.days/self.dt_out)); assert abs(nst*self.dt_out-self.days)<1e-6
        assert nst%nseg==0
        def els():
            out=[]
            for o in s.orbits():
                si=np.sin(o.inc/2)
                out.append([np.log(o.a), o.l, o.e*np.cos(o.pomega), o.e*np.sin(o.pomega), si*np.cos(o.Omega), si*np.sin(o.Omega)])
            return np.array(out)
        res=[E[ni:].copy()]
        if self.gavg>0 and nseg>1:          # dense samples; giants' elements box-averaged over gavg yr (not lambda)
            ds=5; s.integrator.concatenate_steps=ds; nd=nst//ds
            S=[E[ni:].copy()]
            for q in range(nd): s.steps(1); S.append(els())
            S=np.array(S); tt=np.arange(nd+1)*ds*self.dt_out/365.25; Ty=tt[-1]; w=self.gavg
            out=[]
            for q in range(nseg+1):
                tq=Ty*q/nseg; lo=min(max(tq-w/2,0),Ty-w); sel=(tt>=lo)&(tt<=lo+w)
                X=S[int(round(q*nd/nseg))].copy(); X[:,[0,2,3,4,5]]=S[sel][:,:,[0,2,3,4,5]].mean(0); out.append(X)
            out[-1][:,:]=out[-1]; out[-1]=np.where(np.arange(6)[None,:]==1,S[-1],out[-1])
            self._giants_end=S[-1]
            return out
        s.integrator.concatenate_steps=nst//nseg
        for q in range(nseg):
            s.steps(1); res.append(els())
        return res if nseg>1 else res[-1]
    def mean_inner(self, E, n, th=None):
        th=th or self.table(n,E); M=E.copy(); M[:self.ni]=th.mean(E)[:self.ni]; return M
    def osc_inner(self, M, n):
        th=self.table(n,M); Ee=M.copy()
        if self._WS is not None:     # second order: X = Y + W1(Y) + W2(Y), evaluated at the mean elements Y
            Ee[:self.ni]=M[:self.ni]+th.disp(M)[:self.ni]+self._WS.disp(M)[:self.ni]
            Ee[:self.ni,1]=wrap(Ee[:self.ni,1]); return Ee
        for _ in range(3):
            Ee[:self.ni]=M[:self.ni]+th.disp(Ee)[:self.ni]
        return Ee
    def drift(self, M, n, th0=None):
        """Kepler + first-order secular drift of the inner mean elements over one slice"""
        th0=th0 or self.table(n,M); ni=self.ni
        nbar=np.sqrt(mu/np.exp(3*M[:,0]))
        if self.live_sec: R=meanel.Theory(m,M,Ng=self.Ng_sec,inner=ni).secular_rates()
        else: R=th0.secular_rates() if self.sec else np.zeros((8,6))
        M2=M.copy()
        M2[:ni,1]=wrap(M[:ni,1]+(nbar[:ni]+R[:ni,1])*self.T)
        M2[:ni,2:]+=R[:ni,2:]*self.T
        return M2
    def drift_rk(self, M, outs):
        """RK4 over nsub substeps on the inner mean elements; outer elements from the WH samples outs
        (2*nsub+1 samples: substep boundaries and midpoints)"""
        ni=self.ni; h=self.T/self.nsub
        def rate(X, q):
            Y=X.copy(); Y[ni:]=outs[q]
            R=meanel.secular_rates_fast(m,Y,Ng=self.Ng_sec,inner=ni,pair_ng=(self.sec_ng or self.pair_ng) if self.Ng_sec==0 else None)
            D=np.zeros_like(X)
            D[:ni,1]=np.sqrt(mu[:ni]/np.exp(3*X[:ni,0]))+R[:ni,1]
            D[:ni,2:]=R[:ni,2:]
            if self._U2 is not None: D[:ni]+=self._U2[:ni]
            return D
        X=M.copy()
        for q in range(self.nsub):
            k1=rate(X,2*q); k2=rate(X+0.5*h*k1,2*q+1); k3=rate(X+0.5*h*k2,2*q+1); k4=rate(X+h*k3,2*q+2)
            X=X+h/6*(k1+2*k2+2*k3+k4)
        X[:ni,1]=wrap(X[:ni,1])
        return X
    def __call__(self, E, n):
        th=self.table(n,E)
        M=self.mean_inner(E,n,th)
        for _ in range(self.meann):          # tables (divisors, grids) at the mean semi-major axes
            E1=E.copy(); E1[:self.ni,0]=M[:self.ni,0]
            if self.giant_mean:              # giants' mean ln a (first order, giant-giant pairs) for the divisors
                tg=meanel.Theory(m,E,Ng=32,pair_ng=lambda i,j: 32 if i>=self.ni else 0,gr=False)
                E1[self.ni:,0]=tg.mean(E)[self.ni:,0]
            nov=None
            if self.nsec:                    # divisors from the mean motions incl. the secular lambda rates
                nov=np.sqrt(mu/np.exp(3*E1[:,0])); nov[:self.ni]+=th.secular_rates()[:self.ni,1]
            th=meanel.Theory(m,E1,Ng=self.Ng,inner=self.ni,pair_ng=self.pair_ng,n_override=nov); th.comps=self.comps
            M=self.mean_inner(E,n,th)
        self._WS=None; self._U2=None
        if self.w2 is not None and getattr(self,'freeze_U2',0) and getattr(self,'_U2frozen',None) is not None:
            self._U2=self._U2frozen.copy()
        elif self.w2 is not None:                  # second-order transform and secular rates (w2triple)
            import w2triple
            E1=E.copy(); E1[:self.ni,0]=M[:self.ni,0]
            WS=w2triple.W2Set(E1,n=np.sqrt(mu/np.exp(3*E1[:,0])),pairs=self.w2[0],triples=self.w2[1],nu_cut=self.w2[2] if len(self.w2)>2 else 1/500)
            Y=M.copy()
            opt=self.w2[3] if len(self.w2)>3 else ''
            if 'Uonly' in opt:          # X-form first-order transforms kept; only the second-order secular rates
                U=WS.U2()
                if 'Ukh' in opt: U[:,[0,1]]=0
                self._U2=U; WS=None
            if WS is None: pass
            elif 'noW' in opt: WS.disp=lambda X: np.zeros((8,6))
            if 'nolamW' in opt:
                _d=WS.disp; WS.disp=lambda X,_d=_d: (lambda D: (D.__setitem__((slice(None),1),0),D)[1])(_d(X))
            if WS is not None:
                for _ in range(3):
                    Y[:self.ni]=E[:self.ni]-th.disp(Y)[:self.ni]-WS.disp(Y)[:self.ni]
                Y[:self.ni,1]=wrap(Y[:self.ni,1]); M=Y; self._WS=WS; self._U2=WS.U2()
            opt=self.w2[3] if len(self.w2)>3 else ''
            if 'noU' in opt: self._U2=None
            if 'nolamU' in opt and self._U2 is not None: self._U2[:,1]=0
            if getattr(self,'freeze_U2',0): self._U2frozen=None if self._U2 is None else self._U2.copy()
        if getattr(self,'mean_hook',None): M=self.mean_hook(E,th,M)
        if self.nsub:
            outs=self.outerWH(E,2*self.nsub)
            M2=self.drift_rk(M,outs)
            M2[self.ni:]=self._giants_end if (self.gavg>0) else outs[-1]
        else:
            M2=self.drift(M,n,th)
            M2[self.ni:]=self.outerWH(E)
            if self.heun:
                ni=self.ni; th1=self.table(n+1,M2)
                R0=th.secular_rates(); R1=th1.secular_rates()
                nbar=np.sqrt(mu/np.exp(3*M[:,0]))
                M2[:ni,1]=wrap(M[:ni,1]+(nbar[:ni]+0.5*(R0[:ni,1]+R1[:ni,1]))*self.T)
                M2[:ni,2:]=M[:ni,2:]+0.5*(R0[:ni,2:]+R1[:ni,2:])*self.T
                Ee=M2.copy()
                for _ in range(3): Ee[:ni]=M2[:ni]+th1.disp(Ee)[:ni]
                return Ee
        Ee=self.osc_inner(M2,n+1)
        if getattr(self,'end_hook',None): Ee=self.end_hook(Ee,M2)
        return Ee
