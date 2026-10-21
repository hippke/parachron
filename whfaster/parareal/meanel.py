# First-order (in the planetary masses) short-period theory for Jacobi elements, computed numerically.
# For every planet pair (i<j) the perturbing accelerations on both Jacobi vectors are evaluated on a 2-D grid of
# mean longitudes (lambda_i, lambda_j) with all other elements frozen; Gauss' equations dE/dt = (dE/dv) a_pert give
# the element rates, whose 2-D FFT yields the short-period displacements
#     E_sp = sum_{k != 0} f_k exp(i k.lambda) / (i k.nu)       (+ the shear term -1.5 n * lna_sp / (i k.nu) for lambda)
# Mean elements:  M = E - sum_pairs E_sp(E)  (first order). GR (1/r^2 potential, as in WHFast512 gr_potential=1)
# is included as a single-planet term on a 1-D grid.
# Element vector per planet: [ln a, lambda, k=e cos pomega, h=e sin pomega, q=sin(i/2) cos Omega, p=sin(i/2) sin Omega]
import numpy as np
G=1.0
C_GR=10065.32

def el_to_rv(mu, E):
    """E[...,6] -> r[...,3], v[...,3]; mu broadcastable to E[...,0]"""
    a=np.exp(E[...,0]); lam=E[...,1]; k=E[...,2]; h=E[...,3]; q=E[...,4]; p=E[...,5]
    e=np.hypot(k,h); pom=np.arctan2(h,k)
    si=np.hypot(q,p); Om=np.arctan2(p,q); inc=2*np.arcsin(np.minimum(si,1.0))
    M=lam-pom
    Ea=M.copy()
    for _ in range(30):
        Ea=Ea-(Ea-e*np.sin(Ea)-M)/(1-e*np.cos(Ea))
    cE=np.cos(Ea); sE=np.sin(Ea); b=np.sqrt(1-e*e)
    n=np.sqrt(mu/a**3)
    xp=a*(cE-e); yp=a*b*sE
    rr=a*(1-e*cE)
    vxp=-a*n*sE/(1-e*cE); vyp=a*n*b*cE/(1-e*cE)
    om=pom-Om
    co,so=np.cos(om),np.sin(om); cO,sO=np.cos(Om),np.sin(Om); ci,s_i=np.cos(inc),np.sin(inc)
    P=np.stack([cO*co-sO*so*ci, sO*co+cO*so*ci, so*s_i],-1)
    Q=np.stack([-cO*so-sO*co*ci, -sO*so+cO*co*ci, co*s_i],-1)
    r=xp[...,None]*P+yp[...,None]*Q
    v=vxp[...,None]*P+vyp[...,None]*Q
    return r,v

def rv_to_el(mu, r, v):
    rn=np.linalg.norm(r,axis=-1); v2=np.sum(v*v,-1)
    a=1/(2/rn-v2/mu)
    hv=np.cross(r,v); hn=np.linalg.norm(hv,axis=-1)
    inc=np.arccos(np.clip(hv[...,2]/hn,-1,1))
    Om=np.arctan2(hv[...,0],-hv[...,1])
    ev=np.cross(v,hv)/mu[...,None] - r/rn[...,None] if np.ndim(mu) else np.cross(v,hv)/mu - r/rn[...,None]
    e=np.linalg.norm(ev,axis=-1)
    # pomega and mean anomaly via eccentric anomaly
    # node vector
    nv=np.stack([np.cos(Om),np.sin(Om),np.zeros_like(Om)],-1)
    hu=hv/hn[...,None]
    # argument of pericenter: angle from node to e-vector in orbital plane
    om=np.arctan2(np.sum(np.cross(nv,ev)*hu,-1), np.sum(nv*ev,-1))
    pom=Om+om
    rv=np.sum(r*v,-1)
    sE=rv/np.sqrt(mu*a); cE=1-rn/a
    Ea=np.arctan2(sE,cE)
    M=Ea-e*np.sin(Ea)
    lam=M+pom
    si=np.sin(inc/2)
    return np.stack([np.log(a), lam, e*np.cos(pom), e*np.sin(pom), si*np.cos(Om), si*np.sin(Om)],-1)

def dEdv(mu, r, v, rel=1e-6):
    """(...,6,3) partial derivatives of the elements w.r.t. velocity (central differences)"""
    out=[]
    vs=np.linalg.norm(v,axis=-1)[...,None]*rel
    for c in range(3):
        dv=np.zeros_like(v); dv[...,c]=1
        Ep=rv_to_el(mu,r,v+vs*dv); Em=rv_to_el(mu,r,v-vs*dv)
        D=Ep-Em; D[...,1]=(D[...,1]+np.pi)%(2*np.pi)-np.pi
        out.append(D/(2*vs))
    return np.stack(out,-1)

class Theory:
    """First-order short-period theory around given osculating Jacobi elements E0 (8,6)."""
    def __init__(self, masses, E0, Ng=64, gr=True, slow_period=None, inner=None, pair_ng=None):
        """inner: if set, only pairs with i < inner (at least one inner planet) and GR for i < inner"""
        self.m=np.array(masses); m=self.m; Np=len(m)-1; self.Np=Np
        Mint=np.cumsum(m)[:-1]          # interior mass for planet i (index i-1)
        self.mu=G*(Mint+m[1:])
        self.Mint=Mint
        self.E0=E0.copy(); self.Ng=Ng
        self.n=np.sqrt(self.mu/np.exp(3*E0[:,0]))
        self.pairs={}; self.single={}
        lg=2*np.pi*np.arange(Ng)/Ng
        K=np.fft.fftfreq(Ng,1/Ng)
        for i in range(Np):
            for j in range(i+1,Np):
                if inner is not None and i>=inner: continue
                Ngp=Ng if pair_ng is None else pair_ng(i,j)
                if Ngp==0: continue
                lgp=2*np.pi*np.arange(Ngp)/Ngp; Kp=np.fft.fftfreq(Ngp,1/Ngp)
                Ei=np.broadcast_to(E0[i],(Ngp,Ngp,6)).copy(); Ej=np.broadcast_to(E0[j],(Ngp,Ngp,6)).copy()
                Ei[...,1]=lgp[:,None]; Ej[...,1]=lgp[None,:]
                ri,vi=el_to_rv(self.mu[i],Ei); rj,vj=el_to_rv(self.mu[j],Ej)
                d=rj-ri; d3=np.sum(d*d,-1)**1.5; rj2=np.sum(rj*rj,-1); rj3=rj2**1.5
                # on interior i from exterior j
                ai=G*m[j+1]*(d/d3[...,None] - rj/rj3[...,None])
                # on exterior j from interior i
                rdot=np.sum(rj*ri,-1)
                aj=-G*(Mint[j]+m[j+1])/Mint[j]*m[i+1]*(d/d3[...,None] - rj/rj3[...,None] + ri/rj3[...,None]
                                                    - 3*rj*(rdot/rj2/rj3)[...,None])
                fi=np.einsum('...ec,...c->...e',dEdv(self.mu[i],ri,vi),ai)
                fj=np.einsum('...ec,...c->...e',dEdv(self.mu[j],rj,vj),aj)
                Fi=np.fft.fft2(fi,axes=(0,1))/Ngp**2; Fj=np.fft.fft2(fj,axes=(0,1))/Ngp**2
                kk=Kp[:,None]*self.n[i]+Kp[None,:]*self.n[j]
                self.pairs[(i,j)]=(self._sp(Fi,kk,self.n[i],slow_period), self._sp(Fj,kk,self.n[j],slow_period),
                                   Fi[0,0].real, Fj[0,0].real)
        if gr:
            for i in range(Np if inner is None else inner):
                Ei=np.broadcast_to(E0[i],(Ng,6)).copy(); Ei[:,1]=lg
                ri,vi=el_to_rv(self.mu[i],Ei)
                rr=np.linalg.norm(ri,axis=-1)
                agr=(-6*m[0]**2/C_GR**2)*ri/rr[:,None]**4
                fi=np.einsum('...ec,...c->...e',dEdv(self.mu[i],ri,vi),agr)
                Fi=np.fft.fft(fi,axis=0)/Ng
                kk=K*self.n[i]
                self.single[i]=(self._sp(Fi,kk,self.n[i],slow_period), Fi[0].real)
        self.K=K
    def _sp(self, F, kk, n, slow_period):
        """Fourier coefficients of the short-period displacement from those of the rates F[...,6]"""
        with np.errstate(divide='ignore',invalid='ignore'):
            inv=1/(1j*kk)
        zero=(np.abs(kk)<1e-300)
        if slow_period is not None: zero|=(np.abs(kk)<2*np.pi/slow_period)
        inv[zero]=0
        S=F*inv[...,None]
        S[...,1]+= -1.5*n*F[...,0]*inv**2
        return S
    comps=None
    def disp(self, E):
        """short-period displacement (8,6) at osculating elements E (only the lambdas are used: first order)"""
        D=self._disp(E)
        if self.comps is not None: D[:,[c for c in range(6) if c not in self.comps]]=0
        return D
    def _disp(self, E):
        D=np.zeros((self.Np,6),dtype=complex)
        K=self.K
        for (i,j),(Si,Sj,_,_) in self.pairs.items():
            Kp=np.fft.fftfreq(Si.shape[0],1/Si.shape[0])
            ph=np.exp(1j*(Kp[:,None]*E[i,1]+Kp[None,:]*E[j,1]))
            D[i]+=np.einsum('abe,ab->e',Si,ph); D[j]+=np.einsum('abe,ab->e',Sj,ph)
        for i,(Si,_) in self.single.items():
            ph=np.exp(1j*K*E[i,1]); D[i]+=np.einsum('ae,a->e',Si,ph)
        return D.real
    def mean(self, E): return E-self.disp(E)
    def osc(self, M):
        """inverse transform by fixed-point iteration"""
        E=M.copy()
        for _ in range(3): E=M+self.disp(E)
        return E
    def secular_rates(self):
        """first-order secular element rates (8,6) at E0 (lambda rate excludes the Kepler mean motion)"""
        R=np.zeros((self.Np,6))
        for (i,j),(_,_,fi0,fj0) in self.pairs.items(): R[i]+=fi0; R[j]+=fj0
        for i,(_,f0) in self.single.items(): R[i]+=f0
        return R

def secular_rates_fast(masses, E, Ng=24, inner=4, gr=True):
    """first-order secular (averaged) element rates (8,6) at elements E, only pairs with an inner planet;
    rows >= inner are not filled (outer planets are propagated elsewhere)"""
    m=np.asarray(masses); Np=len(m)-1
    Mint=np.cumsum(m)[:-1]; mu=G*(Mint+m[1:])
    R=np.zeros((Np,6))
    lg=2*np.pi*np.arange(Ng)/Ng
    for i in range(min(inner,Np)):
        for j in range(i+1,Np):
            Ei=np.broadcast_to(E[i],(Ng,Ng,6)).copy(); Ej=np.broadcast_to(E[j],(Ng,Ng,6)).copy()
            Ei[...,1]=lg[:,None]; Ej[...,1]=lg[None,:]
            ri,vi=el_to_rv(mu[i],Ei); rj,vj=el_to_rv(mu[j],Ej)
            d=rj-ri; d3=np.sum(d*d,-1)**1.5; rj2=np.sum(rj*rj,-1); rj3=rj2**1.5
            ai=G*m[j+1]*(d/d3[...,None] - rj/rj3[...,None])
            R[i]+=np.einsum('abec,abc->e',dEdv(mu[i],ri,vi),ai)/Ng**2
            if j<inner:
                rdot=np.sum(rj*ri,-1)
                aj=-G*(Mint[j]+m[j+1])/Mint[j]*m[i+1]*(d/d3[...,None] - rj/rj3[...,None] + ri/rj3[...,None]
                                                    - 3*rj*(rdot/rj2/rj3)[...,None])
                R[j]+=np.einsum('abec,abc->e',dEdv(mu[j],rj,vj),aj)/Ng**2
        if gr:
            Ei=np.broadcast_to(E[i],(Ng,6)).copy(); Ei[:,1]=lg
            ri,vi=el_to_rv(mu[i],Ei); rr=np.linalg.norm(ri,axis=-1)
            agr=(-6*m[0]**2/C_GR**2)*ri/rr[:,None]**4
            R[i]+=np.einsum('aec,ac->e',dEdv(mu[i],ri,vi),agr)/Ng
    return R
