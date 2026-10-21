# Round 3, Phase 3: second-order (in the masses) resonant terms from products of first-order pair fields.
# Non-canonical averaging: X = Y + W1(Y) + ..., X' = V0(X) + V1(X); second-order averaged field
#   U2 = < D V1 . W1 - D W1 . V1bar + 1/2 D^2 V0 [W1, W1] >_K   (only the harmonic K is computed here)
# State per planet: [ln a, lambda, k, h, q, p]; V1 = sum over pairs of the first-order Jacobi pair rates (meanel.py);
# W1 = first-order short-period displacement (meanel conventions: osc = mean + W1).
# A real field is sum_k c_k exp(i k.theta); the resonant term is 2 Re(U_K exp(i K.theta)).
import numpy as np
import meanel
from ptlib import mass
m=np.array(mass); Np=8
Mint=np.cumsum(m)[:-1]; MU=Mint+m[1:]

def pair_rates(E, i, j, Ng):
    """rates (Ng,Ng,6) for planets i and j of pair (i<j) on the (lambda_i, lambda_j) grid"""
    lg=2*np.pi*np.arange(Ng)/Ng
    Ei=np.broadcast_to(E[i],(Ng,Ng,6)).copy(); Ej=np.broadcast_to(E[j],(Ng,Ng,6)).copy()
    Ei[...,1]=lg[:,None]; Ej[...,1]=lg[None,:]
    ri,vi=meanel.el_to_rv(MU[i],Ei); rj,vj=meanel.el_to_rv(MU[j],Ej)
    d=rj-ri; d3=np.sum(d*d,-1)**1.5; rj2=np.sum(rj*rj,-1); rj3=rj2**1.5
    ai=m[j+1]*(d/d3[...,None]-rj/rj3[...,None])
    rdot=np.sum(rj*ri,-1)
    aj=-(Mint[j]+m[j+1])/Mint[j]*m[i+1]*(d/d3[...,None]-rj/rj3[...,None]+ri/rj3[...,None]-3*rj*(rdot/rj2/rj3)[...,None])
    fi=np.einsum('...ec,...c->...e',meanel.dEdv(MU[i],ri,vi),ai)
    fj=np.einsum('...ec,...c->...e',meanel.dEdv(MU[j],rj,vj),aj)
    return fi,fj

class Pair:
    """Fourier data of pair (i,j): rate coefficients F[p] (Ng,Ng,6), displacement coefficients W[p],
    slow derivatives dF[p][(q,v)] for q in (i,j), v in (0,2,3,4,5)"""
    def __init__(self, E, i, j, Ng=32, dslow=True, h=1e-4):
        self.i,self.j,self.Ng=i,j,Ng
        n=np.sqrt(MU/np.exp(3*E[:,0])); self.n=n
        K=np.fft.fftfreq(Ng,1/Ng); self.K=K
        fi,fj=pair_rates(E,i,j,Ng)
        self.F={i:np.fft.fft2(fi,axes=(0,1))/Ng**2, j:np.fft.fft2(fj,axes=(0,1))/Ng**2}
        kk=K[:,None]*n[i]+K[None,:]*n[j]
        with np.errstate(divide='ignore',invalid='ignore'): inv=1/(1j*kk)
        inv[np.abs(kk)<1e-300]=0
        self.W={}
        for p in (i,j):
            S=self.F[p]*inv[...,None]
            S[...,1]+=-1.5*n[p]*self.F[p][...,0]*inv**2
            self.W[p]=S
        self.dF={}
        if dslow:
            for q in (i,j):
                for v in (0,2,3,4,5):
                    Ep=E.copy(); Em=E.copy(); Ep[q,v]+=h; Em[q,v]-=h
                    a=pair_rates(Ep,i,j,Ng); b=pair_rates(Em,i,j,Ng)
                    for s,p in enumerate((i,j)):
                        self.dF[(p,q,v)]=np.fft.fft2((a[s]-b[s])/(2*h),axes=(0,1))/Ng**2
    def idx(self, ki, kj):
        """array indices of harmonic (ki, kj) or None"""
        N=self.Ng
        if abs(ki)>=N//2 or abs(kj)>=N//2: return None
        return (ki%N, kj%N)
    def coef(self, arr, kv):
        """coefficient of harmonic kv (dict planet->k) in a (Ng,Ng,...) array"""
        ix=self.idx(kv.get(self.i,0),kv.get(self.j,0))
        return None if ix is None else arr[ix]
    def dlam(self, p, q, kv):
        """d F[p] / d lambda_q coefficient at harmonic kv"""
        c=self.coef(self.F[p],kv); return None if c is None else 1j*kv.get(q,0)*c

def harmonics_split(K, a, b, q, kmax):
    """all splittings K = ka + kb with ka supported on pair a, kb on pair b, shared planet q"""
    out=[]
    for s in range(-kmax,kmax+1):
        ka={}; kb={}
        for p,v in K.items():
            if p==q: continue
            if p in a: ka[p]=v
            elif p in b: kb[p]=v
        ka[q]=s; kb[q]=K.get(q,0)-s
        out.append((ka,kb))
    return out

def resonant_U2(E, K, Ng=32, pairs=None):
    """second-order coefficient U_K[p] (6 complex comps) for the planets of K (dict planet->integer)"""
    planets=sorted(K)
    if pairs is None:
        pairs={}
        for x in range(len(planets)):
            for y in range(x+1,len(planets)):
                i,j=planets[x],planets[y]; pairs[(i,j)]=Pair(E,i,j,Ng)
    U={p:np.zeros(6,dtype=complex) for p in planets}
    n=pairs[next(iter(pairs))].n
    for a,A in pairs.items():
        for b,B in pairs.items():
            if a==b: continue
            shared=set(a)&set(b)
            if len(shared)!=1 or set(a)|set(b)!=set(planets): continue
            q=shared.pop()
            for ka,kb in harmonics_split(K,a,b,q,Ng//2-1):
                Wq=B.coef(B.W[q],kb)                 # displacement of the shared planet (6 comps) at kb
                if Wq is None or not np.any(Wq): continue
                for p in a:                          # rate of planet p from pair a, differentiated w.r.t. planet q
                    dl=A.dlam(p,q,ka)
                    if dl is None: continue
                    tot=dl*Wq[1]
                    for v in (0,2,3,4,5):
                        c=A.coef(A.dF[(p,q,v)],ka); tot=tot+c*Wq[v]
                    U[p]+=tot
    # 1/2 D^2 V0 [W1, W1]: lambda component of the shared planets: 1/2 * 2.25 n_q * (W1_lna,q)^2, cross terms of pairs
    for a,A in pairs.items():
        for b,B in pairs.items():
            if a==b: continue
            shared=set(a)&set(b)
            if len(shared)!=1 or set(a)|set(b)!=set(planets): continue
            q=shared.pop()
            for ka,kb in harmonics_split(K,a,b,q,Ng//2-1):
                wa=A.coef(A.W[q],ka); wb=B.coef(B.W[q],kb)
                if wa is None or wb is None: continue
                U[q][1]+=0.5*2.25*n[q]*wa[0]*wb[0]
    return U
