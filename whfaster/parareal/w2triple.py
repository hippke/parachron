# Round 3, iteration 24: second-order short-period terms (W2) and second-order secular rates within a planet triple
# (e.g. Earth-Mars-Jupiter), on the 3-D torus of their mean longitudes. Non-canonical averaging, X = Y + W1 + W2:
#   f2_q = sum_{pairs P containing q} sum_{r in P} [ dV1^P_q/dlambda_r * W1_{r,lambda} + sum_c dV1^P_q/dx_{r,c} W1_{r,c} ]
#          + 1.125 n_q (W1_{q,lna})^2  (lambda component)
# with W1_r the total first-order displacement of r from all pairs of the triple (same-pair and three-body products).
# W2 = homological solution (with Kepler shear) of f2 - <f2> excluding near-resonant harmonics (|k.nu| < nu_cut, kept
# as resonant dynamics elsewhere); U2 = <f2> (second-order secular rates).
import numpy as np, itertools
import so2
from so2 import MU

def _sp(F, kk, n):
    with np.errstate(divide='ignore', invalid='ignore'): inv=1/(1j*kk)
    inv[np.abs(kk)<1e-300]=0
    S=F*inv[...,None]; S[...,1]+=-1.5*n*F[...,0]*inv**2
    return S

class TripleW2:
    def __init__(self, E, P, Ng=32, n=None, h=1e-4, nu_cut=0.0, planets_out=None, cross_only=False):
        """P: 2 planets (same-pair products) or 3 planets (cross_only=True: only products of different pairs, so that
        triples and pairs can be combined without double counting)"""
        P=tuple(sorted(P)); self.P=P; self.Ng=Ng; D=len(P)
        n=np.sqrt(MU/np.exp(3*E[:,0])) if n is None else np.asarray(n); self.n=n
        K=np.fft.fftfreq(Ng,1/Ng); self.K=K
        ax={p:a for a,p in enumerate(P)}
        def emb(arr2, a, b):
            """(Ng,Ng,...) on axes (ax a, ax b) -> (Ng,)*D + (...)"""
            if D==2: return arr2
            x=np.expand_dims(arr2,[q for q in range(3) if q not in (ax[a],ax[b])][0])
            return np.broadcast_to(x,(Ng,Ng,Ng)+arr2.shape[2:])
        pairs=list(itertools.combinations(P,2))
        # first-order rates, displacements and derivatives per pair
        V={}; W1={p:np.zeros((Ng,)*D+(6,)) for p in P}; W1p={}; dV={}
        for (a,b) in pairs:
            fa,fb=so2.pair_rates(E,a,b,Ng)
            kk2=K[:,None]*n[a]+K[None,:]*n[b]
            for q,f in ((a,fa),(b,fb)):
                F=np.fft.fft2(f,axes=(0,1))/Ng**2
                V[(a,b,q)]=F
                W1p[(a,b,q)]=emb(np.real(np.fft.ifft2(_sp(F,kk2,n[q])*Ng**2,axes=(0,1))),a,b)
                W1[q]=W1[q]+W1p[(a,b,q)]
                # lambda derivatives
                dV[(a,b,q,a,1)]=emb(np.real(np.fft.ifft2(1j*K[:,None,None]*F*Ng**2,axes=(0,1))),a,b)
                dV[(a,b,q,b,1)]=emb(np.real(np.fft.ifft2(1j*K[None,:,None]*F*Ng**2,axes=(0,1))),a,b)
            for r in (a,b):
                for c in (0,2,3,4,5):
                    Ep=E.copy(); Em=E.copy(); Ep[r,c]+=h; Em[r,c]-=h
                    pa,pb=so2.pair_rates(Ep,a,b,Ng); ma,mb=so2.pair_rates(Em,a,b,Ng)
                    dV[(a,b,a,r,c)]=emb((pa-ma)/(2*h),a,b); dV[(a,b,b,r,c)]=emb((pb-mb)/(2*h),a,b)
        self.W1grid=W1
        kk3=(K[:,None,None]*n[P[0]]+K[None,:,None]*n[P[1]]+K[None,None,:]*n[P[2]]) if D==3 else K[:,None]*n[P[0]]+K[None,:]*n[P[1]]
        self.S2={}; self.U2={}
        for q in (planets_out if planets_out is not None else [p for p in P if p<4]):
            f2=np.zeros((Ng,)*D+(6,))
            for (a,b) in pairs:
                if q not in (a,b): continue
                for r in (a,b):
                    if r>=4: continue      # giants' displacements: already in G (osculating giants from their WH run)
                    # displacement of r: all pairs of the group, or (cross_only) only pairs other than (a,b)
                    Wr=W1[r] if not cross_only else W1[r]-W1p[(a,b,r)]
                    for c in range(6):
                        f2=f2+dV[(a,b,q,r,c)]*Wr[...,c:c+1]
            if cross_only:   # 1/2 D^2 V0: cross products of different pairs only
                tot=W1[q][...,0]**2-sum(W1p[(a,b,q)][...,0]**2 for (a,b) in pairs if q in (a,b))
                f2[...,1]+=1.125*n[q]*tot
            else:
                f2[...,1]+=1.125*n[q]*W1[q][...,0]**2
            F2=np.fft.fftn(f2,axes=tuple(range(D)))/Ng**D
            self.U2[q]=np.real(F2[(0,)*D])
            S=_sp(F2,kk3,n[q]); S[(0,)*D]=0
            if nu_cut>0: S[np.abs(kk3)<nu_cut]=0
            # keep only genuinely second-order harmonics involving at least two of the triple's angles is not
            # required: same-pair products are included consistently.
            self.S2[q]=S
    def disp(self, E):
        """W2 (6 comps) of the inner planets of the triple at the angles of E: dict planet -> array(6)"""
        K=self.K; P=self.P
        if len(P)==3:
            ph=np.exp(1j*(K[:,None,None]*E[P[0],1]+K[None,:,None]*E[P[1],1]+K[None,None,:]*E[P[2],1]))
            return {q:np.real(np.einsum('abce,abc->e',S,ph)) for q,S in self.S2.items()}
        ph=np.exp(1j*(K[:,None]*E[P[0],1]+K[None,:]*E[P[1],1]))
        return {q:np.real(np.einsum('abe,ab->e',S,ph)) for q,S in self.S2.items()}

class W2Set:
    """all same-pair groups (pairs with an inner planet) and cross-only triples among planets PL"""
    def __init__(self, E, n=None, PL=(0,1,2,3,4,5), Ng=32, nu_cut=1/500, triples=None, pairs=None):
        pairs=pairs if pairs is not None else [P for P in itertools.combinations(PL,2) if min(P)<4]
        triples=triples if triples is not None else list(itertools.combinations(PL,3))
        self.G=[TripleW2(E,P,Ng=Ng,n=n,nu_cut=nu_cut) for P in pairs]+[TripleW2(E,P,Ng=Ng,n=n,nu_cut=nu_cut,cross_only=True) for P in triples]
    def disp(self, E):
        D=np.zeros((8,6))
        for g in self.G:
            for q,w in g.disp(E).items(): D[q]+=w
        return D
    def U2(self):
        R=np.zeros((8,6))
        for g in self.G:
            for q,u in g.U2.items(): R[q]+=u
        return R
