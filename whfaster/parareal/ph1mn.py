# Round 3 (iterations 21-23): Python coarse variants built from a spec string; parareal evaluation (ph1eval),
# per-slice quality (spec starting with q_) or per-planet k=2 diagnosis (diag=1).
# spec tokens: all (all-component transform) | mn1 (tables at mean a) | dtX (giants' WH step) | nsubX (RK4 secular
# drift with giants sampled) | un (Uranus/Neptune in secular rates) | ngsX (secular grid) | LX (filter, yr) |
# noend | g3 (Earth-Mars-Jupiter 2nd-order term) | gm, ns (rejected variants)
import sys, numpy as np, mgc, ph1eval, gfilt
from ptlib import *
import itertools
W2SETS={'a':([P for P in itertools.combinations(range(6),2) if min(P)<4],[(1,2,4),(2,3,4),(1,2,3)]),
        'p':([P for P in itertools.combinations(range(6),2) if min(P)<4],[])}
def build(spec):
    comps=None if 'all' in spec else [0]
    mn=int(spec.split('mn')[1][0]) if 'mn' in spec else 0
    dto=float(spec.split('dt')[1].split('_')[0]) if '_dt' in spec else 360
    CG=mgc.Coarse(360000,Ng=32,dt_out=dto,live=1,comps=comps,pair_ng=gfilt.pair_rule,heun=1,meann=mn,giant_mean=int('gm' in spec),nsec=int('ns' in spec and 'nsub' not in spec),nsub=int(spec.split('nsub')[1].split('_')[0]) if 'nsub' in spec else 0,Ng_sec=int(spec.split('ngs')[1].split('_')[0]) if 'ngs' in spec else 0,sec_ng=gfilt.sec_rule if 'un' in spec else None,gavg=float(spec.split('ga')[1].split('_')[0]) if '_ga' in spec else 0.,w2=W2SETS.get(spec.split('w2')[1].split('_')[0]) if '_w2' in spec else None)
    G=lambda E: CG(E,0)
    if '_L' in spec:
        L=float(spec.split('_L')[1].split('_')[0]); gf=gfilt.GF2(CG,L,end=int('noend' not in spec))
    if 'g3' in spec:
        import g3; G=g3.G3(base=lambda E: CG(E,0),ode=int(spec.split('ode')[1].split('_')[0]) if 'ode' in spec else 0)
    return G
def calibrated(G, ep, nfit=100):
    """secular-matrix calibration (it. 23): per-slice z = k+ih error of the inner planets modelled as A z_inner,
    A (4x4 complex) fitted from the defects of the nfit slices preceding the window (cal_<ep>.npy)"""
    import seccal
    err=np.load(f'cal_{ep}.npy')[-nfit:]; s0=max(ep-100,0)+100-len(err)
    A=seccal.fit(err,s0,2,cols=4)
    def Gc(E):
        g=G(E); z=E[:4,2]+1j*E[:4,3]; dz=A@z
        g[:4,2]+=dz.real; g[:4,3]+=dz.imag; return g
    return Gc
def diag(G, ep, KM=2):
    D=ph1eval.D; ref=D[f'ref{ep}']; noise=float(D[f'noise{ep}']); E0=ref[0]; N=10
    U=[E0]
    for n in range(N): U.append(G(U[-1]))
    Gold=[None]+U[1:]
    for k in range(KM):
        Fv=[ph1eval.F(U[n]) for n in range(N)]
        Un=[E0]; Gn=[None]
        for n in range(N):
            g=G(Un[n]); Gn.append(g); c=g+ediff(Fv[n],Gold[n+1]); c[:,1]=wrap(c[:,1]); Un.append(c)
        U=Un; Gold=Gn
    pe=np.max([np.max(np.abs(from_el(U[n])[1:,:3]-from_el(ref[n])[1:,:3]),axis=1) for n in range(1,N+1)],axis=0)/noise
    print(f"epoch {ep:5d} k={KM}: "+" ".join(f"{NAMES[p]} {pe[p]:6.2g}" for p in range(8)),flush=True)
if __name__=='__main__':
    EPS=None; DIAG=0
    for spec in sys.argv[1:]:
        if spec.startswith('eps='): EPS=[int(x) for x in spec[4:].split(',')]; continue
        if spec.startswith('diag='): DIAG=int(spec[5:]); continue
        G=build(spec)
        if spec.startswith('q_'):
            import slice_quality; slice_quality.quality(G,1000,40,label=spec)
        elif DIAG:
            print(spec)
            for ep in (EPS or ph1eval.EP): diag(G,ep)
        elif '_cal' in spec:
            R=np.array([ph1eval.evaluate(calibrated(G,ep),spec,eps=[ep],verbose=True)[0] for ep in (EPS or ph1eval.EP)])
            print(f"== {spec} median k1..: "+" ".join(f"{x:.3g}" for x in np.median(R,0))+" | worst: "+" ".join(f"{x:.3g}" for x in R.max(0))+f" | k2<=2 at {np.sum(R[:,1]<=2)}/{len(R)}",flush=True)
        else: ph1eval.evaluate(G, spec, eps=EPS, verbose=True)
