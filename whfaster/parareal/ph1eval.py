# Round 3 evaluation harness: nonlinear parareal (N = 10 x 986 yr) for a coarse propagator at the 10 Phase-0 epochs;
# prints error/noise at k = 1..KM per epoch, and summary (median, worst, #epochs with k=2 <= 2).
import numpy as np, sys, time
from ptlib import *
D=np.load('phase0_jac.npz'); EP=[0,1030,2000,3400,4500,5500,6000,7000,8160,9300]; N=10; days=360000
def F(E): return to_el(prop(from_el(E),days,6))
_Fcache={}
def evaluate(G, label, KM=3, eps=None, verbose=True):
    eps=eps or EP; R=[]; t0=time.time()
    for ep in eps:
        ref=D[f'ref{ep}']; noise=float(D[f'noise{ep}']); E0=ref[0]
        U=[E0]
        for n in range(N): U.append(G(U[-1]))
        Gold=[None]+U[1:]; row=[]
        for k in range(KM):
            Fv=[]
            for n in range(N):
                key=(ep,U[n].tobytes())
                if key not in _Fcache: _Fcache[key]=F(U[n])
                Fv.append(_Fcache[key])
            Un=[E0]; Gn=[None]
            for n in range(N):
                g=G(Un[n]); Gn.append(g); c=g+ediff(Fv[n],Gold[n+1]); c[:,1]=wrap(c[:,1]); Un.append(c)
            U=Un; Gold=Gn
            row.append(max(poserr(from_el(U[n]),from_el(ref[n])) for n in range(1,N+1))/noise)
        R.append(row)
        if verbose: print(f"   {label:24s} epoch {ep:5d}: "+" ".join(f"{x:9.3g}" for x in row)); sys.stdout.flush()
    R=np.array(R)
    print(f"== {label:24s} median k1..: "+" ".join(f"{x:.3g}" for x in np.median(R,0))+" | worst: "+" ".join(f"{x:.3g}" for x in R.max(0))+f" | k2<=2 at {np.sum(R[:,1]<=2)}/{len(eps)} | {time.time()-t0:.0f}s"); sys.stdout.flush()
    return R
if __name__=='__main__':
    import ptc
    for spec in sys.argv[1:]:
        kw={k:float(v) if '.' in v else int(v) for k,v in (x.split('=') for x in spec.split(',') if x)} if spec!='default' else {}
        G=ptc.CoarseC(days,heun=1,**kw)
        evaluate(G,spec,verbose=False)
