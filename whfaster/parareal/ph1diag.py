# per-planet error (max over boundaries, AU) at k=1..KM for one epoch
import numpy as np, sys, mgc, ph1eval, gfilt, g3
from ptlib import *
def run(G, ep, KM=4):
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
        print(f"k={k+1}: "+" ".join(f"{NAMES[p]} {pe[p]:.2g}" for p in range(8)))
CG=mgc.Coarse(360000,Ng=32,dt_out=360,live=1,comps=None,pair_ng=gfilt.pair_rule,heun=1,meann=1)
ep=int(sys.argv[1])
print("all_mn1"); run(lambda E: CG(E,0), ep)
G3=g3.G3(base=lambda E: CG(E,0)); print("all_mn1_g3"); run(G3, ep)
