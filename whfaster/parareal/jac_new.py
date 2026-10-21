# Round 3: slice Jacobians J_G and coarse errors e0 of a Python coarse spec at given epochs (Phase-0 references),
# for the block attribution (jacblock-style).   python3 jac_new.py spec ep [ep...]  -> jacG_<tag>_<ep>.npz
import numpy as np, sys, ph1mn
from ptlib import *
D=np.load('phase0_jac.npz'); N=10
H=np.array([1e-7,1e-4,1e-5,1e-5,1e-5,1e-5])
spec=sys.argv[1]; G=ph1mn.build(spec)
for ep in map(int,sys.argv[2:]):
    ref=D[f'ref{ep}']
    U=[ref[0]]
    for n in range(N): U.append(G(U[-1]))
    e0=np.array([ediff(U[n],ref[n]).ravel() for n in range(N+1)])
    JG=np.zeros((N,48,48))
    for n in range(N):
        for i in range(8):
            for c in range(6):
                Ep=ref[n].copy(); Em=ref[n].copy(); Ep[i,c]+=H[c]; Em[i,c]-=H[c]
                JG[n][:,i*6+c]=ediff(G(Ep),G(Em)).ravel()/(2*H[c])
    np.savez(f'jacG_{ep}.npz',JG=JG,e0=e0,spec=spec); print('done',ep,flush=True)
