# Round 3, it. 23: which mean-longitude combination (<= 3 planets) explains the remaining Mars (k,h) error line?
import numpy as np, itertools
from ptlib import *
d=np.fromfile('../pint/runs/m10_ser2.bin',dtype=np.float64).reshape(-1,50); X=d[:,2:].reshape(-1,8,6)
err=np.load('mars_err_1000.npy'); s0=1000; M=200
L=X[s0:s0+M,:,1]
z=err[:,3,2]+1j*err[:,3,3]; zz=np.sum(abs(z)**2)
K=np.arange(-10,11); out=[]
for tri in itertools.combinations(range(7),3):
    k1,k2,k3=np.meshgrid(K,K,K,indexing='ij'); k1=k1.ravel(); k2=k2.ravel(); k3=k3.ravel()
    sel=(np.abs(k1)+np.abs(k2)+np.abs(k3)<=16)&np.isin(k1+k2+k3,[-1,0,1,2])&((k1!=0)|(k2!=0)|(k3!=0))
    k1,k2,k3=k1[sel],k2[sel],k3[sel]
    ph=np.exp(1j*(np.outer(k1,L[:,tri[0]])+np.outer(k2,L[:,tri[1]])+np.outer(k3,L[:,tri[2]])))
    A=ph.conj()@z/M
    r2=1-np.sum(np.abs(z[None,:]-A[:,None]*ph)**2,1)/zz
    for i in np.argsort(-r2)[:3]: out.append((r2[i],abs(A[i]),{tri[0]:k1[i],tri[1]:k2[i],tri[2]:k3[i]}))
out.sort(key=lambda q:-q[0]); seen=set()
for r2,a,ks in out[:12]:
    key=tuple(sorted((p,k) for p,k in ks.items() if k))
    if key in seen: continue
    seen.add(key); print(f"R2 {r2:.2f} |A| {a:.1e}  "+" ".join(f"{k:+d}{NAMES[p]}" for p,k in key))
