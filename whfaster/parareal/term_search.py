# Round 3, it. 23: harmonic search of per-slice coarse errors (file from slice runs) for one planet/component.
#   python3 term_search.py err.npy s0 planet comp(lna|lam|z) [sums]
import numpy as np, itertools, sys
from ptlib import *
d=np.fromfile('../pint/runs/m10_ser2.bin',dtype=np.float64).reshape(-1,50); X=d[:,2:].reshape(-1,8,6)
err=np.load(sys.argv[1]); s0=int(sys.argv[2]); P=int(sys.argv[3]); comp=sys.argv[4]
sums=[int(x) for x in sys.argv[5].split(',')] if len(sys.argv)>5 else [-2,-1,0,1,2]
M=len(err); L=X[s0:s0+M,:,1]
y={'lna':err[:,P,0],'lam':err[:,P,1],'z':err[:,P,2]+1j*err[:,P,3]}[comp].astype(complex)
cplx=comp=='z'; yy=np.sum(abs(y)**2)
K=np.arange(-10,11); out=[]
for tri in itertools.combinations(range(7),3):
    k1,k2,k3=[a.ravel() for a in np.meshgrid(K,K,K,indexing='ij')]
    sel=(np.abs(k1)+np.abs(k2)+np.abs(k3)<=18)&np.isin(k1+k2+k3,sums)&((k1!=0)|(k2!=0)|(k3!=0))
    if not cplx: sel&=(k1>0)|((k1==0)&(k2>0))|((k1==0)&(k2==0)&(k3>0))
    k1,k2,k3=k1[sel],k2[sel],k3[sel]
    ph=np.exp(1j*(np.outer(k1,L[:,tri[0]])+np.outer(k2,L[:,tri[1]])+np.outer(k3,L[:,tri[2]])))
    if cplx:
        A=ph.conj()@y/M; r2=1-np.sum(np.abs(y[None,:]-A[:,None]*ph)**2,1)/yy
    else:
        c=np.real(ph); s=np.imag(ph); yr=y.real-y.real.mean()
        # 2-parameter LS per row
        cc=np.sum(c*c,1); ss=np.sum(s*s,1); cs=np.sum(c*s,1); cy=c@yr; sy=s@yr
        det=cc*ss-cs**2; a=(ss*cy-cs*sy)/det; b=(cc*sy-cs*cy)/det
        r2=(a*cy+b*sy)/np.sum(yr**2); A=np.hypot(a,b)
    for i in np.argsort(-r2)[:2]: out.append((r2[i],abs(A[i]),tuple((tri[j],kk[i]) for j,kk in enumerate((k1,k2,k3)) if kk[i])))
out.sort(key=lambda q:-q[0]); seen=set(); n=0
print(f"{NAMES[P]} {comp}: rms {np.sqrt(yy/M):.2e}")
for r2,a,key in out:
    if key in seen: continue
    seen.add(key); n+=1; print(f"  R2 {r2:.2f} amp {a:.1e}  "+" ".join(f"{k:+d}{NAMES[p]}" for p,k in key))
    if n>=8: break
