# Round 3: block attribution for the new coarse map (J_G, e0 from jac_new.py; J_F from Phase 0).
#   python3 jacblock2.py ep...
import numpy as np, sys
from ptlib import *
D=np.load('phase0_jac.npz'); N=10
EL={'lna':[0],'lam':[1],'kh':[2,3],'qp':[4,5]}
def idx(planets,els): return [p*6+c for p in planets for e in els for c in EL[e]]
def run(ep, rows=(), cols=(), e0rows=(), kmax=3):
    JF=D[f'JF{ep}']; d=np.load(f'jacG_{ep}.npz'); JG=d['JG'].copy(); e0=d['e0'].copy(); ref=D[f'ref{ep}']; noise=float(D[f'noise{ep}'])
    rows=list(rows); cols=list(cols)
    if rows and cols:
        for n in range(N): JG[n][np.ix_(rows,cols)]=JF[n][np.ix_(rows,cols)]
    e0[:,list(e0rows)]=0
    e=list(e0); res=[]
    for k in range(kmax):
        en=[np.zeros(48)]
        for n in range(N): en.append(JG[n]@en[n]+(JF[n]-JG[n])@e[n])
        e=en; res.append(max(poserr(from_el(ref[n]+e[n].reshape(8,6)),from_el(ref[n])) for n in range(1,N+1))/noise)
    return res
allc=list(range(48))
for ep in map(int,sys.argv[1:]):
    print(f"epoch {ep}: linear model k1..3 {np.round(run(ep),2)}")
    tests=[]
    for p in range(8):
        tests.append((f'{NAMES[p]} rows J+e0 exact', dict(rows=idx([p],EL),cols=allc,e0rows=idx([p],EL))))
    for p in (0,1,2,3):
        for r in EL:
            tests.append((f'{NAMES[p]} {r} rows J+e0', dict(rows=idx([p],[r]),cols=allc,e0rows=idx([p],[r]))))
        tests.append((f'{NAMES[p]} e0 only', dict(e0rows=idx([p],EL))))
        tests.append((f'{NAMES[p]} J only', dict(rows=idx([p],EL),cols=allc)))
    for lab,kw in tests:
        print(f"   {lab:28s} k2 {run(ep,**kw)[1]:8.3g}")
