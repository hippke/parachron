# Which Jacobian blocks must the coarse model get right for k = 2?  (linear model, config A Jacobians)
import numpy as np, pickle, sys, ptc
from ptlib import *
d=pickle.load(open(sys.argv[1] if len(sys.argv)>1 else 'jac2_cA.pkl','rb')); JF=d['JF']; JG=d['JG']; ref=d['ref']; N=len(JF)
G=ptc.CoarseC(360000)
U=[ref[0]]
for n in range(N): U.append(G(U[-1]))
e0=[ediff(U[n],ref[n]).ravel() for n in range(N+1)]
A=np.exp(ref[0][:,0])
def run(JGm, kmax=3):
    e=e0; out=[]
    for k in range(kmax):
        en=[np.zeros(48)]
        for n in range(N): en.append(JGm[n]@en[n]+(JF[n]-JGm[n])@e[n])
        e=en; out.append(max(np.max(np.abs(x[1::6])*A) for x in e))   # ~ position error from lambda
    return out
def idx(planets, comps): return [p*6+c for p in planets for c in comps]
def fix(rows, cols):
    J2=[]
    for n in range(N):
        M=JG[n].copy(); M[np.ix_(rows,cols)]=JF[n][np.ix_(rows,cols)]; J2.append(M)
    return J2
inner=[0,1,2,3]; allp=list(range(8))
cases={'config A':(None,None),
 'lna rows <- lambda cols':(idx(inner,[0]),idx(allp,[1])),
 'lambda rows <- lambda cols':(idx(inner,[1]),idx(allp,[1])),
 'lna,lam rows <- lambda cols':(idx(inner,[0,1]),idx(allp,[1])),
 'lna,lam rows <- k,h cols':(idx(inner,[0,1]),idx(allp,[2,3])),
 'lna,lam rows <- lna cols':(idx(inner,[0,1]),idx(allp,[0])),
 'lna,lam rows <- all cols':(idx(inner,[0,1]),idx(allp,range(6))),
 'k,h rows <- all cols':(idx(inner,[2,3]),idx(allp,range(6))),
 'inner rows <- all cols':(idx(inner,range(6)),idx(allp,range(6))),
 'Mars rows only':(idx([3],range(6)),idx(allp,range(6))),
 'Venus+Earth rows':(idx([1,2],range(6)),idx(allp,range(6))),
}
for name,(r,c) in cases.items():
    J2=JG if r is None else fix(r,c)
    print(f"{name:30s} k=1..3 max a*dlambda [AU]: "+" ".join(f"{x:.1e}" for x in run(J2)))
