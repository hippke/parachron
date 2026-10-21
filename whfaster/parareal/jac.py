# Iteration 1 (round 2): slice-Jacobian test.
# Measures the 1-slice Jacobians of the fine map F (WHFast512, 6 d) and of coarse maps G in Jacobi elements
# (ln a, lambda, k, h, q, p per planet) by central finite differences, at every slice start along the serial fine
# run, and predicts parareal convergence with the linearised error recursion
#     e^{k+1}_{n+1} = G'_n e^{k+1}_n + (F'_n - G'_n) e^k_n ,
# starting from the true coarse errors e^0. Compares with the actual (nonlinear) parareal in elements.
import numpy as np, sys, time, argparse, pickle
from ptlib import *
p=argparse.ArgumentParser()
p.add_argument('--slices',type=int,default=10)
p.add_argument('--slice_days',type=float,default=360000)
p.add_argument('--dtG',type=str,default='8,12,16,24')
p.add_argument('--kmax',type=int,default=6)
p.add_argument('--out',default='jac.pkl')
a=p.parse_args()
N=a.slices; T=a.slice_days
H=np.array([1e-8,1e-4,1e-5,1e-5,1e-5,1e-5])
def mapE(E, dt): return to_el(prop(from_el(E), T, dt))
def jac(E, dt):
    J=np.zeros((48,48))
    for i in range(8):
        for c in range(6):
            Ep=E.copy(); Em=E.copy(); Ep[i,c]+=H[c]; Em[i,c]-=H[c]
            J[:,i*6+c]=ediff(mapE(Ep,dt),mapE(Em,dt)).ravel()/(2*H[c])
    return J
U0=initial()
ref=[to_el(U0)]
for n in range(N): ref.append(to_el(prop(from_el(ref[-1]),T,6)))
res={'ref':ref}
t0=time.time()
JF=[jac(ref[n],6) for n in range(N)]
print(f"JF done {time.time()-t0:.0f}s"); sys.stdout.flush()
res['JF']=JF
def rep(tag, errs):
    # errs: list over slices of (8,6) element errors; report max |dln a| and |dlambda| per planet
    A=np.max([np.abs(e[:,0]) for e in errs],axis=0); L=np.max([np.abs(e[:,1]) for e in errs],axis=0)
    print(f"  {tag:12s} lna: "+" ".join(f"{x:.0e}" for x in A)+" | lam: "+" ".join(f"{x:.0e}" for x in L))
for dtG in [float(x) for x in a.dtG.split(',')]:
    JG=[jac(ref[n],dtG) for n in range(N)]
    res[('JG',dtG)]=JG
    # block norms of mismatch, slice 0
    D=JF[0]-JG[0]
    print(f"\n=== dtG={dtG}: mismatch |F'-G'| (slice 0), rows=lna of planet, max over columns of each planet's (lna,lam,kh,qp)")
    for i in range(8):
        r=D[i*6]
        print(f"  {NAMES[i]} lna-row: "+" ".join(f"{np.max(np.abs(r[j*6:j*6+6])):.0e}" for j in range(8)))
    # actual coarse trajectory errors and actual parareal (element-space correction)
    U=[ref[0]]
    for n in range(N): U.append(mapE(U[-1],dtG))
    Gold=[None]+U[1:]
    e0=[ediff(U[n],ref[n]) for n in range(N+1)]
    # linear prediction
    e=e0
    print(" linear prediction:")
    rep("k=0", e0[1:])
    for k in range(1,a.kmax+1):
        en=[np.zeros((8,6))]
        for n in range(N):
            en.append((JG[n]@en[n].ravel() + (JF[n]-JG[n])@e[n].ravel()).reshape(8,6))
        e=en; rep(f"k={k}", e[1:])
    print(" actual parareal (elements):")
    rep("k=0", e0[1:])
    for k in range(1,a.kmax+1):
        Fv=[mapE(U[n],6) for n in range(N)]
        Un=[ref[0]]; Gn=[None]
        for n in range(N):
            g=mapE(Un[n],dtG); Gn.append(g)
            c=g+ediff(Fv[n],Gold[n+1]); c[:,1]=wrap(c[:,1]); Un.append(c)
        U=Un; Gold=Gn
        rep(f"k={k}", [ediff(U[n],ref[n]) for n in range(1,N+1)])
    sys.stdout.flush()
pickle.dump(res,open(a.out,'wb'))
