# Parareal driver for mgc.Coarse variants (nonlinear iteration, osculating-element correction).
import numpy as np, time, argparse, sys
from ptlib import *
import mgc
p=argparse.ArgumentParser()
p.add_argument('--slices',type=int,default=10)
p.add_argument('--slice_days',type=float,default=360000)
p.add_argument('--kmax',type=int,default=6)
p.add_argument('--Ng',type=int,default=48)
p.add_argument('--dt_out',type=float,default=180)
p.add_argument('--live',type=int,default=0)
p.add_argument('--sec',type=int,default=1)
p.add_argument('--live_sec',type=int,default=0)
p.add_argument('--nsub',type=int,default=0)
p.add_argument('--comps',default='')
p.add_argument('--png',default='', help='pair grid rule: alpha-based tol, e.g. 1e-3; or maxj (drop pairs with j>maxj)')
p.add_argument('--maxj',type=int,default=7)
p.add_argument('--Ng_sec',type=int,default=0)
p.add_argument('--start_days',type=float,default=0, help='start the window after this many days of fine integration')
a,_=p.parse_known_args()
N=a.slices; Td=a.slice_days
E0tmp=to_el(initial()); A=np.exp(E0tmp[:,0]); EC=np.hypot(E0tmp[:,2],E0tmp[:,3])
def png(i,j):
    if j>a.maxj: return 0
    if not a.png: return a.Ng
    if a.png.startswith('io'):   # io<Nout>: inner-inner at Ng, inner-outer at Nout
        return a.Ng if j<4 else int(a.png[2:])
    tol=float(a.png); al=A[i]/A[j]
    K=np.ceil(np.log(tol)/np.log(al))+np.ceil(np.log(tol)/np.log(max(EC[i],EC[j],0.02)))
    n=int(min(a.Ng,max(8,2*K+2))); return n+(n%2)
CG=mgc.Coarse(Td,Ng=a.Ng,dt_out=a.dt_out,live=a.live,sec=a.sec,live_sec=a.live_sec,Ng_sec=a.Ng_sec or None,nsub=a.nsub,comps=[int(x) for x in a.comps.split(',')] if a.comps else None,pair_ng=png)
def F(E): return to_el(prop(from_el(E),Td,6))
def rep(tag, errs, extra=""):
    A=np.max([np.abs(e[:,0]) for e in errs],axis=0); L=np.max([np.abs(e[:,1]) for e in errs],axis=0)
    print(f"  {tag:5s} lna: "+" ".join(f"{x:.0e}" for x in A)+" | lam: "+" ".join(f"{x:.0e}" for x in L)+extra); sys.stdout.flush()
U0=initial()
if a.start_days: U0=prop(U0,a.start_days,6)
E0=to_el(U0); ref=[E0]; ref1=[E0]
t0=time.time()
for n in range(N): ref.append(F(ref[-1]))
tF=(time.time()-t0)/N
for n in range(N): ref1.append(to_el(prop(from_el(ref1[-1]),Td,6,kernel=1)))
noise=[ediff(ref1[n],ref[n]) for n in range(1,N+1)]
pnoise=max(poserr(from_el(ref1[n]),from_el(ref[n])) for n in range(1,N+1))
rep("noise",noise,f"  | pos {pnoise:.1e} AU")
nl=np.max([np.abs(e[:,1]) for e in noise],axis=0); na=np.max([np.abs(e[:,0]) for e in noise],axis=0)
t0=time.time(); U=[E0]
for n in range(N): U.append(CG(U[-1],n))
tG=(time.time()-t0)/N
print(f"slices {N} x {Td/365.25:.0f} yr; python F/slice {tF:.3f}s, G/slice {tG:.3f}s"); sys.stdout.flush()
Gold=[None]+U[1:]
rep("k=0",[ediff(U[n],ref[n]) for n in range(1,N+1)])
kconv=None
for k in range(1,a.kmax+1):
    Fv=[F(U[n]) for n in range(N)]
    Un=[E0]; Gn=[None]
    for n in range(N):
        g=CG(Un[n],n); Gn.append(g)
        c=g+ediff(Fv[n],Gold[n+1]); c[:,1]=wrap(c[:,1]); Un.append(c)
    U=Un; Gold=Gn
    e=[ediff(U[n],ref[n]) for n in range(1,N+1)]
    L=np.max([np.abs(x[:,1]) for x in e],axis=0); A=np.max([np.abs(x[:,0]) for x in e],axis=0)
    pe=max(poserr(from_el(U[n]),from_el(ref[n])) for n in range(1,N+1))
    ok=pe<=3*pnoise
    rep(f"k={k}",e,f"  | pos {pe:.1e} AU"+("  CONVERGED (pos <= 3x noise)" if ok else ""))
    if ok and kconv is None: kconv=k; break
print("k_conv =",kconv)
