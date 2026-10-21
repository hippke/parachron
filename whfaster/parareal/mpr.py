# Mean-element parareal prototype (whfaster round 2, iteration 3).
# F: WHFast512 C kernel, dt = 6 d, GR, corrector 17 (as in pr.py).
# G: osculating Jacobi elements -> mean elements (first-order numerical theory, meanel.py) -> Kepler + first-order
#    secular drift over the slice -> osculating elements. Fourier tables are built once per slice boundary (from the
#    first state seen there) and then frozen, so G is a fixed smooth map.
# Correction: U^{k+1}_{n+1} = G(U^{k+1}_n) + F(U^k_n) - G(U^k_n) in osculating elements (angles wrapped).
import numpy as np, time, argparse, sys
from ptlib import *
import meanel
p=argparse.ArgumentParser()
p.add_argument('--slices',type=int,default=10)
p.add_argument('--slice_days',type=float,default=360000)
p.add_argument('--kmax',type=int,default=6)
p.add_argument('--Ng',type=int,default=64)
p.add_argument('--sec',type=int,default=1, help='include first-order secular rates in G')
p.add_argument('--hybrid',type=int,default=1, help='outer planets by WH in G')
p.add_argument('--dt_out',type=float,default=180)
p.add_argument('--refresh',type=int,default=0, help='rebuild boundary tables every iteration')
a=p.parse_args()
m=np.array(mass); N=a.slices; T=a.slice_days*DAY
Mint=np.cumsum(m)[:-1]; mu=Mint+m[1:]
tables={}
def table(n, E):
    if n not in tables: tables[n]=meanel.Theory(m,E,Ng=a.Ng,inner=4 if a.hybrid else None)
    return tables[n]
def outerWH(E):
    s=rebound.Simulation(); s.add(m=Mint[4])
    for mm,e in zip(m[5:],E[4:]):
        ecc=np.hypot(e[2],e[3]); pom=np.arctan2(e[3],e[2]); si=np.hypot(e[4],e[5]); Om=np.arctan2(e[5],e[4])
        s.add(m=mm, a=np.exp(e[0]), l=e[1], e=ecc, pomega=pom, inc=2*np.arcsin(si), Omega=Om)
    s.move_to_com()
    s.integrator="whfast512"; s.dt=a.dt_out*DAY; s.exact_finish_time=0
    s.integrator.gr_potential=1; s.integrator.corrector=17; s.integrator.kernel=2
    nst=int(round(a.slice_days/a.dt_out)); s.integrator.concatenate_steps=nst; s.steps(1)
    out=[]
    for o in s.orbits():
        si=np.sin(o.inc/2)
        out.append([np.log(o.a), o.l, o.e*np.cos(o.pomega), o.e*np.sin(o.pomega), si*np.cos(o.Omega), si*np.sin(o.Omega)])
    return np.array(out)
def G(E, n):
    th0=table(n,E)
    M=th0.mean(E)
    if a.hybrid: M[4:]=E[4:]
    nbar=np.sqrt(mu/np.exp(3*M[:,0]))
    R=th0.secular_rates() if a.sec else np.zeros((8,6))
    M2=M.copy()
    M2[:,1]=wrap(M[:,1]+(nbar+R[:,1])*T)
    M2[:,2:]+=R[:,2:]*T
    M2[:,0]+=R[:,0]*T
    if not a.hybrid:
        th1=table(n+1,M2)
        return th1.osc(M2)
    Eo=outerWH(E)
    M2[4:]=Eo
    th1=table(n+1,M2)
    Ee=M2.copy()
    for _ in range(3):
        Ee[:4]=M2[:4]+th1.disp(Ee)[:4]
    return Ee
def F(E): return to_el(prop(from_el(E),a.slice_days,6))
def rep(tag, errs, extra=""):
    A=np.max([np.abs(e[:,0]) for e in errs],axis=0); L=np.max([np.abs(e[:,1]) for e in errs],axis=0)
    print(f"  {tag:5s} lna: "+" ".join(f"{x:.0e}" for x in A)+" | lam: "+" ".join(f"{x:.0e}" for x in L)+extra); sys.stdout.flush()
E0=to_el(initial())
t0=time.time(); ref=[E0]
for n in range(N): ref.append(F(ref[-1]))
tF=(time.time()-t0)/N
# round-off noise reference: asm kernel from same ICs
ref1=[E0]
for n in range(N): ref1.append(to_el(prop(from_el(ref1[-1]),a.slice_days,6,kernel=1)))
rep("noise",[ediff(ref1[n],ref[n]) for n in range(1,N+1)])
t0=time.time()
U=[E0]
for n in range(N): U.append(G(U[-1],n))
tG=(time.time()-t0)/N
print(f"slices {N} x {a.slice_days/365.25:.0f} yr; python cost F/slice {tF:.3f}s, G/slice (incl. table builds) {tG:.3f}s")
Gold=[None]+U[1:]
rep("k=0",[ediff(U[n],ref[n]) for n in range(1,N+1)])
for k in range(1,a.kmax+1):
    if a.refresh: tables.clear()
    Fv=[F(U[n]) for n in range(N)]
    Un=[E0]; Gn=[None]
    for n in range(N):
        g=G(Un[n],n); Gn.append(g)
        c=g+ediff(Fv[n],Gold[n+1]); c[:,1]=wrap(c[:,1]); Un.append(c)
    U=Un; Gold=Gn
    e=[ediff(U[n],ref[n]) for n in range(1,N+1)]
    pe=max(poserr(from_el(U[n]),from_el(ref[n])) for n in range(1,N+1))
    rep(f"k={k}",e,f"  | max pos err {pe:.1e} AU")
