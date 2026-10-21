# k needed for given coarse configuration (C coarse via ctypes), Python driver; prints error at k=1..kmax
import numpy as np, sys, argparse, time
from ptlib import *
import ptc
p=argparse.ArgumentParser(); p.add_argument('--N',type=int,default=20); p.add_argument('--kmax',type=int,default=4)
p.add_argument('--days',type=float,default=360000); p.add_argument('--start',type=float,default=0)
p.add_argument('--ng_ii',type=int,default=32); p.add_argument('--ng_io',type=int,default=16); p.add_argument('--maxj',type=int,default=5)
p.add_argument('--dt_out',type=float,default=360); p.add_argument('--heun',type=int,default=1); p.add_argument('--pcut',type=float,default=0)
a=p.parse_args(); N=a.N; Td=a.days
G=ptc.CoarseC(Td,dt_out=a.dt_out,ng_ii=a.ng_ii,ng_io=a.ng_io,maxj=a.maxj,heun=a.heun,pcut=a.pcut)
def F(E,kernel=2): return to_el(prop(from_el(E),Td,6,kernel=kernel))
U0=initial()
if a.start: U0=prop(U0,a.start,6)
E0=to_el(U0); ref=[E0]; ref1=[E0]
for n in range(N): ref.append(F(ref[-1])); ref1.append(F(ref1[-1],1))
noise=max(poserr(from_el(ref1[n]),from_el(ref[n])) for n in range(1,N+1))
t0=time.time(); U=[E0]
for n in range(N): U.append(G(U[-1]))
tG=(time.time()-t0)/N
Gold=[None]+U[1:]; out=[]
for k in range(1,a.kmax+1):
    Fv=[F(U[n]) for n in range(N)]; Un=[E0]; Gn=[None]
    for n in range(N):
        g=G(Un[n]); Gn.append(g); c=g+ediff(Fv[n],Gold[n+1]); c[:,1]=wrap(c[:,1]); Un.append(c)
    U=Un; Gold=Gn
    out.append(max(poserr(from_el(U[n]),from_el(ref[n])) for n in range(1,N+1))/noise)
print(f"{' '.join(sys.argv[1:]):60s} G {1e3*tG:.2f} ms | err/noise at k=1..: "+" ".join(f"{x:.2g}" for x in out))
