# Calibrated coarse propagator: G_cal(U) = G(U) + c, c = mean defect F - G over the previous window's slices.
import numpy as np, sys, argparse
from ptlib import *
import ptc
p=argparse.ArgumentParser(); p.add_argument('--N',type=int,default=10); p.add_argument('--kmax',type=int,default=5)
p.add_argument('--days',type=float,default=360000); p.add_argument('--cal',default='all')
a=p.parse_args(); N=a.N; Td=a.days
G0=ptc.CoarseC(Td)
def F(E,kernel=2): return to_el(prop(from_el(E),Td,6,kernel=kernel))
# window A (previous) and B (tested)
E=to_el(initial()); A=[E]
for n in range(N): A.append(F(A[-1]))
D=np.mean([ediff(A[n+1],G0(A[n])) for n in range(N)],axis=0)
if a.cal=='lam': D[:,[0,2,3,4,5]]=0
if a.cal=='none': D[:]=0
B=[A[-1]]; B1=[A[-1]]
for n in range(N): B.append(F(B[-1])); B1.append(F(B1[-1],1))
noise=max(poserr(from_el(B1[n]),from_el(B[n])) for n in range(1,N+1))
def G(E):
    g=G0(E)+D; g[:,1]=wrap(g[:,1]); return g
U=[B[0]]
for n in range(N): U.append(G(U[-1]))
Gold=[None]+U[1:]
print(f"cal={a.cal}: noise {noise:.2e}; k=0 err {max(poserr(from_el(U[n]),from_el(B[n])) for n in range(1,N+1)):.2e}")
for k in range(1,a.kmax+1):
    Fv=[F(U[n]) for n in range(N)]; Un=[B[0]]; Gn=[None]
    for n in range(N):
        g=G(Un[n]); Gn.append(g); c=g+ediff(Fv[n],Gold[n+1]); c[:,1]=wrap(c[:,1]); Un.append(c)
    U=Un; Gold=Gn
    e=max(poserr(from_el(U[n]),from_el(B[n])) for n in range(1,N+1))
    print(f"  k={k} err {e:.2e}"+("  CONVERGED" if e<=3*noise else "")); sys.stdout.flush()
    if e<=3*noise: break
