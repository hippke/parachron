# Secant-enhanced parareal: U^{k+1}_{n+1} = G(U^{k+1}_n) + D^k_n + A_n (U^{k+1}_n - U^k_n), D = F - G,
# A_n ~ (F' - G') at slice n estimated from the stored iterates of slice n (rank <= k, weighted least squares).
import numpy as np, sys, argparse
from ptlib import *
import ptc
p=argparse.ArgumentParser(); p.add_argument('--N',type=int,default=10); p.add_argument('--kmax',type=int,default=5)
p.add_argument('--days',type=float,default=360000); p.add_argument('--rank',type=int,default=9)
p.add_argument('--start',type=float,default=0); p.add_argument('--neigh',type=int,default=0); p.add_argument('--unw',type=int,default=0); p.add_argument('--damp',type=float,default=1.0)
a=p.parse_args(); N=a.N; Td=a.days
G=ptc.CoarseC(Td)
def F(E,kernel=2): return to_el(prop(from_el(E),Td,6,kernel=kernel))
U0=initial()
if a.start: U0=prop(U0,a.start,6)
E0=to_el(U0); ref=[E0]; ref1=[E0]
for n in range(N): ref.append(F(ref[-1])); ref1.append(F(ref1[-1],1))
noise=max(poserr(from_el(ref1[n]),from_el(ref[n])) for n in range(1,N+1))
# weights: lambda-equivalent (ln a scaled by the shear of the planet over one slice)
shear=1.5*np.sqrt(np.cumsum(mass)[:-1]/np.exp(3*E0[:,0]))*Td*DAY
w=np.ones((8,6)); w[:,0]=shear; w=w.ravel()
if a.unw: w=np.ones(48)
def vec(E): return E.ravel()
X=[[] for n in range(N)]; D=[[] for n in range(N)]   # history per slice: iterate start states and defects
U=[E0]
for n in range(N): U.append(G(U[-1]))
Gold=[None]+U[1:]
print(f"noise {noise:.2e}  k=0 err {max(poserr(from_el(U[n]),from_el(ref[n])) for n in range(1,N+1)):.2e}")
for k in range(1,a.kmax+1):
    Fv=[F(U[n]) for n in range(N)]
    for n in range(N):
        d=ediff(Fv[n],Gold[n+1]); X[n].append(U[n].copy()); D[n].append(d)
    Un=[E0]; Gn=[None]
    for n in range(N):
        g=G(Un[n]); Gn.append(g)
        corr=D[n][-1].copy()
        # secant correction from the history of slice n (and optionally neighbour slices)
        dx=ediff(Un[n],X[n][-1]).ravel()
        pairs=[]
        for m in ([n] + ([n-1] if a.neigh and n>0 else [])):
            for j in range(len(X[m])-1):
                pairs.append((ediff(X[m][j],X[m][-1]).ravel(), (D[m][j]-D[m][-1]).ravel()))
        pairs=pairs[-a.rank:] if a.rank>0 else []
        if pairs and np.any(dx):
            DX=np.array([q[0] for q in pairs]).T*w[:,None]; DD=np.array([q[1] for q in pairs]).T
            coef,*_=np.linalg.lstsq(DX, dx*w, rcond=1e-10)
            corr=corr+a.damp*(DD@coef).reshape(8,6)
        c=g+corr; c[:,1]=wrap(c[:,1]); Un.append(c)
    U=Un; Gold=Gn
    e=max(poserr(from_el(U[n]),from_el(ref[n])) for n in range(1,N+1))
    print(f"  k={k} err {e:.2e}"+("  CONVERGED" if e<=3*noise else "")); sys.stdout.flush()
    if e<=3*noise: break
