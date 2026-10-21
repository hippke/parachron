# Relation between the parareal update (max position change between iterations, known during a run) and the true
# error (vs the serial run), default configuration, 10 x 986 yr windows at 10 epochs.
import numpy as np
from ptlib import *
import ptc
d=np.fromfile('../pint/runs/m10_ser2.bin',dtype=np.float64).reshape(-1,50); SERX=d[:,2:].reshape(-1,8,6)
EP=[0,1030,2000,3400,4500,5500,6000,7000,8160,9300]; N=10; days=360000
G=ptc.CoarseC(days)
def F(E,k=2): return to_el(prop(from_el(E),days,6,kernel=k))
rows=[]
for ep in EP:
    E0=SERX[ep]; ref=[E0]; ref1=[E0]
    for n in range(N): ref.append(F(ref[-1])); ref1.append(F(ref1[-1],1))
    noise=max(poserr(from_el(ref1[n]),from_el(ref[n])) for n in range(1,N+1))
    U=[E0]
    for n in range(N): U.append(G(U[-1]))
    Gold=[None]+U[1:]
    for k in range(1,7):
        Fv=[F(U[n]) for n in range(N)]; Un=[E0]; Gn=[None]
        for n in range(N):
            g=G(Un[n]); Gn.append(g); c=g+ediff(Fv[n],Gold[n+1]); c[:,1]=wrap(c[:,1]); Un.append(c)
        upd=max(poserr(from_el(Un[n]),from_el(U[n])) for n in range(1,N+1))
        U=Un; Gold=Gn
        err=max(poserr(from_el(U[n]),from_el(ref[n])) for n in range(1,N+1))
        rows.append((ep,k,upd,err,noise))
R=np.array(rows); np.save('stoprule.npy',R)
print("epoch k  update[AU]  error[AU]  noise[AU]  err/noise")
for r in R: print(f"{int(r[0]):5d} {int(r[1])} {r[2]:.2e} {r[3]:.2e} {r[4]:.2e} {r[3]/r[4]:.2g}")
for tau in (1e-5,3e-6,1e-6,3e-7,1e-7,3e-8):
    Ks=[]; worst=0
    for ep in EP:
        sub=R[R[:,0]==ep]; ok=np.where(sub[:,2]<=tau)[0]
        k=int(ok[0]) if len(ok) else 5
        Ks.append(k+1); worst=max(worst,sub[k,3]/sub[k,4])
    print(f"stop when update <= {tau:.0e} AU: K per epoch {Ks}, mean {np.mean(Ks):.1f}, worst err/noise {worst:.2g}")
