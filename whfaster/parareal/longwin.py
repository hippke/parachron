# Round 3, Phase 7 feasibility: parareal over one long window (N = 40 slices) — the iterate structure of a sliding /
# continuously pipelined parareal. Error / noise per boundary and iteration (reference: serial run, noise = asm vs C
# kernel difference per slice accumulated as in phase0).
import numpy as np, sys, ptc
from ptlib import *
d=np.fromfile('../pint/runs/m10_ser2.bin',dtype=np.float64).reshape(-1,50); X=d[:,2:].reshape(-1,8,6)
s0=int(sys.argv[1]) if len(sys.argv)>1 else 0; N=int(sys.argv[2]) if len(sys.argv)>2 else 40; KM=int(sys.argv[4]) if len(sys.argv)>4 else 6
G=ptc.CoarseC(360000,heun=1)
if len(sys.argv)>3:
    import ph1mn; G=ph1mn.build(sys.argv[3])
def F(E,k=2): return to_el(prop(from_el(E),360000,6,kernel=k))
# reference: serial fine from X[s0] (same kernel) and noise per boundary (kernel 1 vs 2)
ref=[X[s0]]; ref1=[X[s0]]
for n in range(N): ref.append(F(ref[-1])); ref1.append(F(ref1[-1],1))
noise=np.array([max(poserr(from_el(ref1[n]),from_el(ref[n])),1e-12) for n in range(N+1)])
U=[X[s0]]
for n in range(N): U.append(G(U[-1]))
Gold=[None]+U[1:]
print("k=0 err/noise at n=5,10,20,40:",[f"{poserr(from_el(U[n]),from_el(ref[n]))/noise[n]:.1e}" for n in (5,10,20,N)])
for k in range(1,KM+1):
    Fv=[F(U[n]) for n in range(N)]
    Un=[X[s0]]; Gn=[None]
    for n in range(N):
        g=G(Un[n]); Gn.append(g); c=g+ediff(Fv[n],Gold[n+1]); c[:,1]=wrap(c[:,1]); Un.append(c)
    U=Un; Gold=Gn
    e=np.array([poserr(from_el(U[n]),from_el(ref[n]))/noise[n] for n in range(1,N+1)])
    first_bad=np.argmax(e>2)+1 if np.any(e>2) else None
    print(f"k={k}: err/noise at n=5,10,20,30,40: "+" ".join(f"{e[n-1]:.1e}" for n in (5,10,20,30,N))+f" | max {e.max():.1e} | boundaries <=2x noise up to n={ (first_bad-1) if first_bad else N}",flush=True)
