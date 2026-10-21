# Linearised parareal diagnostics for an arbitrary coarse propagator (default: mgc.Coarse).
# Computes slice Jacobians JF, JG by central differences along the serial fine run, the linear prediction of the
# parareal errors, and the dominant mismatch terms (F'-G')_{rc} e_c (rows in ln a weighted by the planet's shear).
import numpy as np, sys, time, argparse, pickle
from ptlib import *
import mgc
p=argparse.ArgumentParser()
p.add_argument('--slices',type=int,default=10)
p.add_argument('--slice_days',type=float,default=360000)
p.add_argument('--kmax',type=int,default=6)
p.add_argument('--jf',default='')
p.add_argument('--save',default='jac2.pkl')
p.add_argument('--dt_out',type=float,default=180)
p.add_argument('--Ng',type=int,default=64)
p.add_argument('--live',type=int,default=0)
p.add_argument('--cg',default='py')
p.add_argument('--Ngj',type=int,default=0)
a=p.parse_args()
N=a.slices; T=a.slice_days
H=np.array([1e-7,1e-4,1e-5,1e-5,1e-5,1e-5])
CG=mgc.Coarse(T,Ng=a.Ng,dt_out=a.dt_out,live=a.live)
if a.cg=='c':
    import ptc; CG=ptc.CoarseC(T,dt_out=a.dt_out)
def F(E): return to_el(prop(from_el(E),T,6))
def jac(f, E):
    J=np.zeros((48,48))
    for i in range(8):
        for c in range(6):
            Ep=E.copy(); Em=E.copy(); Ep[i,c]+=H[c]; Em[i,c]-=H[c]
            J[:,i*6+c]=ediff(f(Ep),f(Em)).ravel()/(2*H[c])
    return J
E0=to_el(initial()); ref=[E0]
for n in range(N): ref.append(F(ref[-1]))
t0=time.time()
if a.jf: JF=pickle.load(open(a.jf,'rb'))['JF']
else: JF=[jac(F,ref[n]) for n in range(N)]
U=[E0]
for n in range(N): U.append(CG(U[-1],n))      # builds the boundary tables from the coarse prediction
JG=[jac(lambda E: CG(E,n), ref[n]) for n in range(N)]
print(f"Jacobians {time.time()-t0:.0f}s"); sys.stdout.flush()
pickle.dump({'JF':JF,'JG':JG,'ref':ref},open(a.save,'wb'))
shear=np.array([abs(JF[0][i*6+1,i*6]) for i in range(8)])
EL=['lna','lam','k','h','q','p']
def rep(tag, errs):
    A=np.max([np.abs(e[:,0]) for e in errs],axis=0); L=np.max([np.abs(e[:,1]) for e in errs],axis=0)
    print(f"  {tag:5s} lna: "+" ".join(f"{x:.0e}" for x in A)+" | lam: "+" ".join(f"{x:.0e}" for x in L)); sys.stdout.flush()
e0=[ediff(U[n],ref[n]) for n in range(N+1)]
e=[x.ravel() for x in e0]
print("linear prediction:"); rep("k=0",e0[1:])
hist=[e]
for k in range(1,a.kmax+1):
    en=[np.zeros(48)]
    for n in range(N): en.append(JG[n]@en[n]+(JF[n]-JG[n])@e[n])
    e=en; hist.append(e); rep(f"k={k}",[x.reshape(8,6) for x in e[1:]])
# dominant mismatch terms at iterations 0..2
w=np.ones(48);
for i in range(8): w[i*6]=shear[i]
for k in range(3):
    C=np.zeros((48,48))
    for n in range(N): C+=np.abs((JF[n]-JG[n])*hist[k][n][None,:])*w[:,None]
    idx=np.dstack(np.unravel_index(np.argsort(-C.ravel())[:12],C.shape))[0]
    print(f"top mismatch terms acting on e^{k} (summed over slices; lna rows x shear, i.e. in rad of lambda):")
    print("   "+"; ".join(f"{NAMES[r//6]}.{EL[r%6]}<-{NAMES[c//6]}.{EL[c%6]} {C[r,c]:.1e}" for r,c in idx))
