# One-slice mismatch diagnostic for mgc.Coarse variants: J_F, J_G at the window start, ranks (F'-G')_rc e0_c
# (e0 = actual coarse error of the slice), ln a rows weighted by the planet's shear.
import numpy as np, sys, argparse, pickle
from ptlib import *
import mgc
p=argparse.ArgumentParser()
p.add_argument('--slice_days',type=float,default=3600000)
p.add_argument('--live',type=int,default=1); p.add_argument('--nsub',type=int,default=2)
p.add_argument('--Ng',type=int,default=48); p.add_argument('--Ng_sec',type=int,default=32)
p.add_argument('--save',default='jac1.pkl')
a=p.parse_args()
T=a.slice_days
CG=mgc.Coarse(T,Ng=a.Ng,live=a.live,nsub=a.nsub,Ng_sec=a.Ng_sec)
H=np.array([1e-7,1e-4,1e-5,1e-5,1e-5,1e-5])
def F(E): return to_el(prop(from_el(E),T,6))
def jac(f,E):
    J=np.zeros((48,48))
    for i in range(8):
        for c in range(6):
            Ep=E.copy(); Em=E.copy(); Ep[i,c]+=H[c]; Em[i,c]-=H[c]
            J[:,i*6+c]=ediff(f(Ep),f(Em)).ravel()/(2*H[c])
    return J
E0=to_el(initial())
JF=jac(F,E0); JG=jac(lambda E: CG(E,0),E0)
e0=ediff(CG(F(E0),1),F(F(E0))).ravel()      # coarse error of the next slice started from the exact state
pickle.dump({'JF':JF,'JG':JG,'e0':e0},open(a.save,'wb'))
shear=np.array([abs(JF[i*6+1,i*6]) for i in range(8)])
w=np.ones(48)
for i in range(8): w[i*6]=shear[i]
EL=['lna','lam','k','h','q','p']
print("e0 (slice 2 coarse error): lna",' '.join(f'{x:.0e}' for x in e0[0::6]),'| lam',' '.join(f'{x:.0e}' for x in e0[1::6]),
      '| k',' '.join(f'{x:.0e}' for x in e0[2::6]))
C=np.abs((JF-JG)*e0[None,:])*w[:,None]
idx=np.dstack(np.unravel_index(np.argsort(-C.ravel())[:20],C.shape))[0]
print("top (F'-G')_rc e0_c (ln a rows x shear):")
for r,c in idx: print(f"   {NAMES[r//6]}.{EL[r%6]} <- {NAMES[c//6]}.{EL[c%6]}: {C[r,c]:.1e}   (mismatch {JF[r,c]-JG[r,c]:+.2e}, F' {JF[r,c]:+.2e})")
