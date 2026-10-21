# Per-slice coarse prediction error along consecutive slices of the 10 Myr serial run (rms per planet and element).
import numpy as np, sys
from ptlib import *
d=np.fromfile('../pint/runs/m10_ser2.bin',dtype=np.float64).reshape(-1,50); X=d[:,2:].reshape(-1,8,6)
def quality(G, s0=1000, M=200, label=''):
    err=np.array([ediff(G(X[n]),X[n+1]) for n in range(s0,s0+M)])
    out=f"{label:28s}"
    for p,nm in ((1,'V'),(2,'E'),(3,'M')):
        out+=f" | {nm}: lam {err[:,p,1].std():.1e} lna {err[:,p,0].std():.1e} k {err[:,p,2].std():.1e} q {err[:,p,4].std():.1e}"
    print(out); sys.stdout.flush(); return err
if __name__=='__main__':
    import ptc, g3
    s0=int(sys.argv[1]) if len(sys.argv)>1 else 1000
    e0=quality(ptc.CoarseC(360000,heun=1),s0,label='config B')
    e1=quality(g3.G3(),s0,label='G3 (+E-M-J)')
    np.savez('slice_quality_%d.npz'%s0,b=e0,g3=e1)
