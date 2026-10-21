# per-slice coarse errors (F - G as ediff(G,F)) along M slices of the 10 Myr serial run for a spec -> npy
import numpy as np, sys, ph1mn
from ptlib import *
d=np.fromfile('../pint/runs/m10_ser2.bin',dtype=np.float64).reshape(-1,50); X=d[:,2:].reshape(-1,8,6)
spec=sys.argv[1]; s0=int(sys.argv[2]); M=int(sys.argv[3]); out=sys.argv[4]
G=ph1mn.build(spec)
np.save(out,np.array([ediff(G(X[n]),X[n+1]) for n in range(s0,s0+M)]))
