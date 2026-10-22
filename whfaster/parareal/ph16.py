# Iteration 41: does a better coarse step allow K = 2 for 16-slice windows? Nonlinear parareal (N = 16 x 986 yr) from
# the 10 Phase-0 epochs of the 10 Myr serial run (pint/runs/m10_ser2.bin = the exact fine chain), for one coarse
# variant per process. Prints per epoch, for k = 1..3: max true position error at the 16 boundaries [AU], the same
# relative to the asm-vs-C noise of the epoch, and the max boundary update |U^k - U^(k-1)| [AU] (what the stop rule
# of pint_pipe sees; the window stops at k when that is <= tau).
#   python3 ph16.py variant epoch     variant: B (config B, C) | Cmn64 (C, mean-n + 64^2 V-E) | best (round-3 Python G)
#                                     | mix0 (accurate step for the k = 0 prediction only, config B afterwards)
#                                     | any ph1mn.py spec string (e.g. ablations of 'best')
import numpy as np, sys, os, time
from ptlib import *
N=16; days=360000
var=sys.argv[1]; ep=int(sys.argv[2])
if var=='Cmn64':
    os.environ['PINT_MEANN']='1'; os.environ['PINT_NGHI']='64'; os.environ['PINT_ALPHAHI']='0.7'
import ptc
d=np.fromfile('../pint/runs/m10_ser2.bin',dtype=np.float64).reshape(-1,50); X=d[:,2:].reshape(-1,8,6)
ref=X[ep:ep+N+1]
noise=float(np.load('phase0_jac.npz')[f'noise{ep}'])
G0=None
if var in ('B','Cmn64'): G=ptc.CoarseC(days,heun=1)
elif var=='mix0':   # accurate step only for the k = 0 prediction, config B for all corrections
    import ph1mn; G0=ph1mn.build('all_mn1_g3_dt180_nsub1_un'); G=ptc.CoarseC(days,heun=1)
else:   # 'best' = the round-3 configuration; anything else is a ph1mn spec string (ablations)
    import ph1mn; G=ph1mn.build('all_mn1_g3_dt180_nsub2_un_ga120_w2k_frzs' if var=='best' else var)
def F(E): return to_el(prop(from_el(E),days,6))
t0=time.time()
E0=ref[0]; U=[E0]
for n in range(N): U.append((G0 or G)(U[-1]))
Gold=[None]+(U[1:] if G0 is None else [G(U[n]) for n in range(N)])
out=[]
for k in range(1,4):
    Fv=[F(U[n]) for n in range(N)]
    Un=[E0]; Gn=[None]
    for n in range(N):
        g=G(Un[n]); Gn.append(g); c=g+ediff(Fv[n],Gold[n+1]); c[:,1]=wrap(c[:,1]); Un.append(c)
    err=max(poserr(from_el(Un[n]),from_el(ref[n])) for n in range(1,N+1))
    upd=max(poserr(from_el(Un[n]),from_el(U[n])) for n in range(1,N+1))
    out.append((err,err/noise,upd)); U=Un; Gold=Gn
print(f"{var[:42]:42s} epoch {ep:5d}: "+" | ".join(f"k{k+1} err {e:.1e} ({r:.3g}x) upd {u:.1e}" for k,(e,r,u) in enumerate(out))+f" | {time.time()-t0:.0f}s",flush=True)
