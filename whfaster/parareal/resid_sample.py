# Iteration 16 (diagnosis for a second-order theory): sample a fine run every 6 d and record the osculating Jacobi
# elements and the first-order mean ln a (osculating - C transform, live tables) of the inner planets.
#   python3 resid_sample.py out.npz years
import numpy as np, sys, time
from ptlib import *
import ptc
out=sys.argv[1]; yrs=float(sys.argv[2])
G=ptc.CoarseC(360000,ng_hi=int(sys.argv[3]) if len(sys.argv)>3 else 0)
s=sim_from(initial()); s.integrator="whfast512"; s.dt=6*DAY; s.exact_finish_time=0
s.integrator.gr_potential=1; s.integrator.corrector=17; s.integrator.concatenate_steps=1
n=int(yrs*365.25/6)
T=np.zeros(n); EL=np.zeros((n,8,6),dtype=np.float64); ML=np.zeros((n,4))
t0=time.time()
for q in range(n):
    s.steps(1)
    E=to_el(state(s)); T[q]=s.t; EL[q]=E
    ML[q]=E[:4,0]-ptc.lna_disp(G,E)
    if q%50000==0: print(q,n,f"{time.time()-t0:.0f}s"); sys.stdout.flush()
np.savez(out,t=T,el=EL,mlna=ML)
