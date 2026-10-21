# Round 3: dense sampling (every 6 d) of a fine run started from boundary s0 of the 10 Myr serial run; records
# osculating elements and the first-order mean ln a (C transform).   python3 resid_at.py s0 years out.npz [ng_hi]
import numpy as np, sys
from ptlib import *
import ptc
s0=int(sys.argv[1]); yrs=float(sys.argv[2]); out=sys.argv[3]
G=ptc.CoarseC(360000,ng_hi=int(sys.argv[4]) if len(sys.argv)>4 else 0)
d=np.fromfile('../pint/runs/m10_ser2.bin',dtype=np.float64).reshape(-1,50); X=d[:,2:].reshape(-1,8,6)
s=sim_from(from_el(X[s0])); s.integrator="whfast512"; s.dt=6*DAY; s.exact_finish_time=0
s.integrator.gr_potential=1; s.integrator.corrector=17; s.integrator.concatenate_steps=1
n=int(yrs*365.25/6); T=np.zeros(n); EL=np.zeros((n,8,6)); ML=np.zeros((n,4))
for q in range(n):
    s.steps(1); E=to_el(state(s)); T[q]=s.t; EL[q]=E; ML[q]=E[:4,0]-ptc.lna_disp(G,E)
np.savez(out,t=T,el=EL,mlna=ML)
