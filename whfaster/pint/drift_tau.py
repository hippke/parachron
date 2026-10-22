# Iteration 47: energy trend over the first SPAN_MYR (default 100) of the ensemble members: serial (runs/eg/s<q>),
# parareal tau 3e-5 (runs/eg/p<q>) and tau 3e-6 (runs/eg/t3p<q>); mean +- standard error of the per-run linear trends.
#   python3 drift_tau.py [n_members]
import numpy as np, sys, os
SPAN=float(os.environ.get('SPAN_MYR','100'))
def trend(fn):
    d=np.fromfile(fn,dtype=np.float64); d=d[:len(d)//50*50].reshape(-1,50)
    t=d[:,0]/(2*np.pi)/1e9; m=t<=SPAN/1e3; E=d[m,1]; t=t[m]
    return np.polyfit(t[::20],(E[::20]-E[0])/abs(E[0]),1)[0], (E[-1]-E[0])/abs(E[0])
M=int(sys.argv[1]) if len(sys.argv)>1 else 16
G={"serial":"runs/eg/s{q}_ser1.bin","parareal tau 3e-5":"runs/eg/p{q}_pr.bin","parareal tau 3e-6":"runs/eg/t3p{q}_pr.bin"}
R={}
for lbl,pat in G.items():
    v=[trend(pat.format(q=q)) for q in range(1,M+1) if os.path.exists(pat.format(q=q))]
    R[lbl]=np.array([x[0] for x in v])
    print(f"{lbl:18s} n={len(v):2d} trend per Gyr {R[lbl].mean():+.2e} +- {R[lbl].std(ddof=1)/np.sqrt(len(v)):.1e} (sd {R[lbl].std(ddof=1):.1e}); dE/E at {SPAN:.0f} Myr mean {np.mean([x[1] for x in v]):+.1e}")
for lbl in list(G)[1:]:
    d=R[lbl].mean()-R["serial"].mean(); e=np.hypot(R[lbl].std(ddof=1)/np.sqrt(len(R[lbl])),R["serial"].std(ddof=1)/np.sqrt(len(R["serial"])))
    print(f"{lbl} - serial: {d:+.2e} +- {e:.1e} per Gyr ({d/e:+.1f} sigma)")
