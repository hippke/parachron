import numpy as np, sys
for fn in sys.argv[1:]:
    d=np.loadtxt(fn); e1=d[:,2]; eN=d[:,3]
    def st(x):
        x=x[np.isfinite(x)]; return f"median|.| {np.median(abs(x)):.2e} max|.| {abs(x).max():.2e} frac>0 {np.mean(x>0):.3f} mean/sem {x.mean()/(x.std()/np.sqrt(len(x))):+.2f}"
    print(f"{fn:16s} 1 step: {st(e1)}\n{'':16s} 1000 orbits: {st(eN)}   nonfinite {np.sum(~np.isfinite(eN))}")
