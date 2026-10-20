# Summarize bias ensembles (output of ens.c): ensemble mean +- sem of dE/E vs steps, and linear-drift fit
import numpy as np, sys
for a in sys.argv[1:]:
    n,f=a.split('='); d=np.loadtxt(f)
    full=[r for r in np.unique(d[:,0]) if np.any((d[:,0]==r)&(d[:,1]==1e9))]
    d=d[np.isin(d[:,0],full)]
    print(f"{n}: {len(full)} complete runs")
    for N in [1e7,1e8,1e9]:
        v=np.array([d[(d[:,0]==r)&np.isclose(d[:,1],N,rtol=2e-3)][0,2] for r in full])
        print(f"   N={N:.0e}: mean {v.mean():+.2e} +- {v.std(ddof=1)/np.sqrt(len(v)):.1e}   rms {np.sqrt(np.mean(v**2)):.1e}   median|.| {np.median(abs(v)):.1e}")
