import numpy as np, sys
pairs=[a.split('=') for a in sys.argv[1:]]
D={n:np.loadtxt(f) for n,f in pairs}
mask=np.all([np.abs(x[:,2])<1e-15 for x in D.values()],axis=0) & np.all([np.abs(x[:,3])<1e-11 for x in D.values()],axis=0)
print("points:",mask.sum())
for n,x in D.items():
    yN=x[mask,3]; y1=x[mask,2]
    print(f"{n:9s} 1 step frac>0 {np.mean(y1>0):.3f} median {np.median(abs(y1)):.2e} | 1000 orb: frac>0 {np.mean(yN>0):.3f}  median|.| {np.median(abs(yN)):.2e}  mean/sem {yN.mean()/(yN.std()/np.sqrt(len(yN))):+.1f}")
