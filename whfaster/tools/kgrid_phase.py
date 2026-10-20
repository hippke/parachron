import numpy as np, sys
for a in sys.argv[1:]:
    n,f=a.split('='); d=np.loadtxt(f); d=d.reshape(256,16,4)
    ok=np.all(np.abs(d[:,:,2])<1e-15,axis=1) & np.all(np.abs(d[:,:,3])<1e-11,axis=1)
    def chi(col):
        x=d[ok][:,:,col]; m=x.mean(1); s=x.std(1,ddof=1)/4; t=m/s; return np.mean(t**2), np.mean(t>0)
    c1,p1=chi(2); cN,pN=chi(3)
    print(f"{n:9s} points {ok.sum()}: 1 step: <t^2>={c1:.2f} frac(t>0)={p1:.2f} | 1000 orbits: <t^2>={cN:.2f} frac(t>0)={pN:.2f}  median|err| {np.median(np.abs(d[ok][:,:,3])):.2e}")
