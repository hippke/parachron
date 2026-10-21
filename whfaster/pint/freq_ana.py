# Secular-frequency comparison (accuracy contract item 3): NAFF-style frequency analysis of z = k + i h and
# zeta = q + i p (Jacobi elements, one sample per slice boundary) for several runs of the same span.
#   python3 freq_ana.py file1.bin [file2.bin ...]   (records: t, E, 8x6 elements)
import numpy as np, sys
from scipy.optimize import minimize_scalar
NAMES=['Mer','Ven','Ear','Mar','Jup','Sat','Ura','Nep']
AS=180/np.pi*3600   # rad/yr -> arcsec/yr
def load(fn):
    d=np.fromfile(fn,dtype=np.float64).reshape(-1,50); return d[:,0]/(2*np.pi), d[:,2:].reshape(-1,8,6)
def naff(t, z, nf=3):
    z=z-z.mean(); T=t[-1]-t[0]; tc=t-t[0]-T/2
    w=0.5*(1+np.cos(2*np.pi*tc/T))          # Hann
    out=[]
    for _ in range(nf):
        n=len(z); P=8
        F=np.fft.fft(z*w, P*n); fr=np.fft.fftfreq(P*n, d=(t[1]-t[0]))*2*np.pi
        k=np.argmax(np.abs(F)); w0=fr[k]; dw=2*np.pi/(P*n*(t[1]-t[0]))
        f=lambda om: -np.abs(np.sum(z*w*np.exp(-1j*om*tc)))
        om=minimize_scalar(f, bracket=(w0-dw,w0,w0+dw), tol=1e-14).x
        A=np.sum(z*w*np.exp(-1j*om*tc))/np.sum(w)
        out.append((om*AS, abs(A)))
        z=z-A*np.exp(1j*om*tc)
    return out
runs=[(fn,)+load(fn) for fn in sys.argv[1:]]
nmin=min(len(r[1]) for r in runs)
res={}
for fn,t,X in runs:
    t=t[:nmin]; X=X[:nmin]
    res[fn]=[naff(t, X[:,i,2]+1j*X[:,i,3]) for i in range(8)]+[naff(t, X[:,i,4]+1j*X[:,i,5]) for i in range(8)]
names=[f"g({p})" for p in NAMES]+[f"s({p})" for p in NAMES]
ref=sys.argv[1]
print(f"span {runs[0][1][nmin-1]-runs[0][1][0]:.4g} yr, {nmin} samples; leading frequency [arcsec/yr] of each signal, and differences vs {ref}")
print(f"{'signal':8s} {'ref':>12s} "+" ".join(f"{fn.split('/')[-1]:>22s}" for fn in sys.argv[2:]))
for q in range(16):
    f0=res[ref][q][0][0]
    print(f"{names[q]:8s} {f0:12.6f} "+" ".join(f"{res[fn][q][0][0]-f0:+22.2e}" for fn in sys.argv[2:]))
