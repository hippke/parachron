# ensemble of IC realisations: secular-frequency shifts of parareal runs vs their serial reference
import numpy as np, sys
from scipy.optimize import minimize_scalar
AS=180/np.pi*3600
def load(fn):
    d=np.fromfile(fn,dtype=np.float64).reshape(-1,50); return d[:,0]/(2*np.pi), d[:,1], d[:,2:].reshape(-1,8,6)
def naff1(t,z):
    z=z-z.mean(); T=t[-1]-t[0]; tc=t-t[0]-T/2; w=0.5*(1+np.cos(2*np.pi*tc/T)); n=len(z); P=8
    F=np.fft.fft(z*w,P*n); fr=np.fft.fftfreq(P*n,d=(t[1]-t[0]))*2*np.pi; k=np.argmax(np.abs(F)); w0=fr[k]; dw=2*np.pi/(P*n*(t[1]-t[0]))
    f=lambda om: -np.abs(np.sum(z*w*np.exp(-1j*om*tc)))
    return minimize_scalar(f,bracket=(w0-dw,w0,w0+dw),tol=1e-14).x*AS
def freqs(fn):
    t,E,X=load(fn); n=min(len(t),10001); t=t[:n]; X=X[:n]
    return np.array([naff1(t,X[:,i,2]+1j*X[:,i,3]) for i in range(8)]+[naff1(t,X[:,i,4]+1j*X[:,i,5]) for i in range(8)]), E[:n]
names=[f"g{i+1}" for i in range(8)]+[f"s{i+1}" for i in range(8)]
pre=sys.argv[1]; perts=sys.argv[2].split(','); variants=sys.argv[3].split(',')
D={v:[] for v in variants}; dE={v:[] for v in variants}
for p in perts:
    fs,Es=freqs(f"{pre}{p}_{variants[0]}_ser2.bin" if p!='v16' else "runs/p9/v16t3_ser2.bin")
    for v in variants:
        fn=f"{pre}{p}_{v}_pr.bin" if p!='v16' else f"runs/p9/v16{v}_pr.bin"
        fp,Ep=freqs(fn); D[v].append(fp-fs); dE[v].append((Ep[-1]-Es[-1])/abs(Es[0]))
for v in variants:
    A=np.array(D[v]); inner=[0,1,2,3,8,9,10,11]
    print(f"== {v}: {len(A)} realisations; |shift| max inner {np.abs(A[:,inner]).max():.2e}, rms inner {np.sqrt((A[:,inner]**2).mean()):.2e}, rms giants {np.sqrt((A[:,[4,5,6,7,12,13,14,15]]**2).mean()):.2e}; energy at end vs serial: "+" ".join(f"{x:+.1e}" for x in dE[v]))
    print("   mean ± sem per inner frequency: "+"  ".join(f"{names[q]} {A[:,q].mean():+.1e}±{A[:,q].std(ddof=1)/np.sqrt(len(A)):.0e}" for q in inner))
