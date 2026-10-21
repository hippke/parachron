# Long-run comparison (iteration 34+): serial runs vs parareal runs of ~1 Gyr, records (t, E, 8x6 elements) per slice.
# Point-by-point comparison is meaningless after decorrelation (~100 Myr), so compare statistics:
# energy error, secular frequencies in consecutive windows (distributions, KS tests), Mercury's eccentricity.
#   python3 gyr_ana.py ref.bin other.bin [...]   (first file = reference; window length WIN_MYR env, default 20)
import numpy as np, sys, os
from scipy.optimize import minimize_scalar
from scipy.stats import ks_2samp
AS=180/np.pi*3600; WIN=float(os.environ.get('WIN_MYR','20'))
def load(fn):
    d=np.memmap(fn,dtype=np.float64,mode='r'); d=d[:len(d)//50*50].reshape(-1,50); return d
def naff1(t,z):
    z=z-z.mean(); T=t[-1]-t[0]; tc=t-t[0]-T/2; w=0.5*(1+np.cos(2*np.pi*tc/T)); n=len(z); P=8
    F=np.fft.fft(z*w,P*n); fr=np.fft.fftfreq(P*n,d=(t[1]-t[0]))*2*np.pi; k=np.argmax(np.abs(F)); w0=fr[k]; dw=2*np.pi/(P*n*(t[1]-t[0]))
    f=lambda om: -np.abs(np.sum(z*w*np.exp(-1j*om*tc)))
    return minimize_scalar(f,bracket=(w0-dw,w0,w0+dw),tol=1e-14).x*AS
names=['g1','g2','g3','g4','g5','g6','s1','s2','s3','s4']
def analyse(fn, nmax):
    d=load(fn)[:nmax]; t=d[:,0]/(2*np.pi); E=d[:,1]; X=d[:,2:].reshape(-1,8,6)
    dE=(E-E[0])/abs(E[0])
    nw=int(round(WIN*1e6/(t[1]-t[0]))); res=[]
    for a in range(0,len(t)-nw+1,nw):
        s=slice(a,a+nw); tt=t[s]
        row=[naff1(tt,X[s,i,2]+1j*X[s,i,3]) for i in range(6)]+[naff1(tt,X[s,i,4]+1j*X[s,i,5]) for i in range(4)]
        res.append(row)
    eM=np.hypot(X[:,0,2],X[:,0,3])
    return t,dE,np.array(res),eM
files=sys.argv[1:]
nmax=min(len(load(f)) for f in files)
R={f:analyse(f,nmax) for f in files}
t=R[files[0]][0]
print(f"span {t[-1]/1e6:.1f} Myr ({nmax} records), frequency windows of {WIN} Myr: {len(R[files[0]][2])}")
print("energy dE/E (vs own start): final | max | rms of last 10%  ; Mercury e: max | mean")
for f in files:
    _,dE,_,eM=R[f]; n=len(dE)
    print(f"  {os.path.basename(f):22s} {dE[-1]:+.2e} | {np.abs(dE).max():.2e} | {np.sqrt(np.mean(dE[-n//10:]**2)):.2e} ; {eM.max():.4f} | {eM.mean():.4f}")
F0=R[files[0]][2]
print("secular frequencies: mean ± std over windows [\"/yr], KS p-value vs "+os.path.basename(files[0]))
for q,nm in enumerate(names):
    line=f"  {nm:3s} "+" ".join(f"{np.mean(R[f][2][:,q]):+10.5f}±{np.std(R[f][2][:,q]):.4f}" for f in files)
    line+="  p: "+" ".join(f"{ks_2samp(F0[:,q],R[f][2][:,q]).pvalue:.2f}" for f in files[1:])
    print(line)
