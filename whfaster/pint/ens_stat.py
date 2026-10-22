# Iteration 46: statistics of an ensemble of serial vs parareal runs from the same IC realisations
# (runs/eg/s<q>_ser1.bin, runs/eg/p<q>_pr.bin). Per run: secular frequencies in consecutive windows (WIN_MYR, default
# 20), Mercury's maximum eccentricity, energy error at the end and its linear trend. Ensembles compared with KS tests:
# frequencies pooled over runs (windows after the first SKIP_MYR, default 50, i.e. after decorrelation), and the
# per-run maxima of Mercury's eccentricity.
#   python3 ens_stat.py [n_members]
import numpy as np, sys, os
from scipy.optimize import minimize_scalar
from scipy.stats import ks_2samp, mannwhitneyu
AS=180/np.pi*3600; WIN=float(os.environ.get('WIN_MYR','20')); SKIP=float(os.environ.get('SKIP_MYR','50'))
def load(fn): d=np.fromfile(fn,dtype=np.float64); return d[:len(d)//50*50].reshape(-1,50)
def naff1(t,z):
    z=z-z.mean(); T=t[-1]-t[0]; tc=t-t[0]-T/2; w=0.5*(1+np.cos(2*np.pi*tc/T)); n=len(z); P=8
    F=np.fft.fft(z*w,P*n); fr=np.fft.fftfreq(P*n,d=(t[1]-t[0]))*2*np.pi; k=np.argmax(np.abs(F)); w0=fr[k]; dw=2*np.pi/(P*n*(t[1]-t[0]))
    f=lambda om: -np.abs(np.sum(z*w*np.exp(-1j*om*tc)))
    return minimize_scalar(f,bracket=(w0-dw,w0,w0+dw),tol=1e-14).x*AS
names=['g1','g2','g3','g4','s1','s2','s3','s4']
def analyse(fn):
    d=load(fn); t=d[:,0]/(2*np.pi); E=d[:,1]; X=d[:,2:].reshape(-1,8,6)
    dE=(E-E[0])/abs(E[0]); slope=np.polyfit(t[::50]/1e9,dE[::50],1)[0]
    nw=int(round(WIN*1e6/(t[1]-t[0]))); ns=int(round(SKIP*1e6/(t[1]-t[0]))); F=[]
    for a in range(ns,len(t)-nw+1,nw):
        s=slice(a,a+nw); tt=t[s]
        F.append([naff1(tt,X[s,i,2]+1j*X[s,i,3]) for i in range(4)]+[naff1(tt,X[s,i,4]+1j*X[s,i,5]) for i in range(4)])
    return np.array(F), np.hypot(X[:,0,2],X[:,0,3]).max(), dE[-1], slope, t[-1]
M=int(sys.argv[1]) if len(sys.argv)>1 else 16
S=[analyse(f"runs/eg/s{q}_ser1.bin") for q in range(1,M+1)]
P=[analyse(f"runs/eg/p{q}_pr.bin") for q in range(1,M+1)]
print(f"{M} members, span {S[0][4]/1e6:.0f} Myr, frequency windows of {WIN} Myr after {SKIP} Myr")
for lbl,R in (("serial",S),("parareal",P)):
    em=np.array([r[1] for r in R]); de=np.array([r[2] for r in R]); sl=np.array([r[3] for r in R])
    print(f"{lbl:9s} Mercury e_max mean {em.mean():.4f} sd {em.std(ddof=1):.4f} [{em.min():.3f}, {em.max():.3f}]; dE/E end mean {de.mean():+.1e} rms {np.sqrt((de**2).mean()):.1e}; trend per Gyr mean {sl.mean():+.1e} sd {sl.std(ddof=1):.1e}")
print(f"Mercury e_max KS p = {ks_2samp([r[1] for r in S],[r[1] for r in P]).pvalue:.2f}; Mann-Whitney p = {mannwhitneyu([r[1] for r in S],[r[1] for r in P]).pvalue:.2f}")
print("energy trend KS p = %.2f"%ks_2samp([r[3] for r in S],[r[3] for r in P]).pvalue)
FS=np.vstack([r[0] for r in S]); FP=np.vstack([r[0] for r in P])
print("frequencies pooled over runs and windows: serial mean ± sd | parareal mean ± sd | KS p")
for q,nm in enumerate(names):
    print(f"  {nm}: {FS[:,q].mean():+10.5f} ± {FS[:,q].std():.4f} | {FP[:,q].mean():+10.5f} ± {FP[:,q].std():.4f} | {ks_2samp(FS[:,q],FP[:,q]).pvalue:.2f}")
