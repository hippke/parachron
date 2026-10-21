# tau / K_min scan (iteration 37): K histogram, secular-frequency errors vs the dt = 3 d serial run ("truth") and the
# dt = 6 d serial run, energy vs the dt = 6 d serial run. Hard ICs, 10 Myr.
import numpy as np, sys, glob, re, os
from scipy.optimize import minimize_scalar
AS=180/np.pi*3600
def load(fn): d=np.fromfile(fn,dtype=np.float64); d=d[:len(d)//50*50].reshape(-1,50); return d[:,0]/(2*np.pi), d[:,1], d[:,2:].reshape(-1,8,6)
def naff1(t,z):
    z=z-z.mean(); T=t[-1]-t[0]; tc=t-t[0]-T/2; w=0.5*(1+np.cos(2*np.pi*tc/T)); n=len(z); P=8
    F=np.fft.fft(z*w,P*n); fr=np.fft.fftfreq(P*n,d=(t[1]-t[0]))*2*np.pi; k=np.argmax(np.abs(F)); w0=fr[k]; dw=2*np.pi/(P*n*(t[1]-t[0]))
    f=lambda om: -np.abs(np.sum(z*w*np.exp(-1j*om*tc)))
    return minimize_scalar(f,bracket=(w0-dw,w0,w0+dw),tol=1e-14).x*AS
def freqs(fn):
    t,E,X=load(fn); return np.array([naff1(t,X[:,i,2]+1j*X[:,i,3]) for i in range(8)]+[naff1(t,X[:,i,4]+1j*X[:,i,5]) for i in range(8)]), E
idx=[0,1,2,3,5,8,9,11]; nm=['g1','g2','g3','g4','g6','s1','s2','s4']
f3,_=freqs('runs/dt/d3_ser2.bin'); f6,E6=freqs('runs/p9/v16t3_ser2.bin')
print("errors vs dt=3 d [\"/yr]          "+" ".join(f"{n:>9s}" for n in nm)+"   max|err| | max|extra vs dt6| | dE/E vs dt6 | K hist (mean)")
def row(lbl,fn,log=None):
    f,E=freqs(fn); e=(f-f3)[idx]; x=(f-f6)[idx]
    kh=''
    if log and os.path.exists(log):
        m=re.search(r"iterations per window:(.*)\(mean ([0-9.]+)\)",open(log).read()); kh=f"{m.group(1).strip()} ({m.group(2)})" if m else ''
    print(f"{lbl:32s}"+" ".join(f"{v:+9.1e}" for v in e)+f"   {np.abs(e).max():.1e} | {np.abs(x).max():.1e} | {(E[-1]-E6[len(E)-1])/abs(E6[0]):+.1e} | {kh}")
row("serial dt 6 d","runs/p9/v16t3_ser2.bin")
row("serial dt 8 d","runs/dt/d8_ser2.bin")
row("serial dt 12 d","runs/dt/d12_ser2.bin")
row("parareal tau 3e-6 (JR=0)","runs/p9/pp16_3e-6_pr.bin")
row("parareal tau 1e-5 (JR=0)","runs/p9/pp16_1e-5_pr.bin")
for km in (3,2):
    for t in ("3e-5","1e-4","3e-4","1e-3"):
        fn=f"runs/tau/k{km}t{t}_pr.bin"
        if os.path.exists(fn): row(f"parareal Kmin {km} tau {t}",fn,f"runs/tau/k{km}t{t}.txt")
