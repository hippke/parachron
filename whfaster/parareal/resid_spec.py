# Spectrum of the residual (first-order mean) ln a of the inner planets; peaks matched to combinations of mean motions.
#   python3 resid_spec.py samples.npz
import numpy as np, sys, itertools
from ptlib import NAMES
d=np.load(sys.argv[1]); t=d['t']/(2*np.pi); el=d['el']; ml=d['mlna']
dt=t[1]-t[0]
lam=np.unwrap(el[:,:,1],axis=0)
nn=np.array([np.polyfit(t,lam[:,j],1)[0]/(2*np.pi) for j in range(8)])   # mean motions [cycles/yr]
print("mean motions [1/yr]:"," ".join(f"{NAMES[j]} {nn[j]:.6f}" for j in range(8)))
# candidate combinations
cands=[]
for i,j in itertools.combinations(range(8),2):
    for ki in range(-20,21):
        for kj in range(-20,21):
            if ki==0 and kj==0: continue
            if abs(ki)+abs(kj)>24: continue
            f=ki*nn[i]+kj*nn[j]
            if f>0: cands.append((f,f"{ki:+d}{NAMES[i]}{kj:+d}{NAMES[j]}",abs(ki)+abs(kj)))
for i,j,l in itertools.combinations(range(8),3):
    for ks in itertools.product(range(-5,6),repeat=3):
        if 0 in ks: continue
        f=ks[0]*nn[i]+ks[1]*nn[j]+ks[2]*nn[l]
        if f>0: cands.append((f,f"{ks[0]:+d}{NAMES[i]}{ks[1]:+d}{NAMES[j]}{ks[2]:+d}{NAMES[l]}",sum(abs(k) for k in ks)+2))
cf=np.array([c[0] for c in cands])
def match(f, tol):
    idx=np.where(np.abs(cf-f)<tol)[0]
    if not len(idx): return "?"
    best=sorted(idx,key=lambda q:(cands[q][2],abs(cf[q]-f)))[:3]
    return ", ".join(f"{cands[q][1]} (P={1/cf[q]:.3g} yr, ord {cands[q][2]})" for q in best)
def peaks(x, dtt, nk, label, fmin=0, fmax=None):
    x=x-np.polyval(np.polyfit(np.arange(len(x)),x,1),np.arange(len(x)))
    w=np.hanning(len(x)); P=4
    F=np.fft.rfft(x*w,P*len(x)); fr=np.fft.rfftfreq(P*len(x),dtt)
    A=2*np.abs(F)/w.sum()
    sel=(fr>fmin)&((fr<fmax) if fmax else True)
    A2=A.copy(); A2[~sel]=0
    res=f"  {label}: rms {x.std():.2e}\n"
    for _ in range(nk):
        k=np.argmax(A2)
        if A2[k]==0: break
        f=fr[k]; tol=2/(len(x)*dtt)
        res+=f"     {A[k]:.2e} at {f:.5f}/yr (P={1/f:.4g} yr): {match(f,tol)}\n"
        lo=max(0,k-3*P); A2[lo:k+3*P]=0
    print(res,end="")
# short-period part: first 300 yr at full sampling
m300=int(300/dt)
for i in range(4):
    print(f"{NAMES[i]} mean ln a, short-period residual (first 300 yr):")
    peaks(ml[:m300,i], dt, 6, "dense")
# long-period part: 1-yr block averages over the whole run
nb=int(1/dt)
nblk=len(t)//nb
for i in range(4):
    x=ml[:nblk*nb,i].reshape(nblk,nb).mean(1)
    print(f"{NAMES[i]} mean ln a, long-period (1-yr averages, {nblk} yr):")
    peaks(x, nb*dt, 6, "blocks", fmin=1/5000., fmax=0.45)
