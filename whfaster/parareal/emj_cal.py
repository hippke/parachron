# Round 3, it. 24: is the remaining Mars ln a defect an amplitude/phase error of the E-M-J term? Regress the per-slice
# Mars ln a defects on the term's own slice increment basis B = A (e^{i phi1} - e^{i phi0}) (complex factor c - 1).
import numpy as np, sys, so2, g3, ptc
from ptlib import *
d=np.fromfile('../pint/runs/m10_ser2.bin',dtype=np.float64).reshape(-1,50); X=d[:,2:].reshape(-1,8,6)
K=g3.TERMS[0]; G0=ptc.CoarseC(360000,heun=1,dt_out=180); T=360000*DAY
gt=g3.G3()
def basis(E):
    g=G0(E); n=gt.mean_motions(E,g); U=so2.resonant_U2(E,K,Ng=32)
    nu=sum(k*n[p] for p,k in K.items()); p0=sum(k*E[p,1] for p,k in K.items())
    return 2*U[3][0]/(1j*nu)*(np.exp(1j*(p0+nu*T))-np.exp(1j*p0)), nu
for name,s0,f in (('cal_4500',4400,'cal_4500.npy'),('serr_4500',4500,'serr_4500.npy'),('serr_1000',1000,'serr_1000.npy'),('serr_6000',6000,'serr_6000.npy')):
    err=np.load(f)[:100]; Bs=[]; nus=[]
    for q in range(len(err)): b,nu=basis(X[s0+q]); Bs.append(b); nus.append(nu)
    Bs=np.array(Bs); y=err[:,3,0]
    A=np.vstack([Bs.real,-Bs.imag]).T; c,*_=np.linalg.lstsq(A,y,rcond=None); r=y-A@c
    print(f"{name}: Mars lna defect rms {y.std():.1e} -> {r.std():.1e}  (c-1) = {c[0]:+.3f}{c[1]:+.3f}i  |term incr| rms {np.abs(Bs).std():.1e}",flush=True)
