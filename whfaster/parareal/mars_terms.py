# Round 3, iteration 23: identify the long-period terms behind the remaining Mars (k, h) per-slice errors:
# per-slice error z = dk + i dh (200 slices) regressed on single harmonics exp(i phi), phi = sum k_p lambda_p
# with sum k_p = 1 (d'Alembert for an e-type term), frequency < 1/600 yr.
import numpy as np, itertools, sys, ph1mn
from ptlib import *
d=np.fromfile('../pint/runs/m10_ser2.bin',dtype=np.float64).reshape(-1,50); X=d[:,2:].reshape(-1,8,6)
s0=int(sys.argv[1]) if len(sys.argv)>1 else 1000; M=200
try: err=np.load(f'mars_err_{s0}.npy')
except FileNotFoundError:
    G=ph1mn.build('all_mn1_g3_dt180_nsub2_un')
    err=np.array([ediff(G(X[n]),X[n+1]) for n in range(s0,s0+M)]); np.save(f'mars_err_{s0}.npy',err)
L=X[s0:s0+M,:,1]; mu=np.cumsum(mass)[:-1]+np.array(mass[1:]); n=np.sqrt(mu/np.exp(3*X[s0,:,0]))  # rad/(yr/2pi)
T=360000*DAY
pl=[0,1,2,3,4,5]; names=NAMES
for planet in (3,1,2):
    z=err[:,planet,2]+1j*err[:,planet,3]
    res=[]
    for ks in itertools.product(range(-9,10),repeat=len(pl)):
        if sum(ks)!=1 or sum(1 for k in ks if k)>3 or sum(abs(k) for k in ks)>15: continue
        nu=sum(k*n[p] for k,p in zip(ks,pl))
        if abs(nu)*T>2*np.pi*986/600: continue
        ph=np.exp(1j*sum(k*L[:,p] for k,p in zip(ks,pl)))
        A=np.vdot(ph,z)/M; r=z-A*ph
        res.append((1-np.sum(np.abs(r)**2)/np.sum(np.abs(z)**2),abs(A),2*np.pi/abs(nu)/(2*np.pi) if nu else np.inf,ks))
    res.sort(key=lambda q:-q[0])
    print(f"{names[planet]} z-error rms {np.sqrt(np.mean(np.abs(z)**2)):.1e}")
    for r2,a,P,ks in res[:6]:
        print(f"   R2 {r2:.2f} |A| {a:.1e} P {P*2*np.pi/(2*np.pi):8.0f} (yr/2pi units->yr {P:.0f})  "+" ".join(f"{k:+d}{names[p]}" for k,p in zip(ks,pl) if k))
