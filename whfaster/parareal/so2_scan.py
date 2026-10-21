# Round 3: scan of second-order three-planet resonant terms (so2.resonant_U2) for the inner planets at a given state:
# amplitude of the ln a oscillation A = 2|U|/|nu| (and of lambda ~ 1.5 n A / nu) for all K with 3 planets, small nu.
#   python3 so2_scan.py ep [n]
import numpy as np, itertools, sys, so2, ph1eval, ptc
from ptlib import *
ep=int(sys.argv[1]); n0=int(sys.argv[2]) if len(sys.argv)>2 else 0
E=ph1eval.D[f'ref{ep}'][n0]
G0=ptc.CoarseC(360000,heun=1); d,R=ptc.lna_disp_R(G0,E)
n=np.sqrt(so2.MU/np.exp(3*E[:,0])); n[:4]=np.sqrt(so2.MU[:4]/np.exp(3*(E[:4,0]-d)))+R[:,1]
T=360000*DAY
PL=[1,2,3,4,5]
pairs={(i,j):so2.Pair(E,i,j,32) for i,j in itertools.combinations(PL,2)}
out=[]
for tri in itertools.combinations(PL,3):
    if not any(p<4 for p in tri): continue
    for ks in itertools.product(range(-9,10),repeat=3):
        if 0 in ks or sum(map(abs,ks))>16 or abs(sum(ks))>3: continue
        if ks[0]<0: continue
        nu=sum(k*n[p] for k,p in zip(ks,tri))
        if abs(nu)*T>2*np.pi*986/400: continue      # period > 400 yr
        K=dict(zip(tri,ks))
        U=so2.resonant_U2(E,K,Ng=32,pairs={k:v for k,v in pairs.items() if set(k)<=set(tri)})
        for p in tri:
            if p>3: continue
            A=2*abs(U[p][0])/abs(nu); Al=1.5*n[p]*A/abs(nu); Az=2*np.hypot(abs(U[p][2]),abs(U[p][3]))/abs(nu)
            out.append((Al,A,Az,p,K,2*np.pi/abs(nu)/(2*np.pi)))
out.sort(key=lambda q:-q[0])
print(f"epoch {ep}: largest second-order three-planet resonant terms (lambda amplitude)")
for Al,A,Az,p,K,P in out[:20]:
    print(f"  {NAMES[p]}: lam amp {Al:.1e} rad, lna amp {A:.1e}, z amp {Az:.1e}, P {P:7.0f} yr  K="+" ".join(f"{k:+d}{NAMES[q]}" for q,k in K.items()))
