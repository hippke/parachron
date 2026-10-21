# Test: coarse step augmented with the Earth-Mars-Jupiter three-body term psi = 4 lam_E - 8 lam_M + 3 lam_J + pom_S
# in the mean ln a of Mars (and Earth). Amplitudes/phases fitted, for each window, from the mean ln a at the 30 slice
# boundaries before it (as previous windows would provide). Compares parareal convergence with and without it.
import numpy as np, sys
from ptlib import *
import ptc
d=np.fromfile('../pint/runs/m10_ser2.bin',dtype=np.float64).reshape(-1,50); SERX=d[:,2:].reshape(-1,8,6)
EP=[int(x) for x in sys.argv[1].split(',')] if len(sys.argv)>1 else [1030,2000,3400,4500,5500,6000,7000,8160,9300]
N=10; days=360000; T=days*DAY; Tyr=days/365.25
G0=ptc.CoarseC(days)
mu=np.cumsum(mass)[:-1]+np.array(mass[1:])
def psi_of(E): return 4*E[2,1]-8*E[3,1]+3*E[4,1]+np.arctan2(E[5,3],E[5,2])
def fit(ep, nh=30):
    X=SERX[ep-nh:ep+1]
    m=np.array([X[q][[2,3],0]-ptc.lna_disp(G0,X[q])[[2,3]] for q in range(len(X))])
    ps=np.array([psi_of(x) for x in X])
    tt=np.arange(len(X))
    out=[]
    for c in range(2):
        A=np.vstack([np.cos(ps),np.sin(ps),np.ones_like(tt),tt]).T
        cc,*_=np.linalg.lstsq(A,m[:,c],rcond=None); out.append(cc[:2])
    nu=None
    
    return np.array(out),nu
def nu_of(E):
    lm=E[:,0].copy(); lm[:4]-=ptc.lna_disp(G0,E)
    n=np.sqrt(mu/np.exp(3*lm))          # rad per time unit, from the first-order MEAN semi-major axes
    return (4*n[2]-8*n[3]+3*n[4])*T
def make_G(coef,nu):
    def G(E):
        g=G0(E)
        p0=psi_of(E); nu=nu_of(E); p1=p0+nu
        for c,i in enumerate((2,3)):
            a,b=coef[c]
            dl=a*(np.cos(p1)-np.cos(p0))+b*(np.sin(p1)-np.sin(p0))            # change of mean ln a over the slice
            # integrated mean-motion change: int_0^T -1.5 n (lna(t)-lna(0)) dt, lna(t)-lna(0) linear in cos/sin of p0+nu s/T
            ns=np.sqrt(mu[i]/np.exp(3*E[i,0]))*T
            Ic=(np.sin(p1)-np.sin(p0))/nu; Is=-(np.cos(p1)-np.cos(p0))/nu    # (1/T) int cos, sin
            dlam=-1.5*ns*(a*(Ic-np.cos(p0))+b*(Is-np.sin(p0)))
            g[i,0]+=dl; g[i,1]=wrap(g[i,1]+dlam)
        return g
    return G
def F(E,k=2): return to_el(prop(from_el(E),days,6,kernel=k))
def para(G,E0,ref,noise,kmax=5):
    U=[E0]
    for n in range(N): U.append(G(U[-1]))
    Gold=[None]+U[1:]; out=[max(poserr(from_el(U[n]),from_el(ref[n])) for n in range(1,N+1))/noise]
    for k in range(kmax):
        Fv=[F(U[n]) for n in range(N)]; Un=[E0]; Gn=[None]
        for n in range(N):
            g=G(Un[n]); Gn.append(g); c=g+ediff(Fv[n],Gold[n+1]); c[:,1]=wrap(c[:,1]); Un.append(c)
        U=Un; Gold=Gn
        out.append(max(poserr(from_el(U[n]),from_el(ref[n])) for n in range(1,N+1))/noise)
    return out
for ep in EP:
    E0=SERX[ep]; ref=[E0]; ref1=[E0]
    for n in range(N): ref.append(F(ref[-1])); ref1.append(F(ref1[-1],1))
    noise=max(poserr(from_el(ref1[n]),from_el(ref[n])) for n in range(1,N+1))
    coef,nu=fit(ep)
    a=para(G0,E0,ref,noise); b=para(make_G(coef,nu),E0,ref,noise)
    print(f"epoch {ep:5d}: amp Mars {np.hypot(*coef[1]):.1e} Earth {np.hypot(*coef[0]):.1e}, period {2*np.pi/abs(nu_of(E0))*Tyr:.0f} yr")
    print("    default   err/noise k=0..5: "+" ".join(f"{x:.2g}" for x in a))
    print("    + 3-body  err/noise k=0..5: "+" ".join(f"{x:.2g}" for x in b)); sys.stdout.flush()
