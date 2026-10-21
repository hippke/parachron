# Round 3, Phase 0: required accuracy of the coarse map for K = 2.
# J_F, J_G (config B, C coarse) at the 10 slices of windows starting at 10 epochs of the 10 Myr serial run; e0 = actual
# coarse errors. Linearised parareal with J_G(eps) = J_F + eps (J_G - J_F) on selected rows and e0 * eps0.
#   python3 phase0.py compute      -> phase0_jac.npz
#   python3 phase0.py analyse
import numpy as np, sys
from ptlib import *
import ptc
EP=[0,1030,2000,3400,4500,5500,6000,7000,8160,9300]; N=10; days=360000
H=np.array([1e-7,1e-4,1e-5,1e-5,1e-5,1e-5])
d=np.fromfile('../pint/runs/m10_ser2.bin',dtype=np.float64).reshape(-1,50); SERX=d[:,2:].reshape(-1,8,6)
G=ptc.CoarseC(days)
def F(E,k=2): return to_el(prop(from_el(E),days,6,kernel=k))
def jac(f,E):
    J=np.zeros((48,48))
    for i in range(8):
        for c in range(6):
            Ep=E.copy(); Em=E.copy(); Ep[i,c]+=H[c]; Em[i,c]-=H[c]
            J[:,i*6+c]=ediff(f(Ep),f(Em)).ravel()/(2*H[c])
    return J
if sys.argv[1]=='compute':
    out={}
    for ep in EP:
        E0=SERX[ep]; ref=[E0]; ref1=[E0]
        for n in range(N): ref.append(F(ref[-1])); ref1.append(F(ref1[-1],1))
        noise=max(poserr(from_el(ref1[n]),from_el(ref[n])) for n in range(1,N+1))
        U=[E0]
        for n in range(N): U.append(G(U[-1]))
        e0=np.array([ediff(U[n],ref[n]).ravel() for n in range(N+1)])
        JF=np.array([jac(F,ref[n]) for n in range(N)]); JG=np.array([jac(G,ref[n]) for n in range(N)])
        out[f'ref{ep}']=np.array(ref); out[f'JF{ep}']=JF; out[f'JG{ep}']=JG; out[f'e0{ep}']=e0; out[f'noise{ep}']=noise
        print('epoch',ep,'done'); sys.stdout.flush()
    np.savez('phase0_jac.npz',**out)
else:
    D=np.load('phase0_jac.npz')
    def run(ep, eps, eps0, rows, kmax=3):
        JF=D[f'JF{ep}']; JG=D[f'JG{ep}']; e0=D[f'e0{ep}']*eps0; ref=D[f'ref{ep}']; noise=float(D[f'noise{ep}'])
        JGe=JG.copy()
        for n in range(N): JGe[n][rows,:]=JF[n][rows,:]+eps*(JG[n][rows,:]-JF[n][rows,:])
        e=list(e0); res=[]
        for k in range(kmax):
            en=[np.zeros(48)]
            for n in range(N): en.append(JGe[n]@en[n]+(JF[n]-JGe[n])@e[n])
            e=en
            res.append(max(poserr(from_el(ref[n]+e[n].reshape(8,6)),from_el(ref[n])) for n in range(1,N+1))/noise)
        return res
    inner_all=[p*6+c for p in range(4) for c in range(6)]
    inner_al=[p*6+c for p in range(4) for c in (0,1)]
    allrows=list(range(48))
    print("linearised k = 1, 2, 3 error / noise; worst and median over the 10 epochs")
    for rows,lab in ((inner_all,'inner rows (all elements)'),(inner_al,'inner rows (ln a, lambda)'),(allrows,'all rows')):
        for eps,eps0 in ((1,1),(0.3,1),(0.1,1),(0.1,0.1),(0.03,0.1),(0.03,0.03),(0.01,0.1),(0.01,0.01),(3e-3,3e-3),(1e-3,1e-3)):
            R=np.array([run(ep,eps,eps0,rows) for ep in EP])
            print(f"  {lab:27s} eps={eps:<6g} eps0={eps0:<6g}: k=2 worst {R[:,1].max():8.2g} median {np.median(R[:,1]):8.2g} | k=3 worst {R[:,2].max():.2g}  | epochs with k2<=2: {np.sum(R[:,1]<=2)}/10")
