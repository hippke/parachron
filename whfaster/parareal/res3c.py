# Round 3, Phase 3 (empirical variant): the Earth-Mars-Jupiter three-body term in the coarse step, fitted from the
# defects D_n = F(U_n) - G(U_n) of the slices preceding the window (available from converged previous windows).
# Model per planet p in {Earth, Mars}: mean ln a_p(t) = Re(A_p exp(i phi(t))), phi = -4 lam_E + 8 lam_M - 3 lam_J (- pom),
# phi(t) = phi0 + nu t. Over a slice: d ln a = Re(A (e^{i phi1} - e^{i phi0})),
# d lambda = -1.5 n T Re(A [ (e^{i phi1} - e^{i phi0})/(i nu T) - e^{i phi0} ]).
import numpy as np, sys
from ptlib import *
import ptc, ph1eval
d=np.fromfile('../pint/runs/m10_ser2.bin',dtype=np.float64).reshape(-1,50); SERX=d[:,2:].reshape(-1,8,6)
days=360000; T=days*DAY
mu=np.cumsum(mass)[:-1]+np.array(mass[1:])
G0=ptc.CoarseC(days,heun=1)
POM=int(sys.argv[1]) if len(sys.argv)>1 else -1        # planet whose pericentre enters phi (-1: none)
NH=int(sys.argv[2]) if len(sys.argv)>2 else 30
def phi_of(E):
    p=-4*E[2,1]+8*E[3,1]-3*E[4,1]
    if POM>=0: p-=np.arctan2(E[POM,3],E[POM,2])
    return p
def nu_of(E):
    lm=E[:,0].copy(); lm[:4]-=ptc.lna_disp(G0,E)
    n=np.sqrt(mu/np.exp(3*lm)); return (-4*n[2]+8*n[3]-3*n[4])*T, n*T
def basis(E):
    """per planet (Earth, Mars): rows for (d ln a, d lambda) as linear functions of (Re A, Im A)"""
    p0=phi_of(E); nu,nT=nu_of(E); e0=np.exp(1j*p0); e1=np.exp(1j*(p0+nu))
    out={}
    for c,i in enumerate((2,3)):
        ga=e1-e0                                  # d ln a = Re(A ga)
        gl=-1.5*nT[i]*((e1-e0)/(1j*nu)-e0)        # d lambda = Re(A gl)
        # Re(A g) = ReA Re g - ImA Im g
        out[i]=np.array([[ga.real,-ga.imag],[gl.real,-gl.imag]])
    return out
def fit(ep):
    X=SERX[ep-NH:ep+1]; A={}
    rows={2:[],3:[]}; rhs={2:[],3:[]}
    for q in range(NH):
        D=ediff(X[q+1],G0(X[q])); B=basis(X[q])
        for i in (2,3):
            sl=1.0/(1.5*nu_of(X[q])[1][i])           # weight lambda rows like ln a (divide by shear)
            rows[i]+=[B[i][0],B[i][1]*sl]; rhs[i]+=[D[i,0],D[i,1]*sl]
    for i in (2,3):
        c,*_=np.linalg.lstsq(np.array(rows[i]),np.array(rhs[i]),rcond=None); A[i]=c[0]+1j*c[1]
    return A
def make_G(A):
    def G(E):
        g=G0(E); B=basis(E)
        for i in (2,3):
            v=B[i]@np.array([A[i].real,A[i].imag]); g[i,0]+=v[0]; g[i,1]=wrap(g[i,1]+v[1])
        return g
    return G
if __name__=='__main__':
    EPS=[1030,2000,3400,4500,5500,6000,7000,8160,9300]
    R0=[];R1=[]
    for ep in EPS:
        A=fit(ep)
        print(f"epoch {ep}: |A| Earth {abs(A[2]):.2e} Mars {abs(A[3]):.2e}")
        R0.append(ph1eval.evaluate(G0,'config B',eps=[ep],verbose=False)[0])
        R1.append(ph1eval.evaluate(make_G(A),'+ fitted E-M-J term',eps=[ep],verbose=False)[0])
    R0=np.array(R0); R1=np.array(R1)
    print("SUMMARY k=1..3 median:",np.median(R0,0),"->",np.median(R1,0))
    print("        worst:        ",R0.max(0),"->",R1.max(0))
