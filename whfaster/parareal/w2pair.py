# Round 3, iteration 24: same-pair second-order short-period terms (W2) of the inner planet of a pair (i, j), computed
# on the 2-D grid of (lambda_i, lambda_j). Non-canonical averaging, X = Y + W1(Y) + W2(Y):
#   f2 = sum_c dV1_i/dx_{i,c} * W1_{i,c} + 1/2 D^2 V0 [W1, W1]      (lambda: 1.125 n_i W1_lna^2)
#   W2 = homological solution of f2 - <f2> (with the Kepler shear for lambda); U2 = <f2> (second-order secular rates)
# The dependence on the perturber's displacement (of order m_i) is neglected: valid for inner-giant pairs.
# Slow-element derivatives by central differences of the grid rates (step 1e-4).
import numpy as np
import so2
from so2 import MU, m

def _sp(F, kk, n):
    with np.errstate(divide='ignore', invalid='ignore'): inv=1/(1j*kk)
    inv[np.abs(kk)<1e-300]=0
    S=F*inv[...,None]; S[...,1]+=-1.5*n*F[...,0]*inv**2
    return S

class PairW2:
    def __init__(self, E, i, j, Ng=64, n=None, h=1e-4, vbar_i=None, vbar_j=None):
        self.i, self.j, self.Ng = i, j, Ng
        n=np.sqrt(MU/np.exp(3*E[:,0])) if n is None else n
        self.n=n
        K=np.fft.fftfreq(Ng,1/Ng); kk=K[:,None]*n[i]+K[None,:]*n[j]
        fi,_=so2.pair_rates(E,i,j,Ng)                         # V1 of planet i on the grid (Ng,Ng,6)
        F=np.fft.fft2(fi,axes=(0,1))/Ng**2
        S=_sp(F,kk,n[i]); W1=np.real(np.fft.ifft2(S*Ng**2,axes=(0,1)))   # W1_i on the grid
        # derivatives of V1_i: lambdas spectrally, slow elements by central differences
        dlam_i=np.real(np.fft.ifft2(1j*K[:,None,None]*F*Ng**2,axes=(0,1)))
        f2=dlam_i*W1[...,1:2]
        vi=np.real(F[0,0]) if vbar_i is None else vbar_i
        def W1of(f):
            return np.real(np.fft.ifft2(_sp(np.fft.fft2(f,axes=(0,1))/Ng**2,kk,n[i])*Ng**2,axes=(0,1)))
        for c in (0,2,3,4,5):
            Ep=E.copy(); Em=E.copy(); Ep[i,c]+=h; Em[i,c]-=h
            a,_=so2.pair_rates(Ep,i,j,Ng); b,_=so2.pair_rates(Em,i,j,Ng)
            f2=f2+(a-b)/(2*h)*W1[...,c:c+1]
            if c>=2: f2=f2-(W1of(a)-W1of(b))/(2*h)*vi[c]          # - dW1/dx_c * Vbar_c (slow drift of i)
        # - dW1/dlambda_i * Vbar_lambda_i (secular mean-motion correction)
        f2=f2-np.real(np.fft.ifft2(1j*K[:,None,None]*S*Ng**2,axes=(0,1)))*vi[1]
        if vbar_j is not None:                                     # slow drift of the perturber
            for c in (2,3,4,5):
                Ep=E.copy(); Em=E.copy(); Ep[j,c]+=h; Em[j,c]-=h
                a,_=so2.pair_rates(Ep,i,j,Ng); b,_=so2.pair_rates(Em,i,j,Ng)
                f2=f2-(W1of(a)-W1of(b))/(2*h)*vbar_j[c]
        f2=f2.copy(); f2[...,1]+=1.125*n[i]*W1[...,0]**2
        F2=np.fft.fft2(f2,axes=(0,1))/Ng**2
        self.U2=np.real(F2[0,0])
        self.S2=_sp(F2,kk,n[i]); self.S2[0,0]=0
        self.K=K
    def disp(self, E):
        """W2 of planet i (6 comps) at the angles of E"""
        ph=np.exp(1j*(self.K[:,None]*E[self.i,1]+self.K[None,:]*E[self.j,1]))
        return np.real(np.einsum('abe,ab->e',self.S2,ph))
