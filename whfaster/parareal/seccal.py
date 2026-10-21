# Round 3, it. 23: empirical secular-matrix correction. Per-slice z = k+ih and zeta = q+ip errors of the inner planets
# are linear in all planets' z (zeta) vectors: dz_p = sum_q A_pq z_q (R^2 = 1.00 at 4.5 Myr). Fit A from Nfit slices,
# test the residual on other stretches.
import numpy as np, sys
d=np.fromfile('../pint/runs/m10_ser2.bin',dtype=np.float64).reshape(-1,50); X=d[:,2:].reshape(-1,8,6)
def design(s0,M,c):
    return X[s0:s0+M,:,c]+1j*X[s0:s0+M,:,c+1]
def fit(err,s0,c,cols=8):
    Z=design(s0,len(err),c)[:,:cols]; A=np.zeros((4,cols),complex)
    for p in range(4):
        y=err[:,p,c]+1j*err[:,p,c+1]; A[p],*_=np.linalg.lstsq(Z,y,rcond=None)
    return A
def resid(err,s0,c,A):
    Z=design(s0,len(err),c)[:,:A.shape[1]]
    out=[]
    for p in range(4):
        y=err[:,p,c]+1j*err[:,p,c+1]; out.append((np.sqrt(np.mean(abs(y)**2)),np.sqrt(np.mean(abs(y-Z@A[p])**2))))
    return out
if __name__=='__main__':
    E={s:np.load(f'serr_{s}.npy') for s in (1000,4500,6000)}
    for c,lab in ((2,'z'),(4,'zeta')):
        for sf in (1000,4500,6000):
            for nf in (200,30,10):
                A=fit(E[sf][:nf],sf,c)
                row=[]
                for st in (1000,4500,6000):
                    r=resid(E[st],st,c,A); row.append(f"{st}: "+" ".join(f"{a:.1e}->{b:.1e}" for a,b in r[1:]))
                print(f"{lab} fit {sf}[{nf}] | "+" | ".join(row))
