import numpy as np, sys
def load(fn): d=np.fromfile(fn,dtype=np.float64).reshape(-1,50); return d[:,0]/(2*np.pi), d[:,2:].reshape(-1,8,6)
t,S=load('runs/m10_ser2.bin')
for fn in sys.argv[1:]:
    _,X=load(fn); n=min(len(X),len(S))
    out=[]
    for i in (2000,5000,n-1):
        dz=[abs((X[i,p,2]+1j*X[i,p,3])-(S[i,p,2]+1j*S[i,p,3])) for p in range(4)]
        dq=[abs((X[i,p,4]+1j*X[i,p,5])-(S[i,p,4]+1j*S[i,p,5])) for p in range(4)]
        out.append(f"{t[i]/1e6:.0f} Myr |dz| "+" ".join(f"{x:.0e}" for x in dz)+" |dzeta| "+" ".join(f"{x:.0e}" for x in dq))
    print(f"{fn:22s} "+" ; ".join(out))
