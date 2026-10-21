import numpy as np, sys
def load(fn):
    d=np.fromfile(fn,dtype=np.float64).reshape(-1,50); return d[:,0],d[:,1],d[:,2:].reshape(-1,8,6)
pre=sys.argv[1]; N=int(sys.argv[2]) if len(sys.argv)>2 else 10
t,Es,Xs=load(pre+'_ser.bin'); t2,Ep,Xp=load(pre+'_pr.bin')
n=min(len(t),len(t2)); t=t[:n]; Es=Es[:n]; Ep=Ep[:n]; Xs=Xs[:n]; Xp=Xp[:n]
rs=(Es-Es[0])/abs(Es[0]); rp=(Ep-Ep[0])/abs(Ep[0])
yr=t/(2*np.pi)
print(f"span {yr[-1]:.3g} yr, {n-1} boundaries")
print(f"serial   dE/E: final {rs[-1]:+.2e}  max {np.max(np.abs(rs)):.2e}  mean {np.mean(rs):+.2e}")
print(f"parareal dE/E: final {rp[-1]:+.2e}  max {np.max(np.abs(rp)):.2e}  mean {np.mean(rp):+.2e}")
ds=np.diff(rs); dp=np.diff(rp)
win=np.arange(1,n)%N==0
print(f"per-boundary increments: serial std {np.std(ds):.2e}; parareal std {np.std(dp):.2e} (window ends {np.std(dp[win[:len(dp)]]):.2e}, inside {np.std(dp[~win[:len(dp)]]):.2e})")
print(f"  mean increment: serial {np.mean(ds):+.2e}, parareal {np.mean(dp):+.2e}  (sem {np.std(dp)/np.sqrt(len(dp)):.1e})")
# deviation of parareal energy from serial (should be a random walk of the boundary jumps)
dd=rp-rs
print(f"parareal - serial dE/E: final {dd[-1]:+.2e}, rms {np.sqrt(np.mean(dd**2)):.2e}")
# element differences at the end
dl=(Xp[:,:,1]-Xs[:,:,1]+np.pi)%(2*np.pi)-np.pi
for f in [n//10, n//2, n-1]:
    print(f"  t={yr[f]:.3g} yr: |dlambda| "+" ".join(f"{x:.1e}" for x in np.abs(dl[f]))+"  | |d ln a| "+" ".join(f"{x:.1e}" for x in np.abs(Xp[f,:,0]-Xs[f,:,0])))
