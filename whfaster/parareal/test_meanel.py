# Validation of meanel.py: (1) element conversions vs REBOUND, (2) pair decomposition of the Jacobi perturbing
# accelerations vs a direct N-body sum, (3) residual short-period content of the mean elements along a fine run.
import numpy as np, rebound, sys, time
from ptlib import *
import meanel
m=np.array(mass)
U=initial()
E=to_el(U)
Mint=np.cumsum(m)[:-1]; mu=Mint+m[1:]
# (1) Jacobi r,v from elements vs REBOUND
s=sim_from(U)
jac=[]
for i in range(1,9):
    com=s.com(last=i)  # com of particles 0..i-1
    p=s.particles[i]
    jac.append([p.x-com.x,p.y-com.y,p.z-com.z,p.vx-com.vx,p.vy-com.vy,p.vz-com.vz])
jac=np.array(jac)
r,v=meanel.el_to_rv(mu,E)
print("(1) max |r - r_jac|:",np.max(np.abs(r-jac[:,:3])),"  |v - v_jac|:",np.max(np.abs(v-jac[:,3:])))
E2=meanel.rv_to_el(mu,jac[:,:3],jac[:,3:])
D=E2-E; D[:,1]=wrap(D[:,1]); print("    max |E(rv) - E(rebound)|:",np.max(np.abs(D)))
# (2) perturbing accelerations
X=U[:,:3]
acc=np.zeros_like(X)
for i in range(9):
    for j in range(9):
        if i!=j:
            d=X[j]-X[i]; acc[i]+=m[j]*d/np.linalg.norm(d)**3
ajac=[]
for i in range(1,9):
    com_acc=np.sum(m[:i,None]*acc[:i],0)/Mint[i-1]
    ajac.append(acc[i]-com_acc)
ajac=np.array(ajac)
apert=ajac+mu[:,None]*r/np.linalg.norm(r,axis=1)[:,None]**3
asum=np.zeros_like(apert)
for i in range(8):
    for j in range(i+1,8):
        ri,rj=r[i],r[j]; d=rj-ri; d3=np.linalg.norm(d)**3; rj3=np.linalg.norm(rj)**3; rj2=rj@rj
        asum[i]+=m[j+1]*(d/d3-rj/rj3)
        asum[j]+=-(Mint[j]+m[j+1])/Mint[j]*m[i+1]*(d/d3-rj/rj3+ri/rj3-3*rj*(rj@ri)/rj2/rj3)
print("(2) |a_pert| and relative error of the pair sum per planet:")
for i in range(8): print(f"    {NAMES[i]} {np.linalg.norm(apert[i]):.3e}  {np.linalg.norm(asum[i]-apert[i])/np.linalg.norm(apert[i]):.1e}")
# (3) fine run, sampled
Ng=int(sys.argv[1]) if len(sys.argv)>1 else 64
yrs=float(sys.argv[2]) if len(sys.argv)>2 else 30
t0=time.time(); th=meanel.Theory(m,E,Ng=Ng); print(f"(3) theory setup Ng={Ng}: {time.time()-t0:.1f}s")
s=sim_from(U); s.integrator="whfast512"; s.dt=6*DAY; s.integrator.gr_potential=1; s.integrator.corrector=17
s.exact_finish_time=0; s.integrator.concatenate_steps=1
nsamp=int(yrs*365.25/6)
Es=[];Ms=[]; t=[]
for k in range(nsamp):
    s.steps(1)
    Ek=to_el(state(s)); Es.append(Ek); Ms.append(th.mean(Ek)); t.append(s.t)
Es=np.array(Es); Ms=np.array(Ms); t=np.array(t)
def resid(x):
    x=np.unwrap(x) if False else x
    A=np.vstack([np.ones_like(t),t,t*t]).T
    c,*_=np.linalg.lstsq(A,x,rcond=None); return x-A@c
print("    planet   std(ln a): osc  mean   | std(lambda): osc  mean  | std(k): osc mean")
for i in range(8):
    la=np.unwrap(Es[:,i,1]); lm=np.unwrap(Ms[:,i,1])
    print(f"    {NAMES[i]}  {np.std(resid(Es[:,i,0])):.1e} {np.std(resid(Ms[:,i,0])):.1e}   |  {np.std(resid(la)):.1e} {np.std(resid(lm)):.1e}  |  {np.std(resid(Es[:,i,2])):.1e} {np.std(resid(Ms[:,i,2])):.1e}")
np.savez('meanel_test.npz',t=t,Es=Es,Ms=Ms)
