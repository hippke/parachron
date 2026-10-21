# Round 3, it. 24: is the Mars short-period residual a same-pair second-order (Mars-Jupiter x Mars-Jupiter) effect?
# Three-body Sun+Mars+Jupiter (and variants with scaled Jupiter mass): first-order mean ln a residual of Mars.
import numpy as np, rebound, sys, meanel, itertools
from ptlib import *
d=np.load(sys.argv[1]) if len(sys.argv)>1 else None
D=np.load('phase0_jac.npz'); E8=D['ref4500'][0]
def run(mJscale, yrs=200, Ng=64):
    m=np.array([mass[0]+sum(mass[1:4]), mass[4], mass[5]*mJscale])
    E=np.array([E8[3],E8[4]])
    s=rebound.Simulation(); s.add(m=m[0])
    for mm,e in zip(m[1:],E):
        s.add(m=mm,a=np.exp(e[0]),l=e[1],e=np.hypot(e[2],e[3]),pomega=np.arctan2(e[3],e[2]),inc=2*np.arcsin(np.hypot(e[4],e[5])),Omega=np.arctan2(e[5],e[4]))
    s.move_to_com(); s.integrator="ias15"
    def els():
        out=[]
        for o in s.orbits():
            si=np.sin(o.inc/2); out.append([np.log(o.a),o.l,o.e*np.cos(o.pomega),o.e*np.sin(o.pomega),si*np.cos(o.Omega),si*np.sin(o.Omega)])
        return np.array(out)
    X0=els(); th=meanel.Theory(m,X0,Ng=Ng,gr=False); M=th.mean(X0); X1=X0.copy(); X1[:,0]=M[:,0]
    th=meanel.Theory(m,X1,Ng=Ng,gr=False)
    T=[]; ML=[]
    for q in range(int(yrs*365.25/12)):
        s.integrate(s.t+12*DAY); X=els(); T.append(s.t/(2*np.pi)); ML.append(X[0,0]-th.disp(X)[0,0])
    T=np.array(T); ML=np.array(ML); x=ML-np.polyval(np.polyfit(T,ML,2),T)
    return x.std()
for sc in (1.0,0.5,0.25):
    print(f"Jupiter mass x{sc}: Mars first-order mean ln a residual rms {run(sc):.2e}",flush=True)
