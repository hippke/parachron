# Round 3, it. 24: validation of w2pair on the three-body Sun+Mars+Jupiter problem
import numpy as np, rebound, meanel, w2pair
from ptlib import *
D=np.load('phase0_jac.npz'); E8=D['ref4500'][0]
m=np.array([mass[0]+sum(mass[1:4]), mass[4], mass[5]])
s=rebound.Simulation(); s.add(m=m[0])
for mm,e in zip(m[1:],(E8[3],E8[4])):
    s.add(m=mm,a=np.exp(e[0]),l=e[1],e=np.hypot(e[2],e[3]),pomega=np.arctan2(e[3],e[2]),inc=2*np.arcsin(np.hypot(e[4],e[5])),Omega=np.arctan2(e[5],e[4]))
s.move_to_com(); s.integrator="ias15"
def els():
    out=[]
    for o in s.orbits():
        si=np.sin(o.inc/2); out.append([np.log(o.a),o.l,o.e*np.cos(o.pomega),o.e*np.sin(o.pomega),si*np.cos(o.Omega),si*np.sin(o.Omega)])
    return np.array(out)
X0=els(); th=meanel.Theory(m,X0,Ng=64,gr=False); M=th.mean(X0); X1=X0.copy(); X1[:,0]=M[:,0]
import sys
if len(sys.argv)>1 and sys.argv[1]=='allmean':
    for _ in range(3):
        M=th.mean(X0); X1=M.copy(); X1[:,1]=X0[:,1]; th=meanel.Theory(m,X1,Ng=64,gr=False)
th=meanel.Theory(m,X1,Ng=64,gr=False)
E8b=E8.copy(); E8b[3]=X1[0]; E8b[4]=X1[1]
n8=np.sqrt(w2pair.MU/np.exp(3*E8b[:,0]))
w2=w2pair.PairW2(E8b,3,4,Ng=64,n=n8)
print("U2 (second-order secular rates of Mars):",w2.U2)
T=[];A=[];B=[];C=[]
for q in range(int(200*365.25/12)):
    s.integrate(s.t+12*DAY); X=els(); T.append(s.t/(2*np.pi))
    A.append(X[0,0]-th.disp(X)[0,0])
    Y=X.copy()
    for _ in range(4): Y=X-th.disp(Y)
    B.append(Y[0,0])
    Y=X.copy()
    for _ in range(4):
        Y8=E8b.copy(); Y8[3]=Y[0]; Y8[4]=Y[1]
        Y=X-th.disp(Y); Y[0]-=w2.disp(Y8)
    C.append(Y[0,0])
T=np.array(T)
for lab,v in (('W1 at osc',A),('W1 fixed point',B),('W1+W2 fixed point',C)):
    v=np.array(v); x=v-np.polyval(np.polyfit(T,v,2),T); print(f"{lab}: Mars mean ln a residual rms {x.std():.2e}")
