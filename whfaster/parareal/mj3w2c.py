# Round 3, it. 24: validation of w2pair, three-body Sun+Mars+Jupiter, tables rebuilt every 20-yr segment
import numpy as np, rebound, meanel, w2pair, sys
from ptlib import *
D=np.load('phase0_jac.npz'); E8=D['ref4500'][0]
sc=float(sys.argv[1]) if len(sys.argv)>1 else 1.0
m=np.array([mass[0]+sum(mass[1:4]), mass[4], mass[5]*sc])
masses8=np.array(mass); 
s=rebound.Simulation(); s.add(m=m[0])
for mm,e in zip(m[1:],(E8[3],E8[4])):
    s.add(m=mm,a=np.exp(e[0]),l=e[1],e=np.hypot(e[2],e[3]),pomega=np.arctan2(e[3],e[2]),inc=2*np.arcsin(np.hypot(e[4],e[5])),Omega=np.arctan2(e[5],e[4]))
s.move_to_com(); s.integrator="ias15"
def els():
    out=[]
    for o in s.orbits():
        si=np.sin(o.inc/2); out.append([np.log(o.a),o.l,o.e*np.cos(o.pomega),o.e*np.sin(o.pomega),si*np.cos(o.Omega),si*np.sin(o.Omega)])
    return np.array(out)
w2pair.so2.m=np.array(mass); w2pair.so2.m[5]*=sc
res={'W1':[],'W1+W2':[]}
for seg in range(6):
    X0=els(); th=meanel.Theory(m,X0,Ng=64,gr=False); M=th.mean(X0); X1=X0.copy(); X1[:,0]=M[:,0]
    th=meanel.Theory(m,X1,Ng=64,gr=False)
    if len(sys.argv)>2:
        for _ in range(3):
            Y=X0.copy()
            for __ in range(3): Y=X0-th.disp(Y)
            X1=Y.copy(); X1[:,1]=X0[:,1]; th=meanel.Theory(m,X1,Ng=64,gr=False)
    E8b=E8.copy(); E8b[3]=X1[0]; E8b[4]=X1[1]
    w2=w2pair.PairW2(E8b,3,4,Ng=64,n=np.sqrt(w2pair.MU/np.exp(3*E8b[:,0])))
    T=[];A=[];C=[];LM=[];LJ=[]
    for q in range(int(20*365.25/6)):
        s.integrate(s.t+6*DAY); X=els(); T.append(s.t/(2*np.pi)); LM.append(X[0,1]); LJ.append(X[1,1])
        Y=X.copy()
        for _ in range(3): Y=X-th.disp(Y)
        A.append(Y[0,0])
        Y8=E8b.copy(); Y8[3]=Y[0]; Y8[4]=Y[1]; C.append(Y[0,0]-w2.disp(Y8)[0])
    T=np.array(T)
    for lab,v in (('W1',A),('W1+W2',C)):
        v=np.array(v); x=v-np.polyval(np.polyfit(T,v,1),T); res[lab].append(x.std())
    if seg==0:
        x=np.array(C); x=x-np.polyval(np.polyfit(T,x,1),T); LM=np.array(LM); LJ=np.array(LJ); out=[]
        for a in range(-8,9):
            for b in range(-12,13):
                if a<0 or (a==0 and b<=0): continue
                ph=a*LM+b*LJ; c=np.cos(ph); s_=np.sin(ph)
                amp=2*np.hypot(np.mean(x*c),np.mean(x*s_)); out.append((amp,a,b))
        out.sort(reverse=True); print("residual harmonics (kM,kJ):",[(a,b,f"{amp:.1e}") for amp,a,b in out[:10]])
        break
for k,v in res.items(): print(f"Jupiter x{sc}: {k}: residual rms per 20-yr segment",np.array2string(np.array(v),precision=2))
