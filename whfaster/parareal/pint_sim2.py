# Round 3, Phase 7: asynchronous continuously pipelined parareal ("latest-iterate" policy), event-driven simulation
# with the real numerics.
#  * slices n = 0..M-1; boundary 0 is exact; slices enter when n <= front + LOOK (front: last final boundary);
#  * coarse task C(n): new iterate of boundary n+1 = G(U_n latest) + [F(U_n^j) - G(U_n^j)] with j the latest iterate
#    of U_n whose fine solve is done (pure G prediction if none); triggered by a new U_n iterate or a new F result;
#  * fine task F(n): on the latest iterate of U_n (one at a time per slice);
#  * boundary n+1 is final when U_n is final and the fine solve of the final U_n is done (then U_{n+1} = F(U_n)), or
#    when U_n is final and the last update of U_{n+1} was <= tol (AU, positions);
#  * P cores, coarse tasks first, then lowest slice index; costs: fine 1, coarse r.
#   python3 pint_sim2.py spec r P M tol LOOK [s0]
import numpy as np, sys, heapq, time
from ptlib import *
import ptc
d=np.fromfile('../pint/runs/m10_ser2.bin',dtype=np.float64).reshape(-1,50); X=d[:,2:].reshape(-1,8,6)
spec=sys.argv[1]; r=float(sys.argv[2]); P=int(sys.argv[3]); M=int(sys.argv[4]); tol=float(sys.argv[5])
LOOK=int(sys.argv[6]); s0=int(sys.argv[7]) if len(sys.argv)>7 else 0
if spec=='B': Gf=ptc.CoarseC(360000,heun=1)
else:
    import ph1mn; Gf=ph1mn.build(spec)
def F(E): return to_el(prop(from_el(E),360000,6))
def pos(E): return from_el(E)[1:,:3]
Ui={0:[X[s0]]}            # iterates of boundary n
Gi={}                     # Gi[n][j] = G(Ui[n][j])
Fi={}                     # Fi[n] = (j, F(Ui[n][j])) latest completed fine solve
final={0:True}; front=0
lastC={}                  # lastC[n] = (len(Ui[n]), fine j used) at the last coarse update of boundary n+1
lastupd={}
running=set(); busy=[]; t=0.0; CNT=0; nF=0; nG=0; how={'F':0,'tol':0}
def ready():
    out=[]
    for n in sorted(Ui):
        if n>=M: continue
        # coarse update of boundary n+1
        if (n+1) not in final and ('C',n) not in running and n+1<=front+LOOK:
            fj=Fi[n][0] if n in Fi else -1
            state=(len(Ui[n]),fj)
            if lastC.get(n)!=state: out.append((0,n,'C'))
        # fine solve on the latest iterate
        if ('F',n) not in running:
            j=len(Ui[n])-1
            if not (n in Fi and Fi[n][0]==j): out.append((1,n,'F'))
    out.sort()
    return out
def snapshot(typ,n):
    """inputs at launch time"""
    j=len(Ui[n])-1
    if typ=='F': return (j,)
    return (j, Fi.get(n), len(Ui[n]))
def execute(typ,n,snap):
    global nF,nG,front
    if typ=='F':
        j=snap[0]
        if n in Fi and Fi[n][0]>j: return
        Fi[n]=(j,F(Ui[n][j])); nF+=1
        # predecessor final and its fine solve done -> boundary n+1 exact
        if n in final and (n+1) not in final and j==len(Ui[n])-1:
            Ui.setdefault(n+1,[]).append(Fi[n][1]); final[n+1]=True; how['F']+=1
            while front+1 in final: front+=1
    else:
        j,fi,ln=snap
        Gn=Gi.setdefault(n,{})
        if j not in Gn: Gn[j]=Gf(Ui[n][j]); nG+=1
        new=Gn[j].copy()
        if fi is not None:
            jf,f=fi
            if jf not in Gn: Gn[jf]=Gf(Ui[n][jf]); nG+=1
            new=new+ediff(f,Gn[jf]); new[:,1]=wrap(new[:,1])
        lastC[n]=(ln,fi[0] if fi is not None else -1)
        L=Ui.setdefault(n+1,[])
        upd=np.max(np.abs(pos(new)-pos(L[-1]))) if L else np.inf
        L.append(new)
        if n in final and fi is not None and j==len(Ui[n])-1 and np.max(np.abs(pos(Ui[n][fi[0]])-pos(Ui[n][j])))<=tol:
            final[n+1]=True; how['tol']+=1
            while front+1 in final: front+=1
t0=time.time()
while front<M:
    while len(busy)<P:
        rt=[q for q in ready() if (q[2],q[1]) not in running]
        if not rt: break
        q=rt[0]; key=(q[2],q[1]); running.add(key); CNT+=1
        if q[2]=='C': lastC[q[1]]=(len(Ui[q[1]]),Fi[q[1]][0] if q[1] in Fi else -1)
        heapq.heappush(busy,(t+(1.0 if q[2]=='F' else r),CNT,key,snapshot(q[2],q[1])))
    if not busy: raise RuntimeError('deadlock')
    tf,_,key,snap=heapq.heappop(busy); t=tf
    execute(key[0],key[1],snap); running.discard(key)
ref=[X[s0]]
for n in range(M): ref.append(F(ref[-1]))
ref1=[X[s0]]
for n in range(M): ref1.append(to_el(prop(from_el(ref1[-1]),360000,6,kernel=1)))
noise=np.array([max(poserr(from_el(ref1[n]),from_el(ref[n])),1e-13) for n in range(1,M+1)])
err=np.array([poserr(from_el(Ui[n][-1]),from_el(ref[n])) for n in range(1,M+1)])/noise
print(f"spec {spec[:12]} r {r} P {P} M {M} tol {tol} LOOK {LOOK}: wall {t:.1f} -> S = {M/t:.2f}; fine/slice {nF/M:.2f}"
      f" coarse/slice {nG/M:.2f}; final by F {how['F']} by tol {how['tol']}; err/noise max {err.max():.2g} median"
      f" {np.median(err):.2g} ({time.time()-t0:.0f}s)",flush=True)
