# Round 3, Phase 7: continuously pipelined ("sliding window") parareal — event-driven simulation with the real
# numerics. Slices enter the active region when they are within LOOK slices of the converged front; a slice enters
# with a coarse prediction from its predecessor's current best iterate. Iterates follow the parareal update
#   U^k_{n+1} = G(U^k_n) + F(U^{k-1}_n) - G(U^{k-1}_n);
# boundary n+1 is final at iteration k when its predecessor is final and |U^k_{n+1} - U^{k-1}_{n+1}| <= tol (position,
# AU), or automatically one iteration after its predecessor became final (then it equals F of the final state).
# P cores, list scheduling (coarse tasks first, then lowest slice index); task costs: fine 1, coarse r.
# Every task is computed for real (in simulated-time order), so the iterates are those of the schedule.
#   python3 pint_sim.py spec(B|python spec) r P M tol LOOK [s0]
import numpy as np, sys, heapq, time
from ptlib import *
import ptc
d=np.fromfile('../pint/runs/m10_ser2.bin',dtype=np.float64).reshape(-1,50); X=d[:,2:].reshape(-1,8,6)
spec=sys.argv[1]; r=float(sys.argv[2]); P=int(sys.argv[3]); M=int(sys.argv[4]); tol=float(sys.argv[5])
LOOK=int(sys.argv[6]); GATE=int(sys.argv[8]) if len(sys.argv)>8 else 0; s0=int(sys.argv[7]) if len(sys.argv)>7 else 0
if spec=='B': Gf=ptc.CoarseC(360000,heun=1)
else:
    import ph1mn; Gf=ph1mn.build(spec)
def F(E): return to_el(prop(from_el(E),360000,6))
U={0:{0:X[s0]}}; GU={}; FU={}; final={0:0}; kent={0:0}
front=0                       # largest final boundary index
t=0.0; busy=[]                # heap of (finish_time, task)
running=set(); done_tasks=set()
nF=0; nG=0
def pos(E): return from_el(E)[1:,:3]
def ready_tasks():
    out=[]
    for n in sorted(U.keys()):
        if n>=M: continue
        for k in sorted(U[n].keys()):
            # fine task F(n,k): needed unless U_n final at an earlier iteration
            if (n,k,'F') not in done_tasks and (n,k,'F') not in running:
                if not (n in final and final[n]<k):
                    if GATE and LOOK<0 and ((n//(-LOOK))*(-LOOK)) not in final: continue
                    out.append((1,n,k,'F'))
        # coarse update for boundary n+1
        if n+1>M: continue
        if n+1 in final: continue
        if n+1 not in U:
            allowed=(n+1<=front+LOOK) if LOOK>0 else (n+1<((front//(-LOOK))+2)*(-LOOK)+1)
            if allowed and (n,'E') not in running and (n,'E') not in done_tasks:
                out.append((0,n,max(U[n].keys()),'E'))
            continue
        kn=max(U[n+1].keys())+1          # next iterate of boundary n+1
        if (n,kn,'C') in running or (n,kn,'C') in done_tasks: continue
        if n in final and final[n]<=kn-1:   # predecessor final: next iterate = F(final state), exact
            if (n,final[n],'F') in done_tasks: out.append((0,n,kn,'X'))
            continue
        if kn in U[n] and (n,kn-1,'F') in done_tasks and (kn-1) in GU.get(n,{}):
            out.append((0,n,kn,'C'))
    out.sort(key=lambda q:(q[0],q[1],q[2]))
    return out
def execute(task):
    global nF,nG,front
    typ=task[-1]; n=task[1]; k=task[2]
    if typ=='F':
        FU.setdefault(n,{})[k]=F(U[n][k]); nF+=1
    elif typ=='E':
        g=Gf(U[n][k]); GU.setdefault(n,{})[k]=g; nG+=1
        U[n+1]={k:g}; kent[n+1]=k
        if n in final and final[n]<=k:      # predecessor already final: G of the final state is only a prediction
            pass
    elif typ=='X':
        U[n+1][k]=FU[n][final[n]]; final[n+1]=k; how['X']+=1
        while front+1 in final: front+=1
    elif typ=='C':
        g=Gf(U[n][k]); GU.setdefault(n,{})[k]=g; nG+=1
        c=g+ediff(FU[n][k-1],GU[n][k-1]); c[:,1]=wrap(c[:,1])
        U[n+1][k]=c
        if n in final:
            upd=np.max(np.abs(pos(c)-pos(U[n+1][k-1])))
            if final[n]<=k-1 or upd<=tol:
                final[n+1]=k; how['C']+=1
                while front+1 in final: front+=1
def cost(task): return 1.0 if task[-1]=='F' else (0.0 if task[-1]=='X' else r)
t0=time.time(); CNT=[0]; how={'X':0,'C':0}
while front<M:
    # fill idle cores
    while len(busy)<P:
        rt=[q for q in ready_tasks() if (q[1],q[2],q[3]) not in running and not (q[3]=='E' and (q[1],'E') in running)]
        if not rt: break
        q=rt[0]; key=(q[1],q[2],q[3]) if q[3]!='E' else (q[1],'E')
        running.add(key); CNT[0]+=1; heapq.heappush(busy,(t+cost(q),CNT[0],key,q))
    if not busy: raise RuntimeError("deadlock")
    tf,_,key,q=heapq.heappop(busy); t=tf
    execute(q); running.discard(key); done_tasks.add(key)
# accuracy vs serial reference (fine C kernel)
ref=[X[s0]]
for n in range(M): ref.append(F(ref[-1]))
ref1=[X[s0]]
for n in range(M): ref1.append(to_el(prop(from_el(ref1[-1]),360000,6,kernel=1)))
noise=np.array([max(poserr(from_el(ref1[n]),from_el(ref[n])),1e-13) for n in range(1,M+1)])
err=np.array([poserr(from_el(U[n][final[n]]),from_el(ref[n])) for n in range(1,M+1)])/noise
kf=np.array([final[n] for n in range(1,M+1)])
print(f"spec {spec} r {r} P {P} M {M} tol {tol} LOOK {LOOK}: wall {t:.1f} (serial {M}) -> S = {M/t:.2f};"
      f" fine solves/slice {nF/M:.2f}, coarse/slice {nG/M:.2f}; final k mean {kf.mean():.2f} max {kf.max()};"
      f" finalized by update {how['C']} / by predecessor {how['X']};"
      f" err/noise max {err.max():.2g} median {np.median(err):.2g} (real time {time.time()-t0:.0f}s)",flush=True)
