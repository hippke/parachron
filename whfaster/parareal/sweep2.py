# Parameter search for the parareal configuration (round 2, iteration 15).
# For each configuration: error/noise (max over boundaries) after k = 1..kmax iterations, worst over several start
# epochs taken from the 10 Myr serial run; measured coarse cost r = T_G/T_F; modelled 10-core speedup at the
# smallest K with error/noise <= target (default 0.3, the quality of K = 4 in the default configuration, which passed
# the 10 Myr secular-frequency test).
#   python3 sweep2.py "name:key=val,key=val" ...   keys: days, dt_out, ng_ii, ng_io, maxj, heun, N
import numpy as np, sys, time, json
from ptlib import *
import ptc
SER='../pint/runs/m10_ser2.bin'
EPOCHS=[0, 1030, 2000, 3400, 4500, 5500, 6000, 7000, 8160, 9300]
TARGET=float([a for a in sys.argv if a.startswith('--target=')][0].split('=')[1]) if any(a.startswith('--target=') for a in sys.argv) else 1.0
KMAX=6
VERB='--v' in sys.argv
d=np.fromfile(SER,dtype=np.float64).reshape(-1,50); SERX=d[:,2:].reshape(-1,8,6)
def F(E,days,kernel=2): return to_el(prop(from_el(E),days,6,kernel=kernel))
def run(cfg):
    days=cfg.get('days',360000); N=int(cfg.get('N',10))
    G=ptc.CoarseC(days,dt_out=cfg.get('dt_out',360),ng_ii=int(cfg.get('ng_ii',32)),ng_io=int(cfg.get('ng_io',16)),
                  maxj=int(cfg.get('maxj',5)),heun=int(cfg.get('heun',1)))
    worst=np.zeros(KMAX); tF=[]; tG=[]; per={}
    for ep in EPOCHS:
        E0=SERX[ep]; ref=[E0]; ref1=[E0]
        for n in range(N):
            t0=time.time(); ref.append(F(ref[-1],days)); tF.append(time.time()-t0)
            ref1.append(F(ref1[-1],days,1))
        noise=max(poserr(from_el(ref1[n]),from_el(ref[n])) for n in range(1,N+1))
        U=[E0]
        for n in range(N):
            t0=time.time(); U.append(G(U[-1])); tG.append(time.time()-t0)
        Gold=[None]+U[1:]
        for k in range(KMAX):
            Fv=[F(U[n],days) for n in range(N)]; Un=[E0]; Gn=[None]
            for n in range(N):
                g=G(Un[n]); Gn.append(g); c=g+ediff(Fv[n],Gold[n+1]); c[:,1]=wrap(c[:,1]); Un.append(c)
            U=Un; Gold=Gn
            e=max(poserr(from_el(U[n]),from_el(ref[n])) for n in range(1,N+1))/noise
            worst[k]=max(worst[k],e); per.setdefault(ep,[]).append(e)
    r=np.median(tG)/np.median(tF)
    Ks=[]
    for ep in per:
        ok=[k for k,x in enumerate(per[ep]) if x<=TARGET]; Ks.append(ok[0]+1 if ok else KMAX+1)
    K=max(Ks); Km=np.mean(Ks)
    S=(N/(N*r+Km*(1+r)), N/(N*r+K*(1+r)))
    K=(Ks,Km)
    if VERB:
        for ep in per: print(f'    epoch {ep}: '+' '.join(f'{x:.2g}' for x in per[ep]))
    return worst,r,K,S
print(f"target error/noise <= {TARGET}; epochs {EPOCHS}")
for spec in [a for a in sys.argv[1:] if not a.startswith('--')]:
    name,_,kv=spec.partition(':')
    cfg={k:float(v) for k,v in (x.split('=') for x in kv.split(',') if x)}
    worst,r,K,S=run(cfg)
    print(f"{name:26s} r={r:.3f} worst err/noise k=1..: "+" ".join(f"{x:.2g}" for x in worst)+f" | K per epoch {K[0]} mean {K[1]:.1f} -> S(mean K)={S[0]:.2f} S(max K)={S[1]:.2f}"); sys.stdout.flush()
