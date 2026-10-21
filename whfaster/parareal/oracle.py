# Oracle tests with the linearised parareal model (uses jac.pkl from jac.py).
# Replaces blocks of the coarse Jacobian by the fine one (+ noise) to find which accuracy a coarse model must have.
import numpy as np, pickle, sys
from ptlib import NAMES
res=pickle.load(open(sys.argv[1] if len(sys.argv)>1 else 'jac.pkl','rb'))
JF=res['JF']; N=len(JF)
dtG=float(sys.argv[2]) if len(sys.argv)>2 else 16.0
JG=res[('JG',dtG)]
rng=np.random.default_rng(1)
arow=np.arange(0,48,6)
def run(JGm, e0, kmax=6):
    e=e0; out=[]
    for k in range(kmax):
        en=[np.zeros(48)]
        for n in range(N): en.append(JGm[n]@en[n]+(JF[n]-JGm[n])@e[n])
        e=en; out.append(max(np.max(np.abs(x[1::6])) for x in e))   # max |lambda error| over planets/slices
    return out
# e0: synthetic coarse errors: lambda 1e-2, ln a 1e-7 (typical of a cheap coarse model)
e0=[np.zeros(48)]+[np.tile([1e-7,1e-2,1e-5,1e-5,1e-5,1e-5],8)*rng.standard_normal(48) for n in range(N)]
def fmt(v): return " ".join(f"{x:.0e}" for x in v)
print("baseline WH coarse:                 ", fmt(run(JG,e0)))
for eps in [1e-5,1e-6,1e-7,1e-8,1e-9]:
    J2=[]
    for n in range(N):
        M=JF[n].copy()
        M+= eps*rng.standard_normal((48,48))*0   # placeholder
        # coarse = fine everywhere, with random errors: eps in a-rows (absolute), 1e-3 relative elsewhere
        E=1e-3*np.abs(JF[n])*rng.standard_normal((48,48))
        E[arow,:]=eps*rng.standard_normal((8,48))
        J2.append(JF[n]+E)
    print(f"F' + noise: a-rows abs {eps:.0e}, others 1e-3 rel:", fmt(run(J2,e0)))
for rel in [1e-2,1e-1,0.3]:
    J2=[]
    for n in range(N):
        E=rel*np.abs(JF[n])*rng.standard_normal((48,48)); E[arow,:]=1e-8*rng.standard_normal((8,48))
        J2.append(JF[n]+E)
    print(f"F' + noise: a-rows abs 1e-8, others {rel:.0e} rel:", fmt(run(J2,e0)))
# coarse WH with a-rows replaced by fine
J2=[]
for n in range(N):
    M=JG[n].copy(); M[arow,:]=JF[n][arow,:]; J2.append(M)
print("WH coarse with exact a-rows:        ", fmt(run(J2,e0)))
