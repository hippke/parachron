# Round 3: which Jacobian blocks of the coarse map matter? Linearised parareal (Phase 0 data, config B) with
# selected (row-element x column-element) blocks of J_G replaced by J_F, e0 replaced by 0 for the same rows.
import numpy as np, itertools
from ptlib import *
D=np.load('phase0_jac.npz'); EP=[0,1030,2000,3400,4500,5500,6000,7000,8160,9300]; N=10
def run(ep, rows, cols, fix_e0=True, kmax=3):
    JF=D[f'JF{ep}']; JG=D[f'JG{ep}'].copy(); e0=D[f'e0{ep}'].copy(); ref=D[f'ref{ep}']; noise=float(D[f'noise{ep}'])
    for n in range(N): JG[n][np.ix_(rows,cols)]=JF[n][np.ix_(rows,cols)]
    if fix_e0: e0[:,rows]=0
    e=list(e0); res=[]
    for k in range(kmax):
        en=[np.zeros(48)]
        for n in range(N): en.append(JG[n]@en[n]+(JF[n]-JG[n])@e[n])
        e=en; res.append(max(poserr(from_el(ref[n]+e[n].reshape(8,6)),from_el(ref[n])) for n in range(1,N+1))/noise)
    return res
EL={'lna':[0],'lam':[1],'kh':[2,3],'qp':[4,5]}
def idx(planets,els): return [p*6+c for p in planets for e in els for c in EL[e]]
def summ(rows,cols,lab,fix=True):
    R=np.array([run(ep,rows,cols,fix) for ep in EP]); print(f"{lab:45s} k2 median {np.median(R[:,1]):8.3g} worst {R[:,1].max():8.3g}  k3 median {np.median(R[:,2]):6.3g}")
allc=list(range(48)); VEM=[1,2,3]
summ([],[], 'config B (linear)')
summ(idx(VEM,EL),allc,'V,E,M all rows exact (+e0)')
summ(idx(VEM,EL),allc,'V,E,M all rows Jacobian only',False)
summ([],[],'e0 of VEM rows = 0 only') 
R=np.array([run(ep,[],[],True) for ep in EP])
for r in EL:
    summ(idx(VEM,[r]),allc,f'VEM {r} rows exact')
    summ(idx(VEM,[x for x in EL if x!=r]),allc,f'VEM all rows but {r} exact')
for c in EL:
    summ(idx(VEM,EL),idx(range(8),[c]),f'VEM rows, cols {c} (all planets) exact')
for p in range(8):
    summ(idx(VEM,EL),idx([p],EL),f'VEM rows, cols of planet {NAMES[p]} exact')
