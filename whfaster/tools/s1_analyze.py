# Compare s1_acc outputs with a reference. Prints per run: max|dE/E|, and for Mercury (1), Venus(2), Earth(3), Jupiter(5):
# phase error at the end in lambda, slope of pomega/Omega error (arcsec/yr) from linear fit.
import numpy as np, sys, glob, os
def load(f):
    d=np.loadtxt(f); t=d[:,0]; E=d[:,1]; el=d[:,2:].reshape(len(t),8,6); return t,E,el
def wrap(x): return (x+np.pi)%(2*np.pi)-np.pi
ref=sys.argv[1]; files=sys.argv[2:]
tr,Er,er=load(ref)
AS=180/np.pi*3600
print(f"{'run':14s} {'max|dE/E|':>10s} {'dlam1[rad]':>11s} {'dlam5':>10s} {'dg1[as/yr]':>11s} {'dg2':>10s} {'ds1':>10s} {'dg5':>10s} {'da1/a1':>9s}")
for f in files:
    t,E,el=load(f); n=min(len(t),len(tr)); 
    if abs(t[n-1]-tr[n-1])>1e-6*t[n-1]: print(f, "time mismatch"); continue
    d=wrap(el[:n]-er[:n]); tt=t[:n]
    def slope(x): return np.polyfit(tt,np.unwrap(x),1)[0]*AS
    print(f"{os.path.basename(f)[:-4]:14s} {np.max(np.abs(E)):10.2e} {d[-1,0,5]:11.2e} {d[-1,4,5]:10.2e} {slope(d[:,0,4]):11.2e} {slope(d[:,1,4]):10.2e} {slope(d[:,0,3]):10.2e} {slope(d[:,4,4]):10.2e} {np.mean(el[:n,0,0]-er[:n,0,0])/er[0,0,0]:9.1e}")
