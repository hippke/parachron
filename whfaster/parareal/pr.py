# Parareal prototype around WHFast512 (C kernel). F: dt_F (default 6 d), GR, corrector 17.
# G: WHFast512 with dt_G (and optional corrector). Reports, per iteration k, the max position error
# (AU) at slice boundaries vs the serial fine solution, and the reference noise level
# (asm kernel 1 vs C kernel 2, i.e. pure round-off divergence of the fine integrator).
import rebound, numpy as np, sys, time, argparse
DAY = 2*np.pi/365.25
p = argparse.ArgumentParser()
p.add_argument('--slices', type=int, default=10)
p.add_argument('--slice_days', type=float, default=360000)
p.add_argument('--dtF', type=float, default=6)
p.add_argument('--dtG', type=float, default=24)
p.add_argument('--corrG', type=int, default=17)
p.add_argument('--kmax', type=int, default=8)
p.add_argument('--space', default='elem')
a = p.parse_args()
pos=[[-0.008816286905115728,-0.0010954664916791675,0.0002143249385447027],[-0.05942272929227954,-0.46308699693348293,-0.032897989948949075],[-0.7276101005375593,0.006575003332463933,0.041795901908847084],[-0.5789530452882667,-0.8361530119313055,0.0002611520181901174],[-1.45202492400084,0.827519404194876,0.052981833432457694],[4.492983939852296,2.0661626247490354,-0.10909246996001629],[8.4974210980544,-4.8620394993693585,-0.2537835862373596],[12.959111916929283,14.760785302864473,-0.1130656917933948],[29.787987348666505,-2.51460654509393,-0.6347108842010732]]
vel=[[0.00014315746073017681,-0.0004912441820893999,8.127678560998346e-07],[1.2978664284760637,-0.09524541469911743,-0.12677574364801253],[-0.019239782390457125,-1.1813975672919448,-0.01509392594251431],[0.8098712561282222,-0.5682496529341624,2.6169897281383047e-05],[-0.37436417754222295,-0.6365841544564991,-0.004143932260467942],[-0.18818907783656452,0.41919544951404614,0.0024710497024424977],[0.14292308496870448,0.2808676923735748,-0.010574288572728728],[-0.1734971049470612,0.14019515029516152,0.0027683484887051457],[0.014142947617173336,0.18292110872737416,-0.004092845767710294]]
mass=[0.9999999999950272,1.6601208254808336e-07,2.447838287784771e-06,3.0404326489511185e-06,3.2271560828978514e-07,0.0009547919099366768,0.0002858856700231729,4.366249613200406e-05,5.151383772628957e-05]
def initial():
    s=rebound.Simulation()
    for m,x,v in zip(mass,pos,vel): s.add(m=m,x=x[0],y=x[1],z=x[2],vx=v[0],vy=v[1],vz=v[2])
    s.move_to_com()
    return np.array([[q.x,q.y,q.z,q.vx,q.vy,q.vz] for q in s.particles])
def prop(U, days, dt, corr, kernel=2):
    s=rebound.Simulation(); 
    for m,u in zip(mass,U): s.add(m=m,x=u[0],y=u[1],z=u[2],vx=u[3],vy=u[4],vz=u[5])
    s.integrator="whfast512"; s.dt=dt*DAY; s.exact_finish_time=0
    s.integrator.gr_potential=1; s.integrator.corrector=corr; s.integrator.kernel=kernel
    n=int(round(days/dt)); assert abs(n*dt-days)<1e-9
    s.integrator.concatenate_steps=n
    s.steps(1)
    return np.array([[q.x,q.y,q.z,q.vx,q.vy,q.vz] for q in s.particles])
def to_el(U):
    s=rebound.Simulation()
    for m,u in zip(mass,U): s.add(m=m,x=u[0],y=u[1],z=u[2],vx=u[3],vy=u[4],vz=u[5])
    E=[]
    for o in s.orbits():   # Jacobi orbits
        si=np.sin(o.inc/2)
        E.append([o.a, o.l, o.e*np.cos(o.pomega), o.e*np.sin(o.pomega), si*np.cos(o.Omega), si*np.sin(o.Omega)])
    return np.array(E)
def from_el(E, U_like):
    s=rebound.Simulation(); s.add(m=mass[0])
    for m,e in zip(mass[1:],E):
        ecc=np.hypot(e[2],e[3]); pom=np.arctan2(e[3],e[2]); si=np.hypot(e[4],e[5]); Om=np.arctan2(e[5],e[4])
        s.add(m=m, a=e[0], l=e[1], e=ecc, pomega=pom, inc=2*np.arcsin(si), Omega=Om)
    s.move_to_com()
    return np.array([[q.x,q.y,q.z,q.vx,q.vy,q.vz] for q in s.particles])
def wrap(x): return (x+np.pi)%(2*np.pi)-np.pi
def correct(g, f, gold):
    if a.space=='cart': return g+f-gold
    Eg,Ef,Eo=to_el(g),to_el(f),to_el(gold)
    D=Ef-Eo; D[:,1]=wrap(D[:,1])
    return from_el(Eg+D, g)
def err(A,B): return np.max(np.abs(A[1:,:3]-B[1:,:3]))
N=a.slices; T=a.slice_days
U0=initial()
t0=time.time(); ref=[U0]; ref1=[U0]
for n in range(N): ref.append(prop(ref[-1],T,a.dtF,17,2)); ref1.append(prop(ref1[-1],T,a.dtF,17,1))
tF=(time.time()-t0)/2/N
noise=[err(ref[n],ref1[n]) for n in range(N+1)]
t0=time.time(); Gt=prop(U0,T,a.dtG,a.corrG); tG=time.time()-t0
print(f"slices={N} x {T/365.25:.0f} yr, dtF={a.dtF} dtG={a.dtG} corrG={a.corrG}; cost F/slice={tF:.3f}s G/slice={tG:.3f}s r=G/F={tG/tF:.3f}")
print("round-off noise (asm vs C fine) at slice ends:", " ".join(f"{x:.1e}" for x in noise[1:]))
U=[U0]
for n in range(N): U.append(prop(U[-1],T,a.dtG,a.corrG))
Gold=[None]+U[1:]   # Gold[n+1]=G(U[n])
print(f"k=0 maxerr {max(err(U[n],ref[n]) for n in range(N+1)):.2e}  final {err(U[N],ref[N]):.2e}")
for k in range(1,a.kmax+1):
    Fv=[prop(U[n],T,a.dtF,17) for n in range(N)]
    Unew=[U0]; Gnew=[None]
    for n in range(N):
        g=prop(Unew[n],T,a.dtG,a.corrG); Gnew.append(g)
        Unew.append(correct(g,Fv[n],Gold[n+1]))
    U=Unew; Gold=Gnew
    e=[err(U[n],ref[n]) for n in range(N+1)]
    conv=sum(1 for n in range(N+1) if e[n]<=max(noise[n],1e-15))
    print(f"k={k} maxerr {max(e):.2e}  final {e[N]:.2e}  slices within noise: {conv}/{N+1}  per-slice: "+" ".join(f"{x:.0e}" for x in e[1:]))
