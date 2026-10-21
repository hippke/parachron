# Shared helpers for the parallel-in-time experiments (whfaster round 2).
# State representations:
#   U : (9,6) barycentric Cartesian (x,y,z,vx,vy,vz), index 0 = Sun
#   E : (8,6) Jacobi elements per planet: [ln a, lambda, e cos pomega, e sin pomega, sin(i/2) cos Omega, sin(i/2) sin Omega]
import rebound, numpy as np
DAY = 2*np.pi/365.25
YR = 2*np.pi
pos=[[-0.008816286905115728,-0.0010954664916791675,0.0002143249385447027],[-0.05942272929227954,-0.46308699693348293,-0.032897989948949075],[-0.7276101005375593,0.006575003332463933,0.041795901908847084],[-0.5789530452882667,-0.8361530119313055,0.0002611520181901174],[-1.45202492400084,0.827519404194876,0.052981833432457694],[4.492983939852296,2.0661626247490354,-0.10909246996001629],[8.4974210980544,-4.8620394993693585,-0.2537835862373596],[12.959111916929283,14.760785302864473,-0.1130656917933948],[29.787987348666505,-2.51460654509393,-0.6347108842010732]]
vel=[[0.00014315746073017681,-0.0004912441820893999,8.127678560998346e-07],[1.2978664284760637,-0.09524541469911743,-0.12677574364801253],[-0.019239782390457125,-1.1813975672919448,-0.01509392594251431],[0.8098712561282222,-0.5682496529341624,2.6169897281383047e-05],[-0.37436417754222295,-0.6365841544564991,-0.004143932260467942],[-0.18818907783656452,0.41919544951404614,0.0024710497024424977],[0.14292308496870448,0.2808676923735748,-0.010574288572728728],[-0.1734971049470612,0.14019515029516152,0.0027683484887051457],[0.014142947617173336,0.18292110872737416,-0.004092845767710294]]
mass=[0.9999999999950272,1.6601208254808336e-07,2.447838287784771e-06,3.0404326489511185e-06,3.2271560828978514e-07,0.0009547919099366768,0.0002858856700231729,4.366249613200406e-05,5.151383772628957e-05]
NAMES=['Mer','Ven','Ear','Mar','Jup','Sat','Ura','Nep']

def initial():
    s=rebound.Simulation()
    for m,x,v in zip(mass,pos,vel): s.add(m=m,x=x[0],y=x[1],z=x[2],vx=v[0],vy=v[1],vz=v[2])
    s.move_to_com()
    return np.array([[q.x,q.y,q.z,q.vx,q.vy,q.vz] for q in s.particles])

def sim_from(U):
    s=rebound.Simulation()
    for m,u in zip(mass,U): s.add(m=m,x=u[0],y=u[1],z=u[2],vx=u[3],vy=u[4],vz=u[5])
    return s

def state(s): return np.array([[q.x,q.y,q.z,q.vx,q.vy,q.vz] for q in s.particles])

def prop(U, days, dt, corr=17, kernel=2, gr=1):
    s=sim_from(U)
    s.integrator="whfast512"; s.dt=dt*DAY; s.exact_finish_time=0
    s.integrator.gr_potential=gr; s.integrator.corrector=corr; s.integrator.kernel=kernel
    n=int(round(days/dt)); assert abs(n*dt-days)<1e-6, (days,dt)
    s.integrator.concatenate_steps=n
    s.steps(1)
    return state(s)

def to_el(U):
    s=sim_from(U)
    E=[]
    for o in s.orbits():
        si=np.sin(o.inc/2)
        E.append([np.log(o.a), o.l, o.e*np.cos(o.pomega), o.e*np.sin(o.pomega), si*np.cos(o.Omega), si*np.sin(o.Omega)])
    return np.array(E)

def from_el(E):
    s=rebound.Simulation(); s.add(m=mass[0])
    for m,e in zip(mass[1:],E):
        ecc=np.hypot(e[2],e[3]); pom=np.arctan2(e[3],e[2]); si=np.hypot(e[4],e[5]); Om=np.arctan2(e[5],e[4])
        s.add(m=m, a=np.exp(e[0]), l=e[1], e=ecc, pomega=pom, inc=2*np.arcsin(si), Omega=Om)
    s.move_to_com()
    return state(s)

def wrap(x): return (x+np.pi)%(2*np.pi)-np.pi

def ediff(A,B):
    """element difference A-B with wrapped longitudes"""
    D=A-B; D[:,1]=wrap(D[:,1]); return D

def poserr(A,B): return np.max(np.abs(A[1:,:3]-B[1:,:3]))
