# Round 3, iteration 21: first-order transform with tables (small divisors) at the mean instead of the osculating
# semi-major axes; Python coarse (mgc, live, Heun), ln a only or all components.
import sys, mgc, ph1eval, gfilt
EPS=None
for spec in sys.argv[1:]:
    if spec.startswith('eps='): EPS=[int(x) for x in spec[4:].split(',')]; continue
    comps=None if 'all' in spec else [0]
    mn=int(spec.split('mn')[1][0]) if 'mn' in spec else 0
    dto=float(spec.split('dt')[1].split('_')[0]) if '_dt' in spec else 360
    CG=mgc.Coarse(360000,Ng=32,dt_out=dto,live=1,comps=comps,pair_ng=gfilt.pair_rule,heun=1,meann=mn,giant_mean=int('gm' in spec),nsec=int('ns' in spec and 'nsub' not in spec),nsub=int(spec.split('nsub')[1].split('_')[0]) if 'nsub' in spec else 0,Ng_sec=int(spec.split('ngs')[1].split('_')[0]) if 'ngs' in spec else 0,sec_ng=gfilt.sec_rule if 'un' in spec else None)
    if '_L' in spec:
        L=float(spec.split('_L')[1].split('_')[0]); gf=gfilt.GF2(CG,L,end=int('noend' not in spec)); base0=CG
        CG=lambda E,n: base0(E,0)
    if 'g3' in spec:
        import g3; base=CG; G=g3.G3(base=lambda E: base(E,0))
        CG=lambda E,n: G(E)
    if 'q' in spec:
        import slice_quality; slice_quality.quality(lambda E: CG(E,0),1000,40,label=spec)
    else: ph1eval.evaluate(lambda E: CG(E,0), spec, eps=EPS, verbose=True)
