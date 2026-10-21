# Round 3, iteration 21: first-order transform with tables (small divisors) at the mean instead of the osculating
# semi-major axes; Python coarse (mgc, live, Heun), ln a only or all components.
import sys, mgc, ph1eval, gfilt
for spec in sys.argv[1:]:
    comps=None if 'all' in spec else [0]
    mn=int(spec.split('mn')[1][0]) if 'mn' in spec else 0
    CG=mgc.Coarse(360000,Ng=32,dt_out=360,live=1,comps=comps,pair_ng=gfilt.pair_rule,heun=1,meann=mn,giant_mean=int('gm' in spec))
    if 'g3' in spec:
        import g3; base=CG; G=g3.G3(base=lambda E: base(E,0))
        CG=lambda E,n: G(E)
    if 'q' in spec:
        import slice_quality; slice_quality.quality(lambda E: CG(E,0),1000,40,label=spec)
    else: ph1eval.evaluate(lambda E: CG(E,0), spec, verbose=True)
