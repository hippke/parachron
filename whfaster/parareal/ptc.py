# ctypes wrapper for the C coarse propagator (whfaster/pint/libpint.so)
import ctypes, numpy as np, os
from ptlib import mass
_lib=ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','pint','libpint.so'))
class Cfg(ctypes.Structure):
    _fields_=[('masses',ctypes.c_double*9),('slice_days',ctypes.c_double),('dt_out_days',ctypes.c_double),
              ('ng_ii',ctypes.c_int),('ng_io',ctypes.c_int),('maxj',ctypes.c_int),('kernel',ctypes.c_int),('heun',ctypes.c_int),('pcut_yr',ctypes.c_double),('ng_hi',ctypes.c_int),('alpha_hi',ctypes.c_double),('meann',ctypes.c_int),('dual',ctypes.c_int),('gcorr',ctypes.c_int)]
_lib.pint_default_cfg.argtypes=[ctypes.POINTER(Cfg),ctypes.POINTER(ctypes.c_double),ctypes.c_double]
_lib.pint_coarse.argtypes=[ctypes.POINTER(Cfg),ctypes.c_void_p,ctypes.c_void_p]
class CoarseC:
    def __init__(self, slice_days, dt_out=360., ng_ii=32, ng_io=16, maxj=5, kernel=2, heun=0, pcut=0., ng_hi=0, alpha_hi=0.6):
        self.c=Cfg(); m=(ctypes.c_double*9)(*mass)
        _lib.pint_default_cfg(ctypes.byref(self.c),m,slice_days)
        self.c.dt_out_days=dt_out; self.c.ng_ii=ng_ii; self.c.ng_io=ng_io; self.c.maxj=maxj; self.c.kernel=kernel; self.c.heun=heun; self.c.pcut_yr=pcut; self.c.ng_hi=ng_hi; self.c.alpha_hi=alpha_hi
    def __call__(self, E, n=0):
        Ei=np.ascontiguousarray(E,dtype=np.float64); Eo=np.zeros((8,6))
        _lib.pint_coarse(ctypes.byref(self.c),Ei.ctypes.data,Eo.ctypes.data)
        return Eo
_lib.pint_lna_disp.argtypes=[ctypes.POINTER(Cfg),ctypes.c_void_p,ctypes.c_void_p,ctypes.c_void_p]
def lna_disp(cg, E):
    Ei=np.ascontiguousarray(E,dtype=np.float64); d=np.zeros(4)
    _lib.pint_lna_disp(ctypes.byref(cg.c),Ei.ctypes.data,d.ctypes.data,None)
    return d

def lna_disp_R(cg, E):
    Ei=np.ascontiguousarray(E,dtype=np.float64); d=np.zeros(4); R=np.zeros((4,6))
    _lib.pint_lna_disp(ctypes.byref(cg.c),Ei.ctypes.data,d.ctypes.data,R.ctypes.data)
    return d,R
_lib.pint_disp6.argtypes=[ctypes.POINTER(Cfg),ctypes.c_void_p,ctypes.c_void_p,ctypes.c_void_p]
def disp6(cg, E, n_div=None):
    """all-component first-order displacement (4,6) of the inner planets (C); n_div: mean motions for the divisors"""
    Ei=np.ascontiguousarray(E,dtype=np.float64); d=np.zeros((4,6))
    nd=None if n_div is None else np.ascontiguousarray(n_div,dtype=np.float64)
    _lib.pint_disp6(ctypes.byref(cg.c),Ei.ctypes.data,None if nd is None else nd.ctypes.data,d.ctypes.data)
    return d
