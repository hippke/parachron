# ctypes wrapper for the C coarse propagator (whfaster/pint/libpint.so)
import ctypes, numpy as np, os
from ptlib import mass
_lib=ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','pint','libpint.so'))
class Cfg(ctypes.Structure):
    _fields_=[('masses',ctypes.c_double*9),('slice_days',ctypes.c_double),('dt_out_days',ctypes.c_double),
              ('ng_ii',ctypes.c_int),('ng_io',ctypes.c_int),('maxj',ctypes.c_int),('kernel',ctypes.c_int),('heun',ctypes.c_int),('pcut_yr',ctypes.c_double)]
_lib.pint_default_cfg.argtypes=[ctypes.POINTER(Cfg),ctypes.POINTER(ctypes.c_double),ctypes.c_double]
_lib.pint_coarse.argtypes=[ctypes.POINTER(Cfg),ctypes.c_void_p,ctypes.c_void_p]
class CoarseC:
    def __init__(self, slice_days, dt_out=360., ng_ii=32, ng_io=16, maxj=5, kernel=2, heun=0, pcut=0.):
        self.c=Cfg(); m=(ctypes.c_double*9)(*mass)
        _lib.pint_default_cfg(ctypes.byref(self.c),m,slice_days)
        self.c.dt_out_days=dt_out; self.c.ng_ii=ng_ii; self.c.ng_io=ng_io; self.c.maxj=maxj; self.c.kernel=kernel; self.c.heun=heun; self.c.pcut_yr=pcut
    def __call__(self, E, n=0):
        Ei=np.ascontiguousarray(E,dtype=np.float64); Eo=np.zeros((8,6))
        _lib.pint_coarse(ctypes.byref(self.c),Ei.ctypes.data,Eo.ctypes.data)
        return Eo
