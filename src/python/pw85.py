import configparser
import ctypes
import os.path

import numpy as np
import numpy.ctypeslib as npct

from ctypes import c_double

c_double_p = ctypes.POINTER(c_double)

Vector = npct.ndpointer(dtype=np.float64, ndim=1, shape=(3,), flags='C')
Tensor = npct.ndpointer(dtype=np.float64, ndim=1, shape=(6,), flags='C')

# TODO: wrap in a function
path = os.path.expanduser('~/.pw85rc')
if os.path.isfile(path):
    cfg = configparser.ConfigParser()
    cfg.read(path)
    cpw85 = ctypes.cdll.LoadLibrary(cfg['Native']['LibraryPath'])
else:
    # TODO: improve error message
    raise RuntimeError("Configuration file not found.")


cpw85.pw85_spheroid.argtypes = [c_double, c_double, Vector, Tensor]
cpw85.pw85_spheroid.restype = None


def spheroid(a, c, n, q=None):
    # n = Vector()
    # n[0] = 10.
    # n[1] = 20.
    # n[2] = 30.
    if q is None:
        q = np.empty((6,), dtype=np.float64, order='C')
    cpw85.pw85_spheroid(1.0, 0.1, n, q)
    return q
