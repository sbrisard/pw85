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
    """Return the quadratic form associated to a spheroid.

    The spheroid is defined by its equatorial radius `a`, its polar
    radius `c` and the direction of its axis of revolution, `n`
    (array-like of length 3).

    If `q` is not ``None``, then it must be an array of length 6. It
    is modified in place.
    """
    if q is None:
        q = np.empty((6,), dtype=np.float64, order='C')
    cpw85.pw85_spheroid(a, c, n, q)
    return q


cpw85.pw85_det_sym_3x3.argtypes = [Tensor]
cpw85.pw85_det_sym_3x3.restype = c_double
_det_sym_3x3 = cpw85.pw85_det_sym_3x3


cpw85.pw85_det_q_as_poly.argtypes = [Tensor, Tensor, npct.ndpointer(dtype=np.float64, ndim=1, shape=(4,), flags='C')]
cpw85.pw85_det_q_as_poly.restype = None


def det_q_as_poly(q1, q2, b=None):
    if b is None:
        b = np.empty((4,), dtype=np.float64, order='C')
    cpw85.pw85_det_q_as_poly(q1, q2, b)
    return b


_xT_adjA_x = cpw85.pw85__xT_adjA_x
_xT_adjA_x.argtypes = [Vector, Tensor]
_xT_adjA_x.restype = c_double
