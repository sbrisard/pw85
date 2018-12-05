import configparser
import ctypes
import os.path

import numpy as np

from ctypes import c_double

c_double_p = ctypes.POINTER(c_double)


# TODO: wrap in a function
path = os.path.expanduser('~/.pw85rc')
if os.path.isfile(path):
    cfg = configparser.ConfigParser()
    cfg.read(path)
    cpw85 = ctypes.cdll.LoadLibrary(cfg['Native']['LibraryPath'])
else:
    # TODO: improve error message
    raise RuntimeError("Configuration file not found.")


cpw85.pw85__det_sym.argtypes = [c_double_p]
cpw85.pw85__det_sym.restype = c_double

cpw85.pw85__xT_adjA_x.argtypes = [c_double_p, c_double_p]
cpw85.pw85__xT_adjA_x.restype = c_double

cpw85.pw85_spheroid.argtypes = [c_double, c_double, c_double_p, c_double_p]
cpw85.pw85_spheroid.restype = None

cpw85.pw85__rT_adjQ_r_as_poly.argtypes = 5*[c_double_p]
cpw85.pw85__rT_adjQ_r_as_poly.restype = None

cpw85.pw85__detQ_as_poly.argtypes = 5*[c_double_p]
cpw85.pw85__detQ_as_poly.restype = None

cpw85.pw85_contact_function.argtypes = 4*[c_double_p]
cpw85.pw85_contact_function.restype = c_double


def _det_sym(a):
    """Return det(A) (``A``: 3×3, symmetric matrix).

    The symmetric, 3×3 matrix ``A`` is specified trough the array-like
    list `a` of the six coefficients of its triangular upper part,
    listed in row-major order::

          ⎡ a[0] a[1] a[2] ⎤
      A = ⎢      a[3] a[4] ⎥.
          ⎣ sym.      a[5] ⎦

    """
    return cpw85.pw85__det_sym(a.ctypes.data_as(c_double_p))


def _xT_adjA_x(x, a):
    """Return ``xᵀ⋅adj(A)⋅x`` (``A``: 3×3, symmetric matrix; ``x``: 3×1  vector).

    The column vector ``x`` is specified through the array-like list
    of its three coefficients::

          ⎡ x[0] ⎤
      x = ⎢ x[1] ⎥.
          ⎣ x[2] ⎦

    The symmetric, 3×3 matrix ``A`` is specified trough the array-like
    list `a` of the six coefficients of its triangular upper part,
    listed in row-major order::

          ⎡ a[0] a[1] a[2] ⎤
      A = ⎢      a[3] a[4] ⎥.
          ⎣ sym.      a[5] ⎦
    """
    return cpw85.pw85__xT_adjA_x(x.ctypes.data_as(c_double_p),
                                 a.ctypes.data_as(c_double_p))


def _rT_adjQ_r_as_poly(r, q1, q2, q3=None, a=None):
    if q3 is None:
        q3 = 2*q1-q2
    if a is None:
        a = np.empty((3,), dtype=np.float64, order='C')
    args = [arg.ctypes.data_as(c_double_p) for arg in (r, q1, q2, q3, a)]
    cpw85.pw85__rT_adjQ_r_as_poly(*args)
    return a


def _detQ_as_poly(q1, q2, q3=None, q4=None, b=None):
    if q3 is None:
        q3 = 2*q1-q2
    if q4 is None:
        q4 = 0.5*(q1+q2)
    if b is None:
        b = np.empty((4,), dtype=np.float64, order='C')
    args = [arg.ctypes.data_as(c_double_p) for arg in (q1, q2, q3, q4, b)]
    cpw85.pw85__detQ_as_poly(*args)
    return b


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
    cpw85.pw85_spheroid(a, c,
                        n.ctypes.data_as(c_double_p),
                        q.ctypes.data_as(c_double_p))
    return q


def contact_function(r, q1, q2, out=None):
    if out is not None:
        out = out.ctypes.data_as(c_double_p)

    return cpw85.pw85_contact_function(r.ctypes.data_as(c_double_p),
                                       q1.ctypes.data_as(c_double_p),
                                       q2.ctypes.data_as(c_double_p),
                                       out)
