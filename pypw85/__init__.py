import configparser
import ctypes
import os.path

import numpy as np

from ctypes import c_char_p, c_double, c_int

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


def _det_sym(a):
    return cpw85.pw85__det_sym(a.ctypes.data_as(c_double_p))


cpw85.pw85__xT_adjA_x.argtypes = [c_double_p, c_double_p]
cpw85.pw85__xT_adjA_x.restype = c_double


def _xT_adjA_x(x, a):
    return cpw85.pw85__xT_adjA_x(x.ctypes.data_as(c_double_p),
                                 a.ctypes.data_as(c_double_p))


cpw85.pw85__get_flag.argtypes = [c_char_p]
cpw85.pw85__get_flag.restype = c_int


def _get_flag(name):
    if str(name) is name:
        name = name.encode('ascii')
    flag = cpw85.pw85__get_flag(name)
    if flag == -1:
        raise ValueError('unknown flag '+name.decode('ascii'))
    return flag


FLAG_rT_adjQ_r_AS_POLY = _get_flag(b'PW85_FLAG_rT_adjQ_r_AS_POLY')
FLAG_detQ_AS_POLY = _get_flag(b'PW85_FLAG_detQ_AS_POLY')
FLAG_CONTACT_FUNCTION = _get_flag(b'PW85_FLAG_CONTACT_FUNCTION')

cpw85.pw85_spheroid.argtypes = [c_double, c_double, c_double_p, c_double_p]
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
    # TODO Check input
    cpw85.pw85_spheroid(a, c,
                        n.ctypes.data_as(c_double_p),
                        q.ctypes.data_as(c_double_p))
    return q


cpw85.pw85_contact_function.argtypes = [c_double_p, c_double_p, c_double_p,
                                        c_double_p, c_int]
cpw85.pw85_contact_function.restype = c_double


def detQ_as_poly(q1, q2, b=None):
    if b is None:
        b = np.empty((4,), dtype=np.float64, order='C')
    r = np.zeros((3,), dtype=np.float64, order='C')
    cpw85.pw85_contact_function(r.ctypes.data_as(c_double_p),
                                q1.ctypes.data_as(c_double_p),
                                q2.ctypes.data_as(c_double_p),
                                b.ctypes.data_as(c_double_p),
                                FLAG_detQ_AS_POLY)
    return b


def rT_adjQ_r_as_poly(r, q1, q2, a=None):
    if a is None:
        a = np.empty((3,), dtype=np.float64, order='C')
    cpw85.pw85_contact_function(r.ctypes.data_as(c_double_p),
                                q1.ctypes.data_as(c_double_p),
                                q2.ctypes.data_as(c_double_p),
                                a.ctypes.data_as(c_double_p),
                                FLAG_rT_adjQ_r_AS_POLY)
    return a


def contact_function(r, q1, q2, full_output=False):
    out = np.empty((2,), dtype=np.float64, order='C')

    val = cpw85.pw85_contact_function(r.ctypes.data_as(c_double_p),
                                      q1.ctypes.data_as(c_double_p),
                                      q2.ctypes.data_as(c_double_p),
                                      out.ctypes.data_as(c_double_p),
                                      FLAG_CONTACT_FUNCTION)
    if full_output:
        return out
    else:
        return val
