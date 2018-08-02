import configparser
import ctypes
import os.path


from ctypes import c_double

c_double_p = ctypes.POINTER(c_double)

Vector = c_double * 3
Tensor = c_double * 6

# TODO: wrap in a function
path = os.path.expanduser('~/.pw85rc')
if os.path.isfile(path):
    cfg = configparser.ConfigParser()
    cfg.read(path)
    cpw85 = ctypes.cdll.LoadLibrary(cfg['Native']['LibraryPath'])
else:
    # TODO: improve error message
    raise RuntimeError("Configuration file not found.")


cpw85.spheroid.argtypes = [c_double, c_double, Vector, Tensor]
cpw85.spheroid.restype = None


def spheroid(a, c, n, q=None):
    if q is None:
        q = Tensor()
    cpw85.spheroid(1.0, 0.1, Vector(), q)
    return q
