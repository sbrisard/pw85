"""Low-level interface to the PW85 library."""

import ctypes

from ctypes import c_double

c_double_p = ctypes.POINTER(c_double)

# TODO: wrap in a function
# path = os.path.expanduser('~/.pw85rc')
# if os.path.isfile(path):
#     cfg = configparser.ConfigParser()
#     cfg.read(path)
#     cpw85 = ctypes.cdll.LoadLibrary(cfg['Native']['LibraryPath'])
# else:
#     # TODO: improve error message
#     raise RuntimeError("Configuration file not found.")

cpw85 = ctypes.cdll.LoadLibrary('libpw85')

spheroid = cpw85.pw85_spheroid
contact_function = cpw85.pw85_contact_function

_det_sym = cpw85.pw85__det_sym
_xT_adjA_x = cpw85.pw85__xT_adjA_x
_detQ_as_poly = cpw85.pw85__detQ_as_poly
_rT_adjQ_r_as_poly = cpw85.pw85__rT_adjQ_r_as_poly


def init_signature(func, argtypes, restype=None):
    func.argtypes = argtypes
    func.restype = restype

init_signature(spheroid, [c_double, c_double, c_double_p, c_double_p])
init_signature(contact_function, 4*[c_double_p], c_double)

init_signature(_det_sym, [c_double_p], c_double)
init_signature(_xT_adjA_x, [c_double_p, c_double_p], c_double)
init_signature(_detQ_as_poly, 5*[c_double_p])
init_signature(_rT_adjQ_r_as_poly, 5*[c_double_p])
