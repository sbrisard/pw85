"""Low-level interface to the PW85 library."""

import configparser
import ctypes
import pathlib

from ctypes import c_double


def init_signature(func, argtypes, restype=None):
    func.argtypes = argtypes
    func.restype = restype


c_double_p = ctypes.POINTER(c_double)

path = pathlib.Path.home()/'pw85.ini'
if path.is_file():
    cfg = configparser.ConfigParser()
    cfg.read(path)
    cpw85 = ctypes.cdll.LoadLibrary(cfg['pw85']['FullPath'])
else:
    raise RuntimeError('Cannot file configuration file: {}'.format(path))

spheroid = cpw85.pw85_spheroid
contact_function = cpw85.pw85_contact_function

_det_sym = cpw85.pw85__det_sym
_xT_adjA_x = cpw85.pw85__xT_adjA_x
_detQ_as_poly = cpw85.pw85__detQ_as_poly
_rT_adjQ_r_as_poly = cpw85.pw85__rT_adjQ_r_as_poly
_cholesky_decomp = cpw85.pw85__cholesky_decomp
_cholesky_solve = cpw85.pw85__cholesky_solve

init_signature(spheroid, [c_double, c_double, c_double_p, c_double_p])
init_signature(contact_function, 4*[c_double_p], c_double)

init_signature(_det_sym, [c_double_p], c_double)
init_signature(_xT_adjA_x, [c_double_p, c_double_p], c_double)
init_signature(_detQ_as_poly, 5*[c_double_p])
init_signature(_rT_adjQ_r_as_poly, 5*[c_double_p])
init_signature(_cholesky_decomp, 2*[c_double_p])
init_signature(_cholesky_solve, 3*[c_double_p])
