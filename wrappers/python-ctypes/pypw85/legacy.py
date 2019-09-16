import configparser
import ctypes
import pathlib

from ctypes import c_double

import numpy as np

c_double_p = ctypes.POINTER(c_double)


def __load_library():
    path = pathlib.Path.home() / "pw85.ini"
    if path.is_file():
        cfg = configparser.ConfigParser()
        cfg.read(str(path))
        return ctypes.cdll.LoadLibrary(cfg["pw85"]["legacy"])
    else:
        raise RuntimeError("Cannot file configuration file: {}".format(path))

__clib = __load_library()

__clib.pw85_legacy__det_sym.argtypes = [c_double_p]
__clib.pw85_legacy__det_sym.restype = c_double

__clib.pw85_legacy__xT_adjA_x.argtypes = [c_double_p, c_double_p]
__clib.pw85_legacy__xT_adjA_x.restype = c_double

__clib.pw85_legacy__rT_adjQ_r_as_poly.argtypes=5*[c_double_p]

__clib.pw85_legacy__detQ_as_poly.argtypes=5*[c_double_p]

def _det_sym(a):
    """Return ``det(A)``.

    ``A`` is a symmetric matrix represented by the array ``a``.

    """
    return __clib.pw85_legacy__det_sym(a.ctypes.data_as(c_double_p))


def _xT_adjA_x(x, a):
    """Return ``xᵀ⋅adj(A)⋅x``.

    ``x`` is a vector, represented by the array ``x``. ``A`` is a
    symmetric matrix, represented by the array ``a``.

    """
    return __clib.pw85_legacy__xT_adjA_x(x.ctypes.data_as(c_double_p),
                                         a.ctypes.data_as(c_double_p))


def _rT_adjQ_r_as_poly(r, q1, q2, q3=None, a=None):
    """Compute the coefficients of the polynomial ``λ ↦ rᵀ⋅adj[(1-λ)Q₁+λQ₂]⋅r``.

    The symmetric, positive definite, 3×3 matrices ``Q₁`` and ``Q₂`` are
    specified as arrays `q1` and `q2`. If specified, the array `q3` must hold
    the difference ``2Q₁-Q₂``::

      q3[i] = 2*q1[i] - q2[i],

    for ``i = 0, …, 5``. The returned polynomial has degree 2::

      rᵀ⋅adj[(1-λ)Q₁+λQ₂]⋅r = a₀ + a₁λ + a₂λ².

    The coefficients ``aᵢ`` are stored in `a` in *increasing* order:
    ``a[i] = aᵢ``.

    """
    if q3 is None:
        q3 = 2*q1-q2
    if a is None:
        a = np.empty((3,), dtype=np.float64, order='C')
    args = [arg.ctypes.data_as(c_double_p) for arg in (r, q1, q2, q3, a)]
    __clib.pw85_legacy__rT_adjQ_r_as_poly(*args)
    return a


def _detQ_as_poly(q1, q2, q3=None, q4=None, b=None):
    """Compute the coefficients of the polynomial ``λ ↦ det[(1-λ)Q₁+λQ₂]``.

    The symmetric, positive definite, 3×3 matrices ``Q₁`` and ``Q₂`` are
    specified as arrays `q1` and `q2`. If specified, the arrays `q3` and `q4`
    must hold the difference ``2Q₁-Q₂`` and average ``(Q₁+Q₂)/2``,
    respectively::

      q3[i] = 2*q1[i] - q2[i]  and  q4[i] = 0.5*(q1[i] + q2[i]),

    for ``i = 0, …, 5``. The returned polynomial has degree 3::

      det[(1-λ)Q₁+λQ₂] = b₀ + b₁λ + b₂λ² + b₃λ³.

    The coefficients ``bᵢ`` are stored in `b` in *increasing* order:
    ``b[i] = bᵢ``.

    """
    if q3 is None:
        q3 = 2*q1-q2
    if q4 is None:
        q4 = 0.5*(q1+q2)
    if b is None:
        b = np.empty((4,), dtype=np.float64, order='C')
    args = [arg.ctypes.data_as(c_double_p) for arg in (q1, q2, q3, q4, b)]
    __clib.pw85_legacy__detQ_as_poly(*args)
    return b






# def f1(lambda_, r12, q1, q2, out=None):
#     """Return the value of the function ``f`` defined as::

#         f(λ) = λ(1-λ)r₁₂ᵀ⋅Q⁻¹⋅r₁₂,

#     with::

#         Q = (1-λ)Q₁ + λQ₂,

#     where ellipsoids 1 and 2 are defined as the sets of points ``m``
#     (column-vector) such that::

#         (m-cᵢ)⋅Qᵢ⁻¹⋅(m-cᵢ) ≤ 1

#     In the above inequality, ``cᵢ`` is the center; ``r₁₂ = c₂-c₁`` is
#     the center-to-center radius-vector, represented by the
#     ``double[3]`` array `r12`. The symmetric, positive-definite
#     matrices ``Q₁`` and ``Q₂`` are specified through the ``double[6]``
#     arrays `q1` and `q2`.

#     The value of ``λ`` is specified through the parameter `lambda`.

#     This function returns the value of ``f(λ)``. If `out` is not
#     ``None``, then it must be a pre-allocated ``double[3]`` array
#     which is updated with the values of the first and second
#     derivatives:

#     .. code-block:: none

#        out[0] = f(λ),    out[1] = f'(λ)    and    out[2] = f″(λ).

#     This implementation uses :ref:`Cholesky decompositions
#     <implementation-cholesky>`.

#     """
#     return _ll.f1(lambda_,
#                   r12.ctypes.data_as(_ll.c_double_p),
#                   q1.ctypes.data_as(_ll.c_double_p),
#                   q2.ctypes.data_as(_ll.c_double_p),
#                   out if out is None else out.ctypes.data_as(_ll.c_double_p))


# def f2(lambda_, r12, q1, q2, out=None):
#     """Alternative implementation of :py:func:`f`.

#     See :py:func:`f` for the meaning of the parameters `lambda`,
#     `r12`, `q1` and `q2`.

#     This function returns the value of ``f(λ)``. If `out` is not
#     ``None``, then it must be a pre-allocated ``double[1]`` array
#     which is updated with the value of ``f(λ)``.

#     This implementation uses :ref:`rational fractions
#     <implementation-rational-functions>`.

#     """
#     return _ll.f2(lambda_,
#                   r12.ctypes.data_as(_ll.c_double_p),
#                   q1.ctypes.data_as(_ll.c_double_p),
#                   q2.ctypes.data_as(_ll.c_double_p),
#                   out if out is None else out.ctypes.data_as(_ll.c_double_p))
