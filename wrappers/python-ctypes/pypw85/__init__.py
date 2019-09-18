"""Overlap test of two ellipsoids.

This module provides a wrapper around the PW85 C library that implements
the “contact function” defined by Perram and Wertheim (J. Comp. Phys.
58(3), 409–416, DOI:10.1016/0021-9991(85)90171-8) for two ellipsoids.
Given two ellipsoids, this function returns the *square* of the common
factor by which both ellipsoids must be scaled (their centers being
fixed) in order to be tangentially in contact.

This module is released under a BSD 3-Clause License.

Representation of vectors and matrices
--------------------------------------

An ellipsoid is defined from its center ``c`` (a 3×1, column-vector)
and quadratic form ``Q`` (a 3×3, symmetric, positive definite matrix)
as the set of points ``m`` such that::

  (m-c)ᵀ⋅Q⁻¹⋅(m-c) ≤ 1.

In this module, objects referred to as “vectors” are length-3 arrays
of ``double`` coordinates. In other words, the representation of the
vector ``x`` is the ``double[3]`` array ``x`` such that::

      ⎡ x[0] ⎤
  x = ⎢ x[1] ⎥.
      ⎣ x[2] ⎦

Objects referred to as “symmetric matrices” (or “quadratic forms”) are
length-6 arrays of ``double``. Such arrays list in row-major order the
coefficients of the triangular upper part. In other words, the
representation of a the symmetric matrix ``A`` is the array ``a`` such
that::

      ⎡ a[0] a[1] a[2] ⎤
  A = ⎢      a[3] a[4] ⎥.
      ⎣ sym.      a[5] ⎦

The present wrapper around the PW85 C library relies on the NumPy
library.  “``double[n]`` array” should be understood here as “NumPy
array with ``shape == (n,)`` and ``dtype == numpy.float64``.”

Note that true NumPy array *must* be passed (array-likes will *not*
work).

"""
import ctypes

import numpy as np

import pypw85.utils

__version__ = "${version}"
__author__ = "${author}"

__clib = pypw85.utils.load_library("libpw85")

__clib.pw85__cholesky_decomp.argtypes = 2 * [pypw85.utils.c_double_p]
__clib.pw85__cholesky_decomp.restype = None

__clib.pw85__cholesky_solve.argtypes = 3 * [pypw85.utils.c_double_p]
__clib.pw85__cholesky_solve.restype = None

__clib.pw85_spheroid.argtypes = [
    ctypes.c_double,
    ctypes.c_double,
    pypw85.utils.c_double_p,
    pypw85.utils.c_double_p,
]
__clib.pw85_spheroid.restype = None

__clib.pw85_f_neg.argtypes = [ctypes.c_double, pypw85.utils.c_double_p]
__clib.pw85_f_neg.restype = ctypes.c_double

__clib.pw85__residual.argtypes = [ctypes.c_double] + 4 * [pypw85.utils.c_double_p]
__clib.pw85__residual.restype = None

__clib.pw85_contact_function.argtypes = 4 * [pypw85.utils.c_double_p]
__clib.pw85_contact_function.restype = ctypes.c_int


def _cholesky_decomp(a, l=None):
    """Compute the Cholesky decomposition ``A = L⋅Lᵀ`` of a 3×3 matrix.

    ``A`` is a symmetric matrix, ``L`` is a lower matrix, both
    represented by ``double[6]`` arrays.

    This function returns ``l``, suitably updated with the
    coefficients of the Cholesky decomposition. If ``l`` is ``None``,
    then a new array is allocated.

    """
    if l is None:
        l = np.empty((6,), dtype=np.float64, order="C")
    __clib.pw85__cholesky_decomp(
        a.ctypes.data_as(pypw85.utils.c_double_p),
        l.ctypes.data_as(pypw85.utils.c_double_p),
    )
    return l


def _cholesky_solve(l, b, x=None):
    """Compute the solution of the 3×3 linear system ``L⋅Lᵀ⋅x = b.``

    ``L`` is a lower matrix, represented by the ``double[6]`` array
    ``l``; ``x`` and ``b`` are vectors (``double[3]`` arrays).

    This function returns ``x``, suitably updated with the solution to
    the system. If ``x`` is ``None``, then a new array is allocated.

    """
    if x is None:
        x = np.empty((3,), dtype=np.float64, order="C")
    __clib.pw85__cholesky_solve(
        l.ctypes.data_as(pypw85.utils.c_double_p),
        b.ctypes.data_as(pypw85.utils.c_double_p),
        x.ctypes.data_as(pypw85.utils.c_double_p),
    )
    return x


def spheroid(a, c, n, q=None):
    """Return the quadratic form associated to a spheroid.

    The spheroid is defined by its equatorial radius ``a``, its polar
    radius ``c`` and the direction of its axis of revolution, ``n``
    (vector, a.k.a. ``double[3]`` array).

    If ``q`` is not ``None``, then it must be a pre-allocated
    ``double[6]`` array. It is modified in place.

    """
    if q is None:
        q = np.empty((6,), dtype=np.float64, order="C")
    __clib.pw85_spheroid(
        a,
        c,
        n.ctypes.data_as(pypw85.utils.c_double_p),
        q.ctypes.data_as(pypw85.utils.c_double_p),
    )
    return q


def f(lambda_, r12, q1, q2):
    """Return the value of the function ``f`` defined as::

        f(λ) = λ(1-λ)r₁₂ᵀ⋅Q⁻¹⋅r₁₂,

    with::

        Q = (1-λ)Q₁ + λQ₂,

    where ellipsoids 1 and 2 are defined as the sets of points ``m``
    (column-vector) such that::

        (m-cᵢ)⋅Qᵢ⁻¹⋅(m-cᵢ) ≤ 1

    In the above inequality, ``cᵢ`` is the center; ``r₁₂ = c₂-c₁`` is
    the center-to-center radius-vector, represented by the
    ``double[3]`` array ``r12``. The symmetric, positive-definite
    matrices ``Q₁`` and ``Q₂`` are specified through the ``double[6]``
    arrays ``q1`` and ``q2``.

    The value of ``λ`` is specified through the parameter ``lambda_``.

    This function returns the value of ``f(λ)``.

    """
    params = np.empty((15,), dtype=np.float64)
    params[0:3] = r12
    params[3:9] = q1
    params[9:15] = q2
    return -__clib.pw85_f_neg(lambda_, params.ctypes.data_as(pypw85.utils.c_double_p))


def contact_function(r12, q1, q2, out=None):
    """Return the value of the contact function of two ellipsoids.

    See :py:func:`f` for the meaning of the parameters ``r12``, ``q1``
    and ``q2``.

    This function returns the pair ``(μ², λ)``, defined as (see
    :ref:`theory`)::

        μ² = max{ λ(1-λ)r₁₂ᵀ⋅[(1-λ)Q₁ + λQ₂]⁻¹⋅r₁₂, 0 ≤ λ ≤ 1 }

    (the returned value of ``λ`` is the actual maximizer).

    ``μ`` is the common factor by which the two ellipsoids must be
    scaled (their centers being fixed) in order to be tangentially in
    contact.

    If ``out`` is not ``None``, it must be a pre-allocated
    ``double[2]`` array. It is updated with the values of ``μ²``, and
    the maximizer ``λ``::

        out[0] = μ²    and    out[1] = λ.

    """
    if out is None:
        out = np.empty((2,), dtype=np.float64)

    __clib.pw85_contact_function(
        r12.ctypes.data_as(pypw85.utils.c_double_p),
        q1.ctypes.data_as(pypw85.utils.c_double_p),
        q2.ctypes.data_as(pypw85.utils.c_double_p),
        out.ctypes.data_as(pypw85.utils.c_double_p),
    )
    return tuple(out[:2])
