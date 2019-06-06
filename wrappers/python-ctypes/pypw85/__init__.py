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

The present wrapper around the PW85 C library relies on the NumPy library.
“Array of ``double``” should be understood here as “NumPy array with
``dtype == numpy.float64``.”

"""
import numpy as np

from pypw85 import lowlevel as _ll

__version__ = '${version}'
__author__ = '${author}'


def _det_sym(a):
    """Return ``det(A)``.

    ``A`` is a symmetric matrix represented by the array `a`.

    """
    return _ll._det_sym(a.ctypes.data_as(_ll.c_double_p))


def _xT_adjA_x(x, a):
    """Return ``xᵀ⋅adj(A)⋅x``.

    ``x`` is a vector, represented by the array `x`. ``A`` is a
    symmetric matrix, represented by the array `a`.

    """
    return _ll._xT_adjA_x(x.ctypes.data_as(_ll.c_double_p),
                          a.ctypes.data_as(_ll.c_double_p))


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
    args = [arg.ctypes.data_as(_ll.c_double_p) for arg in (r, q1, q2, q3, a)]
    _ll._rT_adjQ_r_as_poly(*args)
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
    args = [arg.ctypes.data_as(_ll.c_double_p) for arg in (q1, q2, q3, q4, b)]
    _ll._detQ_as_poly(*args)
    return b


def spheroid(a, c, n, q=None):
    """Return the quadratic form associated to a spheroid.

    The spheroid is defined by its equatorial radius `a`, its polar
    radius `c` and the direction of its axis of revolution, `n`
    (vector).

    If `q` is not ``None``, then it must the array representation of a
    symmetric matrix. It is modified in place.

    """
    if q is None:
        q = np.empty((6,), dtype=np.float64, order='C')
    _ll.spheroid(a, c, n.ctypes.data_as(_ll.c_double_p),
                 q.ctypes.data_as(_ll.c_double_p))
    return q


def contact_function(r12, q1, q2, out=None):
    """Return the value of the contact function of two ellipsoids.

    Ellipsoids 1 and 2 are defined as the sets of points ``m``
    (column-vector) such that::

      (m-cᵢ)⋅Qᵢ⁻¹⋅(m-cᵢ) ≤ 1

    where ``cᵢ`` is the center; ``r₁₂ = c₂-c₁`` is the
    center-to-center radius-vector, that is represented by the
    ``double[3]`` array `r12`. The symmetric, positive-definite
    matrices ``Q₁`` and ``Q₂`` are specified through the ``double[6]``
    arrays `q1` and `q2`.

    This function returns the value of ``μ²``, defined as (see
    :ref:`theory`)::

      μ² = max{ λ(1-λ)r₁₂ᵀ⋅[(1-λ)Q₁ + λQ₂]⁻¹⋅r₁₂, 0 ≤ λ ≤ 1 }.

    ``μ`` is the common factor by which the two ellipsoids must be
    scaled (their centers being fixed) in order to be tangentially in
    contact.

    If `out` is not ``None``, then a full-output is produced:
    ``out[0]`` is updated with the value of ``μ²``, while ``out[1]``
    is updated with the maximizer ``λ`` .

    """
    if out is not None:
        out = out.ctypes.data_as(_ll.c_double_p)

    return _ll.contact_function(r12.ctypes.data_as(_ll.c_double_p),
                                q1.ctypes.data_as(_ll.c_double_p),
                                q2.ctypes.data_as(_ll.c_double_p),
                                out)



def _cholesky_decomp(a, l=None):
    """Compute the Cholesky decomposition A = L⋅Lᵀ of a 3×3 matrix.

    ``A`` is a symmetric matrix, represented by the array `a`, ``L``
    is a lower matrix, represented by the array `l`.

    This function returns `l`, suitably updated with the coefficients
    of the Cholesky decomposition. If `l` is ``None``, then a new
    array is allocated.

    """
    if l is None:
        l = np.empty((6,), dtype=np.float64, order='C')
    _ll._cholesky_decomp(a.ctypes.data_as(_ll.c_double_p),
                         l.ctypes.data_as(_ll.c_double_p))
    return l


def _cholesky_solve(l, b, x=None):
    """Compute the solution of the 3×3 linear system L⋅Lᵀ⋅x = b.

    ``L`` is a lower matrix, represented by the array `l`. ``x`` and
    ``b`` are vectors, represented by the arrays `x` and `b`.

    This function returns `x`, updated with the solution to the
    system. If `x` is ``None``, then a new array is allocated.

    """
    if x is None:
        x = np.empty((3,), dtype=np.float64, order='C')
    _ll._cholesky_solve(l.ctypes.data_as(_ll.c_double_p),
                        b.ctypes.data_as(_ll.c_double_p),
                        x.ctypes.data_as(_ll.c_double_p))
    return x
