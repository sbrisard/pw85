"""Python wrapper to the legacy API of PW85.

This module offers some implementations of the function ``f`` that
were eventually discarded (for accuracy reasons). These functions are
kept for reference.

The present wrapper around the legacy PW85 C library relies on the
NumPy library.  “``double[n]`` array” should be understood here as
“NumPy array with ``shape == (n,)`` and ``dtype == numpy.float64``.”

Note that true NumPy array *must* be passed (array-likes will *not*
work).

"""
import ctypes

import numpy as np

import pw85.utils


__clib = pw85.utils.load_library("libpw85")

__clib.pw85_legacy__det_sym.argtypes = [pw85.utils.c_double_p]
__clib.pw85_legacy__det_sym.restype = ctypes.c_double

__clib.pw85_legacy__xT_adjA_x.argtypes = [
    pw85.utils.c_double_p,
    pw85.utils.c_double_p,
]
__clib.pw85_legacy__xT_adjA_x.restype = ctypes.c_double

__clib.pw85_legacy__rT_adjQ_r_as_poly.argtypes = 5 * [pw85.utils.c_double_p]

__clib.pw85_legacy__detQ_as_poly.argtypes = 5 * [pw85.utils.c_double_p]

__clib.pw85_legacy_f1.argtypes = [ctypes.c_double] + 4 * [pw85.utils.c_double_p]
__clib.pw85_legacy_f1.restype = ctypes.c_double

__clib.pw85_legacy_f2.argtypes = [ctypes.c_double] + 4 * [pw85.utils.c_double_p]
__clib.pw85_legacy_f2.restype = ctypes.c_double


def _det_sym(a):
    """Return the determinant of ``A``.

    The symmetric matrix ``A`` is specified through the ``double[6]``
    array ``a``.

    """
    return __clib.pw85_legacy__det_sym(a.ctypes.data_as(pw85.utils.c_double_p))


def _xT_adjA_x(x, a):
    """Return the product ``xᵀ⋅adj(A)⋅x``.

    The column vector ``x`` is specified as a ``double[3]`` array. The
    symmetric matrix ``A`` is specified trough the ``double[6]`` array
    ``a``.

    ``adj(A)`` denotes the adjugate matrix of ``A`` (transpose of its
    cofactor matrix), see e.g `Wikipedia
    <https://en.wikipedia.org/wiki/Adjugate_matrix>`_.

    """
    return __clib.pw85_legacy__xT_adjA_x(
        x.ctypes.data_as(pw85.utils.c_double_p),
        a.ctypes.data_as(pw85.utils.c_double_p),
    )


def _rT_adjQ_r_as_poly(r, q1, q2, q3=None, a=None):
    """Compute the coefficients of the polynomial ``λ ↦ rᵀ⋅adj[(1-λ)Q₁+λQ₂]⋅r``.

    The symmetric, positive definite, 3×3 matrices ``Q₁`` and ``Q₂``
    are specified as arrays ``q1`` and ``q2``. If ``q3 is not None``,
    it must hold the difference ``2Q₁-Q₂``::

      q3[i] = 2*q1[i] - q2[i],

    for ``i = 0, …, 5``. The returned polynomial has degree ``2``::

      rᵀ⋅adj[(1-λ)Q₁+λQ₂]⋅r = a₀ + a₁λ + a₂λ².

    If ``a is not None``, it must be a pre-allocated ``double[3]``
    array. It is modified in place with the coefficients ``aᵢ``,
    stored in ``a`` in *increasing* order: ``a[i] = aᵢ``. The function
    returns ``a``.

    If ``a is None``, a new ``double[3]`` array is created and
    returned.

    """
    if q3 is None:
        q3 = 2 * q1 - q2
    if a is None:
        a = np.empty((3,), dtype=np.float64, order="C")
    args = [arg.ctypes.data_as(pw85.utils.c_double_p) for arg in (r, q1, q2, q3, a)]
    __clib.pw85_legacy__rT_adjQ_r_as_poly(*args)
    return a


def _detQ_as_poly(q1, q2, q3=None, q4=None, b=None):
    """Compute the coefficients of the polynomial ``λ ↦ det[(1-λ)Q₁+λQ₂]``.

    The symmetric, positive definite, 3×3 matrices ``Q₁`` and ``Q₂``
    are specified as arrays ``q1`` and ``q2``. If ``q3 is not None``,
    it must hold the difference ``2Q₁-Q₂``; if ``q4 is not None``, it
    must hold the average ``(Q₁+Q₂)/2``::

      q3[i] = 2*q1[i] - q2[i]  and  q4[i] = 0.5*(q1[i] + q2[i]),

    for ``i = 0, …, 5``. The returned polynomial has degree ``3``::

      det[(1-λ)Q₁+λQ₂] = b₀ + b₁λ + b₂λ² + b₃λ³.

    If ``b is not None``, it must be a pre-allocated ``double[4]``
    array. It is modified in place with the coefficients ``bᵢ``,
    stored in *increasing* order: ``b[i] = bᵢ``.

    If ``b is None``, a new ``double[4]`` array is created and
    returned.

    """
    if q3 is None:
        q3 = 2 * q1 - q2
    if q4 is None:
        q4 = 0.5 * (q1 + q2)
    if b is None:
        b = np.empty((4,), dtype=np.float64, order="C")
    args = [arg.ctypes.data_as(pw85.utils.c_double_p) for arg in (q1, q2, q3, q4, b)]
    __clib.pw85_legacy__detQ_as_poly(*args)
    return b


def f1(lambda_, r12, q1, q2, out=None):
    """Return the value of the function ``f`` defined as
    (see :ref:`theory`)::

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

    This function returns the value of ``f(λ)``. If ``out is not
    None``, then it must be a pre-allocated ``double[3]`` array which
    is updated with the values of the first and second derivatives:

    .. code-block:: none

       out[0] = f(λ),    out[1] = f'(λ)    and    out[2] = f″(λ).

    This implementation uses :ref:`Cholesky decompositions
    <implementation-cholesky>`.

    """
    return __clib.pw85_legacy_f1(
        lambda_,
        r12.ctypes.data_as(pw85.utils.c_double_p),
        q1.ctypes.data_as(pw85.utils.c_double_p),
        q2.ctypes.data_as(pw85.utils.c_double_p),
        out if out is None else out.ctypes.data_as(pw85.utils.c_double_p),
    )


def f2(lambda_, r12, q1, q2, out=None):
    """Alternative implementation of :py:func:`f1`.

    See :py:func:`f1` for the meaning of the parameters ``lambda_``,
    ``r12``, ``q1`` and ``q2``.

    This function returns the value of ``f(λ)``. If ``out is not
    None``, then it must be a pre-allocated ``double[1]`` array which
    is updated with the value of ``f(λ)``.

    This implementation uses :ref:`rational fractions
    <implementation-rational-functions>`.

    """
    return __clib.pw85_legacy_f2(
        lambda_,
        r12.ctypes.data_as(pw85.utils.c_double_p),
        q1.ctypes.data_as(pw85.utils.c_double_p),
        q2.ctypes.data_as(pw85.utils.c_double_p),
        out if out is None else out.ctypes.data_as(pw85.utils.c_double_p),
    )