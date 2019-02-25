"""Hi level interface to the PW85 library."""

import numpy as np

from pypw85 import lowlevel as ll


def _det_sym(a):
    """Return ``det(A)``.

    ``A`` is a symmetric matrix represented by the array `a`.

    """
    return ll._det_sym(a.ctypes.data_as(ll.c_double_p))


def _xT_adjA_x(x, a):
    """Return ``xᵀ⋅adj(A)⋅x``.

    ``x`` is a vector, represented by the array `x`. ``A`` is a
    symmetric matrix, represented by the array `a`.

    """
    return ll._xT_adjA_x(x.ctypes.data_as(ll.c_double_p),
                         a.ctypes.data_as(ll.c_double_p))


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
    args = [arg.ctypes.data_as(ll.c_double_p) for arg in (r, q1, q2, q3, a)]
    ll._rT_adjQ_r_as_poly(*args)
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
    args = [arg.ctypes.data_as(ll.c_double_p) for arg in (q1, q2, q3, q4, b)]
    ll._detQ_as_poly(*args)
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
    ll.spheroid(a, c, n.ctypes.data_as(ll.c_double_p),
                q.ctypes.data_as(ll.c_double_p))
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
        out = out.ctypes.data_as(ll.c_double_p)

    return ll.contact_function(r12.ctypes.data_as(ll.c_double_p),
                               q1.ctypes.data_as(ll.c_double_p),
                               q2.ctypes.data_as(ll.c_double_p),
                               out)
