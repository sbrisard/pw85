#########
The C API
#########

.. contents:: Contents
   :local:

.. highlight:: none

Representation of vectors and matrices
======================================

An ellipsoid is defined from its center ``c`` (a 3×1, column-vector)
and quadratic form ``Q`` (a 3×3, symmetric, positive definite matrix)
as the set of points ``m`` such that::

  (m-c)ᵀ⋅Q⁻¹⋅(m-c) ≤ 1.

In this module, objects referred to as “vectors” are ``double[3]``
arrays of coordinates. In other words, the representation of the
vector ``x`` is the ``double[3]`` array ``x`` such that::

      ⎡ x[0] ⎤
  x = ⎢ x[1] ⎥.
      ⎣ x[2] ⎦

Objects referred to as “symmetric matrices” (or “quadratic forms”) are
of type ``double[6]``. Such arrays list in row-major order the
coefficients of the triangular upper part. In other words, the
representation of a the symmetric matrix ``A`` is the ``double[6]``
array ``a`` such that::

      ⎡ a[0] a[1] a[2] ⎤
  A = ⎢      a[3] a[4] ⎥.
      ⎣ sym.      a[5] ⎦


“Public” API
============

These function and macros form the public API of the library.


.. c:macro:: PW85_DIM

  The dimension of the physical space (3).


.. c:macro:: PW85_SYM

  The dimension of the space of symmetric matrices (6).


.. c:function:: void pw85_spheroid(double a, double c, double n[PW85_DIM], double q[PW85_SYM])

  Compute the quadratic form associated to a spheroid.

  The spheroid is defined by its equatorial radius `a`, its polar
  radius `c` and the direction of its axis of revolution, `n`.

  `q` is the representation of a symmetric matrix as a ``double[6]``
  array. It is modified in-place.


.. c:function:: double pw85_f(double lambda, double r12[PW85_DIM], double q1[PW85_SYM], double q2[PW85_SYM], double* out)

  Return the value of the function ``f`` defined as (see
  :ref:`theory`)::

    f(λ) = λ(1-λ)r₁₂ᵀ⋅Q⁻¹⋅r₁₂,

  with::

    Q = (1-λ)Q₁ + λQ₂,

  where ellipsoids 1 and 2 are defined as the sets of points ``m``
  (column-vector) such that::

    (m-cᵢ)⋅Qᵢ⁻¹⋅(m-cᵢ) ≤ 1

  In the above inequality, ``cᵢ`` is the center; ``r₁₂ = c₂-c₁`` is
  the center-to-center radius-vector, represented by the ``double[3]``
  array `r12`. The symmetric, positive-definite matrices ``Q₁`` and
  ``Q₂`` are specified through the ``double[6]`` arrays `q1` and `q2`.

  The value of ``λ`` is specified through the parameter `lambda`.

  This function returns the value of ``f(λ)``. If `out` is not
  ``NULL``, then it must be a pre-allocated ``double[3]`` array which
  is updated with the values of the first and second derivatives::

    out[0] = f(λ),    out[1] = f'(λ)    and    out[2] = f″(λ).

  This implementation uses :ref:`Cholesky decompositions
  <implementation-cholesky>`.


.. c:function:: double pw85_f_alt(double lambda, double r12[PW85_DIM], double q1[PW85_SYM], double q2[PW85_SYM], double* out)

  Alternative implementation of :c:func:`pw85_f`.

  See :c:func:`pw85_f` for the meaning of the parameters `lambda`,
  `r12`, `q1` and `q2`.

  This function returns the value of ``f(λ)``. If `out` is not
  ``NULL``, then it must be a pre-allocated ``double[1]`` array which
  is updated with the value of ``f(λ)``.

  This implementation uses :ref:`rational fractions
  <implementation-rational-functions>`.

.. todo:: This function should also compute the first and second
          derivatives.


.. c:function:: double pw85_contact_function(double r12[PW85_DIM], double q1[PW85_SYM], double q2[PW85_SYM], double *out)

  Return the value of the contact function of two ellipsoids.

  See :c:func:`pw85_f` for the meaning of the parameters `r12`, `q1`
  and `q2`.

  This function returns the value of ``μ²``, defined as (see :ref:`theory`)::

    μ² = max{ λ(1-λ)r₁₂ᵀ⋅[(1-λ)Q₁ + λQ₂]⁻¹⋅r₁₂, 0 ≤ λ ≤ 1 }.

  ``μ`` is the common factor by which the two ellipsoids must be
  scaled (their centers being fixed) in order to be tangentially in
  contact.

  If `out` is not ``NULL``, then a full-output is produced: ``out[0]``
  is updated with the value of ``μ²``, while ``out[1]`` is updated
  with the maximizer ``λ``.


“Private” API
=============

These functions are not really private. They are fully exposed for
testing purposes.  However, they are not really needed for standard
applications of the library.


.. c:function:: double pw85__det_sym(double a[PW85_SYM])

  Return the determinant of ``A``.

  The symmetric matrix ``A`` is specified through the ``double[6]`` array `a`.


.. c:function:: double pw85__xT_adjA_x(double x[PW85_DIM], double a[PW85_SYM])

  Return the product ``xᵀ⋅adj(A)⋅x``.

  The column vector ``x`` is specified through the ``double[3]`` array
  `x`.  The symmetric matrix ``A`` is specified trough the
  ``double[6]`` array `a`.

  ``adj(A)`` denotes the adjugate matrix of ``A`` (transpose of its
  cofactor matrix), see e.g `Wikipedia
  <https://en.wikipedia.org/wiki/Adjugate_matrix>`_.


.. c:function:: void pw85__detQ_as_poly(double q1[PW85_SYM], double q2[PW85_SYM], double q3[PW85_SYM], double q4[PW85_SYM], double b[PW85_DIM+1])

Compute the coefficients of the polynomial ``λ ↦ det[(1-λ)Q₁+λQ₂]``.

The symmetric, positive definite, 3×3 matrices ``Q₁`` and ``Q₂`` are specified
as arrays `q1` and `q2`. The arrays `q3` and `q4` must hold the difference
``2Q₁-Q₂`` and average ``(Q₁+Q₂)/2``, respectively::

  q3[i] = 2*q1[i] - q2[i]  and  q4[i] = 0.5*(q1[i] + q2[i]),

for ``i = 0, …, PW85_SYM-1``. The returned polynomial has degree
:c:macro:`PW85_DIM`::

  det[(1-λ)Q₁+λQ₂] = b₀ + b₁λ + b₂λ² + b₃λ³.

The coefficients ``bᵢ`` are stored in `b` in *increasing* order: ``b[i] = bᵢ``.


.. c:function:: double pw85__rT_adjQ_r_as_poly(double r[PW85_DIM], double q1[PW85_SYM], double q2[PW85_SYM], double q3[PW85_SYM], double a[PW85_DIM])

Compute the coefficients of the polynomial ``λ ↦ rᵀ⋅adj[(1-λ)Q₁+λQ₂]⋅r``.

The symmetric, positive definite, 3×3 matrices ``Q₁`` and ``Q₂`` are specified
as arrays `q1` and `q2`. The array `q3` must hold the difference ``2Q₁-Q₂``::

  q3[i] = 2*q1[i] - q2[i],

for ``i = 0, …, PW85_SYM-1``. The returned polynomial has degree
``PW85_DIM - 1``::

  rᵀ⋅adj[(1-λ)Q₁+λQ₂]⋅r = a₀ + a₁λ + a₂λ².

The coefficients ``aᵢ`` are stored in `a` in *increasing* order: ``a[i] = aᵢ``.
