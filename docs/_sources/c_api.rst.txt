#########
The C API
#########

.. contents:: Contents
   :local:

.. highlight:: none


.. note:: we use the following naming convention

	  - “public” functions are prefixed with ``pw85_`` or ``pw85_legacy_``
            (single underscore),
	  - “private” functions are prefixed with ``pw85__`` or
            ``pw85_legacy__`` (double underscore).

	  Note that “public” and “private” is a matter of convention here,
	  since all functions are exposed (mostly, for testing
	  purposes). However, double underscored functions should not be
	  considered as part of the public API and should not be used, since
	  they are susceptible of incompatible changes (or even removal) in
	  future versions.


Representation of vectors and matrices
======================================

An ellipsoid is defined from its center ``c`` (a 3×1, column-vector) and
quadratic form ``Q`` (a 3×3, symmetric, positive definite matrix) as the set of
points ``m`` such that::

  (m-c)ᵀ⋅Q⁻¹⋅(m-c) ≤ 1.

In this module, objects referred to as “vectors” are ``double[3]`` arrays of
coordinates. In other words, the representation of the vector ``x`` is the
``double[3]`` array ``x`` such that::

      ⎡ x[0] ⎤
  x = ⎢ x[1] ⎥.
      ⎣ x[2] ⎦

Objects referred to as “symmetric matrices” (or “quadratic forms”) are of type
``double[6]``. Such arrays list in row-major order the coefficients of the
triangular upper part. In other words, the representation of a the symmetric
matrix ``A`` is the ``double[6]`` array ``a`` such that::

      ⎡ a[0] a[1] a[2] ⎤
  A = ⎢      a[3] a[4] ⎥.
      ⎣ sym.      a[5] ⎦


The new API
===========

The functions and macros gathered below form the new API that should be invoked
by most users. To use these functions and macros in your code, you must include
the following header:

.. code-block:: C

  #include <pw85.h>

and use the following link directive::

  -lpw85


.. c:macro:: PW85_VERSION

  The current version of the library.


The “legacy” API
================

The functions described below belong to the legacy API. These are functions
that have been superseded by equivalent (more accurate or more efficient)
implementations in the core library. To use these functions in your code, you
must include the following header:

.. code-block:: C

  #include <pw85_legacy.h>

and use the following link directive::

  -lpw85_legacy


.. c:function:: double pw85_legacy__det_sym(double a[PW85_SYM])

  Return the determinant of ``A``.

  The symmetric matrix ``A`` is specified through the ``double[6]`` array ``a``.


.. c:function:: double pw85_legacy__xT_adjA_x(double x[PW85_DIM], double a[PW85_SYM])

  Return the product ``xᵀ⋅adj(A)⋅x``.

  The column vector ``x`` is specified as a ``double[3]`` array. The symmetric
  matrix ``A`` is specified trough the ``double[6]`` array ``a``.

  ``adj(A)`` denotes the adjugate matrix of ``A`` (transpose of its
  cofactor matrix), see e.g `Wikipedia
  <https://en.wikipedia.org/wiki/Adjugate_matrix>`_.


.. c:function:: void pw85_legacy__detQ_as_poly(double q1[PW85_SYM], double q2[PW85_SYM], double q3[PW85_SYM], double q4[PW85_SYM], double b[PW85_DIM+1])

Compute the coefficients of the polynomial ``λ ↦ det[(1-λ)Q₁+λQ₂]``.

The symmetric, positive definite, 3×3 matrices ``Q₁`` and ``Q₂`` are specified
as arrays ``q1`` and ``q2``. The arrays ``q3`` and ``q4`` must hold the
difference ``2Q₁-Q₂`` and average ``(Q₁+Q₂)/2``, respectively::

  q3[i] = 2*q1[i] - q2[i]  and  q4[i] = 0.5*(q1[i] + q2[i]),

for ``i = 0, …, PW85_SYM-1``. The returned polynomial has degree
:c:macro:`PW85_DIM`::

  det[(1-λ)Q₁+λQ₂] = b₀ + b₁λ + b₂λ² + b₃λ³.

The coefficients ``bᵢ`` are stored in ``b`` in *increasing* order: ``b[i] =
bᵢ``.


.. c:function:: double pw85__rT_adjQ_r_as_poly(double r[PW85_DIM], double q1[PW85_SYM], double q2[PW85_SYM], double q3[PW85_SYM], double a[PW85_DIM])

Compute the coefficients of the polynomial ``λ ↦ rᵀ⋅adj[(1-λ)Q₁+λQ₂]⋅r``.

The symmetric, positive definite, 3×3 matrices ``Q₁`` and ``Q₂`` are specified
as arrays ``q1`` and ``q2``. The array ``q3`` must hold the difference
``2Q₁-Q₂``::

  q3[i] = 2*q1[i] - q2[i],

for ``i = 0, …, PW85_SYM-1``. The returned polynomial has degree
``PW85_DIM - 1``::

  rᵀ⋅adj[(1-λ)Q₁+λQ₂]⋅r = a₀ + a₁λ + a₂λ².

The coefficients ``aᵢ`` are stored in ``a`` in *increasing* order: ``a[i] = aᵢ``.


.. c:function:: double pw85_legacy_f1(double lambda, double const r12[PW85_DIM], double const q1[PW85_SYM], double const q2[PW85_SYM], double* out)

  Return the value of the function ``f`` defined as (see :ref:`theory`)::

    f(λ) = λ(1-λ)r₁₂ᵀ⋅Q⁻¹⋅r₁₂,

  with::

    Q = (1-λ)Q₁ + λQ₂,

  where ellipsoids 1 and 2 are defined as the sets of points ``m``
  (column-vector) such that::

    (m-cᵢ)⋅Qᵢ⁻¹⋅(m-cᵢ) ≤ 1

  In the above inequality, ``cᵢ`` is the center; ``r₁₂ = c₂-c₁`` is the
  center-to-center radius-vector, represented by the ``double[3]`` array
  ``r12``. The symmetric, positive-definite matrices ``Q₁`` and ``Q₂`` are
  specified through the ``double[6]`` arrays ``q1`` and ``q2``.

  The value of ``λ`` is specified through the parameter ``lambda``.

  This function returns the value of ``f(λ)``. If ``out`` is not ``NULL``, then
  it must be a pre-allocated ``double[3]`` array which is updated with the
  values of the first and second derivatives::

    out[0] = f(λ),    out[1] = f'(λ)    and    out[2] = f″(λ).

  This implementation uses :ref:`Cholesky decompositions
  <implementation-cholesky>`.


.. c:function:: double pw85_legacy_f2(double lambda, double const r12[PW85_DIM], double const q1[PW85_SYM], double const q2[PW85_SYM], double* out)

  Alternative implementation of :c:func:`pw85_legacy_f1`.

  See :c:func:`pw85_legacy_f1` for the meaning of the parameters ``lambda``,
  ``r12``, ``q1`` and ``q2``.

  This function returns the value of ``f(λ)``. If ``out`` is not ``NULL``, then
  it must be a pre-allocated ``double[1]`` array which is updated with the
  value of ``f(λ)``.

  This implementation uses :ref:`rational fractions
  <implementation-rational-functions>`.

.. todo:: This function should also compute the first and second derivatives.


.. c:function:: int pw85_legacy_contact_function1(double const r12[PW85_DIM], double const q1[PW85_SYM], double const q2[PW85_SYM], double out[2])

  Compute the value of the contact function of two ellipsoids.

  See :c:func:`pw85_contact_function` for the invocation of this function.

  Implementation of this function relies on Newton–Raphson iterations on ``f``;
  it is not robust.

  This function returns ``0``

.. todo:: This function should return an error code.


.. c:function:: int pw85_legacy_contact_function2(double const r12[PW85_DIM], double const q1[PW85_SYM], double const q2[PW85_SYM], double out[2])

  Compute the value of the contact function of two ellipsoids.

  See :c:func:`pw85_contact_function` for the invocation of this
  function.

  This implementation uses the representation of ``f`` as :ref:`rational
  fractions <implementation-rational-functions>`. Finding the maximum of ``f``
  is then equivalent to finding the root of the numerator of the rational
  fraction of ``f'``. For the sake of robustness, bisection is used to compute
  this root.

  This function returns ``0``

.. todo:: This function should return an error code.

Breathe
=======

.. autodoxygenfile:: pw85.hpp
   :project: pw85

.. Local Variables:
.. fill-column: 79
.. End:
