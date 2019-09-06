#########
The C API
#########

.. contents:: Contents
   :local:

.. highlight:: none


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


On “public” and “private” functions
===================================

For testing purposes, all functions are exposed in this library. However,
functions prefixed with two underscores (``pw85__function_name`` or
``pw85_legacy__function_name``) should be considered as private. Allthough they
are exposed and documented, their use is discouraged.


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


.. c:macro:: PW85_DIM

  The dimension of the physical space (3).


.. c:macro:: PW85_SYM

  The dimension of the space of symmetric matrices (6).


.. c:macro:: PW85_LAMBDA_ATOL

  The absolute tolerance for the stopping criterion of Brent’s method (in
  function :c:func:`pw85_contact_function`).


.. c:macro:: PW85_MAX_ITER

  The maximum number of iterations of Brent’s method (in function
  :c:func:`pw85_contact_function`).


.. c:macro:: PW85_NR_ITER

  The total number of iterations of the Newton–Raphson refinement phase (in
  function :c:func:`pw85_contact_function`).


.. c:function:: void pw85__cholesky_decomp(double const a[PW85_SYM], double l[PW85_SYM])

  Compute the Cholesky decomposition of a symmetric, positive matrix.

  Let ``A`` be a symmetric, positive matrix, defined by the ``double[6]`` array
  ``a``. This function computes the lower-triangular matrix ``L``, defined by
  the ``double[6]`` array ``l``, such that ``Lᵀ⋅L = A``.

  The array ``l`` must be pre-allocated; it is modified by this function. Note
  that storage of the coefficients of ``L`` is as follows::

        ⎡ l[0]    0    0 ⎤
    L = ⎢ l[1] l[3]    0 ⎥.
        ⎣ l[2] l[4] l[5] ⎦


.. c:function:: void pw85__cholesky_solve(double const l[PW85_SYM], double const b[PW85_DIM], double x[PW85_DIM])

  Compute the solution to a previously Cholesky decoposed linear system.

  Let ``L`` be a lower-triangular matrix, defined by the ``double[6]`` array
  ``l`` (see :c:func:`pw85__cholesky_decomp` for ordering of the
  coefficients). This function solves (by substitution) the linear system
  ``Lᵀ⋅L⋅x = b``, where the vectors ``x`` and ``b`` are specfied through their
  ``double[3]`` array of coordinates; ``x`` is modified by this function.


.. c:function:: void pw85__residual(double lambda, double const r12[PW85_DIM], double const q1[PW85_SYM], double const q2[PW85_SYM], double out[3])

   Compute the residual ``g(λ) = μ₂² - μ₁²``.

   See :ref:`optimization` for the definition of ``g``. The value of ``λ`` is
   specified through the parameter
   ``lambda``. See :c:func:`pw85_contact_function` for the definition of the
   parameters ``r12``, ``q1`` and ``q2``.

   The preallocated ``double[3]`` array ``out`` is updated with the values of
   ``f(λ)``, ``g(λ)`` and ``g’(λ)``::

     out[0] = f(λ),    out[1] = g(λ)    and    out[2] = g’(λ).

   This function is used in function :c:func:`pw85_contact_function` for the
   final Newton–Raphson refinement step.


.. c:function:: void pw85_spheroid(double a, double c, double n[PW85_DIM], double q[PW85_SYM])

  Compute the quadratic form associated to a spheroid.

  The spheroid is defined by its equatorial radius ``a``, its polar radius
  ``c`` and the direction of its axis of revolution, ``n``.

  ``q`` is the representation of a symmetric matrix as a ``double[6]``
  array. It is modified in-place.


.. c:function:: double pw85_f_neg(double lambda, double cons* params)

  Return the value of the opposite of the function ``f`` defined as (see
  :ref:`theory`)::

    f(λ) = λ(1-λ)r₁₂ᵀ⋅Q⁻¹⋅r₁₂,

  with::

    Q = (1-λ)Q₁ + λQ₂,

  where ellipsoids 1 and 2 are defined as the sets of points ``m``
  (column-vector) such that::

    (m-cᵢ)⋅Qᵢ⁻¹⋅(m-cᵢ) ≤ 1

  In the above inequality, ``cᵢ`` is the center; ``r₁₂ = c₂-c₁`` is the
  center-to-center radius-vector, represented by the first 3 coefficients of
  the array ``params``. The symmetric, positive-definite matrices ``Q₁`` and
  ``Q₂`` are specified through the next 12 coefficients. In other words, if
  ``r12``, ``Q1`` and ``Q2`` were defined as usual by their ``double[3]``,
  ``double[6]`` and ``double[6]`` arrays ``r12``, ``q1`` and ``q2``, then
  ``params`` would be formed as follows::

    double params[] = {r12[0], r12[1], r12[2],
                       q1[0], q1[1], q1[2], q1[3], q1[4], q1[5],
		       q2[0], q2[1], q2[2], q2[3], q2[4], q2[5]};

  The value of ``λ`` is specified through the parameter ``lambda``.

  This function returns the value of ``−f(λ)`` (the “minus” sign comes from the
  fact that we seek the maximum of ``f``, or the minimum of ``−f``).

  This implementation uses :ref:`Cholesky decompositions
  <implementation-cholesky>`. Its somewhat awkward signature is defined in
  accordance with ``gsl_min.h`` from the GNU Scientific Library.


.. c:function:: int pw85_contact_function(double const r12[PW85_DIM], double const q1[PW85_SYM], double const q2[PW85_SYM], double out[2])

  Compute the value of the contact function of two ellipsoids.

  The center-to-center radius-vector is specified by the ``double[3]`` array
  ``r12``. The symmetric, positive-definite matrices ``Q₁`` and ``Q₂`` that
  define the two ellipsoides are specified through the ``double[6]`` arrays
  ``q1`` and ``q2``.

  This function returns the value of ``μ²``, defined as (see :ref:`theory`)::

    μ² = max{ λ(1-λ)r₁₂ᵀ⋅[(1-λ)Q₁ + λQ₂]⁻¹⋅r₁₂, 0 ≤ λ ≤ 1 },

  and the maximizer ``λ``. Both values are stored in the preallocated
  ``double[2]`` array ``out``::

    out[0] = μ²    and    out[1] = λ.

  ``μ`` is the common factor by which the two ellipsoids must be scaled (their
  centers being fixed) in order to be tangentially in contact.

  This function returns ``0``

.. todo:: This function should return an error code.


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

  The column vector ``x`` is specified through the ``double[3]`` array
  ``x``.  The symmetric matrix ``A`` is specified trough the
  ``double[6]`` array ``a``.

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

.. Local Variables:
.. fill-column: 79
.. End:
