#########
The C API
#########

.. highlight:: none

Storage of symmetric matrices
=============================

In this library, the triangular upper part of symmetric, 3×3 matrices is stored
in length-6 arrays, in row-major order. In other words, matrix ``A`` is
reconstructed from the 1d array ``double a[6]`` as follows::

      ⎡ a[0] a[1] a[2] ⎤
  A = ⎢      a[3] a[4] ⎥.
      ⎣ sym.      a[5] ⎦

Macros
======

.. c:macro:: PW85_DIM

  The dimension of the physical space (3).

.. c:macro:: PW85_SYM

  The dimension of the space of symmetric matrices (6).

“Public” functions
==================

These functions form the public API of the library.

.. c:function:: void pw85_spheroid(double a, double c, double n[PW85_DIM], double q[PW85_SYM])

  Return the quadratic form associated to a spheroid.

  The spheroid is defined by its equatorial radius `a`, its polar radius `c` and
  the direction of its axis of revolution, `n` (array of size
  :c:macro:`PW85_DIM`); `q` (array of size :c:macro:`PW85_SYM`) is modified
  in-place with the coefficients of the quadratic form.


.. c:function:: double pw85_contact_function(double* r12, double* q1, double* q2, double* out)

  Return the value of the contact function of two ellipsoids.

  Ellipsoids 1 and 2 are defined as the sets of points ``m`` (column-vector)
  such that::

    (m-cᵢ)⋅Qᵢ⁻¹⋅(m-cᵢ) ≤ 1

  where ``cᵢ`` is the center (column-vector); ``r₁₂ = c₂-c₁`` is the
  center-to-center radius-vector. The symmetric, positive-definite
  matrices ``Q₁`` and ``Q₂`` are specified through the arrays ``q1`` and
  ``q2`` of the coefficients of their upper triangular part (in row-major
  order)::

         ⎡ q1[0] q1[1] q1[2] ⎤            ⎡ q2[0] q2[1] q2[2] ⎤
    Q₁ = ⎢       q1[3] q1[4] ⎥  and  Q₂ = ⎢       q2[3] q2[4] ⎥.
         ⎣ sym.        q1[5] ⎦	          ⎣ sym.        q2[5] ⎦

  This function returns the value of ``μ²``, defined as (see :ref:`theory`)::

    μ² = max{ λ(1-λ)r₁₂ᵀ⋅[(1-λ)Q₁ + λQ₂]⁻¹⋅r₁₂, 0 ≤ λ ≤ 1 }.

  If ``out`` is not null, then a full-output is produced: ``out[0]`` is
  updated with the value of ``μ²``, while ``out[1]`` is updated with the
  maximizer ``λ`` .

“Private” functions
===================

These functions are not really private. They are fully exposed and tested.
However, they are not really needed for standard applications of the library.

.. c:function:: double pw85__det_sym(double a[PW85_SYM])

  Return the determinant of ``A``.

  The symmetric, 3×3 matrix ``A`` is specified trough the coefficients of its
  triangular upper part, listed in row-major order::

        ⎡ a[0] a[1] a[2] ⎤
    A = ⎢      a[3] a[4] ⎥.
        ⎣ sym.      a[5] ⎦


.. c:function:: double pw85__xT_adjA_x(double x[PW85_DIM], double a[PW85_SYM])

  Return the product ``xᵀ⋅adj(A)⋅x``.

  The column vector ``x`` is specified through its coefficients::

        ⎡ x[0] ⎤
    x = ⎢ x[1] ⎥.
        ⎣ x[2] ⎦

  The symmetric, 3×3 matrix ``A`` is specified trough the coefficients of its
  triangular upper part, listed in row-major order::

        ⎡ a[0] a[1] a[2] ⎤
    A = ⎢      a[3] a[4] ⎥.
        ⎣ sym.      a[5] ⎦

  ``adj(A)`` denotes the adjugate matrix of ``A`` (transpose of its cofactor
  matrix), see e.g `Wikipedia <https://en.wikipedia.org/wiki/Adjugate_matrix>`_.


.. c:function:: void pw85__detQ_as_poly(double* q1, double* q2, double* b)

Compute the coefficients of ``det[(1-λ)Q₁+λQ₂]`` as a polynomial of ``λ``.

The symmetric, positive definite, 3×3 matrices ``Q₁`` and ``Q₂`` are specified
as arrays `q1` and `q2` of length :c:macro:`PW85_SYM`. The determinant is a
polynomial of degree :c:macro:`PW85_DIM`::

  det[(1-λ)Q₁+λQ₂] = b₀ + b₁λ + b₂λ² + b₃λ³.

The coefficients ``bᵢ`` are stored in `b` (array of length ``PW85_DIM + 1``) in
*increasing* order: ``b[i] = bᵢ``.


.. c:function:: double pw85__rT_adjQ_r_as_poly(double* r, double* q1, double* q2, double* q3, double* a)

Compute the coefficients of ``rᵀ⋅adj[(1-λ)Q₁+λQ₂]⋅r`` as a polynomial of ``λ``.

The symmetric, positive definite, 3×3 matrices ``Q₁`` and ``Q₂`` are specified as
arrays `q1` and `q2` of length :c:macro:`PW85_SYM`. The array `q3` (also of
length :c:macro:`PW85_SYM`) must hold the difference ``2Q₁-Q₂``::

  q3[i] = 2*q1[i] - q2[i],

for ``i = 0, …, PW85_SYM-1``. The returned polynomial has degree ``PW85_DIM - 1``
::

  rᵀ⋅adj[(1-λ)Q₁+λQ₂]⋅r = a₀ + a₁λ + a₂λ².

The coefficients ``aᵢ`` are stored in `a` (array of length ``PW85_DIM``) in
*increasing* order: ``a[i] = aᵢ``.
