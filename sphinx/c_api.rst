#########
The C API
#########

.. highlight:: none

Storage of symmetric matrices
=============================

In this library, the triangular upper part of symmetric, 3×3 matrices
is stored in length-6 arrays, in row-major order. In other words,
matrix ``A`` is reconstructed from the 1d array ``double a[6]`` as
follows::

      ⎡ a[0] a[1] a[2] ⎤
  A = ⎢      a[3] a[4] ⎥.
      ⎣ sym.      a[5] ⎦

Macros
======

.. c:macro:: PW85_DIM

	     The dimension of the physical space (3).

.. c:macro:: PW85_SYM

	     The dimension of the space of symmetric matrices (6).

Functions
=========

.. c:function:: void pw85_spheroid(double a, double c, double n[PW85_DIM], double q[PW85_SYM])

		Return the quadratic form associated to a spheroid.

		The spheroid is defined by its equatorial radius `a`,
		its polar radius `c` and the direction of its axis of
		revolution, `n` (array of size :c:macro:`PW85_DIM`);
		`q` (array of size :c:macro:`PW85_SYM`) is modified
		in-place with the coefficients of the quadratic form.


.. c:function:: void pw85__axpby(size_t n, double a, double* x, double b, double* y, double* out)

		Compute the linear combination of two vectors.

		This function performs the operation:

		.. code-block:: c

		  for (int i = 0; i < n; i++) {
    		      out[i] = a * x[i] + b * y[i];
		  }

		`n` is the common size of `x`, `y` and `out`.


.. c:function:: double pw85__det_sym(double a0, double a1, double a2, double a3, double a4, double a5)

		Return the determinant of ``A``.

		The symmetric, 3×3 matrix ``A`` is specified trough
		the coefficients of its triangular upper part, listed
		in row-major order::

		      ⎡ a0 a1 a2 ⎤
		  A = ⎢    a3 a4 ⎥.
		      ⎣ sym.  a5 ⎦


.. c:function:: double pw85__xT_adjA_x(double x0, double x1, double x2, double a0, double a1, double a2, double a3, double a4, double a5)

		Return the product ``xᵀ⋅adj(A)⋅x``.

		The column vector ``x`` is specified through its coefficients::

		      ⎡ x0 ⎤
		  x = ⎢ x1 ⎥.
		      ⎣ x2 ⎦

		The symmetric, 3×3 matrix ``A`` is specified trough
		the coefficients of its triangular upper part, listed
		in row-major order::

		      ⎡ a0 a1 a2 ⎤
		  A = ⎢    a3 a4 ⎥.
		      ⎣ sym.  a5 ⎦

		``adj(A)`` denotes the adjugate matrix of ``A``
		(transpose of its cofactor matrix), see e.g `Wikipedia
		<https://en.wikipedia.org/wiki/Adjugate_matrix>`_.


.. c:function:: void pw85_detQ_as_poly(double* q1, double* q2, double* b)

		Compute the coefficients of ``det[(1-λ)Q₁+λQ₂]`` as a polynomial
		of ``λ``.

		The symmetric, positive definite, 3×3 matrices ``Q₁``
		and ``Q₂`` are specified as arrays `q1` and `q2` of
		length :c:macro:`PW85_SYM`. The determinant is a
		polynomial of degree :c:macro:`PW85_DIM`::

		  det[(1-λ)Q₁+λQ₂] = b₀ + b₁λ + b₂λ² + b₃λ³.

		The coefficients ``bᵢ`` are stored in `b` (array of
		length ``PW85_DIM + 1``) in *increasing* order: ``b[i]
		= bᵢ``.


.. c:function:: double pw85_r12T_adjQ_r12_as_poly(double* r12, double* q1, double* q2, double* a)

		Compute the coefficients of
		``r₁₂ᵀ⋅adj[(1-λ)Q₁+λQ₂]⋅r₁₂`` as a polynomial of
		``λ``.

		The symmetric, positive definite, 3×3 matrices ``Q₁``
		and ``Q₂`` are specified as arrays `q1` and `q2` of
		length :c:macro:`PW85_SYM`. The determinant is a
		polynomial of degree ``PW85_DIM - 1``::

		  r₁₂ᵀ⋅adj[(1-λ)Q₁+λQ₂]⋅r₁₂ = a₀ + a₁λ + a₂λ².

		The coefficients ``aᵢ`` are stored in `a` (array of
		length ``PW85_DIM``) in *increasing* order: ``a[i]
		= aᵢ``.
