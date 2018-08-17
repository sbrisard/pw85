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


.. c:function:: double pw85_det_sym_3x3(double a0, double a1, double a2, double a3, double a4, double a5)

		Return the determinant of a 3×3, symmetric matrix.


.. c:function:: void pw85_axpby(size_t n, double a, double* x, double b, double* y, double* out)

		Compute the linear combination of two vectors.

		This function performs the operation::

		  out[i] = a * x[i] + b * y[i] for i = 0, ... n.

		`n` is the common size of `x`, `y` and `out`.


.. c:function:: void pw85_det_q_as_poly(double* q1, double* q2, double* b)

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
