###########
The C++ API
###########

.. highlight:: none


.. note:: functions whose name is prefixed with an underscore should be
	  considered as “private”: these functions are exposed for testing
	  purposes. They should not be used, since they are susceptible of
	  incompatible changes (or even removal) in future versions.


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


API
===

.. autodoxygenfile:: pw85.hpp
   :project: pw85

.. Local Variables:
.. fill-column: 79
.. End:
