"""Overlap test of two ellipsoids.

This module provides an implementation of the “contact function”
defined by Perram and Wertheim (J. Comp. Phys. 58(3), 409–416) for two
ellipsoids. Given two ellipsoids, this function returns the *square*
of the common factor by which both ellipsoids must be scaled (their
centers being fixed) in order to be tangentially in contact.

This module is released under a BSD 3-Clause License.

----

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

from pypw85.hilevel import spheroid, contact_function
from pypw85.hilevel import (_det_sym, _xT_adjA_x, _rT_adjQ_r_as_poly,
                            _detQ_as_poly)
