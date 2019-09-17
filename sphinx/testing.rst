.. _testing:

**************************************************
Testing the implementation of the contact function
**************************************************

This chapter describes how our implementation of the contact function is
tested. The source of the unit tests can be found in the file
``src/test_pw85.c``. Note that the tests described here are repeated over a
large set of tests case, including very flat and very slender sheroids, for
various relative orientations and center-to-center distances.

In the present chapter, we assume that the two ellispoids (their matrices
``Q₁`` and ``Q₂`` are given), as well as their center-to-center radius vector
``r₁₂``. Then, a call to :func:`pw85_contact_function` delivers an estimate of
``λ`` and ``μ²``.

We first assert that ``μ₁²`` and ``μ₂²`` as defined by Eq. :ref:`(3)
<optimization-eq-3>` in chapter :ref:`optimization` are close to the value
returned by :c:func:`pw85_contact_function`. For all the cases considered here,
this is true up to a relative error of ``10⁻¹⁰``.

We also check that ``f’(λ) = 0``, up to an absolute error of ``Δλf”(λ)`` where
``Δλ`` is the absolute tolerance on ``λ`` for the stopping criterion of the
Brent iterations, as defined by the macro :c:macro:`PW85_LAMBDA_ATOL`.

.. note:: The PW85 library offers a “new” API and a “legacy” API. The
          latter is kept for reference. The new API is generally more
          accurate and robust, and should be preferred by most
          users. Both APIs are thoroughly tested; however, we adopted
          two different testing strategies.

	  The legacy API is solely (but fully) tested through its
	  Python wrapper using pytest.

	  The new API is fully tested through pure C tests (using
	  GLib). Then the Python wrapper is also tested (using
	  pytest). However, the python tests do not need to be as
	  thorough, since only the validity of the wrapper itself must
	  be checked, not the validity of the underlying C library.
