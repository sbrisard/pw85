#########
The C API
#########

.. c:macro:: PW85_DIM

	     The dimension of the physical space (3).

.. c:macro:: PW85_SYM

	     The dimension of the space of symmetric matrices (6).


.. c:function:: void pw85_spheroid(double a, double c, double n[PW85_DIM], double q[PW85_SYM])

		Return the quadratic form associated to a spheroid.

		The spheroid is defined by its equatorial radius `a`,
		its polar radius `c` and the direction of its axis of
		revolution, `n`; `q` is modified in-place with the
		coefficients of the quadratic form.
