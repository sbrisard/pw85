.. _tutorial:

********
Tutorial
********

In this tutorial, we consider two ellipsoids, and check wether or not
they overlap.

Ellipsoid ``E₁`` is an `oblate spheroid
<https://en.wikipedia.org/wiki/Spheroid>`_ centered at point ``x₁ = (-0.5, 0.4,
-0.7)``, with equatorial radius ``a₁ = 10``, polar radius ``c₁ = 0.1`` and
polar axis ``(0, 0, 1)``.

Ellipsoid ``E₂`` is a `prolate spheroid
<https://en.wikipedia.org/wiki/Spheroid>`_ centered at point ``(0.2, -0.3,
0.4)``, with equatorial radius ``a₁ = 0.5``, polar radius ``c₁ = 5`` and polar
axis ``(1, 0, 0)``.

To carry out the overlap check, we must first create the representation of
ellipsoids ``Eᵢ`` as quadratic forms ``Qᵢ`` (see :ref:`theory-representation`).
Convenience functions are provided to compute the matrix representation of a
*spheroid*.

.. note:: In principle, the contact function implemented in PW85 applies to
   *any* ellipsoids (with unequal axes). However, at the time of writing this
   tutorial (2019-01-01), convenience functions to compute the matrix
   representation of a general ellipsoid is not yet implemented. Users must
   compute the matrices themselves.

We first check for the overlap of ``E₁`` and ``E₂`` using the Python wrapper of
``pw85``. We will then illustrate the C API.


Python tutorial
===============

The Python module relies on `NumPy <www.numpy.org>`_ for passing arrays to the
underlying C library. We therefore import both modules:

>>> import numpy as np
>>> import pw85

and define the parameters of the simulation:

>>> x1 = np.array([-0.5, 0.4, -0.7])
>>> n1 = np.array([0., 0., 1.])
>>> a1, c1 = 10, 0.1
>>> x2 = np.array([0.2, -0.3, 0.4])
>>> n2 = np.array([1., 0., 0.])
>>> a2, c2 = 0.5, 5.
>>> r12 = x2-x1

where ``r₁₂`` is the vector that joins the center of the first ellipsoid,
``x₁``, to the center of the second ellipsoid, ``x₂``.

We use the function :py:func:`pw85.spheroid` to create the matrix
representations ``q₁`` and ``q₂`` of the two ellipsoids. Note that these arrays
must be *preallocated*:

>>> q1 = np.empty((6,), dtype=np.float64)
>>> pw85.spheroid(a1, c1, n1, q1)
>>> q1
array([ 1.e+02, -0.e+00, -0.e+00,  1.e+02, -0.e+00,  1.e-02])
>>> q2 = np.empty_like(q1)
>>> pw85.spheroid(a2, c2, n2, q2)
>>> q2
array([25.  ,  0.  ,  0.  ,  0.25,  0.  ,  0.25])

We can now compute the value of the contact function — see the documentation of
:py:func:`pw85.contact_function`:

>>> out = np.empty((2,), dtype=np.float64)
>>> pw85.contact_function(r12, q1, q2, out)
>>> mu2, lambda_ = out
>>> print('μ² = {}'.format(mu2))
>>> print('λ = {}'.format(lambda_))
μ² = 3.362706040638343
λ = 0.1668589553405904

We find that ``μ² > 1``, hence ``μ > 1``. In other words, both ellipsoids must
be *swollen* in order to bring them in contact: the ellipsoids do not overlap!


Checking the output
-------------------

The output of this simulation can readily be checked. First, we can check that
``q₁`` and ``q₂`` indeed represent the ellipsoids ``E₁`` and ``E₂``. To do so,
we first construct the symmetric matrices ``Q₁`` and ``Q₂`` from their upper
triangular part

>>> Q1 = np.zeros((3, 3), dtype=np.float64)
>>> i, j = np.triu_indices_from(Q1)
>>> Q1[i, j] = q1
>>> Q1[j, i] = q1
>>> Q1
array([[ 1.e+02, -0.e+00, -0.e+00],
       [-0.e+00,  1.e+02, -0.e+00],
       [-0.e+00, -0.e+00,  1.e-02]])

>>> Q2 = np.zeros_like(Q1)
>>> Q2[i, j] = q2
>>> Q2[j, i] = q2
>>> Q2
array([[25.  ,  0.  ,  0.  ],
       [ 0.  ,  0.25,  0.  ],
       [ 0.  ,  0.  ,  0.25]])

We can now check these matrices for some remarkable points, first for ellipsoid
``E₁``

>>> Q1_inv = np.linalg.inv(Q1)
>>> f1 = lambda x: Q1_inv.dot(x).dot(x)
>>> f1((a1, 0., 0.))
1.0
>>> f1((-a1, 0., 0.))
1.0
>>> f1((0., a1, 0.))
1.0
>>> f1((0., -a1, 0.))
1.0
>>> f1((0., 0., c1))
0.9999999999994884
>>> f1((0., 0., -c1))
0.9999999999994884

then for ellipsoid ``E₂``

>>> Q2_inv = np.linalg.inv(Q2)
>>> f2 = lambda x: Q2_inv.dot(x).dot(x)
>>> f2((c2, 0., 0.))
1.0
>>> f2((-c2, 0., 0.))
1.0
>>> f2((0., a2, 0.))
1.0
>>> f2((0., -a2, 0.))
1.0
>>> f2((0., 0., a2))
1.0
>>> f2((0., 0., -a2))
1.0

Note that in the above tests, we have omitted the centers of ellipsoids ``E₁``
et ``E₂`` (both ellipsoids were translated to the origin).

We will now verify the corectness of the value found for the scaling factor
``μ``. To do so, we will find the coordinates of the contact point of the two
scaled ellipsoids, and check that the normals to the two ellipsoids at that
point coincide.

Although we use formulæ from the :ref:`theory` section to find the coordinates
of the contact point, ``x₀``, it is not essential for the time being to fully
understand this derivation. What really matters is to check that the resulting
point ``x₀`` is indeed the contact point of the two scaled ellipsoids; how the
coordinates of this point were found is irrelevant.

From the value of ``λ`` returned by the function
:py:func:`pw85.contact_function`, we compute ``Q`` defined by Eq. :ref:`(10)
<theory-eq-10>` in section :ref:`theory`, as well as ``x = Q⁻¹⋅r₁₂``

>>> Q = (1-lambda_)*Q1+lambda_*Q2
>>> x = np.linalg.solve(Q, r12)

From which we find ``x₀``, using either Eq. :ref:`(9a) <theory-eq-9>` or
Eq. :ref:`(9b) <theory-eq-9>` (and we can check that both give the same result)

>>> x0a = x1+(1-lambda_)*np.dot(Q1, x)
>>> x0a
array([ 0.16662271, -0.29964969, -0.51687799])
>>> x0b = x2-lambda_*np.dot(Q2, x)
>>> x0b
array([ 0.16662271, -0.29964969, -0.51687799])

We can now check that ``x₀`` belongs to the two scaled ellipsoids, that we
first define, overriding the matrices of the unscaled ellipsoids, that are no
longer needed. We observe that if ellispoid ``Eᵢ`` is scaled by ``μ``, then its
matrix representation ``Qᵢ`` is scaled by ``μ²``, and its inverse ``Qᵢ⁻¹`` is
scaled by ``μ⁻²``.

>>> x0 = x0a
>>> Q1 *= mu2
>>> Q2 *= mu2
>>> Q1_inv /= mu2
>>> Q2_inv /= mu2

>>> x = x0-x1
>>> Q1_inv.dot(x).dot(x)
1.0000000000058238

>>> x = x0-x2
>>> Q2_inv.dot(x).dot(x)
0.9999999999988334

Therefore ``x₀`` indeed belongs to both ellipsoids. We now compute the normal
``nᵢ`` to ellipsoid ``Eᵢ`` at point ``x₀``. Since ellipsoid ``Eᵢ`` is defined
by the level-set: ``(x-xᵢ)ᵀ⋅Qᵢ⁻¹⋅(x-xᵢ) = 1``, the normal to ``Eᵢ`` is given by
``Qᵢ⁻¹⋅(x-xᵢ)`` (which is then suitably normalized)

>>> n1 = Q1_inv.dot(x0-x1)
>>> n1 /= np.linalg.norm(n1)
>>> n1
array([ 3.64031943e-04, -3.82067448e-04,  9.99999861e-01])

>>> n2 = Q2_inv.dot(x0-x2)
>>> n2 /= np.linalg.norm(n2)
>>> n2
array([-3.64031943e-04,  3.82067448e-04, -9.99999861e-01])

We find that ``n₁ = -n₂``. Therefore, ``E₁`` and ``E₂`` are in external
contact. QED

Follow this link to
:download:`download the above Python script<./py_tutorial/tutorial.py>`.

.. _sec20210523205251:

C++ tutorial
============

The Python interface to PW85 has been kept close to the undelying C++ API. The
following C++ program (:download:`download source file
<./cpp_tutorial/tutorial.cpp>`) defines the two ellipsoids, then computes ``μ²``
and ``λ``:

.. literalinclude:: ./cpp_tutorial/tutorial.cpp
   :language: cpp

A ``CMakeLists.txt`` file is provided for the compilation of the tutorial using
CMake_. You can reuse it in one of your own projects (:download:`download
<./cpp_tutorial/CMakeLists.txt>`):

.. literalinclude:: ./cpp_tutorial/CMakeLists.txt
   :language: cmake


``cd`` into the ``cpp_tutorial`` subdirectory. The provided example program
should be compiled and linked against pw85::

  $ mkdir build
  $ cd build
  $ cmake ..
  $ cmake --build . --config Release

An executable called ``tutorial`` should be present in the ``build/Release``
subdirectory. On execution, it prints the following lines to ``stdout``:

.. code-block:: none

  mu^2 = 3.36271
  lambda = 0.166859

.. _CMake: https://cmake.org/

.. Local Variables:
.. fill-column: 79
.. End:
