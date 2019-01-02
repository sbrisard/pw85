.. _tutorial:

********
Tutorial
********

In this tutorial, we consider to ellipsoids, and check wether or not
they overlap.

Ellipsoid ``E₁`` is an `oblate spheroid
<https://en.wikipedia.org/wiki/Spheroid>`_ centered at the origin ``(0, 0,
0)``, with equatorial radius ``a₁ = 10``, polar radius ``c₁ = 0.1`` and polar
axis ``(0, 0, 1)``.

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
>>> import pypw85

and define the parameters of the simulation:

>>> a1, c1, n1 = 10, 0.1, np.array([0., 0., 1.])
>>> a2, c2, n2 = 0.5, 5., np.array([1., 0., 0.])
>>> r12 = np.array([0.2, -0.3, 0.4])

where ``r12`` is the vector that joins the center of the first ellipsoid to the
center of the second ellipsoid.

We use the function :py:func:`pypw85.spheroid` to create the matrix
representations ``q1`` and ``q2`` of the two ellipsoids:

>>> q1 = pypw85.spheroid(a1, c1, n1)
>>> q1
array([ 1.e+02, -0.e+00, -0.e+00,  1.e+02, -0.e+00,  1.e-02])
>>> q2 = pypw85.spheroid(a2, c2, n2)
>>> q2
array([25.  ,  0.  ,  0.  ,  0.25,  0.  ,  0.25])

We can now compute the value of the contact function — see the documentation of
:py:func:`pypw85.contact_function`:

>>> mu2 = pypw85.contact_function(r12, q1, q2)
>>> print('μ² = {}'.format(mu2))
μ² = 0.4446579853323222

We find that ``μ² < 1``, hence ``μ < 1``. In other words, both ellipsoids must
be *shrunk* in order to bring them in contact: the ellipsoids overlap!


Checking the output
-------------------

The output of this simulation can readily be checked. First, we can check that
``q1`` and ``q2`` indeed represent the ellipsoids ``E₁`` and ``E₂``. To do so, we first construct the symmetric matrices ``Q1`` and ``Q2`` from their upper triangular part

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

then for ellipsoid ``E₂`` (translated at the origin)

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


C tutorial
==========
