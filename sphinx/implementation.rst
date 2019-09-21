.. _implementation:

************************************
Implementation of the function ``f``
************************************

.. _implementation-eq-1:

In this chapter, we explain how the contact function is computed. From
Eq.  :ref:`(12) <theory-eq-12>` in chapter :ref:`theory`, the value of
the contact function is found from the solution ``λ`` to equation
``f'(λ) = 0``, where it is recalled that ``f`` is defined as follows::

  (1)    f(λ) = λ(1-λ)r₁₂ᵀ⋅Q⁻¹⋅r₁₂,

with::

  (2)    Q = (1-λ)Q₁ + λQ₂.

In the present chapter, we discuss two implementations for the
evaluation of ``f``.  The :ref:`first implementation
<implementation-cholesky>` uses the Cholesky decomposition of
``Q``. The :ref:`second implementation
<implementation-rational-functions>` uses a representation of ``f`` as
a quotient of two polynomials (rational fraction).


.. _implementation-cholesky:

Implementation #1: using Cholesky decompositions
================================================

Since ``Q`` is a symmetric, positive definite matrix, we can compute
its `Cholesky decomposition
<https://en.wikipedia.org/wiki/Cholesky_decomposition>`_, which reads
as follows::

  (3)    Q = L⋅Lᵀ,

where ``L`` is a lower-triangular matrix. Using this decomposition, it
is straightforward to compute ``s = Q⁻¹⋅r`` (where we write ``r`` as a
shorthand for ``r₁₂``), so that::

  (4)    f(λ) = λ(1-λ)rᵀ⋅s.

In order to solve ``f′(λ) = 0`` numerically, we use a `Newton–Raphson
<https://en.wikipedia.org/wiki/Newton%27s_method>`_ procedure, which
requires the first and second derivatives of ``f``. It is readily
found that::

  (5)    s′ = -Q⁻¹⋅Q′⋅Q⁻¹⋅r = -Q⁻¹⋅u    and    rᵀ⋅s′ = -rᵀ⋅Q⁻¹⋅u = -sᵀ⋅u,

with ``u = Q₁₂⋅s`` and ``Q₁₂ = Q₂-Q₁``. Therefore::

  (6)    f′(λ) = (1-2λ)rᵀ⋅s - λ(1-λ)sᵀ⋅u.

Similarly, introducing ``v = Q⁻¹⋅u``::

  (7)    sᵀ⋅u′ = sᵀ⋅Q₁₂⋅s′ = -sᵀ⋅Q₁₂⋅Q⁻¹⋅u = -uᵀ⋅v,

and::

  (8)    uᵀ⋅s′ = -uᵀ⋅Q⁻¹⋅u = -uᵀ⋅v.

Therefore::

  (9)    f″(λ) = -2rᵀ⋅s - 2(1-2λ)sᵀ⋅u + 2λ(1-λ)uᵀ⋅v.


.. _implementation-rational-functions:

Implementation #2: using rational functions
===========================================

.. _implementation-eq-10:

We observe that ``f(λ)`` is a
rational function [see Eq. :ref:`(1) <implementation-eq-1>`], and we write::

                  λ(1-λ)a(λ)
  (10)    f(λ) =  ──────────,
                    b(λ)

with::

  (11a)    a(λ) = r₁₂ᵀ⋅adj[(1-λ)Q₁+λQ₂]⋅r₁₂ = a₀ + a₁λ + a₂λ²,

  (11b)    b(λ) = det[(1-λ)Q₁+λQ₂] = b₀ + b₁λ + b₂λ² + b₃λ³,

where ``adj(Q)`` denotes the adjugate matrix of ``Q`` (transpose of its cofactor
matrix), see e.g `Wikipedia <https://en.wikipedia.org/wiki/Adjugate_matrix>`_.

The coefficients ``aᵢ`` and ``bᵢ`` are found from the evaluation of ``a(λ)`` and
``b(λ)`` for specific values of ``λ``::

  (12a)    a₀ = a(0),

                a(1) - a(-1)
  (12b)    a₁ = ────────────,
                     2

                a(1) + a(-1)
  (12c)    a₂ = ──────────── - a(0),
		     2

  (12d)    b₀ = b(0),

                8b(½)            b(1)   b(-1)
  (12e)    b₁ = ─────  - 2b(0) - ──── - ─────
                  3                2      6

                b(1) + b(-1)
  (12f)    b₂ = ──────────── - b(0),
         	     2

                 8b(½)                   b(-1)
  (12g)    b₃ = -─────  + 2b(0) + b(1) - ─────.
                   3                       3

This requires the implementation of the determinant and the adjugate matrix of a
3×3, symmetric matrix, see :c:func:`pw85__det_sym` and
:c:func:`pw85__xT_adjA_x`.

Evaluating the derivative of ``f`` with respect to ``λ`` is fairly easy. The following `Sympy <https://www.sympy.org>`_ script will do the job::

  import sympy

  from sympy import Equality, numer, pprint, Symbol

  if __name__ == '__main__':
      sympy.init_printing(use_latex=False, use_unicode=True)
      λ = Symbol('λ')
      a = sum(sympy.Symbol('a{}'.format(i))*λ**i for i in range(3))
      b = sum(sympy.Symbol('b{}'.format(i))*λ**i for i in range(4))
      f = λ*(1-λ)*a/b
      f_prime = f.diff(λ).ratsimp()
      c = numer(f_prime)
      c_dict = c.collect(λ, evaluate=False)
      for i in range(sympy.degree(c, gen=λ)+1):
          pprint(Equality(Symbol('c{}'.format(i)), c_dict[λ**i]))

It is readily found that::

                   c(λ)
  (13)    f′(λ) = ───────,
                   b(λ)²

where ``c(λ)`` is a sixth-order polynomial in λ::

  (14)    c(λ) = c₀ + c₁λ + c₂λ² + c₃λ³ + c₄λ⁴ + c₅λ⁵ + c₆λ⁶,

with::

  (15a)    c₀ = a₀b₀,
  (15b)    c₁ = 2(a₁-a₀)b₀,
  (15c)    c₂ = -a₀(b₁+b₂) + 3b₀(a₂-a₁) + a₁b₁,
  (15d)    c₃ = 2[b₁(a₂-a₁) - a₀b₃] - 4a₂b₀,
  (15e)    c₄ = (a₀-a₁)b₃ + (a₂-a₁)b₂ - 3a₂b₁,
  (15f)    c₅ = -2a₂b₂,
  (15g)    c₆ = -a₂b₃,

Solving ``f'(λ) = 0`` for ``λ`` is therefore equivalent to finding the unique
root of ``c`` in the interval ``0 ≤ λ ≤ 1``. For the sake of robustness, the
`bisection method <https://en.wikipedia.org/wiki/Bisection_method>`_ has been
implemented (more efficient methods will be implemented in future versions).

Once ``λ`` is found, ``μ`` is computed from ``μ² = f(λ)`` using Eq. :ref:`(10)
<implementation-eq-10>`.


Comparison of the two implementations
=====================================

High precision reference data was generated using the `mpmath
<http://mpmath.org/>`_ library. The reference dataset is fully described and
freely downloadable from the `Zenodo platform <http://about.zenodo.org/>`_
(`DOI:10.5281/zenodo.3323683
<https://doi.org/10.5281/zenodo.3323683>`_). Accuracy of both implementations is
then evaluated through the following script (:download:`download source file
<./implementation/f_accuracy/f_accuracy.c>`):

.. literalinclude:: ./implementation/f_accuracy/f_accuracy.c
   :language: C

.. note:: To compute this program, you might need to pass the options
          ``-Dpw85_include=…``, ``-Dpw85_lib=…`` and ``-Dpw85_data=…``
          to ``meson`` (see :ref:`c-tutorial`).

We get the histograms shown in :numref:`implementation-histograms`. These
histograms show that :ref:`Implementation #1 <implementation-cholesky>` is more
accurate than :ref:`Implementation #2 <implementation-rational-functions>`. The
former will therefore be selected as default.

.. _implementation-histograms:
.. figure:: implementation/histograms.*

   Accuracy of the two implementations.
