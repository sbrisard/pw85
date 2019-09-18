.. _optimization:

**********************************
Optimization of the function ``f``
**********************************

It was shown in chapter :ref:`theory` [see Eq. :ref:`(6) <theory-eq-6>`] that
the contact function was defined as the maximum for ``0 ≤ λ ≤ 1`` of the function
``f`` discussed in chapter :ref:`implementation`.

Given that the first and second derivatives of ``f`` can be computed
explicitely (see section :ref:`implementation-cholesky` in chapter
:ref:`implementation`) it would be tempting to use the Newton–Raphson method to
solve ``f’(λ)`` iteratively. However, our experiments show that this method
performs very poorly in the present case, because the variations of ``f`` can
be quite sharp in the neighborhood of ``λ = 0`` or ``λ = 1``. To carry out the
otpimization of ``f``, we therefore proceed in two steps.

In the first step, we use a robust optimization algorithm. We selected here
`Brent's method <https://en.wikipedia.org/wiki/Brent%27s_method>`_, as
implemented in the `GNU Scientific Library (GSL)
<https://www.gnu.org/software/gsl/>`_. However, this method delivers a
relatively low accuracy of the maxmimizer and the maximum.

Therefore, in the second step, we use a few Newton–Raphson iterations to refine
the previously obtained estimates of the minimizer and minimum of ``f``. In the
remainder of this chapter, we describe how these Newton–Raphson iterations are
performed.

.. _optimization-eq-1:

Our starting point is Eqs. :ref:`(9) <theory-eq-9>` and :ref:`(13)
<theory-eq-13>` in chapter :ref:`theory`, from which it results that for a given
value of ``λ`` we can define two values of ``μ²``: one is provided by
Eq. :ref:`(9a) <theory-eq-9>`, the other one is given by Eq. :ref:`(9b)
<theory-eq-9>` (both in chapter :ref:`theory`)::

  (1a)    μ₁² = [x₀(λ₀)-c₁]ᵀ⋅Q₁⁻¹⋅[x₀(λ₀)-c₁] = (1-λ)²sᵀ⋅Q₁⋅s,
  (1b)    μ₂² = [x₀(λ₀)-c₂]ᵀ⋅Q₂⁻¹⋅[x₀(λ₀)-c₂] = λ²sᵀ⋅Q₂⋅s,

.. _optimization-eq-2:

where we have introduced ``s = Q⁻¹⋅r₁₂``. We further define the matrix ``Q₁₂ =
Q₂-Q₁``, so that::

  (2)    Q₁ = Q - λQ₁₂    and    Q₂ = Q + (1-λ)Q₁₂.

.. _optimization-eq-3:

Combining Eqs. :ref:`(1) <optimization-eq-1>` and :ref:`(2) <optimization-eq-2>`
and recalling that ``Q⋅s = r`` then delivers the following expressions::

  (3a)    μ₁² = (1-λ)²rᵀ⋅s - λ(1-λ)²sᵀ⋅u,
  (3b)    μ₂² = λ²rᵀ⋅s + λ²(1-λ)sᵀ⋅u,

where we have introduced ``u = Q₁₂⋅s``.

.. _optimization-eq-4:

The above expressions seem to behave slightly better from a numerical point of
view. Our problem is now to find ``λ`` such that ``μ₁² = μ₂²``. We therefore
define the following residual::

  (4)    g(λ) = μ₂² - μ₁² = (2λ-1)rᵀ⋅s + λ(1-λ)sᵀ⋅u,

.. _optimization-eq-5:

and we need to find ``λ`` such that ``g(λ) = 0``. In order to implement
Newton–Raphson iterations, we need the expression of the derivative of the
residual. Using results that are presented in section
:ref:`implementation-cholesky`, we readily find that::

  (5)    g’(λ) = 2rᵀ⋅s + 2(1-2λ)sᵀ⋅u - 2λ(1-λ)uᵀ⋅v.

Eqs. :ref:`(4) <optimization-eq-4>` and :ref:`(5) <optimization-eq-5>` are then
used for the final, refinement step of determination of ``λ``.
