.. _theory:

******
Theory
******

This chapter provides the theoretical background to the Perram–Wertheim
algorithm [PW85]_. We use matrices rather than tensors: a point/vector is
therefore defined through the 3×1 column-vector of its coordinates. Likewise, a
second-rank tensor is represented by its 3×3 matrix.

Only the global, cartesian frame is considered here, and there is no
ambiguity about the basis to which these column vectors and square matrices
refer.

.. _theory-representation:

Mathematical representation of ellipsoids
=========================================

Ellipsoids are defined from their center ``c`` and a positive-definite quadratic
form ``Q`` as the set of points ``m`` such that::

  (1)    (m-c)ᵀ⋅Q⁻¹⋅(m-c) ≤ 1.

``Q`` is a symmetric, positive-definite matrix::

  (2)    Q = ∑ aᵢ² vᵢ⋅vᵢᵀ,
             ⁱ

where ``a₁``, ``a₂``, ``a₃`` are the lengths of the principal semi-axes and
``v₁``, ``v₂``, ``v₃`` their directions (unit vectors).

In the ``PW85`` library, ``Q`` is represented as a ``double[6]`` array ``q``
which stores the upper triangular part of ``Q`` in row-major order::

             ⎡ q[0] q[1] q[2] ⎤
  (3)    Q = ⎢      q[3] q[4] ⎥.
             ⎣ sym.      q[5] ⎦


The contact function of two ellipsoids
======================================

Let ``E₁`` and ``E₂`` be two ellipsoids, defined by their centers ``c₁`` and
``c₂`` and quadratic forms ``Q₁`` and ``Q₂``, respectively.

.. _theory-eq-4:

For ``0 ≤ λ ≤ 1`` and a point ``x``, we introduce the function::

  (4)    F(x, λ) = λ(x-c₁)ᵀ⋅Q₁⁻¹⋅(x-c₁)+(1-λ)(x-c₂)ᵀ⋅Q₂⁻¹⋅(x-c₂).

.. _theory-eq-5:

For fixed ``λ``, ``F(x, λ)`` has a unique minimum [PW85]_ ``f(λ)``, and we
define::

  (5)    f(λ) = min{ F(x, λ), x ∈ ℝ³ }, 0 ≤ λ ≤ 1.

Now, the function ``f`` has a unique maximum over ``[0, 1]``, and the“contact
function” ``F(r₁₂, Q₁, Q₂)`` of ellipsoids ``E₁`` and ``E₂`` is defined as::

  (6)    F(r₁₂, Q₁, Q₂) = max{ f(λ), 0 ≤ λ ≤ 1 },

where ``r₁₂ = c₂-c₁``. It can be shown that

- if ``F(r₁₂, Q₁, Q₂) < 1`` then ``E₁`` and ``E₂`` overlap,
- if ``F(r₁₂, Q₁, Q₂) = 1`` then ``E₁`` and ``E₂`` are externally tangent,
- if ``F(r₁₂, Q₁, Q₂) > 1`` then ``E₁`` and ``E₂`` do not overlap.

The contact function therefore provides a criterion to check overlap of two
ellipsoids. The ``PW85`` library computes this value.


Geometric interpretation
========================

.. _theory-eq-7:

The scalar ``λ`` being fixed, we introduce the minimizer ``x₀(λ)`` of ``F(x,
λ)``. The stationarity of ``F`` w.r.t to ``x`` reads::

  (7)    ∇F(x₀(λ), λ) = 0,

.. _theory-eq-8:

which leads to::

  (8)    λQ₁⁻¹⋅[x₀(λ)-c₁] + (1-λ)Q₂⁻¹⋅[x₀(λ)-c₂] = 0,

.. _theory-eq-9:

and can be rearranged::

  (9a)    x₀(λ)-c₁ = (1-λ)Q₁⋅Q⁻¹⋅r₁₂,
  (9b)    x₀(λ)-c₂ = -λQ₂⋅Q⁻¹⋅r₁₂,

.. _theory-eq-10:

with::

  (10)    Q = (1-λ)Q₁ + λQ₂.

.. _theory-eq-11:

It results from the above that::

  (11)    f(λ) = F(x₀(λ), λ) = λ(1-λ)r₁₂ᵀ⋅Q⁻¹⋅r₁₂.

.. _theory-eq-12:

Maximization of ``f`` with respect to ``λ`` now delivers the stationarity
condition::

                                            ∂F
  (12)    0 = f'(λ) = ∇F(x₀(λ), λ)⋅x₀'(λ) + ──(x₀(λ), λ).
                                            ∂λ

.. _theory-eq-13:

Using Eqs. :ref:`(4) <theory-eq-4>` and :ref:`(7) <theory-eq-7>`, it is found
that ``f`` is minimum for ``λ = λ₀`` such that::

  (13)    [x₀(λ₀)-c₁]ᵀ⋅Q₁⁻¹⋅[x₀(λ₀)-c₁] = [x₀(λ₀)-c₂]ᵀ⋅Q₂⁻¹⋅[x₀(λ₀)-c₂].

Let ``μ²`` be this common value. It trivially results from Eqs. :ref:`(4)
<theory-eq-4>` and :ref:`(13) <theory-eq-13>` that ``μ² = F(x₀(λ₀), λ₀)``. In
other words, ``μ²`` is the value of the contact function.

We are now in a position to give a geometric interpretation of ``μ``. It results
from Eq. :ref:`(13) <theory-eq-13>` and the definition of ``μ`` that::

  (14a)    [x₀(λ₀)-c₁]ᵀ⋅(μ²Q₁)⁻¹⋅[x₀(λ₀)-c₁] = 1,

and::

  (14b)    [x₀(λ₀)-c₂]ᵀ⋅(μ²Q₂)⁻¹⋅[x₀(λ₀)-c₂] = 1.

The above equations mean that ``x₀(λ₀)`` belongs to both ellipsoids centered at
``cⱼ`` and defined by the symmetric, positive-definite quadratic form ``μ²Qⱼ``
(``j = 1, 2``). These two ellipsoids are nothing but the initial ellipsoids
``E₁`` and ``E₂``, scaled by the *same* factor ``μ``.

Furthermore, Eq. :ref:`(8) <theory-eq-8>` applies for ``λ = λ₀``. Therefore, the
normals to the scaled ellipsoids coincide at ``x₀(λ₀)``: the two scaled
ellipsoids are externally tangent.

To sum up, ``μ`` is the common factor by wich ellipsoids ``E₁`` and ``E₂`` must
be scaled in order for them to be externally tangent at point ``x₀(λ₀)``.


Implementation
==============

.. _theory-eq-15:

In this section, we explain how the contact function is computed. From Eq.
:ref:`(12) <theory-eq-12>`, the value of the contact function is found from
the solution ``λ`` to equation ``f'(λ) = 0``. We observe that ``f(λ)`` is a
rational function [see Eq. :ref:`(11) <theory-eq-11>`], and we write::

                  λ(1-λ)a(λ)
  (15)    f(λ) =  ──────────,
                     b(λ)

with::

  (16a)    a(λ) = r₁₂ᵀ⋅adj[(1-λ)Q₁+λQ₂]⋅r₁₂ = a₀ + a₁λ + a₂λ²,

  (16b)    b(λ) = det[(1-λ)Q₁+λQ₂] = b₀ + b₁λ + b₂λ² + b₃λ³,

where ``adj(Q)`` denotes the adjugate matrix of ``Q`` (transpose of its cofactor
matrix), see e.g `Wikipedia <https://en.wikipedia.org/wiki/Adjugate_matrix>`_.

The coefficients ``aᵢ`` and ``bᵢ`` are found from the evaluation of ``a(λ)`` and
``b(λ)`` for specific values of ``λ``::

  (17a)    a₀ = a(0),

                a(1) - a(-1)
  (17b)    a₁ = ────────────,
		     2

                a(1) + a(-1)
  (17c)    a₂ = ──────────── - a(0),
		     2

  (17d)    b₀ = b(0),

                8b(½)            b(1)   b(-1)
  (17e)    b₁ = ─────  - 2b(0) - ──── - ─────
                  3                2      6

                b(1) + b(-1)
  (17f)    b₂ = ──────────── - b(0),
		     2

                 8b(½)                   b(-1)
  (17g)    b₃ = -─────  + 2b(0) + b(1) - ─────.
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
  (18)    f'(λ) = ───────,
                   b(λ)²

where ``c(λ)`` is a sixth-order polynomial in λ::

  (19)    c(λ) = c₀ + c₁λ + c₂λ² + c₃λ³ + c₄λ⁴ + c₅λ⁵ + c₆λ⁶,

with::

  (20a)    c₀ = a₀b₀,
  (20b)    c₁ = 2(a₁-a₀)b₀,
  (20c)    c₂ = -a₀(b₁+b₂) + 3b₀(a₂-a₁) + a₁b₁,
  (20d)    c₃ = 2[b₁(a₂-a₁) - a₀b₃] - 4a₂b₀,
  (20e)    c₄ = (a₀-a₁)b₃ + (a₂-a₁)b₂ - 3a₂b₁,
  (20f)    c₅ = -2a₂b₂,
  (20g)    c₆ = -a₂b₃,

Solving ``f'(λ) = 0`` for ``λ`` is therefore equivalent to finding the unique
root of ``c`` in the interval ``0 ≤ λ ≤ 1``. For the sake of robustness, the
`bisection method <https://en.wikipedia.org/wiki/Bisection_method>`_ has been
implemented (more efficient methods will be implemented in future versions).

Once ``λ`` is found, ``μ`` is computed from ``μ² = f(λ)`` using Eq. :ref:`(15)
<theory-eq-15>`.

References
==========

.. [PW85] Perram, J. W., & Wertheim, M. S. (1985). Statistical
          mechanics of hard ellipsoids. I. Overlap algorithm and the
          contact function. *Journal of Computational Physics*, 58(3),
          409–416. https://doi.org/10.1016/0021-9991(85)90171-8
