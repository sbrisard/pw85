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

and can be rearranged::

  (9a)    x₀(λ)-c₁ = (1-λ)Q₁⋅Q⁻¹⋅r₁₂,
  (9b)    x₀(λ)-c₂ = λQ₂⋅Q⁻¹⋅r₁₂,

with::

  (10)    Q = (1-λ)Q₁ + λQ₂.

It results from the above that::

  (11)    f(λ) = F(x₀(λ), λ) = λ(1-λ)r₁₂ᵀ⋅Q⁻¹⋅r₁₂.

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


References
==========

.. [PW85] Perram, J. W., & Wertheim, M. S. (1985). Statistical
          mechanics of hard ellipsoids. I. Overlap algorithm and the
          contact function. *Journal of Computational Physics*, 58(3),
          409–416. https://doi.org/10.1016/0021-9991(85)90171-8
