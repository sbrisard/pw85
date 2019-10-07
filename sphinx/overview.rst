.. _overview:

********
Overview
********

Introduction
============

It is quite common in materials science to reason on assemblies of ellipsoids as
model materials. Although simplified upscaling mean-field/effective-field
theories exist for such microstructures, they often fail to capture the finest
details of the microstructure, such as orientation correlations between
anisotropic inclusions, or particle-size distributions. In order to account for
such microstructural details, one must resort to so-called *full-field*
numerical simulations (using dedicated tools such as `Damask
<https://damask.mpie.de/>`_ or `Janus <https://github.com/sbrisard/janus>`_, for
example).

Full-field simulations require *realizations* of the microstructure. For
composites made of ellipsoidal inclusions embedded in a (homogeneous) matrix,
this requires to be able to generate assemblies of (non-overlapping)
ellipsoids. The basic ingredient of such microstructure simulations is of course
the overlap test of two inclusions.

Checking for the overlap (or the absence of it) of two ellipsoids is not as
trivial as checking for the overlap of two spheres. Several criteria can be
found in the literature [@viei1972; @perr1985; @wang2001; @chen2007;
@anou2018]. We propose an implementation of the *contact function* of @perr1985.

The present chapter is organised as follows. We first give a brief description
of the contact function. Then, we discuss two essential features of this
function: robustness with respect to floating-point errors and suitability for
application to Monte-Carlo simulations. Finally, we give a brief description of
the ``pw85`` library.

The contact function of @perr1985
=================================

The origin being fixed, points are represented by the ``3×1`` column-vector of
their coordinates in a global cartesian frame. For ``i = 1, 2``, ``Eᵢ ⊂ ℝ³``
denotes the following ellipsoid::

  (1)    Eᵢ = {m ∈ ℝ³: (m-cᵢ)ᵀ⋅Qᵢ⁻¹⋅(m-cᵢ) ≤ 1},

where ``cᵢ ∈ ℝ³`` is the center of ``Eᵢ``, and ``Qᵢ`` is a positive definite
matrix. @perr1985 define the following function::

 (2)    f(λ; r₁₂, Q₁, Q₂) = λ(1-λ)r₁₂ᵀ⋅Q⁻¹⋅r₁₂,

where ``0 ≤ λ ≤ 1`` is a scalar, ``Q = (1-λ)Q₁ + λQ₂``, and ``r₁₂ = c₂-c₁``
denotes the center-to-center radius-vector. The *contact function* ``μ²(E₁,
E₂)`` of the two ellipsoids is defined as the unique maximum of ``f`` over
``(0, 1)``::

  (3)   μ² = max{f(λ; r₁₂, Q₁, Q₂), 0 ≤ λ ≤ 1}.

It turns out that the contact function has a simple geometric
interpretation. Indeed, ``μ`` is the quantity by which each of the two
ellipsoids ``E₁`` and ``E₂`` must be scaled to bring them in contact. Therefore,
an overlap test could be defined as follows

- ``μ²(E₁, E₂) < 1``: the two ellipsoids overlap,
- ``μ²(E₁, E₂) > 1``: the two ellipsoids do not overlap,
- ``μ²(E₁, E₂) = 1``: the two ellipsoids are tangent.

Despite its apparent complexity, this overlap test has two nice features that
are discussed below.

Features of the overlap test
============================

Robustness with respect to floating-point errors
------------------------------------------------

All overlap tests amount to checking for the sign of a real quantity ``Φ(E₁,
E₂)`` that depends on the two ellipsoids ``E₁`` and ``E₂``. The ellipsoids do
not overlap when ``Φ(E₁, E₂) < 0``; they do overlap when ``Φ(E₁,
E₂) > 0``. Finally, we usually have ``Φ(E₁, E₂) = 0`` when ``E₁`` and ``E₂`` are
in tangent contact (but it is important to note that, depending on the overlap
criterion, the converse is not necessarily true).

In a finite precision setting, we are bound to make wrong decisions about pairs
of ellipsoids that are such that $\Phi$ is small. Indeed, let us consider a pair
of ellipsoids $(\mathcal E_1, \mathcal E_2)$ for which the true value of $\Phi$,
$\Phi_\text{exact}(\mathcal E_1, \mathcal E_2)$, is close to the machine
epsilon. Then, the numerical estimate of $\Phi$, $\Phi_\text{approx}(\mathcal
E_1, \mathcal E_2)$, is also (hopefully) a very small value. However, whether
$\Phi_\text{approx}(\mathcal E_1, \mathcal E_2)$ is the same sign as
$\Phi_\text{exact}(\mathcal E_1, \mathcal E_2)$ (and therefore delivers the
correct answer regarding overlap) is uncertain, owing to accumulation of
round-off errors. Such misclassifications are acceptable provided that they
occur for ellipsoids that are close (nearly in tangent contact). The overlap
criterion will be deemed robust if it is such that $\Phi(\mathcal E_1, \mathcal
E_2)$ is small for nearly tangent ellipsoids only. This is obviously true of the
overlap test based on the contact function of @perr1985. Note that some of the
overlap tests that can be found in the literature do not exhibit such
robustness.

## Application to Monte-Carlo simulations

Generating compact assemblies of hard particles is a notoriously difficult
task. Event-driven simulations [@done2005; @done2005a] are often used, but
require a lot of book-keeping. A comparatively simpler approach [@bris2013] is
similar to atomistic simulations with a non-physical energy. More precisely,
starting from an initial configuration where the $N$ ellipsoids $\mathcal{E}_1,
\ldots, \mathcal{E}_N$ do overlap, a simulated annealing strategy is adopted to
minimize the quantity $U(\mathcal{E}_1,\ldots,\mathcal{E}_N)$ defined as follows

$$U(\mathcal{E}_1,\ldots,\mathcal{E}_N)=\sum_{1\leq i<j\leq N}u(\mathcal{E}_i,
\mathcal{E}_j),$$

where $u(\mathcal{E}_1, \mathcal{E}_2)$ denotes an *ad-hoc* pair-wise
(non-physical) potential, that should vanish when the two ellipsoids do not
overlap, and be “more positive when the overlap is greater” (this sentense being
deliberately kept vague). A possible choice for $u$ is the following

$$u(\mathcal{E}_1, \mathcal{E}_2)=\max\{0, \mu^{-1}(\mathcal{E}_1,
\mathcal{E}_2)\}.$$

Monte-Carlo simulations using previous implementations of the @perr1985 contact
function and the above definition of the energy of the system were successfully
used to produce extremely compact assemblies of ellipsoids [@bris2013].

# Implementation

`pw85` is a C library that implements the contact function of @perr1985. It is
released under a BSD-3 license, and is available at
[https://github.com/sbrisard/pw85](). It is fully documented at
[https://sbrisard.github.io/pw85]().

The core library depends on The [GNU Scientific Library
(GSL)](https://www.gnu.org/software/gsl/) (for its implementation of the Brent
algorithm); the tests also depend on the
[GLib](https://developer.gnome.org/glib/) and
[HDF5](https://portal.hdfgroup.org/) libraries.

The API is extremely simple; in particular it defines no custom objects:
parameters of all functions are either simple types (`size_t`, `double`) or
arrays. Note that all arrays must be pre-allocated and are modified
in-place. This minimizes the risk of creating memory leaks when implementing
wrappers for higher-level (garbage-collected) languages.

A Python wrapper (based on `ctypes`) is also provided. It has the following
(fairly standard) dependencies: [NumPy](https://numpy.org/),
[pytest](https://pytest.org/) and [h5py](https://www.h5py.org/).

Note that when developing the library, several strategies have been tested for
the evaluation of the function $f$ defined above, and its
optimization. Evaluation of $f$ relies on a Cholesky decomposition of
$\mathsf{Q}$; we tested the accuracy of this implementation over a comprehensive
set of large-precision reference values that are available on Zenodo
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3323683.svg)](https://doi.org/10.5281/zenodo.3323683). Optimization
of $f$ first starts with a few iterations of Brent's robust algorithm. Then, the
estimate of the minimizer is refined through a few Newton–Raphson iterations.

# Extensions

Several improvements/extensions are planned for this library:

1. Provide a 2D implementation of the contact function.
2. Allow for early stop of the iterations. If, during the iterations, a value of
   $\lambda$ is found such that $f > 1$, then $\mu^2$ must be greater than 1,
   and the ellipsoids certainly do not overlap, which might be sufficient if the
   user is not interested in the exact value of the contact function.
3. Return error codes when necessary. Note that this would be an extra safety
   net, as the optimization procedure is extremely robust. Indeed, it never
   failed for the thousands of test cases considered (the function to optimize
   has the required convexity over $(0, 1)$).

This project welcomes contributions. We definitely need help for the following
points:

1. Define a “Code of conduct”.
2. Improve the Python wrapper (using Cython or a C extension).
3. Implement wrappers for other languages (Julia, Javascript).

# Acknowledgements

The author would like to thank Prof. Chloé Arson (GeorgiaTech Institute of
Technology, School of Civil and Environmental Engineering) for stimulating
exchanges and research ideas that motivated the exhumation of this project
(which has long been a defunct Java library).

The author would also like to thank Xianda Shen (GeorgiaTech Institute of
Technology, School of Civil and Environmental Engineering) for testing on fruity
operating systems the installation procedure of this and related libraries. His
dedication led him to valiantly fight long battles with `setuptools` and `brew`.

# References

<!-- Local Variables: -->
<!-- fill-column: 80 -->
<!-- End: -->
