/** @file */
#pragma once
#include <cmath>

/** The dimension of the physical space (3).*/
#define PW85_DIM 3

/** The dimension of the space of symmetric matrices (6). */
#define PW85_SYM 6

/*
 * For the Brent algorithm, these two constants should be such that
 *
 * [(3 - sqrt(5)) / 2] ** n < eps
 *
 * where
 *
 *     n = PW85_MAX_ITER
 *     eps = PW85_LAMBDA_ATOL
 */

/**
 * The absolute tolerance for the stopping criterion of Brent’s method (in
 * function contact_function()).
 */
#define PW85_LAMBDA_ATOL 1e-6

/**
 * The maximum number of iterations of Brent’s method (in function
 * contact_function()).
 */
#define PW85_MAX_ITER 25

/**
 * The total number of iterations of the Newton–Raphson refinement phase (in
 * function contact_function()).
 */
#define PW85_NR_ITER 3

#if _WIN32
#define DllExport __declspec(dllexport)
#else
#define DllExport
#endif

namespace pw85 {

/**
 * Compute the Cholesky decomposition of a symmetric, positive matrix.
 *
 * Let `A` be a symmetric, positive matrix, defined by the `double[6]` array
 * `a`. This function computes the lower-triangular matrix `L`, defined by
 * the `double[6]` array `l`, such that `Lᵀ⋅L = A`.
 *
 * The array `l` must be pre-allocated; it is modified by this function. Note
 * that storage of the coefficients of `L` is as follows
 *
 * ```
 *     ⎡ l[0]    0    0 ⎤
 * L = ⎢ l[1] l[3]    0 ⎥.
 *     ⎣ l[2] l[4] l[5] ⎦
 * ```
 */
void _cholesky_decomp(const double *a, double *l) {
  l[0] = sqrt(a[0]);
  l[1] = a[1] / l[0];
  l[2] = a[2] / l[0];
  l[3] = sqrt(a[3] - l[1] * l[1]);
  l[4] = (a[4] - l[1] * l[2]) / l[3];
  l[5] = sqrt(a[5] - l[2] * l[2] - l[4] * l[4]);
}

/**
 * Compute the solution to a previously Cholesky decomposed linear system.
 *
 * Let `L` be a lower-triangular matrix, defined by the `double[6]` array `l`
 * (see pw85::_cholesky_decomp() for ordering of the coefficients). This
 * function solves (by substitution) the linear system `Lᵀ⋅L⋅x = b`, where the
 * vectors `x` and `b` are specified through their `double[3]` array of
 * coordinates; `x` is modified by this function.
 */
DllExport void _cholesky_solve(const double *l, const double *b, double *x);

/**
 * Compute the quadratic form associated to a spheroid.
 *
 * The spheroid is defined by its equatorial radius `a`, its polar radius `c`
 * and the direction of its axis of revolution, `n` (unit-vector, `double[3]`
 * array).
 *
 * `q` is the representation of a symmetric matrix as a `double[6]` array. It is
 * modified in-place.
 */
DllExport void spheroid(double a, double c, const double *n, double *q);

/**
 * Return the value of the opposite of the function ``f`` defined as (see
 * @verbatim embed:rst:inline :ref:`theory`@endverbatim).
 *
 * @f[f(\lambda)=\lambda\bigl(1-\lambda\bigr)r_{12}^{\mathsf{T}}
 *               \cdot Q^{-1}\cdot r_{12},@f]
 *
 * with
 *
 * @f[Q = \bigl(1-\lambda\bigr)Q_1 + \lambda Q₂,@f]
 *
 * where ellipsoids 1 and 2 are defined as the sets of points @f$m@f$
 * (column-vector) such that
 *
 * @f[\bigl(m-c_i\bigr)\cdot Q_i^{-1}\cdot\bigl(m-c_i\bigr)\leq1.@f]
 *
 * In the above inequality, @f$c_i@f$ is the center; @f$r_{12}=c_2-c_1@f$ is the
 * center-to-center radius-vector, represented by the first 3 coefficients of
 * the array `params`. The symmetric, positive-definite matrices @f$Q_1@f$ and
 * @f$Q_2@f$ are specified through the next 12 coefficients. In other words, if
 * @f$r_{12}@f$, @f$Q_1@f$ and @f$Q_2@f$ were defined as usual by their
 * `double[3]`, `double[6]` and `double[6]` arrays `r12`, `q1` and `q2`, then
 * `params` would be formed as follows
 *
 * ```
 * double params[] = {r12[0], r12[1], r12[2],
 *                    q1[0], q1[1], q1[2], q1[3], q1[4], q1[5],
 *                    q2[0], q2[1], q2[2], q2[3], q2[4], q2[5]};
 * ```
 *
 * The value of @f$\lambda@f$ is specified through the parameter ``lambda``.
 *
 * This function returns the value of @f$−f(\lambda)@f$ (the “minus” sign comes
 * from the fact that we seek the maximum of @f$f@f$, or the minimum of
 * @f$−f@f$).
 *
 * This implementation uses
 * @verbatim embed:rst:inline:ref:`Cholesky decompositions <implementation-cholesky>`@endverbatim.
 * Its somewhat awkward signature is defined in accordance with `gsl_min.h` from
 * the GNU Scientific Library.
 */
DllExport double f_neg(double lambda, const double *params);

/**
 * Compute the residual \f$g(\lambda)=\mu_2^2-\mu_1^2\f$.
 *
 * See @verbatim embed:rst:inline :ref:`optimization` @endverbatim for the
 * definition of \f$g\f$. The value of \f$\lambda\f$ is specified through the
 * parameter `lambda`. See contact_function() for the definition of the
 * parameters `r12`, `q1` and `q2`.
 *
 * The preallocated `double[3]` array `out` is updated as follows:
 * `out[0]=`\f$f(\lambda)\f$, `out[1]=`\f$g(\lambda)\f$ and
 * `out[2]=`\f$g'(\lambda)\f$.
 *
 * This function is used in function
 * @verbatim embed:rst:inline :cpp:func:`pw85::contact_function` @endverbatim
 * for the final Newton–Raphson refinement step.
 */
DllExport void _residual(double lambda, const double *r12, const double *q1,
                         const double *q2, double *out);

/**
 * Compute the value of the contact function of two ellipsoids.
 *
 * The center-to-center radius-vector @f$r_{12}@f$ is specified by the
 * `double[3]` array `r12`. The symmetric, positive-definite matrices @f$Q_1@f$
 * and @f$Q_2@f$ that define the two ellipsoids are specified through the
 * `double[6]` arrays `q1` and `q2`.
 *
 * This function computes the value of @f$\mu^2@f$, defined as
 *
 * @f[\mu^2=\max_{0\leq\lambda\leq 1}\bigl\{\lambda\bigl(1-\lambda\bigr)
 * r_{12}^{\mathsf{T}}\cdot\bigl[\bigl(1-\lambda\bigr)Q_1+\lambda Q_2\bigr]^{-1}
 * \cdot r_{12}\bigr\},@f]
 *
 * and the maximizer @f$\lambda@f$, see
 * @verbatim embed:rst:inline:ref:`theory` @endverbatim. Both values are stored
 * in the preallocated `double[2]` array `out`: `out[0] = `@f$\mu^2@f$ and
 * `out[1] = `@f$\lambda@f$.
 *
 * @f$\mu@f$ is the common factor by which the two ellipsoids must be scaled
 * (their centers being fixed) in order to be tangentially in contact.
 *
 * This function returns `0`.
 *
 * @verbatim embed:rst:leading-asterisk
 * .. todo:: This function should return an error code.
 * @endverbatim
 */
DllExport int contact_function(const double *r12, const double *q1,
                               const double *q2, double *out);
}  // namespace pw85