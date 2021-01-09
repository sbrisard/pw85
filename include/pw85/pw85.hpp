#ifndef __PW85_H__
#define __PW85_H__
#include <cmath>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>

#define PW85_DIM 3
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
#define PW85_LAMBDA_ATOL 1e-6
#define PW85_MAX_ITER 25
#define PW85_NR_ITER 3

#if _WIN32
#define DllExport __declspec(dllexport)
#else
#define DllExport
#endif

namespace pw85 {
DllExport void _cholesky_decomp(const double *a, double *l);
DllExport void _cholesky_solve(const double *l, const double *b, double *x);
DllExport void spheroid(double a, double c, const double *n, double *q);
DllExport double f_neg(double lambda, double const *params);
DllExport void _residual(double lambda, const double *r12, const double *q1,
                         const double *q2, double *out);
DllExport int contact_function(const double *r12, const double *q1,
                               const double *q2, double *out);
}  // namespace pw85
#endif
