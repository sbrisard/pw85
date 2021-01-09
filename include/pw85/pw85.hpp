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
DllExport void pw85__cholesky_decomp(double const a[PW85_SYM],
                                     double l[PW85_SYM]);
DllExport void pw85__cholesky_solve(double const l[PW85_SYM],
                                    double const b[PW85_DIM],
                                    double x[PW85_DIM]);
DllExport void pw85_spheroid(double a, double c, double const n[PW85_DIM],
                             double q[PW85_SYM]);
DllExport double pw85_f_neg(double lambda, double const* params);
DllExport void pw85__residual(double lambda, double const r12[PW85_DIM],
                              double const q1[PW85_SYM],
                              double const q2[PW85_SYM], double out[3]);
DllExport int pw85_contact_function(double const r12[PW85_DIM],
                                    double const q1[PW85_SYM],
                                    double const q2[PW85_SYM], double out[2]);
}  // namespace pw85
#endif
