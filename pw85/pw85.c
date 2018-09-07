#include <math.h>
#include <stdio.h>

#define PW85_DIM 3
#define PW85_SYM 6
#define PW85_XX 0
#define PW85_YY 3
#define PW85_ZZ 5
#define PW85_XY 1
#define PW85_XZ 2
#define PW85_YZ 4

#define PW85_FLAG_rT_adjQ_r_AS_POLY 1
#define PW85_FLAG_detQ_AS_POLY 2
#define PW85_FLAG_CONTACT_FUNCTION 4

/* Convenience functions */
/* ===================== */
/* These functions are exported for the sake of testing. */

__declspec(dllexport) void pw85__axpby(size_t n,
                                       double a, double *x,
                                       double b, double *y,
                                       double *out)
{
    for (int i = 0; i < n; i++)
    {
        out[i] = a * x[i] + b * y[i];
    }
}

__declspec(dllexport) double pw85__det_sym(double a0, double a1, double a2,
                                           double a3, double a4, double a5)
{
    return a0 * a3 * a5 + 2 * a1 * a2 * a4 - a0 * a4 * a4 - a3 * a2 * a2 - a5 * a1 * a1;
}

__declspec(dllexport) double pw85__xT_adjA_x(double x0, double x1, double x2,
                                             double a0, double a1, double a2,
                                             double a3, double a4, double a5)
{
    return (x0 * x0 * (a3 * a5 - a4 * a4) + x1 * x1 * (a0 * a5 - a2 * a2) + x2 * x2 * (a0 * a3 - a1 * a1) + 2. * (x0 * x1 * (a2 * a4 - a1 * a5) + x0 * x2 * (a1 * a4 - a2 * a3) + x1 * x2 * (a1 * a2 - a0 * a4)));
}

/* Public API */
/* ========== */

__declspec(dllexport) int pw85__get_flag(char* s) {
    if (strcmp(s, "PW85_FLAG_rT_adjQ_r_AS_POLY") == 0) {
	return PW85_FLAG_rT_adjQ_r_AS_POLY;
    } else if (strcmp(s, "PW85_FLAG_detQ_AS_POLY") == 0) {
	return PW85_FLAG_detQ_AS_POLY;
    } else if (strcmp(s, "PW85_FLAG_CONTACT_FUNCTION") == 0) {
	return PW85_FLAG_CONTACT_FUNCTION;
    } else {
	fprintf(stderr, "unknown flag: %s\n", s); return -1;
    }
}

__declspec(dllexport) void pw85_spheroid(double a, double c, double *n,
                                         double *q)
{
    const double a2 = a * a;
    const double c2_minus_a2 = c * c - a2;
    const double nx = n[0];
    const double ny = n[1];
    const double nz = n[2];
    q[PW85_XX] = nx * nx * c2_minus_a2 + a2;
    q[PW85_YY] = ny * ny * c2_minus_a2 + a2;
    q[PW85_ZZ] = nz * nz * c2_minus_a2 + a2;
    q[PW85_YZ] = ny * nz * c2_minus_a2;
    q[PW85_XZ] = nx * nz * c2_minus_a2;
    q[PW85_XY] = nx * ny * c2_minus_a2;
}


__declspec(dllexport) void pw85_rT_adjQ_r_as_poly(double *r, double *q1,
                                                  double *q2, double *a)
{
    const double r_0 = r[0];
    const double r_1 = r[1];
    const double r_2 = r[2];

    const double q1_0 = q1[0];
    const double q1_1 = q1[1];
    const double q1_2 = q1[2];
    const double q1_3 = q1[3];
    const double q1_4 = q1[4];
    const double q1_5 = q1[5];

    const double q2_0 = q2[0];
    const double q2_1 = q2[1];
    const double q2_2 = q2[2];
    const double q2_3 = q2[3];
    const double q2_4 = q2[4];
    const double q2_5 = q2[5];

    const double f_zero = pw85__xT_adjA_x(r_0, r_1, r_2,
                                          q1_0, q1_1, q1_2, q1_3, q1_4, q1_5);
    const double f_one = pw85__xT_adjA_x(r_0, r_1, r_2,
                                         q2_0, q2_1, q2_2, q2_3, q2_4, q2_5);
    const double q_0 = 2. * q1_0 - q2_0;
    const double q_1 = 2. * q1_1 - q2_1;
    const double q_2 = 2. * q1_2 - q2_2;
    const double q_3 = 2. * q1_3 - q2_3;
    const double q_4 = 2. * q1_4 - q2_4;
    const double q_5 = 2. * q1_5 - q2_5;
    const double f_minus_one = pw85__xT_adjA_x(r_0, r_1, r_2,
                                               q_0, q_1, q_2, q_3, q_4, q_5);
    a[0] = f_zero;
    a[2] = 0.5 * (f_one + f_minus_one) - f_zero;
    a[1] = 0.5 * (f_one - f_minus_one);
}

__declspec(dllexport) double pw85_contact_function(double *r,
						   double *q1,
						   double *q2,
						   double *out,
						   int type) {
    const double r_0 = r[0];
    const double r_1 = r[1];
    const double r_2 = r[2];

    const double q1_0 = q1[0];
    const double q1_1 = q1[1];
    const double q1_2 = q1[2];
    const double q1_3 = q1[3];
    const double q1_4 = q1[4];
    const double q1_5 = q1[5];

    const double q2_0 = q2[0];
    const double q2_1 = q2[1];
    const double q2_2 = q2[2];
    const double q2_3 = q2[3];
    const double q2_4 = q2[4];
    const double q2_5 = q2[5];

    /* Coefficients of (1-λ)Q₁+λQ₂ for λ = -1. */
    const double q_0 = 2. * q1_0 - q2_0;
    const double q_1 = 2. * q1_1 - q2_1;
    const double q_2 = 2. * q1_2 - q2_2;
    const double q_3 = 2. * q1_3 - q2_3;
    const double q_4 = 2. * q1_4 - q2_4;
    const double q_5 = 2. * q1_5 - q2_5;

    size_t out_index = 0;

    /*
     * Compute the coefficients of
     *
     *   a(λ) = rᵀ⋅adj[(1-λ)Q₁+λQ₂]⋅r
     *        = a₀ + a₁λ + a₂λ²
     *
     * as a degree-2 polynomial in λ.
     */
    const double a_zero = pw85__xT_adjA_x(r_0, r_1, r_2,
					  q1_0, q1_1, q1_2, q1_3, q1_4, q1_5);
    const double a_one = pw85__xT_adjA_x(r_0, r_1, r_2,
					 q2_0, q2_1, q2_2, q2_3, q2_4, q2_5);
    const double a_minus_one = pw85__xT_adjA_x(r_0, r_1, r_2,
					       q_0, q_1, q_2, q_3, q_4, q_5);
    const double a0 = a_zero;
    const double a2 = 0.5 * (a_one + a_minus_one) - a_zero;
    const double a1 = 0.5 * (a_one - a_minus_one);
    if (type & PW85_FLAG_rT_adjQ_r_AS_POLY) {
	out[out_index] = a0; ++out_index;
	out[out_index] = a1; ++out_index;
	out[out_index] = a2; ++out_index;
    }

    /*
     * Compute the coefficients of
     *
     *   b(λ) = det[(1-λ)Q₁+λQ₂]⋅r
     *        = b₀ + b₁λ + b₂λ² + b₃λ³
     *
     * as a degree-3 polynomial in λ.
     */
    const double b_zero = pw85__det_sym(q1_0, q1_1, q1_2, q1_3, q1_4, q1_5);
    const double b_one = pw85__det_sym(q2_0, q2_1, q2_2, q2_3, q2_4, q2_5);
    /* Compute det[(1-x)*q1+x*q2] for x = -1. */
    const double b_minus_one = pw85__det_sym(q_0, q_1, q_2, q_3, q_4, q_5);
    /* Compute det[(1-x)*q1+x*q2] for x = 1/2. */
    const double b_one_half = pw85__det_sym(.5 * (q1_0 + q2_0), .5 * (q1_1 + q2_1),
					    .5 * (q1_2 + q2_2), .5 * (q1_3 + q2_3),
					    .5 * (q1_4 + q2_4), .5 * (q1_5 + q2_5));
    const double b0 = b_zero;
    const double b2 = 0.5 * (b_one + b_minus_one) - b_zero;
    const double b1 = (8. * b_one_half - 6. * b_zero - 1.5 * b_one - 0.5 * b_minus_one) / 3.;
    const double b3 = (-8. * b_one_half + 6. * b_zero + 3. * b_one - b_minus_one) / 3.;
    if (type & PW85_FLAG_detQ_AS_POLY) {
	out[out_index] = b0; ++out_index;
	out[out_index] = b1; ++out_index;
	out[out_index] = b2; ++out_index;
	out[out_index] = b3; ++out_index;
    }
    if (type & PW85_FLAG_CONTACT_FUNCTION) {
	const double c0 = a0 * b0;
	const double c1 = 2. * (a1 - a0) * b0;
	const double c2 = -a0 * (b1 + b2) + 3. * b0 * (a2 - a1) + a1 * b1;
	const double c3 = 2. * (b1 * (a2 - a1) - a0 * b3) - 4. * a2 * b0;
	const double c4 = (a0 - a1) * b3 + (a2 - a1) * b2 - 3. * a2 * b1;
	const double c5 = -2. * a2 * b2;
	const double c6 = -a2 * b3;

	double xl = 0.;
	double yl = c0;
	double xr = 1.;
	double yr = c0 + c1 + c2 + c3 + c4 + c5 + c6;
	double x, y;
	while (fabs(xr - xl) > 1E-12) {
	    x = 0.5 * (xl + xr);
	    y = c0 + x * (c1 + x * (c2 + x * (c3 + x * (c4 + x * (c5 + x * c6)))));
	    if (y == 0.0) {
		break;
	    } else {
		if (yl * y < 0) { xr = x; yr = y; }
		else { xl = x; yl = y; }
	    }
	}
	y = x * (1. - x) * (a0 + x * (a1 + x * a2)) / (b0 + x * (b1 + x * (b2 + x * b3)));
        out[out_index] = y; ++out_index;
        out[out_index] = x; ++out_index;
	return y;
    }
    return -INFINITY;
}
