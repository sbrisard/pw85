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

/* Private functions */
/* ================= */
/* These functions are exported for the sake of testing only. */

__declspec(dllexport) double pw85__det_sym(double a[PW85_SYM])
{
    return a[0] * a[3] * a[5] + 2 * a[1] * a[2] * a[4] - a[0] * a[4] * a[4] - a[3] * a[2] * a[2] - a[5] * a[1] * a[1];
}

__declspec(dllexport) double pw85__xT_adjA_x(double x[PW85_DIM],
                                             double a[PW85_SYM])
{
    return (x[0] * x[0] * (a[3] * a[5] - a[4] * a[4]) + x[1] * x[1] * (a[0] * a[5] - a[2] * a[2]) + x[2] * x[2] * (a[0] * a[3] - a[1] * a[1]) + 2. * (x[0] * x[1] * (a[2] * a[4] - a[1] * a[5]) + x[0] * x[2] * (a[1] * a[4] - a[2] * a[3]) + x[1] * x[2] * (a[1] * a[2] - a[0] * a[4])));
}

__declspec(dllexport) void pw85__rT_adjQ_r_as_poly(double *r, double *q1, double *q2, double *q3, double *a)
{
    const double a_zero = pw85__xT_adjA_x(r, q1);
    const double a_one = pw85__xT_adjA_x(r, q2);
    const double a_minus_one = pw85__xT_adjA_x(r, q3);
    a[0] = a_zero;
    a[2] = 0.5 * (a_one + a_minus_one) - a_zero;
    a[1] = 0.5 * (a_one - a_minus_one);
}

__declspec(dllexport) void pw85__detQ_as_poly(double *q1, double *q2, double *q3, double *q4, double *b)
{
    const double b_zero = pw85__det_sym(q1);
    const double b_one = pw85__det_sym(q2);
    /* Compute det[(1-x)*q1+x*q2] for x = -1. */
    const double b_minus_one = pw85__det_sym(q3);
    /* Compute det[(1-x)*q1+x*q2] for x = 1/2. */
    const double b_one_half = pw85__det_sym(q4);
    b[0] = b_zero;
    b[2] = 0.5 * (b_one + b_minus_one) - b_zero;
    b[1] = (8. * b_one_half - 6. * b_zero - 1.5 * b_one - 0.5 * b_minus_one) / 3.;
    b[3] = (-8. * b_one_half + 6. * b_zero + 3. * b_one - b_minus_one) / 3.;
}

/* Public API */
/* ========== */

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


__declspec(dllexport) double pw85_contact_function(double *r, double *q1, double *q2, double *out)
{
    double q3[PW85_SYM]; /* q3 = 2*q1-q2. */
    double q4[PW85_SYM]; /* q4 = (q1+q2)/2. */
    for (int i = 0; i < PW85_SYM; i++)
    {
        q3[i] = 2. * q1[i] - q2[i];
        q4[i] = 0.5 * (q1[i] + q2[i]);
    }

    double a[3];
    pw85__rT_adjQ_r_as_poly(r, q1, q2, q3, a);

    double b[4];
    pw85__detQ_as_poly(q1, q2, q3, q4, b);

    const double c0 = a[0] * b[0];
    const double c1 = 2. * (a[1] - a[0]) * b[0];
    const double c2 = -a[0] * (b[1] + b[2]) + 3. * b[0] * (a[2] - a[1]) + a[1] * b[1];
    const double c3 = 2. * (b[1] * (a[2] - a[1]) - a[0] * b[3]) - 4. * a[2] * b[0];
    const double c4 = (a[0] - a[1]) * b[3] + (a[2] - a[1]) * b[2] - 3. * a[2] * b[1];
    const double c5 = -2. * a[2] * b[2];
    const double c6 = -a[2] * b[3];
    double xl = 0.;
    double yl = c0;
    double xr = 1.;
    double yr = c0 + c1 + c2 + c3 + c4 + c5 + c6;
    double x, y;
    while (fabs(xr - xl) > 1E-12)
    {
        x = 0.5 * (xl + xr);
        y = c0 + x * (c1 + x * (c2 + x * (c3 + x * (c4 + x * (c5 + x * c6)))));
        if (y == 0.0)
        {
            break;
        }
        else
        {
            if (yl * y < 0)
            {
                xr = x;
                yr = y;
            }
            else
            {
                xl = x;
                yl = y;
            }
        }
    }
    y = x * (1. - x) * (a[0] + x * (a[1] + x * a[2])) / (b[0] + x * (b[1] + x * (b[2] + x * b[3])));
    if (out != NULL)
    {
        out[0] = y;
        out[1] = x;
    }
    return y;
}
