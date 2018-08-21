#include <stdio.h>

#define PW85_DIM 3
#define PW85_SYM 6
#define PW85_XX 0
#define PW85_YY 3
#define PW85_ZZ 5
#define PW85_XY 1
#define PW85_XZ 2
#define PW85_YZ 4

__declspec(dllexport) void pw85_spheroid(double a, double c, double* n,
					 double* q)
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


__declspec(dllexport) double pw85_det_sym_3x3(double* a) {
    return a[0]*a[3]*a[5] + 2*a[1]*a[2]*a[4] - a[0]*a[4]*a[4]
	- a[3]*a[2]*a[2] - a[5]*a[1]*a[1];
}


__declspec(dllexport) void pw85_axpby(size_t n,
				      double a, double* x,
				      double b, double* y,
				      double* out) {
    for (int i = 0; i < n; i++) {
	out[i] = a * x[i] + b * y[i];
    }
}


__declspec(dllexport) void pw85_det_q_as_poly(double* q1, double* q2, double* b)
{
    double q[PW85_SYM];
    const double g_zero = pw85_det_sym_3x3(q1);
    const double g_one = pw85_det_sym_3x3(q2);
    pw85_axpby(PW85_SYM, 2., q1, -1.,q2, q);
    const double g_minus_one = pw85_det_sym_3x3(q);
    pw85_axpby(PW85_SYM, .5, q1, .5, q2, q);
    const double g_one_half = pw85_det_sym_3x3(q);

    b[0] = g_zero;
    b[2] = 0.5*(g_one+g_minus_one)-g_zero;
    b[1] = (8.*g_one_half-6.*g_zero-1.5*g_one-0.5*g_minus_one)/3.;
    b[3] = (-8.*g_one_half+6.*g_zero+3.*g_one-g_minus_one)/3.;
}


__declspec(dllexport) double pw85__xT_adjA_x(double* x, double* a) {
    const double a0 = a[0];
    const double a1 = a[1];
    const double a2 = a[2];
    const double a3 = a[3];
    const double a4 = a[4];
    const double a5 = a[5];
    const double x0 = x[0];
    const double x1 = x[1];
    const double x2 = x[2];
    return (x0*x0*(a3*a5-a4*a4)+x1*x1*(a0*a5-a2*a2)+x2*x2*(a0*a3-a1*a1)
	    +2.*(x0*x1*(a2*a4-a1*a5)+x0*x2*(a1*a4-a2*a3)+x1*x2*(a1*a2-a0*a4)));
}


__declspec(dllexport) void pw85_r12T_adjQ_r12_as_poly(double* r12, double* q1,
						      double* q2, double* a)
{
    double q[PW85_SYM];
    const double f_zero = pw85__xT_adjA_x(r12, q1);
    const double f_one = pw85__xT_adjA_x(r12, q2);
    pw85_axpby(PW85_SYM, 2., q1, -1.,q2, q);
    const double f_minus_one = pw85__xT_adjA_x(r12, q);

    a[0] = f_zero;
    a[2] = 0.5*(f_one+f_minus_one)-f_zero;
    a[1] = 0.5*(f_one-f_minus_one);
}
