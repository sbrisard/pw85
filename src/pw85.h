#ifndef __PW85_H__
#define __PW85_H__
#include <math.h>
#include <stdlib.h>

#define PW85_DIM 3
#define PW85_SYM 6

#if _WIN32
#define DllExport __declspec(dllexport)
#else
#define DllExport
#endif

DllExport double pw85__det_sym(double a[PW85_SYM]);
DllExport double pw85__xT_adjA_x(double x[PW85_DIM], double a[PW85_SYM]);
DllExport void pw85__rT_adjQ_r_as_poly(double r[PW85_DIM], double q1[PW85_SYM],
                                       double q2[PW85_SYM], double q3[PW85_SYM],
                                       double a[PW85_DIM]);
DllExport void pw85__detQ_as_poly(double q1[PW85_SYM], double q2[PW85_SYM],
                                  double q3[PW85_SYM], double q4[PW85_SYM],
                                  double b[PW85_DIM + 1]);
DllExport void pw85_spheroid(double a, double c, double n[PW85_DIM],
                             double q[PW85_SYM]);
DllExport double pw85_contact_function(double r12[PW85_DIM],
                                       double q1[PW85_SYM], double q2[PW85_SYM],
                                       double* out);

#endif
