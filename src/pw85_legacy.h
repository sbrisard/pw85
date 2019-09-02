#ifndef __PW85_LEGACY_H__
#define __PW85_LEGACY_H__
#include "pw85.h"

DllExport double pw85_legacy__det_sym(double const a[PW85_SYM]);
DllExport double pw85_legacy__xT_adjA_x(double const x[PW85_DIM],
                                        double const a[PW85_SYM]);
DllExport void pw85_legacy__rT_adjQ_r_as_poly(double const r[PW85_DIM],
                                              double const q1[PW85_SYM],
                                              double const q2[PW85_SYM],
                                              double const q3[PW85_SYM],
                                              double a[PW85_DIM]);
DllExport void pw85_legacy__detQ_as_poly(double const q1[PW85_SYM],
                                         double const q2[PW85_SYM],
                                         double const q3[PW85_SYM],
                                         double const q4[PW85_SYM],
                                         double b[PW85_DIM + 1]);
DllExport double pw85_legacy_f1(double lambda, double const r12[PW85_DIM],
                                double const q1[PW85_SYM],
                                double const q2[PW85_SYM], double* out);
DllExport double pw85_legacy_f2(double lambda, double const r12[PW85_DIM],
                                double const q1[PW85_SYM],
                                double const q2[PW85_SYM], double* out);
DllExport int pw85_legacy_contact_function1(double const r12[PW85_DIM],
                                            double const q1[PW85_SYM],
                                            double const q2[PW85_SYM],
                                            double out[2]);
DllExport int pw85_legacy_contact_function2(double const r12[PW85_DIM],
                                            double const q1[PW85_SYM],
                                            double const q2[PW85_SYM],
                                            double out[2]);
#endif
