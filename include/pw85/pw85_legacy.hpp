#ifndef __PW85_LEGACY_H__
#define __PW85_LEGACY_H__
#include "pw85/pw85.hpp"

namespace pw85_legacy {
DllExport double _det_sym(const double* a);
DllExport double _xT_adjA_x(const double* x, const double* a);
DllExport void _rT_adjQ_r_as_poly(const double* r, const double* q1,
                                  const double* q2, const double* q3,
                                  double* a);
DllExport void _detQ_as_poly(const double* q1, const double* q2,
                             const double* q3, const double* q4, double* b);
DllExport double f1(double lambda, const double* r12, const double* q1,
                    const double* q2, double* out);
DllExport double f2(double lambda, const double* r12, const double* q1,
                    const double* q2, double* out);
DllExport int contact_function1(const double* r12, const double* q1,
                                const double* q2, double* out);
DllExport int contact_function2(const double* r12, const double* q1,
                                const double* q2, double* out);
}  // namespace pw85_legacy
#endif
