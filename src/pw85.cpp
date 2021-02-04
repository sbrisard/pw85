#include "pw85/pw85.hpp"
#include <boost/math/tools/minima.hpp>

namespace pw85 {

void _residual(double lambda, const double *r12, const double *q1,
               const double *q2, double *out) {
  double q[PW85_SYM];
  double q12[PW85_SYM];
  for (size_t i = 0; i < PW85_SYM; i++) {
    q[i] = (1. - lambda) * q1[i] + lambda * q2[i];
    q12[i] = q2[i] - q1[i];
  }
  double l[PW85_SYM];
  _cholesky_decomp(q, l);
  double s[PW85_DIM];
  _cholesky_solve(l, r12, s);
  double u[] = {q12[0] * s[0] + q12[1] * s[1] + q12[2] * s[2],
                q12[1] * s[0] + q12[3] * s[1] + q12[4] * s[2],
                q12[2] * s[0] + q12[4] * s[1] + q12[5] * s[2]};
  double v[PW85_DIM];
  _cholesky_solve(l, u, v);
  double rs = r12[0] * s[0] + r12[1] * s[1] + r12[2] * s[2];
  double su = s[0] * u[0] + s[1] * u[1] + s[2] * u[2];
  double uv = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];

  out[0] = lambda * (1. - lambda) * rs;
  out[1] = (2 * lambda - 1.) * rs + lambda * (1. - lambda) * su;
  out[2] =
      2. * rs + 2. * (1. - 2. * lambda) * su - 2. * lambda * (1. - lambda) * uv;
}

int contact_function(const double *r12, const double *q1, const double *q2,
                     double *out) {
  double params[] = {r12[0], r12[1], r12[2], q1[0], q1[1], q1[2], q1[3], q1[4],
                     q1[5],  q2[0],  q2[1],  q2[2], q2[3], q2[4], q2[5]};

  auto f = [params](double lambda) { return f_neg(lambda, params); };

  auto r = boost::math::tools::brent_find_minima(
      f, 0., 1., std::numeric_limits<double>::digits / 2);
  auto lambda_brent = r.first;
  auto mu2_brent = -r.second;

  double out_res[3];
  _residual(lambda_brent, r12, q1, q2, out_res);
  double res_brent = fabs(out_res[1]);

  /* Try to refine the estimate. */
  double lambda_nr = lambda_brent;
  double mu2_nr = mu2_brent;
  double res_nr = res_brent;

  for (size_t i = 0; i < PW85_NR_ITER; i++) {
    double lambda_trial = lambda_nr - out_res[1] / out_res[2];
    if ((lambda_trial < 0.) || (lambda_trial > 1.)) {
      break;
    }
    _residual(lambda_trial, r12, q1, q2, out_res);
    lambda_nr = lambda_trial;
    mu2_nr = out_res[0];
    res_nr = fabs(out_res[1]);
  }

  if (res_nr < res_brent) {
    out[0] = mu2_nr;
    out[1] = lambda_nr;
  } else {
    out[0] = mu2_brent;
    out[1] = lambda_brent;
  }

  /* TODO: return error code. */
  return 0;
}
}  // namespace pw85