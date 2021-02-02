#include "pw85/pw85.hpp"
#include <boost/math/tools/minima.hpp>

namespace pw85 {
void _cholesky_decomp(const double *a, double *l) {
  l[0] = sqrt(a[0]);
  l[1] = a[1] / l[0];
  l[2] = a[2] / l[0];
  l[3] = sqrt(a[3] - l[1] * l[1]);
  l[4] = (a[4] - l[1] * l[2]) / l[3];
  l[5] = sqrt(a[5] - l[2] * l[2] - l[4] * l[4]);
}

void _cholesky_solve(const double *l, const double *b, double *x) {
  /* Solve L.y = b */
  double const y0 = b[0] / l[0];
  double const y1 = (b[1] - l[1] * y0) / l[3];
  double const y2 = (b[2] - l[2] * y0 - l[4] * y1) / l[5];

  /* Solve L^T.x = y */
  x[2] = y2 / l[5];
  x[1] = (y1 - l[4] * x[2]) / l[3];
  x[0] = (y0 - l[1] * x[1] - l[2] * x[2]) / l[0];
}

void spheroid(double a, double c, const double *n, double *q) {
  double const a2 = a * a;
  double const c2_minus_a2 = c * c - a2;
  double const nx = n[0];
  double const ny = n[1];
  double const nz = n[2];
  q[0] = nx * nx * c2_minus_a2 + a2;
  q[3] = ny * ny * c2_minus_a2 + a2;
  q[5] = nz * nz * c2_minus_a2 + a2;
  q[4] = ny * nz * c2_minus_a2;
  q[2] = nx * nz * c2_minus_a2;
  q[1] = nx * ny * c2_minus_a2;
}

double f_neg(double lambda, const double *params) {
  double const *r12 = params;
  double const *q1_i = params + PW85_DIM;
  double const *q2_i = q1_i + PW85_SYM;
  double q[PW85_SYM];
  double *q_i = q;
  double q12[PW85_SYM];
  double *q12_i = q12;
  for (size_t i = 0; i < PW85_SYM; i++, q1_i++, q2_i++, q_i++, q12_i++) {
    *q_i = (1 - lambda) * (*q1_i) + lambda * (*q2_i);
    *q12_i = (*q2_i) - (*q1_i);
  }
  double l[PW85_SYM];
  _cholesky_decomp(q, l);
  double s[PW85_DIM];
  _cholesky_solve(l, r12, s);
  double const *r_i = r12;
  double *s_i = s;
  double rs = 0.;
  for (size_t i = 0; i < PW85_DIM; i++, r_i++, s_i++) {
    rs += (*r_i) * (*s_i);
  }
  return -lambda * (1. - lambda) * rs;
}

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