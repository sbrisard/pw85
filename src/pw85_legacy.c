#include "pw85_legacy.h"

double pw85_legacy__det_sym(double const a[PW85_SYM]) {
  return a[0] * a[3] * a[5] + 2 * a[1] * a[2] * a[4] - a[0] * a[4] * a[4] -
         a[3] * a[2] * a[2] - a[5] * a[1] * a[1];
}

double pw85_legacy__xT_adjA_x(double const x[PW85_DIM],
                              double const a[PW85_SYM]) {
  return (x[0] * x[0] * (a[3] * a[5] - a[4] * a[4]) +
          x[1] * x[1] * (a[0] * a[5] - a[2] * a[2]) +
          x[2] * x[2] * (a[0] * a[3] - a[1] * a[1]) +
          2. * (x[0] * x[1] * (a[2] * a[4] - a[1] * a[5]) +
                x[0] * x[2] * (a[1] * a[4] - a[2] * a[3]) +
                x[1] * x[2] * (a[1] * a[2] - a[0] * a[4])));
}

void pw85_legacy__rT_adjQ_r_as_poly(double const r[PW85_DIM],
                                    double const q1[PW85_SYM],
                                    double const q2[PW85_SYM],
                                    double const q3[PW85_SYM],
                                    double a[PW85_DIM]) {
  double const a_zero = pw85_legacy__xT_adjA_x(r, q1);
  double const a_one = pw85_legacy__xT_adjA_x(r, q2);
  double const a_minus_one = pw85_legacy__xT_adjA_x(r, q3);
  a[0] = a_zero;
  a[2] = 0.5 * (a_one + a_minus_one) - a_zero;
  a[1] = 0.5 * (a_one - a_minus_one);
}

void pw85_legacy__detQ_as_poly(double const q1[PW85_SYM],
                               double const q2[PW85_SYM],
                               double const q3[PW85_SYM],
                               double const q4[PW85_SYM],
                               double b[PW85_DIM + 1]) {
  double const b_zero = pw85_legacy__det_sym(q1);
  double const b_one = pw85_legacy__det_sym(q2);
  /* Compute det[(1-x)*q1+x*q2] for x = -1. */
  double const b_minus_one = pw85_legacy__det_sym(q3);
  /* Compute det[(1-x)*q1+x*q2] for x = 1/2. */
  double const b_one_half = pw85_legacy__det_sym(q4);
  b[0] = b_zero;
  b[2] = 0.5 * (b_one + b_minus_one) - b_zero;
  b[1] = (8. * b_one_half - 6. * b_zero - 1.5 * b_one - 0.5 * b_minus_one) / 3.;
  b[3] = (-8. * b_one_half + 6. * b_zero + 3. * b_one - b_minus_one) / 3.;
}

double pw85_legacy_f1(double lambda, double const r12[PW85_DIM],
                      double const q1[PW85_SYM], double const q2[PW85_SYM],
                      double *out) {
  double q[PW85_SYM], q12[PW85_SYM];
  double const *q1_i = q1;
  double const *q2_i = q2;
  double *q_i = q;
  double *q12_i = q12;
  for (size_t i = 0; i < PW85_SYM; i++) {
    *q_i = (1 - lambda) * (*q1_i) + lambda * (*q2_i);
    *q12_i = (*q2_i) - (*q1_i);
    ++q1_i;
    ++q2_i;
    ++q_i;
    ++q12_i;
  }

  double l[PW85_SYM];
  pw85__cholesky_decomp(q, l);
  double s[PW85_DIM];
  pw85__cholesky_solve(l, r12, s);
  double u[] = {q12[0] * s[0] + q12[1] * s[1] + q12[2] * s[2],
                q12[1] * s[0] + q12[3] * s[1] + q12[4] * s[2],
                q12[2] * s[0] + q12[4] * s[1] + q12[5] * s[2]};
  double v[PW85_DIM];
  pw85__cholesky_solve(l, u, v);

  double const *r_i = r12;
  double *s_i = s;
  double *u_i = u;
  double *v_i = v;
  double rs = 0., su = 0., uv = 0.;
  for (size_t i = 0; i < PW85_DIM; i++) {
    rs += (*r_i) * (*s_i);
    su += (*s_i) * (*u_i);
    uv += (*u_i) * (*v_i);
    ++r_i;
    ++s_i;
    ++u_i;
    ++v_i;
  }

  double const aux1 = lambda * (1. - lambda);
  double const aux2 = 1. - 2. * lambda;
  if (out) {
    out[0] = aux1 * rs;
    out[1] = aux2 * rs - aux1 * su;
    out[2] = -2. * (rs + aux2 * su - aux1 * uv);
    return out[0];
  } else {
    return aux1 * rs;
  }
}

double pw85_legacy_f2(double lambda, double const r12[PW85_DIM],
                      double const q1[PW85_SYM], double const q2[PW85_SYM],
                      double *out) {
  double q3[PW85_SYM]; /* q3 = 2*q1-q2. */
  double q4[PW85_SYM]; /* q4 = (q1+q2)/2. */
  for (int i = 0; i < PW85_SYM; i++) {
    q3[i] = 2. * q1[i] - q2[i];
    q4[i] = 0.5 * (q1[i] + q2[i]);
  }

  double a[3];
  pw85_legacy__rT_adjQ_r_as_poly(r12, q1, q2, q3, a);

  double b[4];
  pw85_legacy__detQ_as_poly(q1, q2, q3, q4, b);

  double const y = lambda * (1. - lambda) *
                   (a[0] + lambda * (a[1] + lambda * a[2])) /
                   (b[0] + lambda * (b[1] + lambda * (b[2] + lambda * b[3])));
  if (out) out[0] = y;
  return y;
}

double pw85_legacy_contact_function1(double const r12[PW85_DIM],
                                     double const q1[PW85_SYM],
                                     double const q2[PW85_SYM], double *out) {
  double f[3];
  /* Lower-bound for the root. */
  double a = 0., f_prime_a;
  pw85_legacy_f1(a, r12, q1, q2, f);
  f_prime_a = f[1];
  /* Upper-bound for the root. */
  double b = 1., f_prime_b;
  pw85_legacy_f1(b, r12, q1, q2, f);
  f_prime_b = f[1];
  printf("--------------------\n");
  printf("%g\t%g\n", f_prime_a, f_prime_b);
  /* Current estimate of the root f'(x) == 0. */
  double x = (a * f_prime_b - b * f_prime_a) / (f_prime_b - f_prime_a);
  for (int i = 0; i < 100; i++) {
    pw85_legacy_f1(x, r12, q1, q2, f);
    if (f[1] * f_prime_a < 0.) {
      b = x;
      f_prime_b = f[1];
    } else {
      a = x;
      f_prime_a = f[1];
    }
    printf("%d %g %g %g %g\n", i, x, f[0], f[1], f[2]);
    double x_new = x - f[2] / f[1];
    if ((x_new < a) || (x_new > b)) {
      x_new = (a * f_prime_b - b * f_prime_a) / (f_prime_b - f_prime_a);
    }
    x = x_new;
  }
  pw85_legacy_f1(x, r12, q1, q2, f);
  if (out) {
    out[0] = f[0];
    out[1] = x;
  }
  return f[0];
}

double pw85_legacy_contact_function2(double const r12[PW85_DIM],
                                     double const q1[PW85_SYM],
                                     double const q2[PW85_SYM], double *out) {
  double q3[PW85_SYM]; /* q3 = 2*q1-q2. */
  double q4[PW85_SYM]; /* q4 = (q1+q2)/2. */
  for (int i = 0; i < PW85_SYM; i++) {
    q3[i] = 2. * q1[i] - q2[i];
    q4[i] = 0.5 * (q1[i] + q2[i]);
  }

  double a[3];
  pw85_legacy__rT_adjQ_r_as_poly(r12, q1, q2, q3, a);
  printf("a = \n");
  for (size_t i = 0; i < 3; i++) printf("%1.15f\n", a[i]);

  double b[4];
  pw85_legacy__detQ_as_poly(q1, q2, q3, q4, b);
  printf("b = \n");
  for (size_t i = 0; i < 4; i++) printf("%1.15f\n", b[i]);

  double const c0 = a[0] * b[0];
  double const c1 = 2. * (a[1] - a[0]) * b[0];
  double const c2 =
      -a[0] * (b[1] + b[2]) + 3. * b[0] * (a[2] - a[1]) + a[1] * b[1];
  double const c3 =
      2. * (b[1] * (a[2] - a[1]) - a[0] * b[3]) - 4. * a[2] * b[0];
  double const c4 =
      (a[0] - a[1]) * b[3] + (a[2] - a[1]) * b[2] - 3. * a[2] * b[1];
  double const c5 = -2. * a[2] * b[2];
  double const c6 = -a[2] * b[3];
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
      if (yl * y < 0) {
        xr = x;
        yr = y;
      } else {
        xl = x;
        yl = y;
      }
    }
  }
  y = x * (1. - x) * (a[0] + x * (a[1] + x * a[2])) /
      (b[0] + x * (b[1] + x * (b[2] + x * b[3])));
  if (out != NULL) {
    out[0] = y;
    out[1] = x;
  }
  return y;
}
