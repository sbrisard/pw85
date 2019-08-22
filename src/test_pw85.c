#include <glib.h>
#include <math.h>
#include <pw85.h>
#include <stdio.h>

void gen_directions(double *out) {
  const double u_abs = sqrt(2. / (5. + sqrt(5.)));
  const double v_abs = sqrt((3 + sqrt(5.)) / (5. + sqrt(5.)));
  const double u_values[] = {-u_abs, u_abs};
  const double v_values[] = {-v_abs, v_abs};
  double *n = out;
  for (size_t i = 0; i < 2; i++) {
    const double u = u_values[i];
    for (size_t j = 0; j < 2; j++) {
      const double v = v_values[j];
      n[0] = 0.;
      n[1] = u;
      n[2] = v;
      n += 3;
      n[0] = v;
      n[1] = 0.;
      n[2] = u;
      n += 3;
      n[0] = u;
      n[1] = v;
      n[2] = 0.;
      n += 3;
    }
  }
}

void test_pw85_spheroid(const double *data) {
  /*
   * Relative and absolute tolerance on the coefficients of the matrix
   * q to be computed and tested.
   */
  const double rtol = 1e-15;
  const double atol = 1e-15;

  double exp, act, tol;

  const double a = data[0];
  const double c = data[1];
  const double a2 = a * a;
  const double c2 = c * c;

  const double *n = data + 2;
  double abs_n[PW85_DIM];
  for (size_t i = 0; i < PW85_DIM; i++) {
    abs_n[i] = fabs(n[i]);
  }

  double q[PW85_SYM];
  pw85_spheroid(a, c, n, q);
  double abs_q[PW85_SYM], delta_q[PW85_SYM];
  for (size_t i = 0; i < PW85_SYM; i++) {
    abs_q[i] = fabs(q[i]);
    delta_q[i] = rtol * fabs(q[i]) + atol;
  }

  /* Check that n is an eigenvector. */
  double qn[] = {q[0] * n[0] + q[1] * n[1] + q[2] * n[2],
                 q[1] * n[0] + q[3] * n[1] + q[4] * n[2],
                 q[2] * n[0] + q[4] * n[1] + q[5] * n[2]};
  double delta_qn[] = {
      delta_q[0] * abs_n[0] + delta_q[1] * abs_n[1] + delta_q[2] * abs_n[2],
      delta_q[1] * abs_n[0] + delta_q[3] * abs_n[1] + delta_q[4] * abs_n[2],
      delta_q[2] * abs_n[0] + delta_q[4] * abs_n[1] + delta_q[5] * abs_n[2]};
  for (size_t i = 0; i < PW85_DIM; i++) {
    g_assert_cmpfloat(fabs(qn[i] - c2 * n[i]), <=, delta_qn[i]);
  }

  /*
   * Check that the eigenvalues are [a, a, c]. It is sufficient to
   * compute the characteristic polynomial of the matrix q.
   */

  /* Check tr(q). */
  exp = 2. * a2 + c2;
  act = q[0] + q[3] + q[5];
  tol = delta_q[0] + delta_q[3] + delta_q[5];
  g_assert_cmpfloat(fabs(act - exp), <=, tol);

  /* Check det(q). */
  exp = a2 * a2 * c2;
  act = q[0] * q[3] * q[5] + 2. * q[1] * q[2] * q[4] - q[0] * q[4] * q[4] -
        q[3] * q[2] * q[2] - q[5] * q[1] * q[1];
  tol =
      delta_q[0] * abs_q[3] * abs_q[5] + 2. * delta_q[1] * abs_q[2] * abs_q[4] +
      delta_q[0] * abs_q[4] * abs_q[4] + delta_q[3] * abs_q[2] * abs_q[2] +
      delta_q[5] * abs_q[1] * abs_q[1] + abs_q[0] * delta_q[3] * abs_q[5] +
      +abs_q[0] * abs_q[3] * delta_q[5] +
      2. *
          (abs_q[1] * delta_q[2] * abs_q[4] + abs_q[0] * delta_q[4] * abs_q[4] +
           abs_q[3] * delta_q[2] * abs_q[2] + abs_q[5] * delta_q[1] * abs_q[1] +
           abs_q[1] * abs_q[2] * delta_q[4]);
  g_assert_cmpfloat(fabs(act - exp), <=, tol);

  /* Check [tr(q)^2 - tr(q^2)]/2 */
  exp = a2 * (a2 + 2. * c2);
  act = q[0] * q[3] + q[3] * q[5] + q[5] * q[0] - q[1] * q[1] - q[2] * q[2] -
        q[4] * q[4];
  tol = delta_q[0] * abs_q[3] + abs_q[0] * delta_q[3] + delta_q[3] * abs_q[5] +
        abs_q[3] * delta_q[5] + delta_q[5] * abs_q[0] + abs_q[5] * delta_q[0] +
        2. * delta_q[1] * abs_q[1] + 2. * delta_q[2] * abs_q[2] +
        2. * delta_q[4] * abs_q[4];
  g_assert_cmpfloat(fabs(act - exp), <=, tol);
}

int main(int argc, char **argv) {
  g_test_init(&argc, &argv, NULL);

  const size_t num_directions = 12;
  double *directions =
      (double *)malloc(num_directions * PW85_DIM * sizeof(double));
  gen_directions(directions);

  const size_t num_radii = 3;
  const double radii[] = {0.0199, 1.999, 9.999};

  for (size_t i = 0; i < num_radii; i++) {
    const double a = radii[i];
    for (size_t j = 0; j < num_radii; j++) {
      const double c = radii[j];
      for (size_t k = 0; k < num_directions; k++) {
        const double *n = directions + (k * PW85_DIM);
        double *data = g_new(double, 2 + PW85_SYM);
        data[0] = a;
        data[1] = c;
        data[2] = n[0];
        data[3] = n[1];
        data[4] = n[2];
        char path[255];
        sprintf(path, "/pw85/spheroid/a=%g,c=%g,n=[%g,%g,%g]", a, c, n[0], n[1],
                n[2]);
        g_test_add_data_func_full(path, data, test_pw85_spheroid, g_free);
      }
    }
  }

  return g_test_run();
}
