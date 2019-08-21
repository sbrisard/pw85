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
  const double a = data[0];
  const double c = data[1];
  const double *n = data + 2;
  double q[PW85_SYM];
  pw85_spheroid(a, c, n, q);

  /* Check that n is an eigenvector. */
  double qn[] = {q[0] * n[0] + q[1] * n[1] + q[2] * n[2],
                 q[1] * n[0] + q[3] * n[1] + q[4] * n[2],
                 q[2] * n[0] + q[4] * n[1] + q[5] * n[2]};
  for (size_t i = 0; i < PW85_DIM; i++) {
    const double act = qn[i];
    const double exp = c * c * n[i];
    const double rtol = 1e-10;
    const double atol = 1e-15;
    g_assert_cmpfloat(fabs(act - exp), <=, rtol * fabs(exp) + atol);
  }
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
