#include <math.h>
#include <stdio.h>

#include <glib.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <pw85.h>

/* #define TEST_PW85_NUM_DIRECTIONS 12 */

/* double *test_pw85_gen_directions() { */
/*   double *directions = g_new(double, TEST_PW85_NUM_DIRECTIONS *PW85_DIM); */
/*   double u_abs = sqrt(2. / (5. + sqrt(5.))); */
/*   double v_abs = sqrt((3 + sqrt(5.)) / (5. + sqrt(5.))); */
/*   double u_values[] = {-u_abs, u_abs}; */
/*   double v_values[] = {-v_abs, v_abs}; */
/*   double *n = directions; */
/*   for (size_t i = 0; i < 2; i++) { */
/*     double u = u_values[i]; */
/*     for (size_t j = 0; j < 2; j++) { */
/*       double v = v_values[j]; */
/*       n[0] = 0.; */
/*       n[1] = u; */
/*       n[2] = v; */
/*       n += PW85_DIM; */
/*       n[0] = v; */
/*       n[1] = 0.; */
/*       n[2] = u; */
/*       n += PW85_DIM; */
/*       n[0] = u; */
/*       n[1] = v; */
/*       n[2] = 0.; */
/*       n += PW85_DIM; */
/*     } */
/*   } */
/*   return directions; */
/* } */

#define TEST_PW85_NUM_DIRECTIONS 3

double *test_pw85_gen_directions() {
  double *directions = g_new(double, TEST_PW85_NUM_DIRECTIONS *PW85_DIM);
  double u = sqrt(2. / (5. + sqrt(5.)));
  double v = sqrt((3 + sqrt(5.)) / (5. + sqrt(5.)));
  directions[0] = 0.;
  directions[1] = -u;
  directions[2] = -v;
  directions[3] = -v;
  directions[4] = 0.;
  directions[5] = u;
  directions[6] = u;
  directions[7] = -v;
  directions[8] = 0.;
  return directions;
}

double *test_pw85_gen_radius_vectors(size_t num_distances, double *distances,
                                     size_t num_directions,
                                     double *directions) {
  size_t num_radius_vectors = num_distances * num_directions;
  size_t num_doubles = num_radius_vectors * PW85_DIM;
  double *radius_vectors = g_new(double, num_doubles);
  double *r12 = radius_vectors;
  for (size_t i = 0; i < num_distances; i++) {
    double r = distances[i];
    double *n = directions;
    for (size_t j = 0; j < num_directions; j++) {
      for (size_t k = 0; k < PW85_DIM; k++) {
        r12[k] = r * n[k];
      }
      n += PW85_DIM;
      r12 += PW85_DIM;
    }
  }
  return radius_vectors;
}

double *test_pw85_gen_spheroids(size_t num_radii, double *radii,
                                size_t num_directions, double *directions) {
  size_t num_spheroids = num_radii * num_radii * num_directions;
  size_t num_doubles = num_spheroids * PW85_SYM;
  double *spheroids = g_new(double, num_doubles);
  double *q = spheroids;
  for (size_t i = 0; i < num_radii; i++) {
    double a = radii[i];
    for (size_t j = 0; j < num_radii; j++) {
      double c = radii[j];
      double *n = directions;
      for (size_t k = 0; k < num_directions; k++) {
        pw85_spheroid(a, c, n, q);
        n += PW85_DIM;
        q += PW85_SYM;
      }
    }
  }
  return spheroids;
}

void test_pw85_spheroid_test(double const *data) {
  /*
   * Relative and absolute tolerance on the coefficients of the matrix
   * q to be computed and tested.
   */
  double rtol = 1e-15;
  double atol = 1e-15;

  double exp, act, tol;

  double a = data[0];
  double c = data[1];
  double a2 = a * a;
  double c2 = c * c;

  double const *n = data + 2;
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

void test_pw85_contact_function_test(double const *data) {
  double atol = 1e-15;
  double rtol = 1e-10;

  double const *r12 = data;
  double const *q1 = data + PW85_DIM;
  double const *q2 = data + PW85_DIM + PW85_SYM;
  double out[2];
  double mu2 = pw85_contact_function(r12, q1, q2, out);
  g_assert_cmpfloat(mu2, ==, out[0]);
  double lambda = out[1];

  gsl_vector *r = gsl_vector_calloc(PW85_DIM);
  gsl_vector *s = gsl_vector_calloc(PW85_DIM);
  gsl_vector *u = gsl_vector_calloc(PW85_DIM);
  gsl_vector *v = gsl_vector_calloc(PW85_DIM);
  gsl_vector *dx0 = gsl_vector_calloc(PW85_DIM);
  gsl_vector *c1x0 = gsl_vector_calloc(PW85_DIM);
  gsl_vector *c2x0 = gsl_vector_calloc(PW85_DIM);

  gsl_matrix *Q1 = gsl_matrix_alloc(PW85_DIM, PW85_DIM);
  gsl_matrix *Q2 = gsl_matrix_alloc(PW85_DIM, PW85_DIM);
  gsl_matrix *Q = gsl_matrix_alloc(PW85_DIM, PW85_DIM);

  for (size_t i = 0, ij = 0; i < PW85_DIM; i++) {
    gsl_vector_set(r, i, r12[i]);
    for (size_t j = i; j < PW85_DIM; j++, ij++) {
      gsl_matrix_set(Q1, i, j, q1[ij]);
      gsl_matrix_set(Q1, j, i, q1[ij]);
      gsl_matrix_set(Q2, i, j, q2[ij]);
      gsl_matrix_set(Q2, j, i, q2[ij]);
      gsl_matrix_set(Q, i, j, (1 - lambda) * q1[ij] + lambda * q2[ij]);
      gsl_matrix_set(Q, j, i, (1 - lambda) * q1[ij] + lambda * q2[ij]);
    }
  }

  /*
   * Let x0 be the contact point. It is found from the following
   * identities
   *
   *     c1x0 = (1-lambda) Q1.s,
   *     c2x0 = -lambda Q2.s,
   *
   * where r = Q.s.
   *
   * We comupte dx0, which is the difference between the coordinates
   * of x0 found from the above two formulas:
   *
   *     dx0 = r12 - c1x0 + c2x0
   *         = r12 - (1-lambda) Q1.s - lambda Q2.s.
   */

  gsl_linalg_cholesky_decomp1(Q);
  gsl_linalg_cholesky_solve(Q, r, s);

  gsl_blas_dgemv(CblasNoTrans, 1. - lambda, Q1, s, 0., c1x0);
  gsl_blas_dgemv(CblasNoTrans, -lambda, Q2, s, 0., c2x0);

  gsl_vector_memcpy(dx0, r);
  gsl_blas_dgemv(CblasNoTrans, -1. + lambda, Q1, s, 1., dx0);
  gsl_blas_dgemv(CblasNoTrans, -lambda, Q2, s, 1., dx0);

  /*
   * We first assert that dx0 is small compared to r12.
   */
  g_assert_cmpfloat(gsl_blas_dnrm2(dx0), <=, rtol * gsl_blas_dnrm2(r) + atol);

  /*
   * We could also check that x0 belong to both (scaled) ellipsoids:
   *
   *     mu^2 = c1x0.Q1^(-1).c1x0 = (1-lambda) c1x0.s
   *
   * and
   *
   *     mu^2 = c2x0.Q2^(-1).c2x0 = -lambda c2x0.s.
   *
   * The code is provided below. However, the accuracy seems to be
   * quite poor.
   */
  double rtol_mu2 = 1e-4;
  double atol_mu2 = 1e-15;

  double mu2_1;
  gsl_blas_ddot(c1x0, s, &mu2_1);
  mu2_1 *= 1. - lambda;

  g_assert_cmpfloat(fabs(mu2_1 - mu2), <=, rtol_mu2 * mu2 + atol_mu2);

  double mu2_2;
  gsl_blas_ddot(c2x0, s, &mu2_2);
  mu2_2 *= -lambda;
  g_assert_cmpfloat(fabs(mu2_2 - mu2), <=, rtol_mu2 * mu2 + atol_mu2);

  /*
   * Finally, we check that f'(lambda) = 0.
   *
   * The tolerance is given by the tolerance eps on lambda (defined in
   * the call to gsl_min_test_interval in the function
   * pw85_contact_function.
   *
   * The test therefore reads
   *
   *     |f'(lambda)| <= eps * f"(lambda).
   *
   * We first compute f' and f".
   */

  gsl_blas_dgemv(CblasNoTrans, 1., Q2, s, 0., u);
  gsl_blas_dgemv(CblasNoTrans, -1., Q1, s, 1., u);
  gsl_linalg_cholesky_solve(Q, u, v);

  double rs, su, uv;
  gsl_blas_ddot(r, s, &rs);
  gsl_blas_ddot(s, u, &su);
  gsl_blas_ddot(u, v, &uv);

  double f1 = (1. - 2. * lambda) * rs - lambda * (1. - lambda) * su;
  double f2 = -2. * rs - 2. * (1. - 2. * lambda) * su +
              2. * lambda * (1. - lambda) * uv;

  g_assert_cmpfloat(fabs(f1), <=, 1000. * PW85_LAMBDA_ATOL * fabs(f2));

  gsl_vector_free(r);
  gsl_vector_free(s);
  gsl_vector_free(u);
  gsl_vector_free(v);
  gsl_vector_free(dx0);
  gsl_vector_free(c1x0);
  gsl_vector_free(c2x0);
  gsl_matrix_free(Q1);
  gsl_matrix_free(Q2);
  gsl_matrix_free(Q);
}

int main(int argc, char **argv) {
  g_test_init(&argc, &argv, NULL);

  size_t num_radii = 3;
  double radii[] = {0.0199, 1.999, 9.999};

  size_t num_distances = 3;
  double distances[] = {0.15, 1.1, 11.};

  size_t num_directions = TEST_PW85_NUM_DIRECTIONS;
  double *directions = test_pw85_gen_directions();

  size_t num_radius_vectors = num_distances * num_directions;
  double *radius_vectors = test_pw85_gen_radius_vectors(
      num_distances, distances, num_directions, directions);

  size_t num_spheroids = num_radii * num_radii * num_directions;
  double *spheroids =
      test_pw85_gen_spheroids(num_radii, radii, num_directions, directions);

  for (size_t i = 0; i < num_radii; i++) {
    double a = radii[i];
    for (size_t j = 0; j < num_radii; j++) {
      double c = radii[j];
      double *n = directions;
      for (size_t k = 0; k < num_directions; k++) {
        double *data = g_new(double, 2 + PW85_SYM);
        data[0] = a;
        data[1] = c;
        data[2] = n[0];
        data[3] = n[1];
        data[4] = n[2];
        char path[255];
        sprintf(path, "/pw85/spheroid/a=%g,c=%g,n=[%g,%g,%g]", a, c, n[0], n[1],
                n[2]);
        g_test_add_data_func_full(path, data, test_pw85_spheroid_test, g_free);

        n += PW85_DIM;
      }
    }
  }

  double *r12 = radius_vectors;
  for (size_t i = 0; i < num_radius_vectors; i++) {
    double *q1 = spheroids;
    for (size_t j1 = 0; j1 < num_spheroids; j1++) {
      double *q2 = spheroids;
      for (size_t j2 = 0; j2 < num_spheroids; j2++) {
        double *data = g_new(double, PW85_DIM + 2 * PW85_SYM);
        for (size_t k = 0; k < PW85_DIM; k++) {
          data[k] = r12[k];
        }
        for (size_t k = 0; k < PW85_SYM; k++) {
          data[PW85_DIM + k] = q1[k];
          data[PW85_DIM + PW85_SYM + k] = q2[k];
        }

        char path[255];
        sprintf(path, "/pw85/contact_function(r12[%d],q1[%d],q2[%d])", (int)i,
                (int)j1, (int)j2);
        g_test_add_data_func_full(path, data, test_pw85_contact_function_test,
                                  g_free);

        q2 += PW85_SYM;
      }
      q1 += PW85_SYM;
    }
    r12 += PW85_DIM;
  }

  g_free(directions);
  g_free(radius_vectors);
  g_free(spheroids);

  return g_test_run();
}
