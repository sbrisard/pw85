#include <math.h>
#include <stdio.h>

#include <glib.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <hdf5_hl.h>

#include <pw85.h>

void test_pw85_read_dataset_double(hid_t const hid, char const *dset_name,
                                   size_t *size, double **buffer) {
  int ndims;
  H5LTget_dataset_ndims(hid, dset_name, &ndims);
  hsize_t *dim = g_new(hsize_t, ndims);
  H5LTget_dataset_info(hid, dset_name, dim, NULL, NULL);
  *size = 1;
  for (size_t i = 0; i < ndims; i++) {
    *size *= dim[i];
  }
  *buffer = g_new(double, *size);
  H5LTread_dataset_double(hid, dset_name, *buffer);
  g_free(dim);
}

/* This is a global variable that holds the parameters used for most
 * test-cases. */

struct {
  size_t num_radii;
  double *radii;
  size_t num_directions;
  double *directions;
  size_t num_spheroids;
  double *spheroids;
  size_t num_lambdas;
  double *lambdas;
  size_t num_f;
  double *f;
  size_t num_distances;
  double *distances;
} test_pw85_context;

void test_pw85_init_context(hid_t const hid) {
  test_pw85_read_dataset_double(hid, "/directions",
                                &test_pw85_context.num_directions,
                                &test_pw85_context.directions);
  test_pw85_context.num_directions /= PW85_DIM;

  test_pw85_read_dataset_double(hid, "/radii", &test_pw85_context.num_radii,
                                &test_pw85_context.radii);

  test_pw85_read_dataset_double(hid, "/spheroids",
                                &test_pw85_context.num_spheroids,
                                &test_pw85_context.spheroids);
  test_pw85_context.num_spheroids /= PW85_SYM;

  test_pw85_read_dataset_double(hid, "/lambdas", &test_pw85_context.num_lambdas,
                                &test_pw85_context.lambdas);

  test_pw85_read_dataset_double(hid, "/F", &test_pw85_context.num_f,
                                &test_pw85_context.f);

  test_pw85_context.num_distances = 3;
  test_pw85_context.distances = g_new(double, test_pw85_context.num_distances);
  test_pw85_context.distances[0] = 0.15;
  test_pw85_context.distances[1] = 1.1;
  test_pw85_context.distances[2] = 11.;
}

void test_pw85_free_context() {
  g_free(test_pw85_context.directions);
  g_free(test_pw85_context.radii);
  g_free(test_pw85_context.spheroids);
  g_free(test_pw85_context.lambdas);
  g_free(test_pw85_context.f);
  g_free(test_pw85_context.distances);
}

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

void test_pw85_cholesky_decomp_test(double const *data) {
  double const *a = data;
  double const *exp = data + PW85_SYM;
  double const rtol = exp[PW85_SYM];

  double act[PW85_SYM];
  pw85__cholesky_decomp(data, act);

  double *exp_i = exp;
  double *act_i = act;
  for (size_t i = 0; i < PW85_SYM; ++i, ++exp_i, ++act_i) {
    g_assert_cmpfloat(fabs(*exp_i - *act_i), <=, rtol * fabs(*act_i));
  }
}

void pw85_test_add_cholesky_decomp_test() {
  char path_template[] =
      "/pw85/cholesky_decomp/a=[%1.3e,%1.3e,%1.3e,%1.3e,%1.3e,%1.3e]";
  char path[92];
  double *data1 = g_new(double, 2 * PW85_SYM + 1);
  data1[0] = 4;
  data1[1] = 2;
  data1[2] = 6;
  data1[3] = 17;
  data1[4] = 23;
  data1[5] = 70;
  data1[6] = 2;
  data1[7] = 1;
  data1[8] = 3;
  data1[9] = 4;
  data1[10] = 5;
  data1[11] = 6;
  data1[12] = 1e-15;

  sprintf(path, path_template, data1[0], data1[1], data1[2], data1[3], data1[4],
          data1[5]);
  g_test_add_data_func_full(path, data1, test_pw85_cholesky_decomp_test,
                            g_free);

  double *data2 = g_new(double, 2 * PW85_SYM + 1);
  data2[0] = 4;
  data2[1] = -2;
  data2[2] = 6;
  data2[3] = 17;
  data2[4] = -23;
  data2[5] = 70;
  data2[6] = 2;
  data2[7] = -1;
  data2[8] = 3;
  data2[9] = 4;
  data2[10] = -5;
  data2[11] = 6;
  data2[12] = 1e-15;

  sprintf(path, path_template, data2[0], data2[1], data2[2], data2[3], data2[4],
          data2[5]);
  g_test_add_data_func_full(path, data2, test_pw85_cholesky_decomp_test,
                            g_free);

  double *data3 = g_new(double, 2 * PW85_SYM + 1);
  data3[0] = 1e10;
  data3[1] = -2;
  data3[2] = -3;
  data3[3] = 16 + 1. / 25e8;
  data3[4] = -0.02 + 3. / 5e9;
  data3[5] = 29. / 8e3 - 9e-10;
  data3[6] = 1e5;
  data3[7] = -2e-5;
  data3[8] = -3e-5;
  data3[9] = 4;
  data3[10] = -5e-3;
  data3[11] = 6e-2;
  data3[12] = 1e-6;
  sprintf(path, path_template, data3[0], data3[1], data3[2], data3[3], data3[4],
          data3[5]);
  g_test_add_data_func_full(path, data3, test_pw85_cholesky_decomp_test,
                            g_free);
}

void test_pw85_cholesky_solve_test(double const *data) {
  double const *l = data;
  double const *b = data + PW85_SYM;
  double const *exp = b + PW85_DIM;
  double const rtol = exp[PW85_DIM];

  double const act[PW85_DIM];
  pw85__cholesky_solve(l, b, act);

  double *exp_i = exp;
  double *act_i = act;
  for (size_t i = 0; i < PW85_DIM; ++i, ++exp_i, ++act_i) {
    g_assert_cmpfloat(fabs(*exp_i - *act_i), <=, rtol * fabs(*act_i));
  }
}

void pw85_test_add_cholesky_solve_test() {
  char path_template[] =
      "/pw85/cholesky_solve/l=[%g,%g,%g,%g,%g,%g],b=[%g,%g,%g]";
  char path[255];
  double *data1 = g_new(double, PW85_SYM + 2 * PW85_DIM + 1);
  /* L */
  data1[0] = 1;
  data1[1] = 2;
  data1[2] = 3;
  data1[3] = 4;
  data1[4] = 5;
  data1[5] = 6;
  /* b */
  data1[6] = 11.5;
  data1[7] = 82.6;
  data1[8] = 314.2;
  /* x */
  data1[9] = 1.2;
  data1[10] = -3.4;
  data1[11] = 5.7;
  /* rtol */
  data1[12] = 1e-15;

  sprintf(path, path_template, data1[0], data1[1], data1[2], data1[3], data1[4],
          data1[5], data1[6], data1[7], data1[8]);
  g_test_add_data_func_full(path, data1, test_pw85_cholesky_solve_test, g_free);

  double *data2 = g_new(double, PW85_SYM + 2 * PW85_DIM + 1);
  /* L */
  data2[0] = 1;
  data2[1] = -2;
  data2[2] = -3;
  data2[3] = 4;
  data2[4] = -5;
  data2[5] = 6;
  /* b */
  data2[6] = -9.1;
  data2[7] = -150.2;
  data2[8] = 443;
  /* x */
  data2[9] = 1.2;
  data2[10] = -3.4;
  data2[11] = 5.7;
  /* rtol */
  data2[12] = 4e-15;

  sprintf(path, path_template, data2[0], data2[1], data2[2], data2[3], data2[4],
          data2[5], data2[6], data2[7], data2[8]);
  g_test_add_data_func_full(path, data2, test_pw85_cholesky_solve_test, g_free);
}

void test_pw85_spheroid_elementary_test(double a, double c,
                                        double const n[PW85_DIM]) {
  /*
   * Relative and absolute tolerance on the coefficients of the matrix
   * q to be computed and tested.
   */
  double rtol = 1e-15;
  double atol = 1e-15;

  double exp, act, tol;

  double a2 = a * a;
  double c2 = c * c;

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

void test_pw85_spheroid_test() {
  for (size_t i = 0; i < test_pw85_context.num_radii; i++) {
    double const a = test_pw85_context.radii[i];
    for (size_t j = 0; j < test_pw85_context.num_radii; j++) {
      double const c = test_pw85_context.radii[j];
      for (size_t k = 0; k < test_pw85_context.num_directions; k++) {
        double *const n = test_pw85_context.directions + PW85_DIM * k;
        test_pw85_spheroid_elementary_test(a, c, n);
      }
    }
  }
}

void test_pw85_contact_function_elementary_test(double r12[PW85_DIM],
                                                double q1[PW85_SYM],
                                                double q2[PW85_SYM]) {
  double atol = 1e-15;
  double rtol = 1e-10;

  double out[2];
  pw85_contact_function(r12, q1, q2, out);
  double mu2 = out[0];
  double lambda = out[1];
  double lambda1 = 1. - lambda;

  gsl_vector *r = gsl_vector_calloc(PW85_DIM);

  gsl_matrix *Q12 = gsl_matrix_alloc(PW85_DIM, PW85_DIM);
  gsl_matrix *Q = gsl_matrix_alloc(PW85_DIM, PW85_DIM);

  for (size_t i = 0, ij = 0; i < PW85_DIM; i++) {
    gsl_vector_set(r, i, r12[i]);
    for (size_t j = i; j < PW85_DIM; j++, ij++) {
      gsl_matrix_set(Q12, i, j, q2[ij] - q1[ij]);
      gsl_matrix_set(Q12, j, i, q2[ij] - q1[ij]);
      gsl_matrix_set(Q, i, j, lambda1 * q1[ij] + lambda * q2[ij]);
      gsl_matrix_set(Q, j, i, lambda1 * q1[ij] + lambda * q2[ij]);
    }
  }

  gsl_linalg_cholesky_decomp1(Q);

  gsl_vector *s = gsl_vector_calloc(PW85_DIM);
  gsl_linalg_cholesky_solve(Q, r, s);

  gsl_vector *u = gsl_vector_calloc(PW85_DIM);
  gsl_blas_dgemv(CblasNoTrans, 1., Q12, s, 0., u);

  gsl_vector *v = gsl_vector_calloc(PW85_DIM);
  gsl_linalg_cholesky_solve(Q, u, v);

  double rs, su, uv;
  gsl_blas_ddot(r, s, &rs);
  gsl_blas_ddot(s, u, &su);
  gsl_blas_ddot(u, v, &uv);

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

  double mu2_1 = lambda1 * lambda1 * (rs - lambda * su);
  double mu2_2 = lambda * lambda * (rs + lambda1 * su);

  g_assert_cmpfloat(fabs(mu2_1 - mu2), <=, rtol * mu2 + atol);
  g_assert_cmpfloat(fabs(mu2_2 - mu2), <=, rtol * mu2 + atol);

  /*
   * Finally, we check that f'(lambda) = 0.
   *
   * The tolerance is given by the tolerance eps on lambda (defined in
   * the call to gsl_min_test_interval in the function
   * pw85_contact_function.
   *
   * The test therefore reads
   *
   *     |f'(lambda)| <= eps * |f"(lambda)|.
   *
   * We first compute f' and f".
   */

  double lambda2 = 1. - 2. * lambda;
  double f1 = lambda2 * rs - lambda * lambda1 * su;
  double f2 = -2. * rs - 2. * lambda2 * su + 2. * lambda * lambda1 * uv;
  g_assert_cmpfloat(fabs(f1), <=, PW85_LAMBDA_ATOL * fabs(f2));

  gsl_vector_free(r);
  gsl_vector_free(s);
  gsl_vector_free(u);
  gsl_vector_free(v);
  gsl_matrix_free(Q12);
  gsl_matrix_free(Q);
}

void test_pw85_contact_function_test() {
  for (size_t i = 0; i < test_pw85_context.num_distances; i++) {
    double const r = test_pw85_context.distances[i];
    for (size_t j = 0; j < test_pw85_context.num_directions; j++) {
      double const *n = test_pw85_context.directions + PW85_DIM * j;
      double const r12[] = {r * n[0], r * n[1], r * n[2]};
      for (size_t k1 = 0; k1 < test_pw85_context.num_spheroids; k1++) {
        double const *q1 = test_pw85_context.spheroids + PW85_SYM * k1;
        for (size_t k2 = 0; k2 < test_pw85_context.num_spheroids; k2++) {
          double const *q2 = test_pw85_context.spheroids + PW85_SYM * k2;
          test_pw85_contact_function_elementary_test(r12, q1, q2);
        }
      }
    }
  }
}

void test_pw85_f_neg_test() {
  double params[2 * PW85_SYM + PW85_DIM];
  double rtol = 1e-10;
  double *exp = test_pw85_context.f;
  for (size_t i1 = 0; i1 < test_pw85_context.num_spheroids; i1++) {
    memcpy(params + PW85_DIM, test_pw85_context.spheroids + PW85_SYM * i1,
           PW85_SYM * sizeof(double));
    for (size_t i2 = 0; i2 < test_pw85_context.num_spheroids; i2++) {
      memcpy(params + PW85_DIM + PW85_SYM,
             test_pw85_context.spheroids + PW85_SYM * i2,
             PW85_SYM * sizeof(double));
      for (size_t i = 0; i < test_pw85_context.num_directions; i++) {
        memcpy(params, test_pw85_context.directions + PW85_DIM * i,
               PW85_DIM * sizeof(double));
        for (size_t j = 0; j < test_pw85_context.num_lambdas; j++, exp++) {
          double const act = -pw85_f_neg(test_pw85_context.lambdas[j], params);
          g_assert_cmpfloat(fabs(act - *exp), <=, rtol * (*exp));
        }
      }
    }
  }
}

int main(int argc, char **argv) {
  g_test_init(&argc, &argv, NULL);

  hid_t const hid = H5Fopen(PW85_REF_DATA_PATH, H5F_ACC_RDONLY, H5P_DEFAULT);
  test_pw85_init_context(hid);
  H5Fclose(hid);

  pw85_test_add_cholesky_decomp_test();
  pw85_test_add_cholesky_solve_test();

  g_test_add_func("/pw85/spheroid", test_pw85_spheroid_test);
  g_test_add_func("/pw85/f_neg", test_pw85_f_neg_test);
  g_test_add_func("/pw85/contact_function", test_pw85_contact_function_test);

  int out = g_test_run();

  test_pw85_free_context();
  return out;
}
