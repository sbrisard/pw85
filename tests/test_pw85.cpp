#include <algorithm>
#include <cmath>
#include <cstring>

#include "catch2/catch.hpp"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <hdf5_hl.h>

#include "pw85/pw85.hpp"

using Vec = std::array<double, 3>;
using Sym = std::array<double, 6>;

template <typename E, typename A>
void assert_cmp_doubles(E first_expected, E last_expected, A first_actual,
                        double rtol, double atol) {
  auto act = first_actual;
  for (auto exp = first_expected; exp != last_expected; ++exp, ++act) {
    REQUIRE(*act == Catch::Detail::Approx(*exp).epsilon(rtol).margin(atol));
  }
}

void test_pw85_read_dataset_double(hid_t const hid, char const *dset_name,
                                   size_t *size, double **buffer) {
  int ndims;
  H5LTget_dataset_ndims(hid, dset_name, &ndims);
  hsize_t *dim = new hsize_t[ndims];
  H5LTget_dataset_info(hid, dset_name, dim, NULL, NULL);
  *size = 1;
  for (size_t i = 0; i < ndims; i++) {
    *size *= dim[i];
  }
  *buffer = new double[*size];
  H5LTread_dataset_double(hid, dset_name, *buffer);
  delete[] dim;
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
  if (hid > 0) {
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

    test_pw85_read_dataset_double(hid, "/lambdas",
                                  &test_pw85_context.num_lambdas,
                                  &test_pw85_context.lambdas);

    test_pw85_read_dataset_double(hid, "/F", &test_pw85_context.num_f,
                                  &test_pw85_context.f);
  } else {
    test_pw85_context.num_radii = 0;
    test_pw85_context.radii = NULL;
    test_pw85_context.num_directions = 0;
    test_pw85_context.directions = NULL;
    test_pw85_context.num_spheroids = 0;
    test_pw85_context.spheroids = NULL;
    test_pw85_context.num_lambdas = 0;
    test_pw85_context.lambdas = NULL;
    test_pw85_context.num_f = 0;
    test_pw85_context.f = NULL;
  }

  test_pw85_context.num_distances = 3;
  test_pw85_context.distances = static_cast<double *>(
      malloc(sizeof(double) * test_pw85_context.num_distances));
  test_pw85_context.distances[0] = 0.15;
  test_pw85_context.distances[1] = 1.1;
  test_pw85_context.distances[2] = 11.;
}

void test_pw85_free_context() {
  free(test_pw85_context.directions);
  free(test_pw85_context.radii);
  free(test_pw85_context.spheroids);
  free(test_pw85_context.lambdas);
  free(test_pw85_context.f);
  free(test_pw85_context.distances);
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

std::vector<std::array<double, PW85_DIM>> test_pw85_gen_directions() {
  double u = sqrt(2. / (5. + sqrt(5.)));
  double v = sqrt((3 + sqrt(5.)) / (5. + sqrt(5.)));
  return {{0., -u, -v}, {-v, 0., u}, {u, -v, 0.}};
}

std::vector<Vec> test_pw85_gen_radius_vectors(std::vector<double> distances,
                                              std::vector<Vec> directions) {
  std::vector<Vec> radius_vectors{};
  for (const auto r : distances) {
    for (const auto n : directions) {
      radius_vectors.push_back({r * n[0], r * n[1], r * n[2]});
    }
  }
  return radius_vectors;
}

std::vector<Sym> test_pw85_gen_spheroids(std::vector<double> radii,
                                         std::vector<Vec> directions) {
  std::vector<Sym> spheroids{};
  for (const auto a : radii) {
    for (const auto c : radii) {
      for (const auto n : directions) {
        Sym q{};
        pw85::spheroid(a, c, n.data(), q.data());
        spheroids.push_back(q);
      }
    }
  }
  return spheroids;
}

void test_cholesky_decomp(Sym const a, Sym const exp, double const rtol) {
  Sym act;
  pw85::_cholesky_decomp(a.data(), act.data());
  assert_cmp_doubles(exp.cbegin(), exp.cend(), act.cbegin(), rtol, 0.);
}

void test_cholesky_solve(Sym const l, Vec const b, Vec const exp, double rtol) {
  Vec act;
  pw85::_cholesky_solve(l.data(), b.data(), act.data());
  assert_cmp_doubles(exp.cbegin(), exp.cend(), act.cbegin(), rtol, 0.0);
}

void test_spheroid(double a, double c, Vec n) {
  /*
   * Relative and absolute tolerance on the coefficients of the matrix
   * q to be computed and tested.
   */
  double rtol = 1e-15;
  double atol = 1e-15;

  double exp, act, tol;

  double a2 = a * a;
  double c2 = c * c;

  auto abs_f = [](double x) { return fabs(x); };
  auto tol_f = [rtol, atol](double x) { return rtol * fabs(x) + atol; };

  Vec abs_n;
  std::transform(n.cbegin(), n.cend(), abs_n.begin(), abs_f);

  Sym q, abs_q, delta_q;
  pw85::spheroid(a, c, n.data(), q.data());
  std::transform(q.cbegin(), q.cend(), abs_q.begin(), abs_f);
  std::transform(q.cbegin(), q.cend(), delta_q.begin(), tol_f);

  /* Check that n is an eigenvector. */
  Vec qn{q[0] * n[0] + q[1] * n[1] + q[2] * n[2],
         q[1] * n[0] + q[3] * n[1] + q[4] * n[2],
         q[2] * n[0] + q[4] * n[1] + q[5] * n[2]};
  Vec delta_qn{
      delta_q[0] * abs_n[0] + delta_q[1] * abs_n[1] + delta_q[2] * abs_n[2],
      delta_q[1] * abs_n[0] + delta_q[3] * abs_n[1] + delta_q[4] * abs_n[2],
      delta_q[2] * abs_n[0] + delta_q[4] * abs_n[1] + delta_q[5] * abs_n[2]};
  for (size_t i = 0; i < PW85_DIM; i++) {
    REQUIRE(qn[i] == Catch::Detail::Approx(c2 * n[i]).margin(delta_qn[i]));
  }

  /*
   * Check that the eigenvalues are [a, a, c]. It is sufficient to
   * compute the characteristic polynomial of the matrix q.
   */

  /* Check tr(q). */
  exp = 2. * a2 + c2;
  act = q[0] + q[3] + q[5];
  tol = delta_q[0] + delta_q[3] + delta_q[5];
  REQUIRE(act == Catch::Detail::Approx(exp).margin(tol));

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
  REQUIRE(act == Catch::Detail::Approx(exp).margin(tol));

  /* Check [tr(q)^2 - tr(q^2)]/2 */
  exp = a2 * (a2 + 2. * c2);
  act = q[0] * q[3] + q[3] * q[5] + q[5] * q[0] - q[1] * q[1] - q[2] * q[2] -
        q[4] * q[4];
  tol = delta_q[0] * abs_q[3] + abs_q[0] * delta_q[3] + delta_q[3] * abs_q[5] +
        abs_q[3] * delta_q[5] + delta_q[5] * abs_q[0] + abs_q[5] * delta_q[0] +
        2. * delta_q[1] * abs_q[1] + 2. * delta_q[2] * abs_q[2] +
        2. * delta_q[4] * abs_q[4];
  REQUIRE(act == Catch::Detail::Approx(exp).margin(tol));
}

void test_pw85_contact_function_test(double const *r12, double const *q1,
                                     double const *q2) {
  double atol = 1e-15;
  double rtol = 1e-10;

  double out[2];
  pw85::contact_function(r12, q1, q2, out);
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

  REQUIRE(mu2 == Catch::Detail::Approx(mu2_1).scale(rtol).margin(atol));
  REQUIRE(mu2 == Catch::Detail::Approx(mu2_2).scale(rtol).margin(atol));

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
  REQUIRE(fabs(f1) <= PW85_LAMBDA_ATOL * fabs(f2));

  gsl_vector_free(r);
  gsl_vector_free(s);
  gsl_vector_free(u);
  gsl_vector_free(v);
  gsl_matrix_free(Q12);
  gsl_matrix_free(Q);
}

void test_pw85_f_neg_test(double lambda, double const r12[PW85_DIM],
                          double const q1[PW85_SYM], double const q2[PW85_SYM],
                          double exp, double rtol, double atol) {
  double params[2 * PW85_SYM + PW85_DIM];
  memcpy(params, r12, PW85_DIM * sizeof(double));
  memcpy(params + PW85_DIM, q1, PW85_SYM * sizeof(double));
  memcpy(params + PW85_DIM + PW85_SYM, q2, PW85_SYM * sizeof(double));
  double act = -pw85::f_neg(lambda, params);
  REQUIRE(act == Catch::Detail::Approx(exp).scale(rtol).margin(atol));
}

TEST_CASE("pw85") {
  hid_t const hid = H5Fopen(PW85_REF_DATA_PATH, H5F_ACC_RDONLY, H5P_DEFAULT);
  test_pw85_init_context(hid);
  H5Fclose(hid);

  SECTION("cholesky_decomp") {
    test_cholesky_decomp({4, 2, 6, 17, 23, 70}, {2, 1, 3, 4, 5, 6}, 1e-15);

    test_cholesky_decomp({4, -2, 6, 17, -23, 70}, {2, -1, 3, 4, -5, 6}, 1e-15);

    test_cholesky_decomp(
        {1e10, -2, -3, 16 + 1. / 25e8, -0.02 + 3. / 5e9, 29. / 8e3 - 9e-10},
        {1e5, -2e-5, -3e-5, 4, -5e-3, 6e-2}, 1e-6);
  }

  SECTION("cholesky_solve") {
    test_cholesky_solve({1, 2, 3, 4, 5, 6}, {11.5, 82.6, 314.2},
                        {1.2, -3.4, 5.7}, 1e-15);

    test_cholesky_solve({1, -2, -3, 4, -5, 6}, {-9.1, -150.2, 443},
                        {1.2, -3.4, 5.7}, 4e-15);
  }

  SECTION("spheroid") {
    for (size_t i = 0; i < test_pw85_context.num_radii; i++) {
      double const a = test_pw85_context.radii[i];
      for (size_t j = 0; j < test_pw85_context.num_radii; j++) {
        double const c = test_pw85_context.radii[j];
        for (size_t k = 0; k < test_pw85_context.num_directions; k++) {
          // TODO: remove raw pointer
          double *const n = test_pw85_context.directions + PW85_DIM * k;
          test_spheroid(a, c, {n[0], n[1], n[2]});
        }
      }
    }
  }

  SECTION("f_neg") {
    double rtol = 1e-10;
    double *q_begin = test_pw85_context.spheroids;
    double *q_end = q_begin + test_pw85_context.num_spheroids * PW85_SYM;
    double *n_begin = test_pw85_context.directions;
    double *n_end = n_begin + test_pw85_context.num_directions * PW85_DIM;
    double *lambda_begin = test_pw85_context.lambdas;
    double *lambda_end = lambda_begin + test_pw85_context.num_lambdas;
    double *exp = test_pw85_context.f;

    for (double *q1 = q_begin; q1 < q_end; q1 += PW85_SYM) {
      for (double *q2 = q_begin; q2 < q_end; q2 += PW85_SYM) {
        for (double *n = n_begin; n < n_end; n += PW85_DIM) {
          for (double *lambda = lambda_begin; lambda < lambda_end;
               lambda++, exp++) {
            test_pw85_f_neg_test(*lambda, n, q1, q2, *exp, rtol, 0.0);
          }
        }
      }
    }
  }

  SECTION("contact_function") {
    // TODO: these assertions should be removed.
    assert(test_pw85_context.distances != nullptr);
    assert(test_pw85_context.directions != nullptr);
    assert(test_pw85_context.spheroids != nullptr);

    double *q_begin = test_pw85_context.spheroids;
    double *q_end = q_begin + test_pw85_context.num_spheroids * PW85_SYM;
    double *n_begin = test_pw85_context.directions;
    double *n_end = n_begin + test_pw85_context.num_directions * PW85_DIM;

    for (size_t i = 0; i < test_pw85_context.num_distances; i++) {
      double const r = test_pw85_context.distances[i];
      for (double *n = n_begin; n < n_end; n += PW85_DIM) {
        double const r12[] = {r * n[0], r * n[1], r * n[2]};
        for (double *q1 = q_begin; q1 < q_end; q1 += PW85_SYM) {
          for (double *q2 = q_begin; q2 < q_end; q2 += PW85_SYM) {
            test_pw85_contact_function_test(r12, q1, q2);
          }
        }
      }
    }
  }

  test_pw85_free_context();
}