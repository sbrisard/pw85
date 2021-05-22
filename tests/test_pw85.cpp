#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstring>
#include <iostream>

#include <Eigen/Dense>

#include <hdf5_hl.h>

#include "pw85/pw85.hpp"

using Vec = std::array<double, 3>;
using Sym = std::array<double, 6>;

void assert_cmp_double(double exp, double act, double rtol, double atol) {
  double err = fabs(act - exp);
  double tol = rtol * fabs(exp) + atol;
  if (err > tol) {
    std::cerr << "exp = " << exp << std::endl;
    std::cerr << "act = " << act << std::endl;
  }
  assert(err <= tol);
}

template <typename E, typename A>
void assert_cmp_doubles(E first_expected, E last_expected, A first_actual,
                        double rtol, double atol) {
  auto act = first_actual;
  for (auto exp = first_expected; exp != last_expected; ++exp, ++act) {
    assert_cmp_double(*exp, *act, rtol, atol);
  }
}

std::vector<double> test_pw85_read_dataset_double(hid_t const hid,
                                                  char const *dset_name) {
  // TODO Use C++ HDF interface
  int ndims;
  H5LTget_dataset_ndims(hid, dset_name, &ndims);
  hsize_t *dim = new hsize_t[ndims];
  H5LTget_dataset_info(hid, dset_name, dim, NULL, NULL);
  hsize_t size = 1;
  for (size_t i = 0; i < ndims; i++) {
    size *= dim[i];
  }
  std::vector<double> dataset(size);
  H5LTread_dataset_double(hid, dset_name, dataset.data());
  delete[] dim;
  return dataset;
}

/* This is a global variable that holds the parameters used for most
 * test-cases. */

struct {
  std::vector<double> radii;
  std::vector<Vec> directions;
  std::vector<Sym> spheroids;
  std::vector<double> lambdas;
  std::vector<double> f;
  std::vector<double> distances;
} test_pw85_context;

void test_pw85_init_context(hid_t const hid) {
  if (hid > 0) {
    // TODO: this is ugly
    auto directions = test_pw85_read_dataset_double(hid, "/directions");
    for (auto n = directions.cbegin(); n != directions.cend(); n += pw85::dim) {
      test_pw85_context.directions.push_back({n[0], n[1], n[2]});
    }

    test_pw85_context.radii = test_pw85_read_dataset_double(hid, "/radii");

    // TODO: this is ugly
    auto spheroids = test_pw85_read_dataset_double(hid, "/spheroids");
    for (auto q = spheroids.cbegin(); q != spheroids.cend(); q += pw85::sym) {
      test_pw85_context.spheroids.push_back(
          {q[0], q[1], q[2], q[3], q[4], q[5]});
    }

    test_pw85_context.lambdas = test_pw85_read_dataset_double(hid, "/lambdas");

    test_pw85_context.f = test_pw85_read_dataset_double(hid, "/F");
  } else {
    exit(-1);
  }

  test_pw85_context.distances = {0.15, 1.1, 11.};
}

/* #define TEST_PW85_NUM_DIRECTIONS 12 */

/* double *test_pw85_gen_directions() { */
/*   double *directions = g_new(double, TEST_PW85_NUM_DIRECTIONS * pw85::dim); */
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
/*       n += pw85::dim; */
/*       n[0] = v; */
/*       n[1] = 0.; */
/*       n[2] = u; */
/*       n += pw85::dim; */
/*       n[0] = u; */
/*       n[1] = v; */
/*       n[2] = 0.; */
/*       n += pw85::dim; */
/*     } */
/*   } */
/*   return directions; */
/* } */

std::vector<std::array<double, pw85::dim>> test_pw85_gen_directions() {
  double u = sqrt(2. / (5. + sqrt(5.)));
  double v = sqrt((3 + sqrt(5.)) / (5. + sqrt(5.)));
  return {{0., -u, -v}, {-v, 0., u}, {u, -v, 0.}};
}

std::vector<Vec> test_pw85_gen_radius_vectors(
    const std::vector<double> &distances, const std::vector<Vec> &directions) {
  std::vector<Vec> radius_vectors{};
  for (const auto r : distances) {
    for (const auto n : directions) {
      radius_vectors.push_back({r * n[0], r * n[1], r * n[2]});
    }
  }
  return radius_vectors;
}

std::vector<Sym> test_pw85_gen_spheroids(const std::vector<double> &radii,
                                         const std::vector<Vec> &directions) {
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
  for (size_t i = 0; i < pw85::dim; i++) {
    assert_cmp_double(c2 * n[i], qn[i], 0., delta_qn[i]);
  }

  /*
   * Check that the eigenvalues are [a, a, c]. It is sufficient to
   * compute the characteristic polynomial of the matrix q.
   */

  /* Check tr(q). */
  exp = 2. * a2 + c2;
  act = q[0] + q[3] + q[5];
  tol = delta_q[0] + delta_q[3] + delta_q[5];
  assert_cmp_double(exp, act, 0., tol);

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
  assert_cmp_double(exp, act, 0., tol);

  /* Check [tr(q)^2 - tr(q^2)]/2 */
  exp = a2 * (a2 + 2. * c2);
  act = q[0] * q[3] + q[3] * q[5] + q[5] * q[0] - q[1] * q[1] - q[2] * q[2] -
        q[4] * q[4];
  tol = delta_q[0] * abs_q[3] + abs_q[0] * delta_q[3] + delta_q[3] * abs_q[5] +
        abs_q[3] * delta_q[5] + delta_q[5] * abs_q[0] + abs_q[5] * delta_q[0] +
        2. * delta_q[1] * abs_q[1] + 2. * delta_q[2] * abs_q[2] +
        2. * delta_q[4] * abs_q[4];
  assert_cmp_double(exp, act, 0., tol);
}

void test_pw85_contact_function_test(Vec const r12, Sym const q1,
                                     Sym const q2) {
  double atol = 1e-15;
  double rtol = 1e-10;

  double out[2];
  pw85::contact_function(r12.data(), q1.data(), q2.data(), out);
  double mu2 = out[0];
  double lambda = out[1];
  double lambda1 = 1. - lambda;

  Eigen::Vector3d r{r12[0], r12[1], r12[2]};
  Eigen::Matrix3d Q1, Q2;
  // clang-format off
  Q1 << q1[0], q1[1], q1[2],
        q1[1], q1[3], q1[4],
        q1[2], q1[4], q1[5];
  Q2 << q2[0], q2[1], q2[2],
        q2[1], q2[3], q2[4],
        q2[2], q2[4], q2[5];
  // clang-format on
  auto Q12 = Q2 - Q1;
  auto Q = (1 - lambda) * Q1 + lambda * Q2;
  auto cholesky_Q = Q.llt();
  auto s = cholesky_Q.solve(r);

  auto u = Q12 * s;
  auto v = cholesky_Q.solve(u);

  double rs = r.dot(s);
  double su = s.dot(u);
  double uv = u.dot(v);

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

  assert_cmp_double(mu2_1, mu2, rtol, atol);
  assert_cmp_double(mu2_2, mu2, rtol, atol);

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
  assert_cmp_double(0., f1, 0., pw85::lambda_atol * fabs(f2));
}

void test_pw85_f_neg_test(double lambda, const Vec r12, const Sym q1,
                          const Sym q2, double exp, double rtol, double atol) {
  double params[2 * pw85::sym + pw85::dim];
  // TODO: this is ugly
  memcpy(params, r12.data(), pw85::dim * sizeof(double));
  memcpy(params + pw85::dim, q1.data(), pw85::sym * sizeof(double));
  memcpy(params + pw85::dim + pw85::sym, q2.data(), pw85::sym * sizeof(double));
  double act = -pw85::f_neg(lambda, params);
  assert_cmp_double(exp, act, rtol, atol);
}

int main() {
  hid_t const hid = H5Fopen(PW85_REF_DATA_PATH, H5F_ACC_RDONLY, H5P_DEFAULT);
  test_pw85_init_context(hid);
  H5Fclose(hid);


  for (const auto a : test_pw85_context.radii) {
    for (const auto c : test_pw85_context.radii) {
      for (const auto n : test_pw85_context.directions) {
        test_spheroid(a, c, n);
      }
    }
  }

  auto exp = test_pw85_context.f.cbegin();
  for (const auto q1 : test_pw85_context.spheroids) {
    for (const auto q2 : test_pw85_context.spheroids) {
      for (const auto n : test_pw85_context.directions) {
        for (const auto lambda : test_pw85_context.lambdas) {
          test_pw85_f_neg_test(lambda, n, q1, q2, *exp, 1e-10, 0.0);
          ++exp;
        }
      }
    }
  }

  for (auto const r : test_pw85_context.distances) {
    for (auto const n : test_pw85_context.directions) {
      Vec r12{r * n[0], r * n[1], r * n[2]};
      for (auto const q1 : test_pw85_context.spheroids) {
        for (auto const q2 : test_pw85_context.spheroids) {
          test_pw85_contact_function_test(r12, q1, q2);
        }
      }
    }
  }
}
