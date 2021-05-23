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

int main() {
  hid_t const hid = H5Fopen(PW85_REF_DATA_PATH, H5F_ACC_RDONLY, H5P_DEFAULT);
  test_pw85_init_context(hid);
  H5Fclose(hid);

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
