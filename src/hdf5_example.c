#include <glib.h>
#include <hdf5.h>
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

int main() {
  hid_t const hid = H5Fopen("../pw85_ref_data.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

  size_t num_directions;
  double *directions;
  test_pw85_read_dataset_double(hid, "/directions", &num_directions,
                                &directions);
  num_directions /= PW85_DIM;
  printf("Total number of directions = %d\n", num_directions);
  for (size_t i = 0; i < num_directions; i++) {
    double const *n = directions + PW85_DIM * i;
    printf("n[%d] = [%g, %g, %g]\n", i, n[0], n[1], n[2]);
  }

  size_t num_lambdas;
  double *lambdas;
  test_pw85_read_dataset_double(hid, "/lambdas", &num_lambdas, &lambdas);
  printf("Total number of lambdas = %d\n", num_lambdas);
  for (size_t i = 0; i < num_lambdas; i++) {
    printf("lambdas[%d] = %g\n", i, lambdas[i]);
  }

  size_t num_radii;
  double *radii;
  test_pw85_read_dataset_double(hid, "/radii", &num_radii, &radii);
  printf("Total number of radii = %d\n", num_radii);
  for (size_t i = 0; i < num_radii; i++) {
    printf("radii[%d] = %g\n", i, radii[i]);
  }

  size_t num_spheroids;
  double *spheroids;
  test_pw85_read_dataset_double(hid, "/spheroids", &num_spheroids, &spheroids);
  num_spheroids /= PW85_SYM;
  printf("Total number of spheroids = %d\n", num_spheroids);
  for (size_t i = 0; i < num_spheroids; i++) {
    double *q = spheroids + PW85_SYM * i;
    printf("spheroid[%d] = [%g, %g, %g, %g, %g, %g]\n", i, q[0], q[1], q[2],
           q[3], q[4], q[5]);
  }

  size_t num_expecteds;
  double *expecteds;
  test_pw85_read_dataset_double(hid, "/F", &num_expecteds, &expecteds);
  printf("Total number of expecteds = %d\n", num_expecteds);

  double *exp = expecteds;
  double params[2 * PW85_SYM + PW85_DIM];
  for (size_t i1 = 0; i1 < num_spheroids; i1++) {
    memcpy(params + PW85_DIM, spheroids + PW85_SYM * i1,
           PW85_SYM * sizeof(double));
    for (size_t i2 = 0; i2 < num_spheroids; i2++) {
      memcpy(params + PW85_DIM + PW85_SYM, spheroids + PW85_SYM * i2,
             PW85_SYM * sizeof(double));
      for (size_t i = 0; i < num_directions; i++) {
        memcpy(params, directions + PW85_DIM * i, PW85_DIM * sizeof(double));
        for (size_t j = 0; j < num_lambdas; j++, exp++) {
          double const act = pw85_f_neg(lambdas[j], params);
          printf("%g, %g\n", act, *exp);
        }
      }
    }
  }

  g_free(spheroids);
  g_free(radii);
  g_free(lambdas);
  g_free(directions);
  H5Fclose(hid);

  return 0;
}
