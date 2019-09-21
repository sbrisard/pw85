#include <glib.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <pw85.h>
#include <pw85_legacy.h>

void read_dataset_double(hid_t const hid, char const *dset_name, size_t *size,
                         double **buffer) {
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

void update_histogram(double act, double exp, size_t num_bins, size_t *hist) {
  double const err = fabs((act - exp) / exp);
  int prec;
  if (err == 0.0) {
    prec = num_bins - 1;
  } else {
    prec = (int)(floor(-log10(err)));
    if (prec <= 0) {
      prec = 0;
    }
    if (prec >= num_bins) {
      prec = num_bins - 1;
    }
  }
  ++hist[prec];
}

int main() {
  hid_t const hid = H5Fopen(PW85_REF_DATA_PATH, H5F_ACC_RDONLY, H5P_DEFAULT);

  size_t num_directions;
  double *directions;
  read_dataset_double(hid, "/directions", &num_directions, &directions);
  num_directions /= PW85_DIM;

  size_t num_lambdas;
  double *lambdas;
  read_dataset_double(hid, "/lambdas", &num_lambdas, &lambdas);

  size_t num_radii;
  double *radii;
  read_dataset_double(hid, "/radii", &num_radii, &radii);

  size_t num_spheroids;
  double *spheroids;
  read_dataset_double(hid, "/spheroids", &num_spheroids, &spheroids);
  num_spheroids /= PW85_SYM;

  size_t num_expecteds;
  double *expecteds;
  read_dataset_double(hid, "/F", &num_expecteds, &expecteds);

  double *exp = expecteds;
  double params[2 * PW85_SYM + PW85_DIM];

  size_t num_bins = 16;
  size_t hist1[num_bins];
  size_t hist2[num_bins];
  for (size_t i = 0; i < num_bins; i++) {
    hist1[i] = 0;
    hist2[i] = 0;
  }
  for (size_t i1 = 0; i1 < num_spheroids; i1++) {
    memcpy(params + PW85_DIM, spheroids + PW85_SYM * i1,
           PW85_SYM * sizeof(double));
    for (size_t i2 = 0; i2 < num_spheroids; i2++) {
      memcpy(params + PW85_DIM + PW85_SYM, spheroids + PW85_SYM * i2,
             PW85_SYM * sizeof(double));
      for (size_t i = 0; i < num_directions; i++) {
        memcpy(params, directions + PW85_DIM * i, PW85_DIM * sizeof(double));
        for (size_t j = 0; j < num_lambdas; j++, exp++) {
          double const act1 = -pw85_f_neg(lambdas[j], params);
          update_histogram(act1, *exp, num_bins, hist1);
	  double out[2];
	  pw85_legacy_f2(lambdas[j], params, params + PW85_DIM,
			 params + PW85_DIM + PW85_SYM, out);
	  double const act2 = out[0];
          update_histogram(act2, *exp, num_bins, hist2);
        }
      }
    }
  }

  FILE *f = fopen(HISTOGRAM_PATH, "w");
  for (size_t i = 0; i < num_bins; i++) {
    fprintf(f, "%d,%g,%g\n", (int)i,
	    100. * ((double)hist1[i]) / ((double)num_expecteds),
	    100. * ((double)hist2[i]) / ((double)num_expecteds));
  }
  fclose(f);

  g_free(spheroids);
  g_free(radii);
  g_free(lambdas);
  g_free(directions);
  H5Fclose(hid);

  return 0;
}
