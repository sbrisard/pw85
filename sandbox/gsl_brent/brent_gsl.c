#include <float.h>
#include <math.h>
#include <stdio.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>

#include <pw85.h>

double f(double x, void *params) {
  const double golden_ratio = (1. + sqrt(5.)) / 2.;
  const double norm = sqrt(1 + golden_ratio * golden_ratio);
  const double u_abs = 1. / norm;
  const double v_abs = golden_ratio / norm;

  double r12[] = {0., u_abs, v_abs};
  double n1[] = {0., u_abs, v_abs};
  double n2[] = {0., u_abs, v_abs};
  const double a1 = 1.999;
  const double c1 = 0.199;
  const double a2 = 0.1999;
  const double c2 = 2.0199;

  double q1[PW85_SYM];
  double q2[PW85_SYM];
  pw85_spheroid(a1, c1, n1, q1);
  pw85_spheroid(a2, c2, n2, q2);

  return -pw85_f(x, r12, q1, q2, NULL);
}

int main() {
  int status, iter = 0, max_iter = 100;
  double a = 0.;
  double b = 1.;
  double m = 0.5;
  double fm = f(m, NULL);

  const gsl_function F = {.function = &f, .params = NULL};
  gsl_min_fminimizer *s = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
  gsl_min_fminimizer_set(s, &F, m, a, b);

  printf("using %s method\n", gsl_min_fminimizer_name(s));
  printf("%5s %9s [%9s, %9s] %9s %10s\n", "iter", "value", "lower", "upper",
         "min", "err(est)");
  printf("%5d %.7f [%.7f, %.7f] %.7f %+.7f\n", iter, fm, a, b, m, b - a);

  do {
    iter++;
    status = gsl_min_fminimizer_iterate(s);

    m = gsl_min_fminimizer_x_minimum(s);
    a = gsl_min_fminimizer_x_lower(s);
    b = gsl_min_fminimizer_x_upper(s);
    fm = gsl_min_fminimizer_f_minimum(s);

    status = gsl_min_test_interval(a, b, 1e-8, 0.0);

    if (status == GSL_SUCCESS)
      printf("Converged:\n");

    printf("%5d %.7f [%.7f, %.7f] %.7f %.7f\n", iter, fm, a, b, m, b - a);
  } while (status == GSL_CONTINUE && iter < max_iter);

  gsl_min_fminimizer_free(s);
}
