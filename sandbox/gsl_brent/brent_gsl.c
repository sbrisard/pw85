#include <math.h>
#include <stdio.h>

#include <pw85.h>

int main() {
  const double golden_ratio = (1. + sqrt(5.)) / 2.;
  const double norm = sqrt(1 + golden_ratio * golden_ratio);
  const double u_abs = 1. / norm;
  const double v_abs = golden_ratio / norm;

  double r[] = {0., u_abs, v_abs};
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

  for (size_t i = 0; i < PW85_SYM; i++)
    printf("%g, ", q1[i]);
  printf("\n");
  for (size_t i = 0; i < PW85_SYM; i++)
    printf("%g, ", q2[i]);
  printf("\n");

  return 0;
}
