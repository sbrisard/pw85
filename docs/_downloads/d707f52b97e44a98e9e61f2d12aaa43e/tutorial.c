#include <pw85.h>
#include <stdio.h>
#include <stdlib.h>

#define DIM 3
#define SYM 6

int main() {
  double x1[] = {-0.5, 0.4, -0.7};
  double n1[] = {0., 0., 1.};
  double a1 = 10.;
  double c1 = 0.1;

  double x2[] = {0.2, -0.3, 0.4};
  double n2[] = {1., 0., 0.};
  double a2 = 0.5;
  double c2 = 5.;

  double r12[DIM];
  for (int i = 0; i < DIM; i++) r12[i] = x2[i] - x1[i];

  double q1[SYM], q2[SYM];
  pw85_spheroid(a1, c1, n1, q1);
  pw85_spheroid(a2, c2, n2, q2);

  double mu2 = pw85_contact_function(r12, q1, q2, NULL);
  printf("Simplified output\n");
  printf("-----------------\n");
  printf("mu^2 = %g\n", mu2);
  printf("\n");

  double out[2];
  pw85_contact_function(r12, q1, q2, out);
  printf("Full output\n");
  printf("-----------\n");
  printf("mu^2 = %g\n", out[0]);
  printf("lambda = %g\n", out[1]);
  printf("\n");
}
