#include <stdio.h>

#define DIM 3
#define SYM 6

void spheroid(double a, double c, double n[DIM], double q[SYM]) {
    for (int i = 0; i < SYM; i++) {
	q[i] = i;
    }
}

int main(void) {
    double q[SYM];
    spheroid(1.0, 0.1, NULL, q);
    for (int i = 0; i < SYM; i++) {
	printf("q[%d] = %g\n", i, q[i]);
    }
}
