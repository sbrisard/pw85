#include <stdio.h>

#define PW85_DIM 3
#define PW85_SYM 6
#define PW85_XX 0
#define PW85_YY 3
#define PW85_ZZ 5
#define PW85_XY 1
#define PW85_XZ 2
#define PW85_YZ 4

__declspec(dllexport) void pw85_spheroid(double a, double c, double n[PW85_DIM], double q[PW85_SYM])
{
    const double inv_a2 = 1.0 / (a * a);
    const double inv_c2_minus_inv_a2 = 1.0 / (c * c) - inv_a2;
    const double nx = n[0];
    const double ny = n[1];
    const double nz = n[2];
    printf("%g\t%g\t%g\n", nx, ny, nz);
    const double inv_n2 = 1.0 / (nx * nx + ny * ny + nz * nz);
    q[PW85_XX] = inv_a2 + inv_n2 * nx * nx * inv_c2_minus_inv_a2;
    q[PW85_YY] = inv_a2 + inv_n2 * ny * ny * inv_c2_minus_inv_a2;
    q[PW85_ZZ] = inv_a2 + inv_n2 * nz * nz * inv_c2_minus_inv_a2;
    q[PW85_YZ] = inv_n2 * ny * nz * inv_c2_minus_inv_a2;
    q[PW85_XZ] = inv_n2 * nx * nz * inv_c2_minus_inv_a2;
    q[PW85_XY] = inv_n2 * nx * ny * inv_c2_minus_inv_a2;
    for (int i = 0; i < PW85_SYM; i++) {
        printf("%g\n", q[i]);
    }
}