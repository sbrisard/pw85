import itertools

import h5py
import mpmath
import numpy as np


def spheroid_mp(a, c, n):
    a2 = a * a
    c2_minus_a2 = c * c - a2
    return (c * c - a2) * n * n.T + a2 * mpmath.mp.eye(3)


def f_mp(lambda_, r12, q1, q2):
    one_minus_lambda = 1 - lambda_
    q = one_minus_lambda * q1 + lambda_ * q2
    s12, _ = mpmath.mp.qr_solve(q, r12)
    return lambda_ * one_minus_lambda * (r12.T * s12)[0, 0]


def gen_directions():
    golden_ratio = (mpmath.mp.mpf(1) + mpmath.mp.sqrt(mpmath.mp.mpf(5))) / (
        mpmath.mp.mpf(2)
    )
    u_abs = 1.0 / mpmath.mp.sqrt(1 + golden_ratio ** 2)
    v_abs = golden_ratio * u_abs
    directions = []
    for u in (-u_abs, u_abs):
        for v in (-v_abs, v_abs):
            directions += [
                mpmath.mp.matrix((0.0, u, v)),
                mpmath.mp.matrix((v, 0.0, u)),
                mpmath.mp.matrix((u, v, 0.0)),
            ]
    return directions


def gen_spheroids(radii, directions):
    return [
        spheroid_mp(ai, ci, ni)
        for ai, ci, ni in itertools.product(radii, radii, directions)
    ]


if __name__ == "__main__":
    mpmath.mp.dps = 50

    radii = [mpmath.mp.mpf(r) for r in ["0.01999", "1.999", "9.999"]]
    num_radii = len(radii)

    directions = gen_directions()
    num_directions = len(directions)

    spheroids = gen_spheroids(radii, directions)
    num_spheroids = len(spheroids)

    num_lambdas = 11
    lambdas = mpmath.mp.linspace(0, 1, num_lambdas)

    num_blocks = num_radii * num_spheroids ** 2

    with h5py.File("pw85_reference_data.h5", "w") as f:
        f["radii"] = [float(r) for r in radii]
        f["directions"] = [(float(x), float(y), float(z)) for x, y, z in directions]
        f["lambdas"] = [float(lambda_) for lambda_ in lambdas]
        dset = f.create_dataset("spheroids", shape=(num_spheroids, 6), dtype="d")
        for i, q in enumerate(spheroids):
            dset[i, 0] = float(q[0, 0])
            dset[i, 1] = float(q[0, 1])
            dset[i, 2] = float(q[0, 2])
            dset[i, 3] = float(q[1, 1])
            dset[i, 4] = float(q[1, 2])
            dset[i, 5] = float(q[2, 2])

        shape = (num_spheroids, num_spheroids, num_directions, num_radii, num_lambdas)

        dset = f.create_dataset("f", shape=shape, dtype="d")
        block = np.empty((num_directions, num_radii, num_lambdas), dtype=np.float64)
        for i1, q1 in enumerate(spheroids):
            print(f"{i1+1}/{num_spheroids}")
            for i2, q2 in enumerate(spheroids):
                for i, r12_dir in enumerate(directions):
                    for j, r12_norm in enumerate(radii):
                        r12 = r12_norm * r12_dir
                        for k, lambda_ in enumerate(lambdas):
                            block[i, j, k] = float(f_mp(lambda_, r12, q1, q2))
                dset[i1, i2, :, :, :] = block
                f.flush()
