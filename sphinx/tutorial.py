import numpy as np
import pypw85

if __name__ == '__main__':
    a1, c1, n1 = 10, 0.1, np.array([0., 0., 1.])
    a2, c2, n2 = 0.5, 5., np.array([1., 0., 0.])
    r12 = np.array([0.2, -0.3, 0.4])

    q1 = pypw85.spheroid(a1, c1, n1)
    q2 = pypw85.spheroid(a2, c2, n2)

    print('q1 = {}'.format(q1))
    print('q2 = {}'.format(q2))

    mu2 = pypw85.contact_function(r12, q1, q2)
    print('μ² = {}'.format(mu2))

    Q1 = np.zeros((3, 3), dtype=np.float64)
    i, j = np.triu_indices_from(Q1)

    Q1[i, j] = q1
    Q1[j, i] = q1

    Q2 = np.zeros_like(Q1)
    Q2[i, j] = q2
    Q2[j, i] = q2

    Q1_inv = np.linalg.inv(Q1)
    Q2_inv = np.linalg.inv(Q2)

    f1 = lambda x: Q1_inv.dot(x).dot(x)

    print(f1((0., 0., c1)))

    f2 = lambda x: Q2_inv.dot(x).dot(x)

    print(f2((c2, 0., 0.)))
    print(f2((-c2, 0., 0.)))
    print(f2((0., a2, 0.)))
    print(f2((0., -a2, 0.)))
    print(f2((0., 0., a2)))
    print(f2((0., 0., -a2)))
