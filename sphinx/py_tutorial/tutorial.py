import numpy as np
import pypw85

if __name__ == "__main__":
    x1 = np.array([-0.5, 0.4, -0.7])
    n1 = np.array([0.0, 0.0, 1.0])
    a1, c1 = 10, 0.1
    x2 = np.array([0.2, -0.3, 0.4])
    n2 = np.array([1.0, 0.0, 0.0])
    a2, c2 = 0.5, 5.0
    r12 = x2 - x1

    q1 = np.empty((6,), dtype=np.float64)
    pypw85.spheroid(a1, c1, n1, q1)
    print(repr(q1))

    q2 = np.empty_like(q1)
    pypw85.spheroid(a2, c2, n2, q2)
    print(repr(q2))

    out = np.empty((2,), dtype=np.float64)
    pypw85.contact_function(r12, q1, q2, out)
    mu2, lambda_ = out
    print("μ² = {}".format(mu2))
    print("λ = {}".format(lambda_))

    Q1 = np.zeros((3, 3), dtype=np.float64)
    i, j = np.triu_indices_from(Q1)

    Q1[i, j] = q1
    Q1[j, i] = q1
    print(repr(Q1))

    Q2 = np.zeros_like(Q1)
    Q2[i, j] = q2
    Q2[j, i] = q2
    print(repr(Q2))

    Q1_inv = np.linalg.inv(Q1)
    Q2_inv = np.linalg.inv(Q2)

    f1 = lambda x: Q1_inv.dot(x).dot(x)
    print(f1((-a1, 0.0, 0.0)))
    print(f1((a1, 0.0, 0.0)))
    print(f1((0.0, -a1, 0.0)))
    print(f1((0.0, a1, 0.0)))
    print(f1((0.0, 0.0, -c1)))
    print(f1((0.0, 0.0, c1)))

    f2 = lambda x: Q2_inv.dot(x).dot(x)
    print(f2((c2, 0.0, 0.0)))
    print(f2((-c2, 0.0, 0.0)))
    print(f2((0.0, a2, 0.0)))
    print(f2((0.0, -a2, 0.0)))
    print(f2((0.0, 0.0, a2)))
    print(f2((0.0, 0.0, -a2)))

    Q = (1 - lambda_) * Q1 + lambda_ * Q2
    x = np.linalg.solve(Q, r12)

    x0a = x1 + (1 - lambda_) * np.dot(Q1, x)
    print(repr(x0a))

    x0b = x2 - lambda_ * np.dot(Q2, x)
    print(repr(x0b))

    x0 = x0a
    Q1 *= mu2
    Q2 *= mu2
    Q1_inv /= mu2
    Q2_inv /= mu2

    x = x0 - x1
    print(Q1_inv.dot(x).dot(x))
    x = x0 - x2
    print(Q2_inv.dot(x).dot(x))

    n1 = Q1_inv.dot(x0 - x1)
    n1 /= np.linalg.norm(n1)
    print(repr(n1))

    n2 = Q2_inv.dot(x0 - x2)
    n2 /= np.linalg.norm(n2)
    print(repr(n2))
