import numpy as np

import pypw85

def to_array_2d(a):
    return np.array([[a[0], a[1], a[2]],
                     [a[1], a[3], a[4]],
                     [a[2], a[4], a[5]]])

if __name__ == '__main__':
    golden_ratio = (1.+np.sqrt(5.))/2.
    norm = np.sqrt(1+golden_ratio**2)
    u_abs = 1./norm
    v_abs = golden_ratio/norm

    r = np.array([0., u_abs, v_abs])
    n1 = np.array([0., u_abs, v_abs])
    n2 = np.array([0., u_abs, v_abs])
    a1 = c1 = 1.999
    a2 = c2 = 0.0199

    q1 = pypw85.spheroid(a1, c1, n1)
    q2 = pypw85.spheroid(a2, c2, n2)

    out = np.empty((2,), dtype=np.float64)
    pypw85.contact_function(r, q1, q2, out)
    μ2, λ = out

    Q1, Q2 = to_array_2d(q1), to_array_2d(q2)
    Q = (1-λ)*Q1 + λ*Q2

    y = np.linalg.solve(Q, r)

    c1_x0 = np.dot((1-λ)*Q1, y)
    c2_x0 = -np.dot(λ*Q2, y)
    r_actual = c1_x0-c2_x0
    print(r_actual, r)

    μ2_1 = np.linalg.solve(Q1, c1_x0).dot(c1_x0)
    μ2_2 = np.linalg.solve(Q2, c2_x0).dot(c2_x0)

    print(μ2, μ2_1, μ2_2)
