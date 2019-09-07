import os.path

import h5py
import numpy as np
import pytest
import scipy.linalg
import scipy.optimize

import pypw85

from numpy.testing import assert_allclose


def setup_module():
    np.random.seed(20180813)


def gen_directions():
    # Use some vertices of icosahedron
    golden_ratio = (1.+np.sqrt(5.))/2.
    norm = np.sqrt(1+golden_ratio**2)
    u_abs = 1./norm
    v_abs = golden_ratio/norm
    # out = []
    # for u in (-u_abs, u_abs):
    #     for v in (-v_abs, v_abs):
    #         out += [np.array((0., u, v), dtype=np.float64),
    #                 np.array((v, 0., u), dtype=np.float64),
    #                 np.array((u, v, 0.), dtype=np.float64)]
    out = [[0., u_abs, v_abs],
           [v_abs, 0., u_abs],
           [u_abs, v_abs, 0.]]
    return np.array(out, dtype=np.float64)


RADII = np.array([0.0199, 1.999, 9.999], dtype=np.float64)
DIRECTIONS = gen_directions()


def to_array_2d(a):
    return np.array([[a[0], a[1], a[2]],
                     [a[1], a[3], a[4]],
                     [a[2], a[4], a[5]]])


@pytest.mark.parametrize('l', [[1., 2., 3., 4., 5., 6.],
                               [1., -2., -3., 4., -5., 6.],
                               [1e5, -2e-5, -3e-5, 4., -5e-5, 6.e-5]])
def test__cholesky_decomp(l, rtol=1e-15, atol=1e-15):
    l = np.asarray(l)
    l_mat = np.array([[l[0], 0., 0.],
                      [l[1], l[3], 0.],
                      [l[2], l[4], l[5]]], dtype=np.float64)
    a_mat = l_mat.dot(l_mat.T)
    a = a_mat[[0, 0, 0, 1, 1, 2],
              [0, 1, 2, 1, 2, 2]]
    assert_allclose(pypw85._cholesky_decomp(a), l,
                    rtol=rtol, atol=atol)


@pytest.mark.parametrize('l', [[1., 2., 3., 4., 5., 6.],
                               [1., -2., -3., 4., -5., 6.]])
@pytest.mark.parametrize('x', [[1.2, -3.4, 5.7]])
def test__cholesky_solve(l, x, rtol=3e-15, atol=3e-15):
    l = np.asarray(l)
    x = np.asarray(x)
    l_mat = np.array([[l[0], 0., 0.],
                      [l[1], l[3], 0.],
                      [l[2], l[4], l[5]]], dtype=np.float64)
    a_mat = l_mat.dot(l_mat.T)
    b = a_mat.dot(x)
    assert_allclose(pypw85._cholesky_solve(l, b), x,
                    rtol=rtol, atol=atol)


@pytest.mark.parametrize('a', RADII)
@pytest.mark.parametrize('c', RADII)
@pytest.mark.parametrize('n', DIRECTIONS)
@pytest.mark.parametrize('in_place', [False, True])
def test_spheroid(a, c, n, in_place, rtol=1E-10, atol=1E-12):
    a2 = a**2
    c2 = c**2
    n = np.asarray(n)
    out = np.empty((6,), dtype=np.float64) if in_place else None
    q = pypw85.spheroid(a, c, n, out)
    if in_place:
        assert q is out
    q_mat = np.zeros((3, 3), dtype=np.float64)
    # scipy.linalg.eigh only needs the lower-part of the matrix to be filled
    q_mat[:, 0] = q[:3]
    q_mat[1:, 1] = q[3:5]
    q_mat[2, 2] = q[5]
    w_act, v_act = scipy.linalg.eigh(q_mat)
    if a2 <= c2:
        w_exp = (a2, a2, c2)
        n_act = v_act[:, -1]
    else:
        w_exp = (c2, a2, a2)
        n_act = v_act[:, 0]
    assert_allclose(w_act, w_exp, rtol, atol)
    if a != c:
        s = np.sign(np.dot(n_act, n))
        assert_allclose(s*n_act, n, rtol, atol)


def _test_contact_function(r, a1, c1, n1, a2, c2, n2, rtol=1E-9, atol=1E-11):
    # Test contact function with full output (μ² and λ)
    q1 = pypw85.spheroid(a1, c1, n1)
    q2 = pypw85.spheroid(a2, c2, n2)

    out = np.empty((2,), dtype=np.float64)
    pypw85.contact_function(r, q1, q2, out)
    μ2, λ = out

    Q1, Q2 = to_array_2d(q1), to_array_2d(q2)
    Q = (1-λ)*Q1 + λ*Q2

    # Check that stationarity condition
    #     λQ₁⁻¹⋅[x₀(λ)-c₁] + (1-λ)Q₂⁻¹⋅[x₀(λ)-c₂] = 0
    # is satisfied.
    y = np.linalg.solve(Q, r)
    c1_x0 = np.dot((1-λ)*Q1, y)
    c2_x0 = -np.dot(λ*Q2, y)
    r_actual = c1_x0-c2_x0
    assert_allclose(r_actual, r, rtol, atol)

    # Check that point x0 belongs to both scaled ellipsoids
    assert_allclose(μ2, np.linalg.solve(Q1, c1_x0).dot(c1_x0), 2e-5, 1e-6)
    assert_allclose(μ2, np.linalg.solve(Q2, c2_x0).dot(c2_x0), 2e-5, 1e-6)

    # Check that full and partial outputs are consistent
    assert_allclose(μ2, pypw85.contact_function(r, q1, q2), rtol, atol)


@pytest.mark.parametrize('r', DIRECTIONS)
@pytest.mark.parametrize('a1', RADII)
@pytest.mark.parametrize('c1', RADII)
@pytest.mark.parametrize('n1', DIRECTIONS)
@pytest.mark.parametrize('a2', RADII)
@pytest.mark.parametrize('c2', RADII)
@pytest.mark.parametrize('n2', DIRECTIONS)
def test_contact_function_fixed_cc_distance(r, a1, c1, n1, a2, c2, n2):
    _test_contact_function(r, a1, c1, n1, a2, c2, n2)


@pytest.mark.parametrize('r', DIRECTIONS[0, :]*RADII[:, None])
@pytest.mark.parametrize('a1', RADII)
@pytest.mark.parametrize('c1', RADII)
@pytest.mark.parametrize('n1', [DIRECTIONS[1]])
@pytest.mark.parametrize('a2', RADII)
@pytest.mark.parametrize('c2', RADII)
@pytest.mark.parametrize('n2', [DIRECTIONS[2]])
def test_contact_function_variable_cc_distance(r, a1, c1, n1, a2, c2, n2):
    _test_contact_function(r, a1, c1, n1, a2, c2, n2)


@pytest.mark.parametrize('r12_dir', DIRECTIONS)
@pytest.mark.parametrize('a1', RADII)
@pytest.mark.parametrize('c1', RADII)
@pytest.mark.parametrize('n1', DIRECTIONS)
@pytest.mark.parametrize('a2', RADII)
@pytest.mark.parametrize('c2', RADII)
@pytest.mark.parametrize('n2', DIRECTIONS)
def test_f(r12_dir, a1, c1, n1, a2, c2, n2,
           num_lambdas=11, rtol=1e-10, atol=1e-10):
    q1 = pypw85.spheroid(a1, c1, n1)
    q2 = pypw85.spheroid(a2, c2, n2)
    q1_mat = to_array_2d(q1)
    q2_mat = to_array_2d(q2)
    expected = []
    actual = []
    for lambda_i in np.linspace(0., 1., num=num_lambdas):
        q_mat = (1.-lambda_i)*q1_mat+lambda_i*q2_mat
        exp_ref = (lambda_i*(1.-lambda_i)
                   *np.dot(np.linalg.solve(q_mat, r12_dir), r12_dir))
        for r12_norm_j in RADII:
            expected.append(r12_norm_j**2*exp_ref)
            actual.append(pypw85.f(lambda_i, r12_norm_j*r12_dir, q1, q2))
    assert_allclose(actual, expected, rtol=rtol, atol=atol)
