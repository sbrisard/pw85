import os.path

import h5py
import numpy as np
import pytest
import scipy.linalg
import scipy.optimize

import pypw85.legacy

from numpy.testing import assert_allclose


def setup_module():
    np.random.seed(20180813)


def gen_directions():
    # Use some vertices of icosahedron
    golden_ratio = (1.0 + np.sqrt(5.0)) / 2.0
    norm = np.sqrt(1 + golden_ratio ** 2)
    u_abs = 1.0 / norm
    v_abs = golden_ratio / norm
    # out = []
    # for u in (-u_abs, u_abs):
    #     for v in (-v_abs, v_abs):
    #         out += [np.array((0., u, v), dtype=np.float64),
    #                 np.array((v, 0., u), dtype=np.float64),
    #                 np.array((u, v, 0.), dtype=np.float64)]
    out = [[0.0, u_abs, v_abs], [v_abs, 0.0, u_abs], [u_abs, v_abs, 0.0]]
    return np.array(out, dtype=np.float64)


RADII = np.array([0.0199, 1.999, 9.999], dtype=np.float64)
DIRECTIONS = gen_directions()


def to_array_2d(a):
    return np.array([[a[0], a[1], a[2]], [a[1], a[3], a[4]], [a[2], a[4], a[5]]])


@pytest.mark.parametrize("a", np.random.rand(100, 6))
def test__det_sym(a, rtol=1e-12, atol=1e-14):

    expected = np.linalg.det(to_array_2d(a))
    actual = pypw85.legacy._det_sym(a)
    assert_allclose(actual, expected, rtol, atol)

def minor(A, i, j):
    assert A.ndim >= 2
    assert A.shape[-2] == A.shape[-1]
    n = A.shape[-1]
    rows = list(range(n))
    rows.remove(i)
    cols = list(range(n))
    cols.remove(j)
    rows, cols = np.meshgrid(rows, cols, indexing="ij")
    return A[..., rows, cols]

def adjugate(A):
    adjA = np.empty_like(A)
    for i in range(3):
        for j in range(3):
            adjA[..., i, j] = (-1) ** (i + j) * np.linalg.det(minor(A, i, j))
    return adjA


@pytest.mark.parametrize("x", np.random.rand(5, 3))
@pytest.mark.parametrize("a", np.random.rand(5, 6))
def test__xT_adjA_x(x, a, rtol=1e-12, atol=1e-14):
    actual = pypw85.legacy._xT_adjA_x(x, a)
    expected = np.dot(x, np.dot(adjugate(to_array_2d(a)), x))
    assert_allclose(actual, expected, rtol, atol)


@pytest.mark.parametrize("a1", RADII)
@pytest.mark.parametrize("c1", RADII)
@pytest.mark.parametrize("n1", DIRECTIONS)
@pytest.mark.parametrize("a2", RADII)
@pytest.mark.parametrize("c2", RADII)
@pytest.mark.parametrize("n2", DIRECTIONS)
def test__detQ_as_poly(a1, c1, n1, a2, c2, n2, rtol=1e-10, atol=1e-8):
    q1 = pypw85.spheroid(a1, c1, n1)
    Q1 = to_array_2d(q1)
    q2 = pypw85.spheroid(a2, c2, n2)
    Q2 = to_array_2d(q2)
    x = np.linspace(0.0, 1.0, num=11)
    b = np.poly1d(pypw85.legacy._detQ_as_poly(q1, q2)[::-1])
    actual = b(x)

    x = x[:, None, None]
    expected = np.linalg.det((1.0 - x) * Q1 + x * Q2)
    assert_allclose(actual, expected, rtol, atol)






# def _test__rT_adjQ_r_as_poly(r, a1, c1, n1, a2, c2, n2, rtol=1e-10, atol=1e-8):
#     q1 = pypw85.spheroid(a1, c1, n1)
#     q2 = pypw85.spheroid(a2, c2, n2)
#     x = np.linspace(0.0, 1.0, num=11)
#     actual = np.poly1d(pypw85._rT_adjQ_r_as_poly(r, q1, q2)[::-1])(x)
#     x = x[:, None, None]
#     Q = (1 - x) * to_array_2d(q1) + x * to_array_2d(q2)
#     expected = np.dot(np.dot(adjugate(Q), r), r)
#     assert_allclose(actual, expected, rtol, atol)


# @pytest.mark.parametrize("r", DIRECTIONS)
# @pytest.mark.parametrize("a1", RADII)
# @pytest.mark.parametrize("c1", RADII)
# @pytest.mark.parametrize("n1", DIRECTIONS)
# @pytest.mark.parametrize("a2", RADII)
# @pytest.mark.parametrize("c2", RADII)
# @pytest.mark.parametrize("n2", DIRECTIONS)
# def test__rT_adjQ_r_as_poly_fixed_cc_distance(r, a1, c1, n1, a2, c2, n2):
#     _test__rT_adjQ_r_as_poly(r, a1, c1, n1, a2, c2, n2)


# @pytest.mark.parametrize("r", DIRECTIONS[0, :] * RADII[:, None])
# @pytest.mark.parametrize("a1", RADII)
# @pytest.mark.parametrize("c1", RADII)
# @pytest.mark.parametrize("n1", [DIRECTIONS[1]])
# @pytest.mark.parametrize("a2", RADII)
# @pytest.mark.parametrize("c2", RADII)
# @pytest.mark.parametrize("n2", [DIRECTIONS[2]])
# def test__rT_adjQ_r_as_poly_variable_cc_distance(r, a1, c1, n1, a2, c2, n2):
#     _test__rT_adjQ_r_as_poly(r, a1, c1, n1, a2, c2, n2)
