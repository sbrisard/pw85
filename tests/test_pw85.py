import numpy as np
import pytest
import scipy.linalg

import pypw85

from numpy.testing import assert_allclose


def setup_module():
    np.random.seed(20180813)


def generate_directions(num_theta, num_phi):
    cos_theta = np.linspace(-1.0, 1.0, num=num_theta)
    sin_theta = np.sqrt(1.0-cos_theta**2)
    out = []
    for cos_theta in np.linspace(-1.0, 1.0, num=num_theta):
        sin_theta = np.sqrt(1.0-cos_theta**2)
        for phi in np.linspace(0., 2.*np.pi, num=num_phi, endpoint=False):
            out.append(np.array((sin_theta*np.cos(phi),
                                 sin_theta*np.sin(phi),
                                 cos_theta)))
    return out


RADII = [0.0199, 1.999, 9.999]
DIRECTIONS = generate_directions(7, 7)
BOOLEANS = [False, True]


@pytest.mark.parametrize('a', RADII)
@pytest.mark.parametrize('c', RADII)
@pytest.mark.parametrize('n', DIRECTIONS)
@pytest.mark.parametrize('in_place', BOOLEANS)
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


def to_array_2d(a):
    return np.array([[a[0], a[1], a[2]],
                     [a[1], a[3], a[4]],
                     [a[2], a[4], a[5]]])


@pytest.mark.parametrize('a', np.random.rand(100, 6))
def test__det_sym(a, rtol=1E-12, atol=1E-14):

    expected = np.linalg.det(to_array_2d(a))
    actual = pypw85._det_sym(*a)
    assert_allclose(actual, expected, rtol, atol)


@pytest.mark.parametrize('a1', [2.0])
@pytest.mark.parametrize('c1', [3.0])
@pytest.mark.parametrize('n1', DIRECTIONS[:1])
@pytest.mark.parametrize('a2', [0.04])
@pytest.mark.parametrize('c2', [5.0])
@pytest.mark.parametrize('n2', DIRECTIONS[:1])
def test_detQ_as_poly(a1, c1, n1, a2, c2, n2, rtol=1E-12, atol=1E-14):
    q1 = pypw85.spheroid(a1, c1, n1)
    Q1 = to_array_2d(q1)
    q2 = pypw85.spheroid(a2, c2, n2)
    Q2 = to_array_2d(q2)
    x = np.linspace(0., 1., num=11)
    b = np.poly1d(pypw85.detQ_as_poly(q1, q2)[::-1])
    actual = b(x)

    x = x[:, None, None]
    expected = np.linalg.det((1.-x)*Q1+x*Q2)
    assert_allclose(actual, expected, rtol, atol)


def minor(A, i, j):
    assert A.ndim >= 2
    assert A.shape[-2] == A.shape[-1]
    n = A.shape[-1]
    rows = list(range(n))
    rows.remove(i)
    cols = list(range(n))
    cols.remove(j)
    rows, cols = np.meshgrid(rows, cols, indexing='ij')
    return A[..., rows, cols]


def adjugate(A):
    adjA = np.empty_like(A)
    for i in range(3):
        for j in range(3):
            adjA[..., i, j] = (-1)**(i+j)*np.linalg.det(minor(A, i, j))
    return adjA


@pytest.mark.parametrize('x', np.random.rand(5, 3))
@pytest.mark.parametrize('a', np.random.rand(5, 6))
def test__xT_adjA_x(x, a, rtol=1E-12, atol=1E-14):
    actual = pypw85._xT_adjA_x(*x, *a)
    expected = np.dot(x, np.dot(adjugate(to_array_2d(a)), x))
    assert_allclose(actual, expected, rtol, atol)


@pytest.mark.parametrize('r', [np.array([3.0, 4.0, 5.0])])
@pytest.mark.parametrize('a1', [2.0])
@pytest.mark.parametrize('c1', [3.0])
@pytest.mark.parametrize('n1', DIRECTIONS[:1])
@pytest.mark.parametrize('a2', [0.04])
@pytest.mark.parametrize('c2', [5.0])
@pytest.mark.parametrize('n2', DIRECTIONS[:1])
def test_rT_adjQ_r_as_poly(r, a1, c1, n1, a2, c2, n2, rtol=1E-12, atol=1E-14):
    q1 = pypw85.spheroid(a1, c1, n1)
    q2 = pypw85.spheroid(a2, c2, n2)
    x = np.linspace(0., 1., num=11)
    actual = np.poly1d(pypw85.rT_adjQ_r_as_poly(r, q1, q2)[::-1])(x)
    x = x[:, None, None]
    Q = (1-x)*to_array_2d(q1)+x*to_array_2d(q2)
    expected = np.dot(np.dot(adjugate(Q), r), r)
    assert_allclose(actual, expected, rtol, atol)
