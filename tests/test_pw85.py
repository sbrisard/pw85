import itertools

import numpy as np
import pytest
import scipy.linalg

import pypw85

from numpy.testing import assert_allclose


def generate_directions(num_theta, num_phi):
    cos_theta = np.linspace(-1.0, 1.0, num=num_theta)
    sin_theta = np.sqrt(1.0-cos_theta**2)
    out = []
    for cos_theta in np.linspace(-1.0, 1.0, num=num_theta):
        sin_theta = np.sqrt(1.0-cos_theta**2)
        for phi in np.linspace(0., 2.*np.pi, num=num_phi, endpoint=False):
            out.append((sin_theta*np.cos(phi), sin_theta*np.sin(phi), cos_theta))
    return out


RADII = [0.0199, 1.999, 9.999]
DIRECTIONS = generate_directions(7, 7)
BOOLEANS = [False, True]


@pytest.mark.parametrize('a, c, n, in_place',
                         itertools.product(RADII, RADII, DIRECTIONS, BOOLEANS))
def test_spheroid(a, c, n, in_place, rtol=1E-7, atol=0):
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
        assert_allclose(s*n_act, n)
