import os.path
import shutil

import h5py
import numpy as np
import pytest
import requests

import pypw85.legacy

from numpy.testing import assert_allclose


PROXIES = {"http": "http://proxy.enpc.fr:3128", "https": "https://proxy.enpc.fr:3128"}

REF_DATA_URL = "https://zenodo.org/record/3323683/files/pw85_ref_data-20190712.h5"
REF_DATA_PATH = "pw85_ref_data.h5"


def setup_module():
    np.random.seed(20180813)
    if not os.path.exists(REF_DATA_PATH):
        r = requests.get(REF_DATA_URL, proxies=PROXIES, stream=True)
        if r.status_code == 200:
            with open(REF_DATA_PATH, "wb") as f:
                shutil.copyfileobj(r.raw, f)
        else:
            raise RuntimeError("could not retrieve reference data")


@pytest.fixture(scope="module")
def directions():
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


@pytest.fixture(scope="module")
def radii():
    return np.array([0.0199, 1.999, 9.999], dtype=np.float64)


@pytest.fixture(scope="module")
def distances():
    return np.array([0.15, 1.1, 11.0], dtype=np.float64)


@pytest.fixture(scope="module")
def spheroids(radii, directions):
    num_radii = radii.shape[0]
    num_directions = directions.shape[0]
    num_spheroids = num_radii * num_radii * num_directions
    out = np.empty((num_spheroids, 6), dtype=np.float64)
    i = 0
    for a in radii:
        for c in radii:
            for n in directions:
                out[i, :] = pypw85.spheroid(a, c, n)
    return out


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


def test__detQ_as_poly(spheroids, rtol=1e-10, atol=1e-8):
    for q1 in spheroids:
        Q1 = to_array_2d(q1)
        for q2 in spheroids:
            Q2 = to_array_2d(q2)
            x = np.linspace(0.0, 1.0, num=11)
            b = np.poly1d(pypw85.legacy._detQ_as_poly(q1, q2)[::-1])
            actual = b(x)

            x = x[:, None, None]
            expected = np.linalg.det((1.0 - x) * Q1 + x * Q2)
            assert_allclose(actual, expected, rtol, atol)


def test__rT_adjQ_r_as_poly(distances, directions, spheroids, rtol=1e-14, atol=1e-15):
    for r12 in distances:
        for n12 in directions:
            r12_vec = r12 * n12
            for q1 in spheroids:
                for q2 in spheroids:
                    x = np.linspace(0.0, 1.0, num=11)
                    actual = np.poly1d(
                        pypw85.legacy._rT_adjQ_r_as_poly(r12_vec, q1, q2)[::-1]
                    )(x)
                    x = x[:, None, None]
                    Q = (1 - x) * to_array_2d(q1) + x * to_array_2d(q2)
                    expected = np.dot(np.dot(adjugate(Q), r12_vec), r12_vec)
                    assert_allclose(actual, expected, rtol, atol)


def test_f1():
    rtol = 1e-10
    atol = 1e-15
    with h5py.File(REF_DATA_PATH, "r") as f:
        spheroids = np.array(f["spheroids"])
        directions = np.array(f["directions"])
        lambdas = np.array(f["lambdas"])
        expecteds = np.array(f["F"])
    for i1, q1 in enumerate(spheroids):
        for i2, q2 in enumerate(spheroids):
            for i, r12 in enumerate(directions):
                for j, lambda_ in enumerate(lambdas):
                    exp = expecteds[i1, i2, i, j]
                    act = pypw85.legacy.f1(lambda_, r12, q1, q2)
                    assert np.abs(act - exp) <= rtol * exp + atol
