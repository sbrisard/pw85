import h5py
import numpy as np
import pathlib
import pypw85
import pytest
from numpy.testing import assert_allclose

# Initialization
REF_DATA = {}
with h5py.File(str(pathlib.Path.cwd() / ".." / "data" / "pw85_ref_data.h5")) as f:
    REF_DATA["radii"] = np.asarray(f["radii"])
    REF_DATA["directions"] = np.asarray(f["directions"])
    REF_DATA["spheroids"] = np.asarray(f["spheroids"])
    REF_DATA["lambdas"] = np.asarray(f["lambdas"])
    REF_DATA["F"] = np.asarray(f["F"])

REF_DATA["distances"] = np.array([0.15, 1.1, 11.0], dtype=np.float64)


def empty_vec():
    return np.empty((3,), dtype=np.float64)


def empty_sym():
    return np.empty((6,), dtype=np.float64)


@pytest.mark.parametrize(
    "a, exp, rtol",
    [
        (
            [4.0, 2.0, 6.0, 17.0, 23.0, 70.0],
            [2.0, 1.0, 3.0, 4.0, 5.0, 6.0],
            1e-15,
        ),
        (
            [4.0, -2.0, 6.0, 17.0, -23.0, 70.0],
            [2.0, -1.0, 3.0, 4.0, -5.0, 6.0],
            1e-15,
        ),
        (
            [1e10, -2, -3, 16 + 1.0 / 25e8, -0.02 + 3.0 / 5e9, 29.0 / 8e3 - 9e-10],
            [1e5, -2e-5, -3e-5, 4, -5e-3, 6e-2],
            1e-6,
        ),
    ],
)
def test_cholesky_decomp(a, exp, rtol):
    act = empty_sym()
    pypw85._cholesky_decomp(a, act)
    assert_allclose(act, exp, rtol=rtol, atol=0.0)


@pytest.mark.parametrize(
    "l, b, exp, rtol",
    [
        ([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], [11.5, 82.6, 314.2], [1.2, -3.4, 5.7], 1e-15),
        ([1, -2, -3, 4, -5, 6], [-9.1, -150.2, 443], [1.2, -3.4, 5.7], 4e-15),
    ],
)
def test_cholesky_solve(l, b, exp, rtol):
    l = np.asarray(l)
    b = np.asarray(b)
    exp = np.asarray(exp)
    act = empty_vec()
    pypw85._cholesky_solve(l, b, act)
    assert_allclose(act, exp, rtol=rtol, atol=0.0)


@pytest.mark.parametrize("a", REF_DATA["radii"])
@pytest.mark.parametrize("c", REF_DATA["radii"])
@pytest.mark.parametrize("n", REF_DATA["directions"])
def test_spheroid(a, c, n):
    # Relative and absolute tolerance on the coefficients of the matrix
    # q to be computed and tested.
    rtol = atol = 1e-15

    a2, c2 = a ** 2, c ** 2

    q = empty_sym()
    pypw85.spheroid(a, c, n, q)

    # Check that n is an eigenvector.
    qn = np.array(
        [
            q[0] * n[0] + q[1] * n[1] + q[2] * n[2],
            q[1] * n[0] + q[3] * n[1] + q[4] * n[2],
            q[2] * n[0] + q[4] * n[1] + q[5] * n[2],
        ],
        dtype=np.float64,
    )
    abs_q = np.abs(q)
    Δq = rtol * abs_q + atol
    abs_n = np.abs(n)
    Δqn = np.array(
        [
            Δq[0] * abs_n[0] + Δq[1] * abs_n[1] + Δq[2] * abs_n[2],
            Δq[1] * abs_n[0] + Δq[3] * abs_n[1] + Δq[4] * abs_n[2],
            Δq[2] * abs_n[0] + Δq[4] * abs_n[1] + Δq[5] * abs_n[2],
        ],
        dtype=np.float64,
    )

    assert np.all(np.abs(qn - c2 * n) <= Δqn)

    # Check that the eigenvalues are [a, a, c]. It is sufficient to
    # compute the characteristic polynomial of the matrix q.

    # Check tr(q)
    exp = 2.0 * a2 + c2
    act = q[0] + q[3] + q[5]
    tol = Δq[0] + Δq[3] + Δq[5]
    assert np.abs(act - exp) <= tol

    # Check det(q)
    exp = a2 * a2 * c2
    act = (
        q[0] * q[3] * q[5]
        + 2.0 * q[1] * q[2] * q[4]
        - q[0] * q[4] * q[4]
        - q[3] * q[2] * q[2]
        - q[5] * q[1] * q[1]
    )
    tol = (
        Δq[0] * abs_q[3] * abs_q[5]
        + 2.0 * Δq[1] * abs_q[2] * abs_q[4]
        + Δq[0] * abs_q[4] * abs_q[4]
        + Δq[3] * abs_q[2] * abs_q[2]
        + Δq[5] * abs_q[1] * abs_q[1]
        + abs_q[0] * Δq[3] * abs_q[5]
        + +abs_q[0] * abs_q[3] * Δq[5]
        + 2.0
        * (
            abs_q[1] * Δq[2] * abs_q[4]
            + abs_q[0] * Δq[4] * abs_q[4]
            + abs_q[3] * Δq[2] * abs_q[2]
            + abs_q[5] * Δq[1] * abs_q[1]
            + abs_q[1] * abs_q[2] * Δq[4]
        )
    )
    assert np.abs(act - exp) <= tol

    # Check [tr(q)^2 - tr(q^2)]/2
    exp = a2 * (a2 + 2.0 * c2)
    act = (
        q[0] * q[3]
        + q[3] * q[5]
        + q[5] * q[0]
        - q[1] * q[1]
        - q[2] * q[2]
        - q[4] * q[4]
    )
    tol = (
        Δq[0] * abs_q[3]
        + abs_q[0] * Δq[3]
        + Δq[3] * abs_q[5]
        + abs_q[3] * Δq[5]
        + Δq[5] * abs_q[0]
        + abs_q[5] * Δq[0]
        + 2.0 * Δq[1] * abs_q[1]
        + 2.0 * Δq[2] * abs_q[2]
        + 2.0 * Δq[4] * abs_q[4]
    )
    assert np.abs(act - exp) <= tol


def test_f_neg():
    rtol = 1e-10
    for i1, q1 in enumerate(REF_DATA["spheroids"]):
        for i2, q2 in enumerate(REF_DATA["spheroids"]):
            for j, n in enumerate(REF_DATA["directions"]):
                for k, lambda_ in enumerate(REF_DATA["lambdas"]):
                    exp = REF_DATA["F"][i1, i2, j, k]
                    act = -pypw85.f_neg(lambda_, n, q1, q2)
                    assert np.abs(act - exp) <= rtol * np.abs(exp)


def _test_contact_function(r12, q1, q2):
    atol = 1e-15
    rtol = 1e-10

    out = np.empty((2,), dtype=np.float64)
    pypw85.contact_function(r12, q1, q2, out)
    mu2 = out[0]
    lambda_ = out[1]
    lambda1 = 1.0 - lambda_

    Q1 = np.array(
        [[q1[0], q1[1], q1[2]], [q1[1], q1[3], q1[4]], [q1[2], q1[4], q1[5]]],
        dtype=np.float64,
    )

    Q2 = np.array(
        [[q2[0], q2[1], q2[2]], [q2[1], q2[3], q2[4]], [q2[2], q2[4], q2[5]]],
        dtype=np.float64,
    )
    Q12 = Q2 - Q1
    Q = (1 - lambda_) * Q1 + lambda_ * Q2
    s = np.linalg.solve(Q, r12)
    u = Q12 @ s
    v = np.linalg.solve(Q, u)

    rs = r12.dot(s)
    su = s.dot(u)
    uv = u.dot(v)
    # We could also check that x0 belong to both (scaled) ellipsoids:
    #
    # mu^2 = c1x0.Q1^(-1).c1x0 = (1-lambda) c1x0.s
    #
    # and
    #
    # mu^2 = c2x0.Q2^(-1).c2x0 = -lambda c2x0.s.
    #
    # The code is provided below. However, the accuracy seems to be
    # quite poor.
    mu2_1 = lambda1 * lambda1 * (rs - lambda_ * su)
    mu2_2 = lambda_ * lambda_ * (rs + lambda1 * su)

    assert np.abs(mu2 - mu2_1) <= rtol * mu2_1 + atol
    assert np.abs(mu2 - mu2_2) <= rtol * mu2_2 + atol

    # Finally, we check that f'(lambda) = 0.
    # The tolerance is given by the tolerance eps on lambda (defined in
    # contact_function). The test therefore reads
    #
    #     |f'(lambda)| <= eps * |f"(lambda)|.
    #
    #  We first compute f' and f".
    lambda2 = 1.0 - 2.0 * lambda_
    f1 = lambda2 * rs - lambda_ * lambda1 * su
    f2 = -2.0 * rs - 2.0 * lambda2 * su + 2.0 * lambda_ * lambda1 * uv

    assert np.abs(f1) <= pypw85.lambda_atol * np.abs(f2)


def test_contact_function():
    for r in REF_DATA["distances"]:
        for n in REF_DATA["directions"]:
            r12 = r * n
            for q1 in REF_DATA["spheroids"]:
                for q2 in REF_DATA["spheroids"]:
                    _test_contact_function(r12, q1, q2)
