import pathlib

import h5py
import numpy as np
import pytest

import pypw85

from numpy.testing import assert_allclose, assert_equal

# Initialization
with h5py.File(str(pathlib.Path.cwd() / ".." / "data" / "pw85_ref_data.h5")) as f:
    radii = np.asarray(f["radii"])
    directions = np.asarray(f["directions"])


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


@pytest.mark.parametrize("a", radii)
@pytest.mark.parametrize("c", radii)
@pytest.mark.parametrize("n", directions)
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


# @pytest.mark.parametrize(
#     "lambda_, r12, q1, q2, expected, rtol",
#     [
#         (
#             0.1234,
#             np.array([1, -2, 3], dtype=np.float64),
#             np.array(
#                 [
#                     9.9800100000000000e001,
#                     0.0000000000000000e000,
#                     0.0000000000000000e000,
#                     7.2503590263748606e001,
#                     -4.4166680527497185e001,
#                     2.8336909736251414e001,
#                 ]
#             ),
#             np.array(
#                 [
#                     7.5286760167683042e001,
#                     -0.0000000000000000e000,
#                     -4.6523719335366117e001,
#                     9.8010000000000007e-003,
#                     0.0000000000000000e000,
#                     2.8763040832316932e001,
#                 ]
#             ),
#             9.0045499998758230e-002,
#             1e-15,
#         )
#     ],
# )
# def test_f(lambda_, r12, q1, q2, expected, rtol, atol=1e-15):
#     actual = pw85.f(lambda_, r12, q1, q2)
#     assert np.abs(actual - expected) <= rtol * np.abs(expected) + atol
#
#
# @pytest.mark.parametrize(
#     "r12, q1, q2, expected, rtol",
#     [
#         (
#             np.array([1, -2, 3], dtype=np.float64),
#             np.array(
#                 [
#                     9.9800100000000000e001,
#                     0.0000000000000000e000,
#                     0.0000000000000000e000,
#                     7.2503590263748606e001,
#                     -4.4166680527497185e001,
#                     2.8336909736251414e001,
#                 ]
#             ),
#             np.array(
#                 [
#                     7.5286760167683042e001,
#                     -0.0000000000000000e000,
#                     -4.6523719335366117e001,
#                     9.8010000000000007e-003,
#                     0.0000000000000000e000,
#                     2.8763040832316932e001,
#                 ]
#             ),
#             np.array([1.9630383437733556e-001, 9.7798461290800798e-001]),
#             1e-15,
#         )
#     ],
# )
# @pytest.mark.parametrize("in_place", [False, True])
# def test_contact_function(r12, q1, q2, in_place, expected, rtol, atol=1e-15):
#     out = np.empty((2,), dtype=np.float64) if in_place else None
#     actual = np.array(pw85.contact_function(r12, q1, q2, out))
#     assert_allclose(actual, expected, rtol=rtol, atol=atol)
#     if in_place:
#         assert_equal(out, actual)
