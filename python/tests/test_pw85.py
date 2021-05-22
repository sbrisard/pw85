import numpy as np
import pytest

import pypw85

from numpy.testing import assert_allclose, assert_equal


def empty_vec():
    return np.empty((3,), dtype=np.float64)


def empty_sym():
    return np.empty((6,), dtype=np.float64)


@pytest.mark.parametrize(
    "a, exp, rtol",
    [
        (
            np.array([4, 2, 6, 17, 23, 70], dtype=np.float64),
            np.array([2, 1, 3, 4, 5, 6], dtype=np.float64),
            1e-15,
        ),
        (
            np.array([4, -2, 6, 17, -23, 70], dtype=np.float64),
            np.array([2, -1, 3, 4, -5, 6], dtype=np.float64),
            1e-15,
        ),
        (
            np.array(
                [1e10, -2, -3, 16 + 1.0 / 25e8, -0.02 + 3.0 / 5e9, 29.0 / 8e3 - 9e-10],
                dtype=np.float64,
            ),
            np.array([1e5, -2e-5, -3e-5, 4, -5e-3, 6e-2], dtype=np.float64),
            1e-6,
        ),
    ],
)
def test_cholesky_decomp(a, exp, rtol):
    a = np.asarray(a)
    exp = np.asarray(exp)
    act = empty_sym()
    pypw85._cholesky_decomp(a, act)
    assert_allclose(act, exp, rtol=rtol, atol=0.0)


# @pytest.mark.parametrize(
#     "a, expected, rtol",
#     [
#         (
#             np.array([4, 2, 6, 17, 23, 70], dtype=np.float64),
#             np.array([2, 1, 3, 4, 5, 6], dtype=np.float64),
#             1e-15,
#         ),
#         (
#             np.array([4, -2, 6, 17, -23, 70], dtype=np.float64),
#             np.array([2, -1, 3, 4, -5, 6], dtype=np.float64),
#             1e-15,
#         ),
#         (
#             np.array([1e10, -2, -3, 16 + 1 / 25e8, -0.02 + 3 / 5e9, 29 / 8000 - 9e-10]),
#             np.array([1e5, -2e-5, -3e-5, 4, -5e-3, 6e-2]),
#             1e-6,
#         ),
#     ],
# )
# def test__cholesky_decomp(a, expected, rtol, atol=1e-15):
#     actual = pw85._cholesky_decomp(a)
#     assert_allclose(actual, expected, rtol=rtol, atol=atol)
#
#
# @pytest.mark.parametrize(
#     "l, b, expected, rtol",
#     [
#         (
#             np.array([1, 2, 3, 4, 5, 6], dtype=np.float64),
#             np.array([11.5, 82.6, 314.2]),
#             np.array([1.2, -3.4, 5.7]),
#             1e-15,
#         ),
#         (
#             np.array([1, -2, -3, 4, -5, 6], dtype=np.float64),
#             np.array([-9.1, -150.2, 443]),
#             np.array([1.2, -3.4, 5.7]),
#             4e-15,
#         ),
#     ],
# )
# def test__cholesky_solve(l, b, expected, rtol, atol=1e-15):
#     actual = pw85._cholesky_solve(l, b)
#     assert_allclose(actual, expected, rtol=rtol, atol=atol)
#
#
# @pytest.mark.parametrize(
#     "a, c, n, expected, rtol",
#     [
#         (
#             9.9900000000000002e000,
#             1.0200000000000000e000,
#             np.array(
#                 [
#                     0.0000000000000000e000,
#                     -5.2573111211913359e-001,
#                     -8.5065080835203988e-001,
#                 ]
#             ),
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
#             1e-15,
#         )
#     ],
# )
# @pytest.mark.parametrize("in_place", [False, True])
# def test_spheroid(a, c, n, in_place, expected, rtol, atol=1e-15):
#     out = np.empty((6,), dtype=np.float64) if in_place else None
#     actual = pw85.spheroid(a, c, n, out)
#     if in_place:
#         assert actual is out
#     assert_allclose(actual, expected, rtol=rtol, atol=atol)
#
#
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
