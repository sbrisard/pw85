import numpy as np
import scipy.linalg

import pw85


def generate_directions(num_theta, num_phi):
    cos_theta = np.linspace(-1.0, 1.0, num=num_theta)
    sin_theta = np.sqrt(1.0-cos_theta**2)
    phi = np.linspace(0., 2.*np.pi, num=num_phi, endpoint=False)
    out = np.empty((num_theta, num_phi, 3), dtype=np._float64)
    out[:, :, 0] = sin_theta*np.cos(phi)
    out[:, :, 1] = sin_theta*np.sin(phi)
    out[:, :, 2] = cos_theta
    return out


def test_spheroid():
    n = np.array([0., 0., 1.], dtype=np.float64)
    q = pw85.spheroid(1.0, 0.1, n, None)
    q = np.asarray(q)
    q_mat = np.empty((3, 3), dtype=np.float64)
    # scipy.linalg.eigh only needs the lower-part of the matrix to be filled
    q_mat[:, 0] = q[:3]
    q_mat[1:, 1] = q[3:5]
    q_mat[2, 2] = q[5]
    w, v = scipy.linalg.eigh(q_mat)
    print(w)
