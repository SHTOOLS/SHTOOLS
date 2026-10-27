"""
This script tests the rotation of spherical harmonic coefficients
"""
import numpy as np
import pyshtools as pysh

pysh.utils.figstyle()


def main():
    test_SHRotations()
    test_x_rotation_degrees()


def test_SHRotations():
    # ---- input parameters ----
    lmax = 3
    alpha, beta, gamma = 20., 90., 90.

    # ---- derived parameters ----
    mask = np.zeros((2, lmax + 1, lmax + 1), dtype=bool)
    for l in np.arange(lmax + 1):
        mask[:, l, :l + 1] = True
    mask[1, :, 0] = False

    angles = np.radians([alpha, beta, gamma])

    print('\n---- testing djpi2 ----')
    print('computing rotation matrix for Euler angles: ' +
          '({:2.2f},{:2.2f},{:2.2f})'
          .format(alpha, beta, gamma))
    dj_matrix = pysh.rotate.djpi2(lmax)

    print('\n---- testing SHRotateRealCoef ----')
    print('generating normal distributed complex coefficients with ' +
          'variance 1...')
    rcoeffs = np.random.normal(size=(2, lmax + 1, lmax + 1))
    rcoeffs[np.invert(mask)] = 0.
    rcoeffs_rot = pysh.rotate.SHRotateRealCoef(rcoeffs, angles, dj_matrix)
    print(rcoeffs_rot)

    print('\n---- testing SHRotateCoef ----')
    print('generating normal distributed complex coefficients with ' +
          'variance 1...')
    ccoeffs = np.random.normal(loc=0., scale=1.,
                               size=(2, (lmax + 1) * (lmax + 2) // 2))
    ccoeffs_rot = pysh.rotate.SHRotateCoef(angles, ccoeffs, dj_matrix)
    print(ccoeffs_rot)


def test_x_rotation_degrees():
    coeffs = np.zeros((2, 4, 4))
    coeffs[0, 1, 0] = 1.
    coeffs[0, 2, 1] = 0.2
    coeffs[1, 3, 2] = -0.4
    models = (
        pysh.SHCoeffs.from_array(coeffs),
        pysh.SHGravCoeffs.from_array(coeffs, gm=3.986e14, r0=6371200.),
        pysh.SHMagCoeffs.from_array(coeffs, r0=6371200.))
    angles = np.array([10., 30., 20.])
    for model in models:
        for body in (False, True):
            in_degrees = model.rotate(*angles, convention='x', body=body,
                                      backend='shtools')
            in_radians = model.rotate(*np.radians(angles), convention='x',
                                      degrees=False, body=body,
                                      backend='shtools')
            in_y = model.rotate(angles[0] - 90., angles[1], angles[2] + 90.,
                                convention='y', body=body, backend='shtools')
            np.testing.assert_allclose(in_degrees.coeffs, in_radians.coeffs,
                                       rtol=1.e-13, atol=1.e-13)
            np.testing.assert_allclose(in_degrees.coeffs, in_y.coeffs,
                                       rtol=1.e-13, atol=1.e-13)


# ==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
