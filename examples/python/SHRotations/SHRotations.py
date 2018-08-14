#!/usr/bin/env python
"""
This script tests the rotation of spherical harmonic coefficients
"""
from __future__ import absolute_import, division, print_function

import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), "../../.."))
import pyshtools
from pyshtools import rotate

pyshtools.utils.figstyle()

def main():
    test_SHRotations()


def test_SHRotations():
    # ---- input parameters ----
    lmax = 3
    alpha, beta, gamma = 20., 90., 90.

    # ---- derived parameters ----
    mask = np.zeros((2, lmax + 1, lmax + 1), dtype=np.bool)
    for l in np.arange(lmax + 1):
        mask[:, l, :l + 1] = True
    mask[1, :, 0] = False

    angles = np.radians([alpha, beta, gamma])

    print('\n---- testing djpi2 ----')
    print('computing rotation matrix for Euler angles: ' +
          '({:2.2f},{:2.2f},{:2.2f})'
          .format(alpha, beta, gamma))
    dj_matrix = rotate.djpi2(lmax)

    print('\n---- testing SHRotateRealCoef ----')
    print('generating normal distributed complex coefficients with ' +
          'variance 1...')
    rcoeffs = np.random.normal(size=(2, lmax + 1, lmax + 1))
    rcoeffs[np.invert(mask)] = 0.
    rcoeffs_rot = rotate.SHRotateRealCoef(rcoeffs, angles, dj_matrix)
    print(rcoeffs_rot)

    print('\n---- testing SHRotateCoef ----')
    print('generating normal distributed complex coefficients with ' +
          'variance 1...')
    ccoeffs = np.random.normal(loc=0., scale=1.,
                               size=(2, (lmax + 1) * (lmax + 2) // 2))
    ccoeffs_rot = rotate.SHRotateCoef(angles, ccoeffs, dj_matrix)
    print(ccoeffs_rot)

# ==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
