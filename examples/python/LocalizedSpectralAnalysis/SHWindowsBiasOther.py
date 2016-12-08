#!/usr/bin/env python
"""
This script tests other routines related to localized spectral analyses
"""
from __future__ import absolute_import, division, print_function

import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), "../../.."))
from pyshtools import shtools

sys.path.append(os.path.join(os.path.dirname(__file__), "../Common"))
from FigStyle import style_shtools

# set shtools plot style:
mpl.rcParams.update(style_shtools)


def main():
    test_LocalizationWindows()
    # test_LocalizationBias()
    # test_Other()


def test_LocalizationWindows():

    print('\n---- testing SphericalCapCoef ----')
    lmax = 15
    theta = 50.
    print('generating {:2.1f} degrees cap:'.format(theta))
    coeffsm0 = shtools.SphericalCapCoef(np.radians(theta), lmax)
    print(coeffsm0)

    print('\n---- testing SHBias ----')
    winpower = coeffsm0**2
    ldata = 20
    globalpower = np.random.rand(ldata)
    localpower = shtools.SHBias(winpower, globalpower)
    print(localpower[:min(ldata, 20)])

    print('\n---- testing Curve2Mask ----')
    # defines lat/lon square (bug!?)
    nlat = 100
    dlat = 180. / nlat
    npoints = 4
    points = np.empty((npoints, 2))
    points[0] = [40., 20.]
    points[1] = [40., 50.]
    points[2] = [80., 50.]
    points[3] = [80., 20.]
    hasnorthpole = False
    dhmask = shtools.Curve2Mask(nlat, points, hasnorthpole)
    # compute covered area as a check
    thetas = np.linspace(0 + dlat / 2., 180. - dlat / 2., nlat)
    weights = 2 * np.sin(np.radians(thetas))
    maskarea = np.sum(dhmask * weights[:, None] * dlat**2)
    globearea = 4 * np.pi * (180 / np.pi)**2
    print('mask covers {:2.2f}%% of the globe'
          .format(100 * maskarea / globearea))
    fig = plt.figure()
    plt.imshow(dhmask)
    fig.savefig('mask.png')

    print('\n---- testing ComputeDMap ----')
    nlat = 180
    lmax = 20
    dlat = 180. / nlat
    lats = np.linspace(90. - dlat / 2., -90. + dlat / 2., nlat)
    lons = np.linspace(0. + dlat, 360. - dlat / 2., nlat)
    latgrid, longrid = np.meshgrid(lats, lons, indexing='ij')
    dh_mask = np.logical_and(5. < latgrid, latgrid < 20.)
    print('dij matrix[0,:lmax={:d}]:'.format(lmax))
    dij_matrix = shtools.ComputeDMap(dh_mask, lmax)
    print(dij_matrix[0, :lmax])

    print('\n---- testing SHReturnTapersMap ----')
    tapers, evalues = shtools.SHReturnTapersMap(dh_mask, lmax, ntapers=1)
    print('best taper concentration: {:2.2f}'.format(evalues[0]))

# ==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
