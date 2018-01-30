#!/usr/bin/env python
"""
This script tests the python class interface
"""

from __future__ import absolute_import, division, print_function

import os
import sys

import numpy as np
import matplotlib as mpl

sys.path.append(os.path.join(os.path.dirname(__file__), "../../.."))
import pyshtools

sys.path.append(os.path.join(os.path.dirname(__file__), "../Common"))
from FigStyle import style_shtools

# set shtools plot style:
mpl.rcParams.update(style_shtools)


# ==== MAIN FUNCTION ====


def main():
    example1()


def example1():
    # generate random spherical harmonics coefficients with a given
    # power spectrum and plot its bandspectrum
    degrees = np.arange(201)
    scale = 10
    power = 1. / (1. + (degrees / scale)**2)**2

    coeffs = pyshtools.SHCoeffs.from_random(power)
    coeffs.plot_spectrum(show=False, fname='power.png')

    # expand coefficients into a spatial grid and plot it
    grid1 = coeffs.expand(grid='DH1')
    grid1.plot(show=False, fname='DHGrid_unrotated.png')

    grid2 = coeffs.expand(grid='GLQ')
    grid2.plot(show=False, fname='GLQGrid.png')

    # rotate coefficients, expand to grid and plot again
    coeffs_rot = coeffs.rotate(40., 0., 0., degrees=True)
    grid3 = coeffs_rot.expand(grid='DH1')
    grid3.plot(show=False, fname='DHGrid_rotated.png')


# ==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
