#!/usr/bin/env python
"""
This script tests the python class interface
"""

from __future__ import absolute_import, division, print_function

# standard imports:
import os
import sys
import numpy as np
import matplotlib as mpl

# import shtools:
sys.path.append(os.path.join(os.path.dirname(__file__), "../../.."))
import pyshtools as shtools

# set shtools plot style:
sys.path.append(os.path.join(os.path.dirname(__file__), "../Common"))
from FigStyle import style_shtools
mpl.rcParams.update(style_shtools)

#==== MAIN FUNCTION ====


def main():
    example1()

def example1():
    # generate random spherical harmonics coefficients with a given
    # power spectrum and plot its bandspectrum
    degrees = np.arange(201)
    scale = 10
    power = 1. / (1. + (degrees / scale)**2)**2

    coeffs = shtools.SHCoeffs.from_random(power)
    coeffs.plot_powerperband(show=False, fname='power.png')

    # expand coefficients into a spatial grid and plot it
    grid1 = coeffs.expand(grid='DH1')
    grid1.plot_rawdata(show=False, fname='DHGrid_unrotated.png')
    
    grid2 = coeffs.expand(grid='GLQ')
    grid2.plot_rawdata(show=False, fname='GLQGrid.png')

    # rotate coefficients, expand to grid and plot again
    coeffs.rotate(40., 0., 0., degrees=True)
    grid3 = coeffs.expand(grid='DH1')
    grid3.plot_rawdata(show=False, fname='DHGrid_rotated.png')


#==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
