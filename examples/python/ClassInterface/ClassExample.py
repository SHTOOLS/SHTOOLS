#!/usr/bin/env python
"""
This script tests the python class interface
"""

from __future__ import division
from __future__ import print_function

#standard imports:
import os, sys
import numpy as np
import matplotlib as mpl

#import shtools:
sys.path.append(os.path.join(os.path.dirname(__file__), "../../.."))
from pyshtools import SHCoeffs

#set shtools plot style:
sys.path.append(os.path.join(os.path.dirname(__file__), "../Common"))
from FigStyle import style_shtools
mpl.rcParams.update(style_shtools)

#==== MAIN FUNCTION ====
def main():
    # generate random spherical harmonics coefficients with a given
    # power spectrum and plot its bandspectrum
    degrees = np.arange(201)
    scale   = 10
    power   = 1./(1.+(degrees/scale)**2)**2

    coeffs = SHCoeffs.from_random(power)
    coeffs.plot_powerperband(show=False)

    #expand coefficients into a spatial grid and plot it
    grid1 = coeffs.expand(kind='DH1')
    grid1.plot_rawdata(show=False)

    #rotate coefficients, expand to grid and plot again
    coeffs.rotate(40.,0.,0.,degrees=True)
    grid2 = coeffs.expand(kind='DH1')
    grid2.plot_rawdata()

#==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
