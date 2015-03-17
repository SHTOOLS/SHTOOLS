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
import matplotlib.pyplot as plt

#import shtools:
sys.path.append(os.path.join(os.path.dirname(__file__), "../../.."))
import pyshtools as shtools

#set shtools plot style:
sys.path.append(os.path.join(os.path.dirname(__file__), "../Common"))
from FigStyle import style_shtools
mpl.rcParams.update(style_shtools)

#==== MAIN FUNCTION ====
def main():
    example()

#==== PLOT POWER SPECTRA ====
def example():
    lmax   = 200
    scale  = 10
    power  = 1./(1.+(np.arange(lmax+1)/scale)**2)**2

    coeffs = shtools.SHCoefficients.from_random(power)
    coeffs.plot_bandpower(show=False)

    grid = coeffs.expand(kind='DH1')
    grid.plot_rawdata(show=False)

    coeffs.rotate(40.,0.,0.,degrees=True)
    grid2 = coeffs.expand(kind='DH1')
    grid2.plot_rawdata()

#==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
