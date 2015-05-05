#!/usr/bin/env python
"""
This script tests the conversions between real and complex spherical harmonics
coefficients
"""

# standard imports:
import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# import shtools:
sys.path.append(os.path.join(os.path.dirname(__file__), "../../.."))
import pyshtools as shtools

# set shtools plot style:
sys.path.append(os.path.join(os.path.dirname(__file__), "../Common"))
from FigStyle import style_shtools
mpl.rcParams.update(style_shtools)


def main():
    test_LocalizationWindows()
    # test_LocalizationBias()
    # test_Other()


def test_LocalizationWindows():
    print '\n---- testing SHMTCouplingMatrix ----'
    lmax  = 30
    lwin = 10
    nwins = 1

    sqrt_taper_power = np.zeros( (lwin+1,nwins) )
    sqrt_taper_power[:,0] = np.hanning(lwin+1)
    Mmt = shtools.SHMTCouplingMatrix(lmax,sqrt_taper_power)
    print Mmt

#==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
