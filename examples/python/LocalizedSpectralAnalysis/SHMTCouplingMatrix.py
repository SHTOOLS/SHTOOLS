#!/usr/bin/env python
"""
This script tests the SHMTCouplingMatrix routine
"""
from __future__ import absolute_import, division, print_function

import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from pyshtools import localizedspectralanalysis
from FigStyle import style_shtools

sys.path.append(os.path.join(os.path.dirname(__file__), "../../.."))

# set shtools plot style:
sys.path.append(os.path.join(os.path.dirname(__file__), "../Common"))
mpl.rcParams.update(style_shtools)


def main():
    test_CouplingMatrix()


def test_CouplingMatrix():
    print('\n---- testing SHMTCouplingMatrix ----')
    lmax = 30
    lwin = 10
    nwins = 1

    sqrt_taper_power = np.zeros((lwin+1, nwins))
    sqrt_taper_power[:, 0] = np.hanning(lwin+1)
    Mmt = localizedspectralanalysis.SHMTCouplingMatrix(lmax, sqrt_taper_power)
    print(Mmt)

# ==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
