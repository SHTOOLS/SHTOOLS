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
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), "../../.."))
import pyshtools

sys.path.append(os.path.join(os.path.dirname(__file__), "../Common"))
from FigStyle import style_shtools

# set shtools plot style:
mpl.rcParams.update(style_shtools)


# ==== MAIN FUNCTION ====

def main():
    example1()
    example2()


# ==== EXAMPLES ====
def example1():
    # generate cap window
    lmax = 20
    nwin = 20
    theta = 25.
    cap = pyshtools.SHWindow.from_cap(theta, lmax, nwin=nwin)
    cap.info()
    cap.plot_windows(20, show=False, fname='cap_tapers.png')
    cap.plot_coupling_matrix(30, nwin=5, show=False, fname='cap_coupling.png')


# ==== EXAMPLES ====
def example2():
    # generate cap window
    lmax = 15
    nwins = 15

    coeffs = pyshtools.SHCoeffs.from_file('../../ExampleDataFiles/srtmp300.msl')
    topo = coeffs.expand(grid='DH2')
    dh_mask = topo.data > 0.
    print(dh_mask.shape)
    region = pyshtools.SHWindow.from_mask(dh_mask, lmax, nwins)
    region.info()
    region.plot_windows(nwins, show=False, fname='continent_tapers.png')
    region.plot_coupling_matrix(30, 5, show=False,
                                fname='continent_coupling.png')

# ==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
