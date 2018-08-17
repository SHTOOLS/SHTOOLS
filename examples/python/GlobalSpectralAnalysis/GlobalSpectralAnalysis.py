#!/usr/bin/env python
"""
This script tests the different Spherical Harmonics Transforms on the Mars
topography data set
"""
from __future__ import absolute_import, division, print_function

import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), "../../.."))
import pyshtools
from pyshtools import spectralanalysis
from pyshtools import shio
from pyshtools import expand

pyshtools.utils.figstyle()

# ==== MAIN FUNCTION ====

def main():
    example()


# ==== PLOT POWER SPECTRA ====

def example():
    """
    example that plots the power spectrum of Mars topography data
    """
    # --- input data filename ---
    infile = os.path.join(os.path.dirname(__file__),
                          '../../ExampleDataFiles/MarsTopo719.shape')
    coeffs, lmax = shio.shread(infile)

    # --- plot grid ---
    grid = expand.MakeGridDH(coeffs, csphase=-1)
    fig_map = plt.figure()
    plt.imshow(grid)

    # ---- compute spectrum ----
    ls = np.arange(lmax + 1)
    pspectrum = spectralanalysis.spectrum(coeffs, unit='per_l')
    pdensity = spectralanalysis.spectrum(coeffs, unit='per_lm')

    # ---- plot spectrum ----
    fig_spectrum, ax = plt.subplots(1, 1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('degree l')
    ax.grid(True, which='both')

    ax.plot(ls[1:], pspectrum[1:], label='power per degree l')
    ax.plot(ls[1:], pdensity[1:], label='power per degree l and order m')

    ax.legend()

    fig_map.savefig('SHRtopography_mars.png')
    fig_spectrum.savefig('SHRspectrum_mars.png')
    print('mars topography and spectrum saved')

    # plt.show()

# ==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
