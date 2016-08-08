#!/usr/bin/env python
"""
This script tests functions to compute different Power Spectra from
complex SH coefficients
"""
from __future__ import absolute_import, division, print_function

import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), "../../.."))
from pyshtools import shtools

# set shtools plot style:
sys.path.append(os.path.join(os.path.dirname(__file__), "../Common"))
from FigStyle import style_shtools
mpl.rcParams.update(style_shtools)


# ==== MAIN FUNCTION ====

def main():
    test_ComplexSpectralAnalysis()
    example()


def test_ComplexSpectralAnalysis():
    # ---- input parameters ----
    lmax = 5
    ls = np.arange(lmax + 1)
    mask = np.zeros((2, lmax + 1, lmax + 1), dtype=np.bool)
    for l in np.arange(lmax + 1):
        mask[:, l, :l + 1] = True
    mask[1, :, 0] = False

    print('creating {:d} random coefficients'.format(2 * (lmax + 1) *
                                                     (lmax + 1)))

    print('\n---- testing SHPower/DensityLC, SHPowerSpectrum/DensityC ----')
    print('generating normal distributed complex coefficients with ' +
          'variance 1...')
    coeffs1 = (np.random.normal(loc=0., scale=1.,
                                size=(2, lmax + 1, lmax + 1)) +
               1j * np.random.normal(loc=0., scale=1.,
                                     size=(2, lmax + 1, lmax + 1)))
    coeffs1[np.invert(mask)] = 0.

    spec1 = np.array([shtools.SHPowerLC(coeffs1, l) for l in ls])
    spec2 = shtools.SHPowerSpectrumC(coeffs1)
    print('tot power computed with SHPowerL={:2.2f}'.format(np.sum(spec1)))
    print('tot power computed with SHPowerSpectrum={:2.2f}'
          .format(np.sum(spec2)))

    spec1 = np.array([shtools.SHPowerDensityLC(coeffs1, l) for l in ls])
    spec2 = shtools.SHPowerSpectrumDensityC(coeffs1)
    print('tot power computed with SHPowerDensityL={:2.2f}'
          .format(np.sum(spec1 * (2 * ls + 1))))
    print('tot power computed with SHPowerSpectrumDensity={:2.2f}'
          .format(np.sum(spec2 * (2 * ls + 1))))

    print('\n---- testing SHCrossPower/DensityLC, ' +
          'SHCrossCrossPowerSpectrum/DensityC ----')
    print('generating two sets of normal distributed complex coefficients ' +
          'with variance 1...')
    coeffs2 = (np.random.normal(loc=0., scale=1.,
                                size=(2, lmax + 1, lmax + 1)) +
               1j * np.random.normal(loc=0., scale=1.,
                                     size=(2, lmax + 1, lmax + 1)))
    coeffs2[np.invert(mask)] = 0.

    spec1 = np.array([shtools.SHCrossPowerLC(coeffs1, coeffs2, l) for l in ls])
    spec2 = shtools.SHCrossPowerSpectrumC(coeffs1, coeffs2)
    print('tot cpower computed with SHCrossPowerL={:2.2f}'
          .format(np.sum(spec1)))
    print('tot cpower computed with SHCrossPowerSpectrum={:2.2f}'
          .format(np.sum(spec2)))

    spec1 = np.array([shtools.SHCrossPowerDensityLC(coeffs1, coeffs2, l)
                      for l in ls])
    spec2 = shtools.SHCrossPowerSpectrumDensityC(coeffs1, coeffs2)
    print('tot cpower computed with SHCrossPowerDensityL={:2.2f}'
          .format(np.sum(spec1 * (2 * ls + 1))))
    print('tot cpower computed with SHCrossPowerSpectrumDensity={:2.2f}'
          .format(np.sum(spec2 * (2 * ls + 1))))


# ==== PLOT POWER SPECTRA ====

def example():
    """
    example that plots the power spectrum of Mars topography data
    """
    # --- input data filename ---
    infile = '../../ExampleDataFiles/MarsTopo719.shape'
    coeffs, lmax = shtools.SHRead(infile, 719)
    lmax = coeffs.shape[1] - 1

    # --- plot grid ---
    grid = shtools.MakeGridDH(coeffs, lmax, csphase=-1)
    fig_map = plt.figure()
    plt.imshow(grid)

    # ---- compute spectrum ----
    ls = np.arange(lmax + 1)
    pspectrum = shtools.SHPowerSpectrum(coeffs)
    pdensity = shtools.SHPowerSpectrumDensity(coeffs)

    # ---- plot spectrum ----
    fig_spectrum, ax = plt.subplots(1, 1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('degree l')
    ax.grid(True, which='both')

    ax.plot(ls[1:], pspectrum[1:], label='power per degree l')
    ax.plot(ls[1:], pdensity[1:], label='power per degree l and order m')

    ax.legend()

    fig_map.savefig('SHCtopography_mars.png')
    fig_spectrum.savefig('SHCspectrum_mars.png')
    print('mars topography and spectrum saved')

    # plt.show()

# ==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
