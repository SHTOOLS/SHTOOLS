"""
This script tests the different Spherical Harmonics Transforms on the Mars
topography data set
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pyshtools as pysh

pysh.utils.figstyle()


# ==== MAIN FUNCTION ====
def main():
    example()


# ==== PLOT POWER SPECTRA ====

def example():
    """
    example that plots the power spectrum of Mars topography data
    """
    # --- input data filename ---
    topofile = 'MarsTopo719.shape'
    if len(sys.argv) > 1:
        topofile = os.path.join(sys.argv[1], topofile)
    else:
        topofile = os.path.join('../../ExampleDataFiles', topofile)

    coeffs, lmax = pysh.shio.shread(topofile)

    # --- plot grid ---
    grid = pysh.expand.MakeGridDH(coeffs, csphase=-1)
    fig_map = plt.figure()
    plt.imshow(grid)

    # ---- compute spectrum ----
    ls = np.arange(lmax + 1)
    pspectrum = pysh.spectralanalysis.spectrum(coeffs, unit='per_l')
    pdensity = pysh.spectralanalysis.spectrum(coeffs, unit='per_lm')

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
