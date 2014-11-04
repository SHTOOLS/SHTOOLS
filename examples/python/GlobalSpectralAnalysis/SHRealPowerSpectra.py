#!/usr/bin/env python
"""
This script tests the different Spherical Harmonics Transforms on the Mars topography data set
"""

import os, sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), "../../.."))
import pyshtools as shtools

#==== MAIN FUNCTION ====
def main():
    #--- input data filename ---
    infile = '../../ExampleDataFiles/MarsTopo719.shape'
    coeffs,lmax = shtools.SHRead(infile,719)

    #--- plot grid ---
    grid,nlat = shtools.MakeGridDH(coeffs,lmax,csphase=-1)
    plt.figure()
    plt.imshow(grid)

    #--- compute and plot power spectrum ---
    compute_plot_power(coeffs)

    plt.show()

#==== PLOT POWER SPECTRA ====
def compute_plot_power(coeffs):
    ls = np.arange(coeffs.shape[1])
    pspectrum = shtools.SHPowerSpectrum(coeffs)
    pdensity  = shtools.SHPowerSpectrumDensity(coeffs)

    fig,ax = plt.subplots(1,1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('degree l')
    ax.grid(True,which='both')

    ax.plot(ls[1:],pspectrum[1:],label='power per degree l')
    ax.plot(ls[1:],pdensity[1:],label ='power per degree l and order m')

    ax.legend()

#==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()

