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
    coeffs1, lmax = shtools.SHRead(infile,719)
    coeffs1 = coeffs1[:,:lmax+1,:lmax+1]

    #--- convert to complex coefficients ---
    coeffs2 = np.empty( (2,lmax+1,lmax+1), dtype=np.complex )
    coeffs2_buf = shtools.SHrtoc(coeffs1,convention=1,switchcs=0)
    coeffs2[0,:,:].real = coeffs2_buf[0,:,:]
    coeffs2[0,:,:].imag = coeffs2_buf[1,:,:]
    coeffs2[1] = coeffs2[0].conjugate()*((-1)**np.arange(lmax+1)).reshape(1,1,lmax+1)

    #--- plot grid ---
    grid1,nlat = shtools.MakeGridDH(coeffs1,lmax,csphase=-1)
    grid2,nlat = shtools.MakeGridDHC(coeffs2,lmax,csphase=-1)

    plt.figure()
    plt.imshow(grid1)

    plt.figure()
    plt.imshow(np.real(grid2))

    plt.show()

#==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()

