#!/usr/bin/env python
"""
This script should ultimately test all the wrapped python functions.
"""

import numpy as np
import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), "../.."))
import pyshtools as shtools
import matplotlib.pyplot as plt

print 'parameters are set to:'

#==== MAIN FUNCTION ====
def main():
    test_shtransforms()
    plt.show()

#==== TEST FUNCTIONS ====
def test_shtransforms():
    lmax = 2
    #create random coefficients
    coeffs    = np.random.normal( loc=0., scale=1.,size=2*lmax*lmax ).reshape(2,lmax,lmax)
    print coeffs
    grid,nlat = shtools.MakeGridDH(coeffs,norm=1,csphase=1,sampling=2)
    #coeffs2   = shtools.SHExpandDH(grid,sampling=2)

    print grid.shape
    print nlat
    plt.imshow(grid)

#==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
