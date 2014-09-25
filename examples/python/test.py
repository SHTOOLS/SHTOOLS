#!/usr/bin/env python
"""
This script should ultimately test all the wrapped python functions.
"""

import numpy as np
import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), "../../lib"))
import pyshtools as shtools
import matplotlib.pyplot as plt

print 'parameters are set to:'
print 'norm:',    shtools.params.norm
print 'csphase:', shtools.params.csphase

#==== MAIN FUNCTION ====
def main():
    test_legendre()
    test_shtransforms()
    plt.show()

#==== TEST FUNCTIONS ====
def test_legendre():
    """
    tests the different Legendre polynomial functions
    """
    shtools.pyplmbar()

def test_shtransforms():
    lmax = 200
    #create random coefficients
    coeffs    = np.random.normal( loc=0., scale=1.,size=2*lmax*lmax ).reshape(2,lmax,lmax)
    grid      = shtools.pymakegriddh(coeffs,sampling=2)
    coeffs2   = shtools.pyshexpanddh(grid,sampling=2)
    plt.imshow(grid)

#==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
