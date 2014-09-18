#!/usr/bin/env python
"""
This script should ultimately test all the wrapped python functions.
"""

import numpy as np
import os, sys
sys.path = [os.path.join(os.path.dirname(__file__), "../../lib")] + sys.path
import pyshtools as shtools
import matplotlib.pyplot as plt

print shtools.params.norm

#==== MAIN FUNCTION ====
def main():
    test_shtransforms()

#==== TEST FUNCTIONS ====
def test_shtransforms():
    lmax = 200
    coeffs    = np.random.normal( loc=0., scale=1.,size=2*lmax*lmax ).reshape(2,lmax,lmax)

    grid      = shtools.pymakegriddh(coeffs,sampling=2)
    coeffs2   = shtools.pyshexpanddh(grid,sampling=2)
    print coeffs[1,3,2]
    print coeffs2[1,3,2]

    plt.imshow(grid)
    plt.show()

#==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
