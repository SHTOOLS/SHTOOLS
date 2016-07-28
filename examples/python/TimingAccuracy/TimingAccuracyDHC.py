#!/usr/bin/env python
"""
This script is a python version of TimingAccuracyDHC. We use numpy functions to
simplify the creation of random coefficients.
"""
from __future__ import absolute_import, division, print_function

import os
import sys
import time
import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), "../../.."))
import pyshtools.shtools as shtools

#==== MAIN FUNCTION ====
def main():
    TimingAccuracyDHC(2)

#==== TEST FUNCTIONS ====
def TimingAccuracyDHC(sampling=1):
    #---- input parameters ----
    maxdeg = 2800
    ls = np.arange(maxdeg + 1)
    beta = -1.5
    print('Driscoll-Healy (complex), sampling =', sampling)
    
    #---- create mask to filter out m<=l ----
    mask = np.zeros((2, maxdeg + 1, maxdeg + 1), dtype=np.bool)
    mask[0, 0, 0] = True
    for l in ls:
        mask[:, l, :l + 1] = True
    mask[1, :, 0] = False

    #---- create Gaussian powerlaw coefficients ----
    print('creating {:d} random coefficients'.format(2 * (maxdeg + 1) * (maxdeg + 1)))
    cilm = np.zeros((2, (maxdeg + 1), (maxdeg + 1)), dtype=np.complex)
    cilm.imag = np.random.normal(loc=0., scale=1.,
                                 size=(2, maxdeg + 1, maxdeg + 1))
    cilm.real = np.random.normal(loc=0., scale=1.,
                                 size=(2, maxdeg + 1, maxdeg + 1))

    cilm[:, 1:, :] *= np.sqrt((ls[1:]**beta) / (2. * ls[1:] + 1.))[None, :, None]

    #---- time spherical harmonics transform for lmax set to increasing powers of 2 ----
    lmax = 2
    print('lmax    maxerror    rms         tinverse    tforward')
    while lmax <= maxdeg:
        # trim coefficients to lmax
        cilm_trim = cilm[:, :lmax + 1, :lmax + 1]
        mask_trim = mask[:, :lmax + 1, :lmax + 1]

        #synthesis / inverse
        tstart = time.time()
        grid = shtools.MakeGridDHC(cilm_trim, sampling=sampling)
        tend = time.time()
        tinverse = tend - tstart

        #analysis / forward
        tstart = time.time()
        cilm2_trim = shtools.SHExpandDHC(grid, sampling=sampling)
        tend = time.time()
        tforward = tend - tstart

        # compute error
        err = np.abs(cilm_trim[mask_trim] - cilm2_trim[mask_trim]) / np.abs(cilm_trim[mask_trim])
        maxerr = err.max()
        rmserr = np.mean(err**2)

        print('{:4d}    {:1.2e}    {:1.2e}    {:1.1e}s    {:1.1e}s'.format(
            lmax, maxerr, rmserr, tinverse, tforward))
        lmax = lmax * 2

#==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
