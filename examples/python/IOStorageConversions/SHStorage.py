#!/usr/bin/env python
"""
This script tests the conversions between real and complex spherical harmonics
coefficients
"""
from __future__ import absolute_import, division, print_function

# standard imports:
import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# import shtools:
sys.path.append(os.path.join(os.path.dirname(__file__), "../../.."))
import pyshtools.shtools as shtools

# set shtools plot style:
sys.path.append(os.path.join(os.path.dirname(__file__), "../Common"))
from FigStyle import style_shtools
mpl.rcParams.update(style_shtools)


def main():
    test_SHStorage()


def test_SHStorage():
    #---- input parameters ----
    lmax = 1
    mask = np.zeros((2, lmax + 1, lmax + 1), dtype=np.bool)
    for l in np.arange(lmax + 1):
        mask[:, l, :l + 1] = True
    mask[1, :, 0] = False
    #---- creating random coefficients ----
    coeffs = np.random.normal(size=(2, lmax + 1, lmax + 1))
    coeffs[np.invert(mask)] = 0.

    print('\n---- testing SHCilmToCindex and SHCindexToCilm ----')
    coeffs_indexed = shtools.SHCilmToCindex(coeffs)
    coeffs_recomp = shtools.SHCindexToCilm(coeffs_indexed)
    print('input coeffs (l={:d}):'.format(lmax))
    print(coeffs)
    print('indexed coeffs:')
    print(coeffs_indexed)
    print('recomputed coeffs:')
    print(coeffs_recomp)

    print('\n---- testing SHCilmToVector and SHVectorToCilm ----')
    coeffs_indexed = shtools.SHCilmToVector(coeffs)
    coeffs_recomp = shtools.SHVectorToCilm(coeffs_indexed)
    print('\ninput coeffs (l={:d}):'.format(lmax))
    print(coeffs)
    print('\nindexed coeffs:')
    print(coeffs_indexed)
    print('\nrecomputed coeffs:')
    print(coeffs_recomp)

#==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
