#!/usr/bin/env python
"""
This script tests the conversions between real and complex spherical harmonics
coefficients
"""
from __future__ import absolute_import, division, print_function

import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), "../../.."))
from pyshtools import shio
from pyshtools import expand

sys.path.append(os.path.join(os.path.dirname(__file__), "../Common"))
from FigStyle import style_shtools

# set shtools plot style:
mpl.rcParams.update(style_shtools)


def main():
    test_SHConversions()
    example()


def test_SHConversions():
    print('---- testing SHrtoc and SHctor ----')
    lmax = 10
    coeffs1 = np.random.normal(loc=0., scale=1., size=(2, lmax + 1, lmax + 1))
    mask = np.zeros((2, lmax + 1, lmax + 1), dtype=np.bool)
    for l in np.arange(lmax + 1):
        mask[:, l, :l + 1] = True
    mask[1, :, 0] = False
    coeffs2 = shio.SHrtoc(coeffs1)
    coeffs3 = shio.SHctor(coeffs2)
    error = np.sqrt(np.sum((coeffs3[mask] - coeffs1[mask])**2))
    print('error after real to complex to real conversion: ', error)


def example():
    print('---- SHrtoc example ----')
    # --- input data filename ---
    infile = '../../ExampleDataFiles/MarsTopo719.shape'
    coeffs1, lmax = shio.shread(infile)
    coeffs1 = coeffs1[:, :lmax + 1, :lmax + 1]

    # --- convert to complex coefficients, fill negative order coefficients ---
    coeffs2 = np.empty((2, lmax + 1, lmax + 1), dtype=np.complex)
    coeffs2_buf = shio.SHrtoc(coeffs1, convention=1, switchcs=0)
    coeffs2[0, :, :].real = coeffs2_buf[0, :, :]
    coeffs2[0, :, :].imag = coeffs2_buf[1, :, :]
    coeffs2[1] = (coeffs2[0].conjugate() *
                  ((-1)**np.arange(lmax + 1))[np.newaxis, :])

    # --- compute and plot grid ---
    grid1 = expand.MakeGridDH(coeffs1, lmax, csphase=-1)
    grid2 = expand.MakeGridDHC(coeffs2, lmax, csphase=-1)

    gridmin = min(grid1.min(), grid2.real.min())
    gridmax = max(grid1.max(), grid2.real.max())
    norm = plt.Normalize(gridmin, gridmax)

    fig, axes = plt.subplots(1, 2)
    im1 = axes[0].imshow(grid1, norm=norm)
    axes[0].set_title('from real coefficients')

    im2 = axes[1].imshow(grid2.real, norm=norm)
    axes[1].set_title('from complex coefficients')

    fig.tight_layout(pad=1)
    fig.savefig('topography_mars.png')
    print('mars topography plotted and saved to file')

    # plt.show()

# ==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
