#!/usr/bin/env python
"""
This script tests and plots all Geodesy normalized Legendre functions.
Parameters can be changed in the main function.
"""
from __future__ import absolute_import, division, print_function

import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), "../../.."))
import pyshtools
from pyshtools import legendre

pyshtools.utils.figstyle()


# ==== MAIN FUNCTION ====

def main():
    # --- input parameters (change here) ---
    normalization = ''
    # normalization should be one of ['Bar','Schmidt','ON','']
    lmax = 40  # maximum degree
    mplot = min(lmax, 10)
    # maximum plotting order (all degrees are plotted)

    # --- run tests ---
    test_associatedlegendre(lmax, mplot, normalization)
    test_legendre(lmax, normalization)


# ==== TEST LEGENDRE FUNCTIONS ====

def test_legendre(lmax, normalization):
    print('testing Pl{0} and Pl{0}_d1...'.format(normalization))
    # --- import function from shtools ---
    if normalization == '':
        Pl = legendre.PLegendre
        Pl_d1 = legendre.PLegendre_d1
    else:
        Pl = getattr(legendre, 'Pl' + normalization)
        Pl_d1 = getattr(legendre, 'Pl' + normalization + '_d1')

    # --- derived parameters ---
    npoints = 5 * lmax
    ls = np.arange(lmax)
    cost = np.cos(np.linspace(np.pi / npoints, np.pi - np.pi / npoints,
                              npoints))

    # --- create arrays to store Legendre functions of degrees l and orders
    # --- m at all points cost
    Pl1 = np.zeros((npoints, lmax))
    Pl2 = np.zeros((npoints, lmax))
    dPl2 = np.zeros((npoints, lmax))

    for iz, z in enumerate(cost):
        Pl1_buf = Pl(lmax, z)
        Pl2_buf, dPl2_buf = Pl_d1(lmax, z)
        for l in ls:
            Pl1[iz, l] = Pl1_buf[l]
            Pl2[iz, l] = Pl2_buf[l]
            dPl2[iz, l] = dPl2_buf[l]

    # ---- check if both subroutines computed the same Legendre functions ---
    if not np.allclose(Pl1, Pl2, rtol=1e-10):
        raise Exception('Legendre functions from PlmON and PlmON_d1 are ' +
                        'different (rtol>1e-10)')

    # ---- plot the legendre functions and derivatives up to maximum
    # ---- order mplot
    fig, ax = plt.subplots(1, 2, sharey=True, figsize=(15, 6))
    fig.suptitle('orthonormalized Legendre functions (col1) and ' +
                 'derivatives (col2)')
    ax[0].imshow(Pl1[:, :], extent=(0., lmax, 0., np.pi), aspect='auto')
    ax[1].imshow(dPl2[:, :], extent=(0., lmax, 0., np.pi), aspect='auto')
    ax[1].set_xlabel('l')
    fig.savefig('legendre.png')


# ==== TEST ASSOCIATED LEGENDRE FUNCTIONS ====

def test_associatedlegendre(lmax, mplot, normalization):
    print('testing Plm{0} and Plm{0}_d1...'.format(normalization))
    # --- import function from shtools ---
    if normalization == '':
        Plm = legendre.PLegendreA
        Plm_d1 = legendre.PLegendreA_d1
    else:
        Plm = getattr(legendre, 'Plm' + normalization)
        Plm_d1 = getattr(legendre, 'Plm' + normalization + '_d1')

    # --- derived parameters ---
    npoints = 5 * lmax
    ls = np.arange(lmax)
    cost = np.cos(np.linspace(np.pi / npoints, np.pi - np.pi / npoints,
                              npoints))

    # --- create arrays to store Legendre functions of degrees l and orders
    # ----m at all points cost
    Plm1 = np.zeros((npoints, lmax, lmax))
    Plm2 = np.zeros((npoints, lmax, lmax))
    dPlm2 = np.zeros((npoints, lmax, lmax))

    for iz, z in enumerate(cost):
        Plm1_buf = Plm(lmax, z)
        Plm2_buf, dPlm2_buf = Plm_d1(lmax, z)
        for l in ls:
            for m in np.arange(l):
                ind = legendre.PlmIndex(l, m) - 1  # Fortran indexing
                Plm1[iz, l, m] = Plm1_buf[ind]
                Plm2[iz, l, m] = Plm2_buf[ind]
                dPlm2[iz, l, m] = dPlm2_buf[ind]

    # ---- check if both subroutines computed the same Legendre functions ---
    if not np.allclose(Plm1_buf, Plm2_buf, rtol=1e-10):
        raise Exception('Legendre functions from PlmON and PlmON_d1 are ' +
                        'different (rtol>1e-10)')

    # ---- plot the legendre functions and derivatives up to maximum
    # --- order mplot
    fig, ax = plt.subplots(2, mplot, sharey=True, sharex=True, figsize=(15, 6))
    fig.suptitle('orthonormalized associated Legendre functions (row1) and ' +
                 'derivatives (row2)')
    for m in range(mplot):
        ax[0, m].imshow(Plm1[:, :, m], extent=(0., lmax, 0., np.pi),
                        aspect='auto')
        ax[0, m].set_title('m=%d' % m)
        ax[1, m].imshow(dPlm2[:, :, m], extent=(0., lmax, 0., np.pi),
                        aspect='auto')
        ax[1, m].set_xlabel('l')
    fig.savefig('associatedlegendre.png')

# ==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
