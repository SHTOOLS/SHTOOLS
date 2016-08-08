#!/usr/bin/env python
"""
This script tests the gravity and magnetics routines. 
"""
from __future__ import absolute_import, division, print_function

import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), "../../.."))
from pyshtools import shtools

# set shtools plot style:
sys.path.append(os.path.join(os.path.dirname(__file__), "../Common"))
from FigStyle import style_shtools
mpl.rcParams.update(style_shtools)


# ==== MAIN FUNCTION ====

def main():
    TestCircle()
    TestEig()
    TestWigner()
    TestRandom()


# ==== TEST FUNCTIONS ====

def TestCircle():
    coord = shtools.MakeCircleCoord(30, 10, 30)
    lat = coord[:, 0]
    lon = coord[:, 1]
    fig = plt.figure()
    plt.axis([-180, 180, -90, 90])
    plt.plot(lon, lat, 'r-', 10, 30, 'r.')
    coord = shtools.MakeCircleCoord(-75, -45, 10)
    plt.plot(coord[:, 1], coord[:, 0], 'b-', -45, -75, 'b.')
    coord = shtools.MakeEllipseCoord(0, 45, 20, 30, 10)
    plt.plot(coord[:, 1], coord[:, 0], 'g-', 45, 0, 'g.')

    fig.savefig('Circles.png')


def TestEig():
    a = np.array([[1, -1., 0],
                  [-1, 2., -1],
                  [0, -1, 1]])
    eig, vec = shtools.EigValVecSym(a)
    eig2 = shtools.EigValSym(a)
    eigtri, vectri = shtools.EigValVecSymTri(a)
    print("Symmetric matrix")
    print(a)
    print("Eigenvalues (three different routines)")
    print(eig)
    print(eig2)
    print(eigtri)
    print("Eigenvectors (two different routines)")
    print(vec)
    print(vectri)


def TestWigner():
    w3j, jmin, jmax = shtools.Wigner3j(4, 2, 0, 0, 0)
    print("< J, 4, 2 / 0, 0, 0 >")
    print("jmin = ", jmin)
    print("jmax = ", jmax)
    print(w3j)
    w3j, jmin, jmax = shtools.Wigner3j(10, 14, -1, -4, 5)
    print("< J, 10, 14 / -1, -4, 5 >")
    print("jmin = ", jmin)
    print("jmax = ", jmax)
    print(w3j)


def TestRandom():
    seed = -1232
    first, nseed = shtools.RandomN(seed)
    list = np.array([shtools.RandomN(nseed)[0] for x in range(0, 10000)])
    fig, axes = plt.subplots(1, 2)
    axes[0].hist(list, bins=50, range=(0, 1))

    seed = -1232
    first, nseed = shtools.RandomGaussian(seed)
    list = np.array([shtools.RandomGaussian(nseed)[0]
                     for x in range(0, 10000)])
    axes[1].hist(list, bins=100)
    axes[1].axis([-5., 5, 0, 400])

    fig.savefig('Histogram.png')

# ==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
