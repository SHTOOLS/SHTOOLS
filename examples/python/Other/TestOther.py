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
from pyshtools import utils

sys.path.append(os.path.join(os.path.dirname(__file__), "../Common"))
from FigStyle import style_shtools

# set shtools plot style:
mpl.rcParams.update(style_shtools)


# ==== MAIN FUNCTION ====

def main():
    TestCircle()
    TestWigner()


# ==== TEST FUNCTIONS ====

def TestCircle():
    coord = utils.MakeCircleCoord(30, 10, 30)
    lat = coord[:, 0]
    lon = coord[:, 1]
    fig = plt.figure()
    plt.axis([-180, 180, -90, 90])
    plt.plot(lon, lat, 'r-', 10, 30, 'r.')
    coord = utils.MakeCircleCoord(-75, -45, 10)
    plt.plot(coord[:, 1], coord[:, 0], 'b-', -45, -75, 'b.')
    coord = utils.MakeEllipseCoord(0, 45, 20, 30, 10)
    plt.plot(coord[:, 1], coord[:, 0], 'g-', 45, 0, 'g.')

    fig.savefig('Circles.png')


def TestWigner():
    w3j, jmin, jmax = utils.Wigner3j(4, 2, 0, 0, 0)
    print("< J, 4, 2 / 0, 0, 0 >")
    print("jmin = ", jmin)
    print("jmax = ", jmax)
    print(w3j)
    w3j, jmin, jmax = utils.Wigner3j(10, 14, -1, -4, 5)
    print("< J, 10, 14 / -1, -4, 5 >")
    print("jmin = ", jmin)
    print("jmax = ", jmax)
    print(w3j)


# ==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
