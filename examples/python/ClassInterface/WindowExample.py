#!/usr/bin/env python
"""
This script tests the python class interface
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
from pyshtools import SHWindow

# set shtools plot style:
sys.path.append(os.path.join(os.path.dirname(__file__), "../Common"))
from FigStyle import style_shtools
mpl.rcParams.update(style_shtools)

#==== MAIN FUNCTION ====

def main():
    example1()
    example2()

#==== EXAMPLES ====
def example1():
    #generate cap window
    lmax  = 20
    nwins = 20
    theta = 25.
    cap = SHWindow.from_cap(lmax,nwins,theta)
    cap.info()
    cap.plot(20,show=False,fname='cap_tapers.png')
    cap.plot_couplingmatrix(30,5,show=False,fname='cap_coupling.png')

#==== EXAMPLES ====
def example2():
    #generate cap window
    lmax  = 15
    nwins = 15

    topo = np.loadtxt('topo.dat.gz')
    dh_mask = topo > 0.
    print(dh_mask.shape)
    region = SHWindow.from_mask(lmax, nwins, dh_mask, sampling=2)
    region.info()
    region.plot(nwins,show=False,fname='continent_tapers.png')
    region.plot_couplingmatrix(30,5,show=False,fname='continent_coupling.png')

#==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
