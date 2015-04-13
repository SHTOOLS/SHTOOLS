#!/usr/bin/env python
"""
This script tests the python class interface
"""

from __future__ import division
from __future__ import print_function

# standard imports:
import os
import sys
import numpy as np
import matplotlib as mpl

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

#==== EXAMPLES ====
def example1():
    #generate cap window
    lmax  = 20
    nwins = 20
    theta = 30.
    cap = SHWindow.cap(lmax,nwins,theta)
    cap.info()
    cap.plot(20)

#==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
