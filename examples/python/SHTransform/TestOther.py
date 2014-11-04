#!/usr/bin/env python
"""
This script is a python version of TimingAccuracy. We use some numpy functions
to simplify the creation of random coefficients.
"""

import os, sys, time
import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), "../../.."))
import pyshtools as shtools

#==== MAIN FUNCTION ====
def main():
    test_DHaj()
    test_GLQGridCoord()
    test_PreGLQ()
    test_NGLQ()

#==== TEST FUNCTIONS ====
def test_DHaj():
	print '--- Testing DHaj ---'
	n=20
	aj = shtools.DHaj(n)
	print 'n = ', n
	print 'Driscoll-Healy weights = ', aj
def test_GLQGridCoord():
	print '--- Testing GLQGridCoord ---'
	lmax = 10
	print 'lmax = ', lmax
	lat, lon, nlat, nlong = shtools.GLQGridCoord(lmax)
	print '(nlat,nlong) = ', (nlat,nlong)
	print 'latitude = ', lat
	print 'longitude = ', lon	
def test_PreGLQ():
	print '--- Testing PreGLQ ---'
	lower, upper, n = (-1.0, 1.0, 20)
	zero, w = shtools.PreGLQ(lower,upper,n)
	print 'lower limit = ', lower
	print 'upper limit = ', upper
	print 'N = ', n
	print 'zeros = ', zero
	print 'weights = ', w
def test_NGLQ():
	print '--- Testing NGLQ, NGLQSH, NGLQSHN ---'
	lmax = 20
	n = 8
	print 'lmax = ', lmax
	print 'n = ', n
	print 'NGLQ = ', shtools.NGLQ(lmax)
	print 'NGLQSH = ', shtools.NGLQSH(lmax)
	print 'NGLQNSH = ', shtools.NGLQSHN(lmax, n)
#==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()

