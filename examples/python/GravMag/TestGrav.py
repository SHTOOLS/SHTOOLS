#!/usr/bin/env python
"""
This script tests the gravity and magnetics routines. 
"""

#standard imports:
import os, sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

#import shtools:
sys.path.append(os.path.join(os.path.dirname(__file__), "../../.."))
import pyshtools as shtools

#set shtools plot style:
sys.path.append(os.path.join(os.path.dirname(__file__), "../Common"))
from FigStyle import style_shtools
mpl.rcParams.update(style_shtools)

#==== MAIN FUNCTION ====
def main():
	TestNormalGravity()
	TestMakeGeoidGrid()
	TestGravGrad()

#==== TEST FUNCTIONS ====
def TestNormalGravity():
	gm = shtools.constants.gm_mars
	omega = shtools.constants.omega_mars
	a = shtools.constants.a_mars
	b = shtools.constants.b_mars
	lat = np.arange(-90.,90.,1.)
	ng = np.array([shtools.NormalGravity(x,gm,omega,a,b) for x in lat])
	fig = plt.figure()
	plt.plot(lat, ng, '-')
	fig.savefig('Mars_normalgravity.png')
		
def TestMakeGeoidGrid():
	infile = '../../ExampleDataFiles/jgmro_110b_sha.tab'
	header = np.zeros(2,dtype=float)
	clm,lmax, header = shtools.SHReadH(infile,110,header)
	gm = header[1] * 1.e9
	r0 = header[0] * 1.e3
	clm[0,0,0] = 1.0
	geoid = shtools.MakeGeoidGridDH(clm,r0,gm,shtools.constants.w0_mars,a=shtools.constants.a_mars,f=shtools.constants.f_mars,omega=shtools.constants.omega_mars)
	geoid = geoid / 1.e3 # convert to meters
	fig_map = plt.figure()
	plt.imshow(geoid)
	fig_map.savefig('MarsGeoid.png')

def TestGravGrad():
	#---- input parameters ----
	lmax = 100
	clm = np.zeros((2,lmax+1,lmax+1),dtype=float)
	clm[0,2,2] = 1.0
	gm = 1.0
	r0 = 1.0
	a = 1.0
	f = 0.0
	
	vxx,vyy,vzz,vxy,vxz,vyz = shtools.MakeGravGradGridDH(clm,gm,r0,a,f)
	
	print "Maximum Trace(Vxx+Vyy+Vzz) = ", np.max(vxx+vyy+vzz)
	print "Minimum Trace(Vxx+Vyy+Vzz) = ", np.min(vxx+vyy+vzz)

	fig, axes = plt.subplots(2,3)
	fig.suptitle("Gravity gradient tensor", fontsize=10)
		
	for num, vv, s in ((0,vxx,"$V_{xx}$"), (1,vyy,"$V_{yy}$"), (2,vzz,"$V_{zz}$"), (3,vxy,"$V_{xy}$"), \
		(4,vxz,"$V_{xz}$"), (5,vyz,"$V_{yz}$")):
		axes.flat[num].imshow(vv, vmin=-5,vmax=5)
		axes.flat[num].set_title(s)
		axes.flat[num].set_xticks(())
		axes.flat[num].set_yticks(())
	
	fig.savefig('GravGrad_C22.png')
			
#==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()

