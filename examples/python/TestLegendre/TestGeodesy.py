#!/usr/bin/env python
"""
This script tests and plots all Geodesy normalized Legendre functions
"""

import os, sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), "../../.."))
import pyshtools as shtools

#==== MAIN FUNCTION ====
def main():
    #--- input parameters (change here) ---
    lmax  = 40            #maximum degree
    mplot = min(lmax,10)  #maximum plotting order (all degrees are plotted)

    #--- derived parameters ---
    npoints = 5*lmax
    ls = np.arange(lmax)
    cost  = np.cos(np.linspace(np.pi/npoints,np.pi-np.pi/npoints,npoints))

    #--- create arrays to store Legendre functions of degrees l and orders m at all points cost ---
    Plm1  = np.zeros( (npoints,lmax,lmax) )
    Plm2  = np.zeros( (npoints,lmax,lmax) )
    dPlm2 = np.zeros( (npoints,lmax,lmax) )

    for iz,z in enumerate(cost):
        Plm1_buf = shtools.PlmBar(lmax,z)
        Plm2_buf,dPlm2_buf = shtools.PlmBar_d1(lmax,z)
        for l in ls:
            for m in np.arange(l):
                ind = shtools.PlmIndex(l,m)-1 #Fortran indexing
                Plm1[iz,l,m] = Plm1_buf[ind]
                Plm2[iz,l,m] = Plm2_buf[ind]
                dPlm2[iz,l,m] = dPlm2_buf[ind]

    #---- check if both subroutines computed the same Legendre functions ---
    if not np.allclose(Plm1_buf,Plm2_buf,rtol=1e-10):
        raise Exception('Legendre functions from PlmBar and PlmBar_d1 are different (rtol>1e-10)')

    #---- plot the legendre functions and derivatives up to maximum order mmax ---
    fig,ax = plt.subplots(2,mplot,sharey=True,sharex=True)
    fig.suptitle(r'$4\pi$ normalized associated Legendre functions (row1) and derivatives (row2)')
    for m in range(mplot):
        ax[0,m].imshow(Plm1[:,:,m],extent=(0.,lmax,0.,np.pi),aspect='auto')
        ax[0,m].set_title('m=%d'%m)
        ax[1,m].imshow(dPlm2[:,:,m],extent=(0.,lmax,0.,np.pi),aspect='auto')
        ax[1,m].set_xlabel('l')
    plt.show()


#==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()

