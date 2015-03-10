"""
This file contains some classes that facilitate in particular the interactive
examination of spherical data.

Matthias Meschede and Mark Wieczorek, 2015
"""

import numpy as np
import matplotlib.pyplot as plt

from _SHTOOLS import *

#=========== COEFFICIENT CLASSES ===============================================
class SHCoefficients(object):
    """
    abstract class that defines the interface for either complex or real
    spherical harmonics coefficients
    """
    def __init__(self):
        print('use one of the following methods to initialize sh-coefficients:\n\n'+
              '>> SHCoefficients.from_array(...)\n'+
              '>> SHCoefficients.from_random(...)\n'+
              '>> SHCoefficients.from_file(...)')

    #==== factory methods: ====
    @classmethod
    def from_array(self, array, kind='real'):
        for cls in self.__subclasses__():
            if cls.istype(kind):
                return cls(array)


    @classmethod
    def from_random(self, power, kind='real'):
        lmax = len(power)-1
        for cls in self.__subclasses__():
            if cls.istype(kind):
                if   kind=='real':
                    coeffs = np.random.normal(size=2*(lmax+1)*(lmax+1) ).reshape(2,lmax+1,lmax+1)
                elif kind=='complex':
                    coeffs = np.random.normal( loc=0., scale=1.,size=2*(lmax+1)*(lmax+1) ) + \
                           1j*np.random.normal( loc=0., scale=1.,size=2*(lmax+1)*(lmax+1) )
                    coeffs = coeffs1.reshape(2,lmax+1,lmax+1)
                else:
                    raise ValueError("kind='{:s}' should be 'real' or 'complex'".format(str(kind)))
                return cls(coeffs)


    @classmethod
    def from_file(self, fname, lmax, format='shtools'):
        """
        reads coefficients from a spherical harmonics file.
        """
        coeffs,lmax = SHRead(fname,lmax)
        lmax = coeffs.shape[1]-1

        for cls in self.__subclasses__():
            if format=='shtools' and cls.istype('real'):
                return cls(coeffs)

    #---- extracting data ----
    def get_powerperdegree(self):
        return self._powerperdegree()

    def get_degrees(self):
        return np.arange(self.lmax+1)

    #---- plotting routines ----
    def plot_powerperdegree(self,loglog=True):
        power = self.get_powerperdegree()
        ls    = self.get_degrees()

        fig,ax = plt.subplots(1,1)
        ax.set_xlabel('degree l')
        ax.set_ylabel('power per degree')
        if loglog:
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.grid(True,which='both')
            ax.plot(ls[1:],power[1:],label='power per degree l')
        plt.show()

    def plot_bandpower(self,bandwidth=2):
        ls    = self.get_degrees()
        power = self.get_powerperdegree()*ls*np.log(bandwidth)

        fig,ax = plt.subplots(1,1)
        ax.set_xlabel('degree l')
        ax.set_ylabel('bandpower')
        ax.set_xscale('log',basex=bandwidth)
        ax.set_yscale('log',basey=bandwidth)
        ax.grid(True,which='both')
        ax.plot(ls[1:],power[1:],label='power per degree l')
        plt.show()

    #---- rotation ----
    def rotate():
        self._rotate(angles)

    #---- expansion ----
    def expand():
        raise NotImplementedError('not implemented by subclass %s'%s.__class__)

#---- implementation for real spherical harmonics ----
class SHRealCoeffients(SHCoefficients):
    """
    Real Spherical Harmonics Coefficients class.
    """
    @staticmethod
    def istype(kind):
        return kind=='real'

    def __init__(self,coeffs):
        self.coeffs = coeffs
        self.lmax = coeffs.shape[1]-1

    def _powerperdegree(self):
        return SHPowerSpectrum(self.coeffs)

    def _rotate(self,angles,dj_matrix=None):
        if djmatrix is None:
            dj_matrix = shtools.djpi2(self.lmax)
        self.coeffs = shtools.SHRotateRealCoef(self.coeffs,angles,dj_matrix)

    def _expandDH(self):
        return shtools.MakeGridDH(cilm_trim,sampling=sampling) 

    def _expandGLQ(self):
        zeros, weights = shtools.PreCompute(lmax)
        return shtools.MakeGridGLQ(cilm_trim,zeros)

#---- implementation for complex spherical harmonics ----
class SHComplexCoeffients(SHCoefficients):
    """
    Complex Spherical Harmonics Coefficients class.
    """
    @staticmethod
    def istype(kind):
        return kind=='complex'

    def __init__(self,coeffs):
        self.coeffs = coeffs
        self.lmax = coeffs.shape[1]-1

    def _powerperdegree():
        return SHCPowerSpectrum(self.coeffs)

    def _rotate():
        raise NotImplementedError('not implemented by subclass %s'%s.__class__)

    def _expandDH():
        raise NotImplementedError('not implemented by subclass %s'%s.__class__)

    def _expandGLQ():
        raise NotImplementedError('not implemented by subclass %s'%s.__class__)

#=========== GRID CLASSES ===============================================
