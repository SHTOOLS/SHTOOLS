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
    Spherical Harmonics Coefficient super class. Coefficients can be initialized
    using one of the constructor methods:

    >> SHCoefficients.from_array( np.zeros(2*(lmax+1)*(lmax+1)) )
    >> SHCoefficients.from_random( np.exp(-ls**2) )
    >> SHCoefficients.from_file( 'fname.dat' )
    """
    def __init__(self):
        print('use one of the following methods to initialize sh-coefficients:\n\n'+
              '>> SHCoefficients.from_array(...)\n'+
              '>> SHCoefficients.from_random(...)\n'+
              '>> SHCoefficients.from_file(...)')

    #---- factory methods:
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
                    coeffs  = np.random.normal(size=2*(lmax+1)*(lmax+1) ).reshape(2,lmax+1,lmax+1)
                    coeffs *= np.sqrt(power.reshape(1,lmax+1,1))
                elif kind=='complex':
                    coeffs = np.random.normal( loc=0., scale=1.,size=2*(lmax+1)*(lmax+1) ) + \
                           1j*np.random.normal( loc=0., scale=1.,size=2*(lmax+1)*(lmax+1) )
                    coeffs = coeffs1.reshape(2,lmax+1,lmax+1)
                    coeffs *= np.sqrt(power.reshape(1,lmax+1,1))
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

    #---- extracting data:
    def get_powerperdegree(self):
        return self._powerperdegree()

    def get_degrees(self):
        return np.arange(self.lmax+1)

    #---- rotation:
    def rotate(self,alpha,beta,gamma,degrees=True):
        if degrees:
            angles = np.radians([alpha,beta,gamma])
        else:
            angles = np.array([alpha,beta,gamma])
        self._rotate(angles)

    #---- expansion:
    def expand(self,kind='DH1'):
        if kind == 'DH1':
            grid = self._expandDH(sampling=1)
        else:
            raise NotImplementedError('grid type {:s} not found/implemented'.format(kind))
        return grid

    #---- plotting routines:
    def plot_powerperdegree(self,loglog=True,show=True):
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
        if show: plt.show()

    def plot_bandpower(self,bandwidth=2,show=True):
        ls    = self.get_degrees()
        power = self.get_powerperdegree()*ls*np.log(bandwidth)

        fig,ax = plt.subplots(1,1)
        ax.set_xlabel('degree l')
        ax.set_ylabel('bandpower')
        ax.set_xscale('log',basex=bandwidth)
        ax.set_yscale('log',basey=bandwidth)
        ax.grid(True,which='both')
        ax.plot(ls[1:],power[1:],label='power per degree l')
        if show: plt.show()

#---- implementation for real spherical harmonics ----
class SHRealCoeffients(SHCoefficients):
    """
    Real Spherical Harmonics Coefficients class.
    """
    @staticmethod
    def istype(kind):
        return kind=='real'

    def __init__(self,coeffs):
        #---- create mask to filter out m<=l ----
        lmax = coeffs.shape[1]-1
        mask = np.zeros(2*(lmax+1)*(lmax+1),dtype=np.bool).reshape(2,lmax+1,lmax+1)
        mask[0,0,0] = True
        for l in np.arange(lmax+1):
            mask[:,l,:l+1] = True
        mask[1,:,0] = False

        self.lmax   = lmax
        self.coeffs = np.copy(coeffs)
        self.coeffs[np.invert(mask)] = 0.

    def _powerperdegree(self):
        return SHPowerSpectrum(self.coeffs)

    def _rotate(self,angles,dj_matrix=None):
        if dj_matrix is None:
            dj_matrix = djpi2(self.lmax+1)
        self.coeffs = SHRotateRealCoef(self.coeffs, angles, dj_matrix)

    def _expandDH(self,sampling):
        data = MakeGridDH(self.coeffs,sampling=sampling)
        grid = SphericalGrid.from_array(data)
        return grid

    def _expandGLQ(self):
        zeros, weights = PreCompute(lmax)
        return MakeGridGLQ(cilm_trim,zeros)

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
class SphericalGrid(object):
    """
    Spherical Grid Class that can deal with spatial data on the sphere that is
    defined on different grids. Can be constructed from:

    >> SphericalGrid.from_array(...)
    >> SphericalGrid.from_file(...)
    """
    def __init__():
        print('use one of the following methods to initialize the grid:\n\n'+
              '>> SphericalGrid.from_array(...)\n'+
              '>> SphericalGrid.from_file(...)')

    #---- constructors:
    @classmethod
    def from_array(self, array, kind='DH'):
        for cls in self.__subclasses__():
            if cls.istype(kind):
                return cls(array)

    #---- plotting routines:
    def plot_rawdata(self,show=True):
        self._plot_rawdata()
        if show: plt.show()

#---- implementation of the Driscoll and Healy Grid class ----
class DHGrid(SphericalGrid):
    """
    Driscoll and Healy Grid (publication?)
    """
    @staticmethod
    def istype(kind):
        return kind == 'DH'

    def __init__(self, array):
        nlat,nlon = array.shape
        if nlat == 2*nlon:
            self.sampling = 2
        if nlat == nlon:
            self.sampling = 1
        else:
            raise ValueError('input array with shape (nlat={:d},nlon={:d})\n'+
                             'it needs nlat=nlon or nlat=2*nlon'.format(nlat,nlon))
        self.data = array

    def _plot_rawdata(self):
        fig,ax = plt.subplots(1,1)
        ax.imshow(self.data,origin='top',extent=(0.,360.,-90.,90.))
        ax.set_title('Driscoll Healy Grid')
        ax.set_xlabel('longitude')
        ax.set_ylabel('latitude')
        fig.tight_layout(pad=0.5)

#---- implementation of the Gauss-Legendre Grid class ----
class GLQGrid(SphericalGrid):
    """
    Gauss Legendre Class
    """
    @staticmethod
    def istype(kind):
        return kind == 'DH'

    def __init__(self, array):
        nlat,nlon = array.shape
        if nlat == 2*nlon:
            self.sampling = 2
        if nlat == nlon:
            self.sampling = 1
        else:
            raise ValueError('input array with shape (nlat={:d},nlon={:d})\n'+
                             'it needs nlat=nlon or nlat=2*nlon'.format(nlat,nlon))
        self.data = array

    def _plot_rawdata(self):
        fig,ax = plt.subplots(1,1)
        ax.imshow(self.data,origin='top',extent=(0.,360.,-90.,90.))
        ax.set_title('Driscoll Healy Grid')
        ax.set_xlabel('longitude')
        ax.set_ylabel('latitude')
        fig.tight_layout(pad=0.5)
