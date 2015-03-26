"""
This file contains some classes that facilitate in particular the interactive
examination of spherical data.

Matthias Meschede and Mark Wieczorek, 2015
"""

import numpy as np
import matplotlib.pyplot as plt

from _SHTOOLS import *







#===============================================================================
#=========== COEFFICIENT CLASSES ===============================================
#===============================================================================



class SHCoeffs(object):
    """
    Spherical Harmonics Coefficient class. Coefficients can be initialized
    using one of the constructor methods:

    >> SHCoeffs.from_array( np.zeros(2*(lmax+1)*(lmax+1)) )
    >> SHCoeffs.from_random( np.exp(-ls**2) )
    >> SHCoeffs.from_file( 'fname.dat' )
    """
    def __init__(self):
        print('use one of the following methods to initialize sh-coefficients:\n\n'+
              '>> SHCoeffs.from_array(...)\n'+
              '>> SHCoeffs.from_random(...)\n'+
              '>> SHCoeffs.from_file(...)')

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

    #---- extracting data ----
    def get_degrees(self):
        """returns the array [0,...,lmax] that contains the degrees l"""
        return np.arange(self.lmax+1)

    def get_powerperdegree(self):
        """returns the power per degree l spectrum"""
        return self._powerperdegree()

    def get_powerperband(self,bandwidth):
        """returns the power per log_{bandwidth} l spectrum"""
        ls = self.get_degrees()
        return self._powerperdegree()*ls*np.log(bandwidth)

    #---- conversions ----
    def get_coeffs(self,normalization='4pi'):
        raise NotImplementedError('Not yet implemented')

    #---- rotation ----
    def rotate(self,alpha,beta,gamma,degrees=True):
        """
        rotates the spherical harmonics coefficients by
        alpha, beta, gamma
        """
        if degrees:
            angles = np.radians([alpha,beta,gamma])
        else:
            angles = np.array([alpha,beta,gamma])
        self._rotate(angles)

    #---- expansion ----
    def expand(self,kind='DH1'):
        """
        expands the coefficients to a spherical grid. 
        Available grid types:
           kind ='DH1': equidistant lat/lon grid with   nlat=nlon
           kind ='DH2': equidistant lat/lon grid with 2*nlat=nlon
           kind ='GLQ': Gauss Legendre Grid
        """
        if   kind == 'DH1':
            grid = self._expandDH(sampling=1)
        elif kind == 'DH2':
            grid = self._expandDH(sampling=2)
        else:
            raise NotImplementedError('grid type {:s} not found/implemented'.format(kind))
        return grid

    #---- plotting routines ----
    def plot_powerperdegree(self,loglog=True,show=True):
        """
        plots the power per degree spectrum. This is in particular useful to
        analyze global isotropic power at a certain wavelength.
        """
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

    def plot_powerperband(self,bandwidth=2,show=True):
        """
        plots the power per log_{bandwidth}(degree) spectrum. This is in
        particular useful to analyze local heterogeneity strength.
        """
        power = self.get_powerperband(bandwidth)
        ls    = self.get_degrees()

        fig,ax = plt.subplots(1,1)
        ax.set_xlabel('degree l')
        ax.set_ylabel('bandpower')
        ax.set_xscale('log',basex=bandwidth)
        ax.set_yscale('log',basey=bandwidth)
        ax.grid(True,which='both')
        ax.plot(ls[1:],power[1:],label='power per degree l')
        if show: plt.show()


#================== REAL SPHERICAL HARMONICS ================

class SHRealCoefficients(SHCoeffs):
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

    def make_complex(self,convention=1,switchcs=0):
        """converts the real coefficient class to the complex harmonic coefficient class"""
        complex_coeffs = SHrtoc(self.coeffs,convention=convention,switchcs=switchcs)
        return SHCoeffs.from_array(complex_coeffs,kind='complex')

    def _powerperdegree(self):
        """-> use powerperdegree instead of _powerperdegree"""
        return SHPowerSpectrum(self.coeffs)

    def _rotate(self,angles,dj_matrix=None):
        """-> use rotate instead of _rotate"""
        if dj_matrix is None:
            dj_matrix = djpi2(self.lmax+1)
        self.coeffs = SHRotateRealCoef(self.coeffs, angles, dj_matrix)

    def _expandDH(self,sampling):
        """-> use expand(kind='DH1') instead of _expandDH"""
        data = MakeGridDH(self.coeffs,sampling=sampling)
        grid = SHGrid.from_array(data)
        return grid

    def _expandGLQ(self,zeros=None):
        """-> use expand(kind='GLQ') instead of _expandGLQ"""
        if zeros is None:
            zeros, weights = PreCompute(lmax)
        return MakeGridGLQ(cilm_trim,zeros)



#=============== COMPLEX SPHERICAL HARMONICS ================


class SHComplexCoefficients(SHCoeffs):
    """
    Complex Spherical Harmonics Coefficients class.
    """
    @staticmethod
    def istype(kind):
        return kind=='complex'

    def __init__(self,coeffs):
        self.coeffs = coeffs
        self.lmax = coeffs.shape[1]-1

    def make_real(self,convention=1,switchcs=0):
        """converts the complex coefficient class to the real harmonic coefficient class"""
        complex_coeffs = SHctor(self.coeffs,convention=convention,switchcs=switchcs)
        return SHCoeffs.from_array(complex_coeffs,kind='real')

    def _powerperdegree():
        """-> use powerperdegree instead of _powerperdegree"""
        return SHCPowerSpectrum(self.coeffs)

    def _rotate(self,angles,dj_matrix=None):
        """-> use rotate instead of _rotate"""
        if dj_matrix is None:
            dj_matrix = djpi2(self.lmax+1)
        self.coeffs = SHRotateRealCoef(self.coeffs, angles, dj_matrix)

    def _expandDH(self,sammpling):
        """-> use expand(kind='DH1') instead of _expandDH"""
        data = MakeGridDHC(self.coeffs,sampling=sampling)
        grid = SHGrid.from_array(data)
        return grid

    def _expandGLQ(self, zeros=None):
        """-> use expand(kind='GLQ') instead of _expandGLQ"""
        if zeros is None:
            zeros, weights = PreCompute(lmax)
        return MakeGridGLQC(cilm_trim,zeros)






#========================================================================
#======      GRID CLASSES      ==========================================
#========================================================================



class SHGrid(object):
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

    def expand(self):
        return self._expand()

#---- implementation of the Driscoll and Healy Grid class ----
class DHGrid(SHGrid):
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

    def _expand(self):
        """-> use expand instead of _expand"""
        cilm   = SHExpandDH(self.data)
        coeffs = SHCoeffs.from_array(cilm)
        return coeffs

    def _plot_rawdata(self):
        """-> use plot_rawdata instead of _plot_rawdata"""
        fig,ax = plt.subplots(1,1)
        ax.imshow(self.data,origin='top',extent=(0.,360.,-90.,90.))
        ax.set_title('Driscoll Healy Grid')
        ax.set_xlabel('longitude')
        ax.set_ylabel('latitude')
        fig.tight_layout(pad=0.5)

#---- implementation of the Gauss-Legendre Grid class ----
class GLQGrid(SHGrid):
    """
    Gauss Legendre Class
    """
    @staticmethod
    def istype(kind):
        return kind == 'DH'

    def __init__(self, array):
        raise NotImplementedError('GLQ grid not yet implemented')

    def _plot_rawdata(self):
        fig,ax = plt.subplots(1,1)
        ax.imshow(self.data,origin='top',extent=(0.,360.,-90.,90.))
        ax.set_title('Driscoll Healy Grid')
        ax.set_xlabel('longitude')
        ax.set_ylabel('latitude')
        fig.tight_layout(pad=0.5)
