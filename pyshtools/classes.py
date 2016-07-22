"""
pyshtools defines several classes that facilitate the interactive
examination of geographical gridded data and their associated 
spherical harmonic coefficients. Subclasses are used to handle different 
internal data types and superclasses are used to implement interface
functions and the documentation.

The following classes and subclasses are defined:

    SHCoeffs
        SHRealCoefficients
        SHComplexCoefficients

    SHGrid
        DHGrid
        GLQGrid

    SHWindow
        SymmetricWindow
        AsymmetricWindow
"""

from __future__ import absolute_import, division, print_function

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from ._SHTOOLS import *

#===============================================================================
#=========    COEFFICIENT CLASSES    ===========================================
#===============================================================================

class SHCoeffs(object):
    """
    Spherical Harmonics Coefficient class. Coefficients can be initialized
    using one of the three constructor methods:

    >> x = SHCoeffs.from_array( np.zeros( (2, lmax+1, lmax+1) ) )
    >> x = SHCoeffs.from_random( np.exp(-ls**2) )
    >> x = SHCoeffs.from_file( 'fname.dat' )
    """

    def __init__(self):
        pass

    #---- factory methods:
    @classmethod
    def from_array(self, coeffs, kind='real', normalization='4pi', csphase=False):
        """"
        Initialize the spherical harmonic coefficients of the object
        using an input array dimensioned as (2, lmax+1, lmax+1).
        """
        for cls in self.__subclasses__():
            if cls.istype(kind):
                return cls(coeffs, normalization=normalization, csphase=csphase)

    @classmethod
    def from_random(self, power, kind='real', normalization='4pi', csphase=False):
        """"
        Initialize the spherical harmonic coefficients of the object
        using Gaussian random variables with a given input power spectrum.
        """
        nl = len(power)
        for cls in self.__subclasses__():
            if cls.istype(kind):
                if kind == 'real':
                    coeffs = np.random.normal(size=(2, nl, nl))
                    coeffs *= np.sqrt(power)[np.newaxis, :, np.newaxis]
                elif kind == 'complex':
                    coeffs = (np.random.normal(loc=0., scale=1.,
                                               size=(2, nl, nl)) +
                              1j * np.random.normal(loc=0., scale=1.,
                                                    size=(2, nl, nl)))
                    coeffs *= np.sqrt(power)[np.newaxis, :, np.newaxis]
                else:
                    raise ValueError("kind='{:s}': Should be 'real' or 'complex'".format(str(kind)))
                return cls(coeffs, normalization=normalization, csphase=csphase)

    @classmethod
    def from_file(self, fname, lmax, format='shtools', kind='real', normalization='4pi', csphase=False):
        """
        Initialize the spherical harmonic coefficients of the object
        by reading the coefficients from a specified file name.
        """
        if format == 'shtools' and kind == 'real':
            coeffs, lmax = SHRead(fname, lmax)
        else:
            raise NotImplementedError('Not yet implemented')
            
        lmax = coeffs.shape[1] - 1

        for cls in self.__subclasses__():
            if cls.istype(kind):
                return cls(coeffs, normalization=normalization, csphase=csphase)

    #---- extracting data ----
    def get_degrees(self):
        """Return an array listing the spherical harmonic degrees from 0 to lmax."""
        return np.arange(self.lmax + 1)

    def get_powerperdegree(self):
        """Return the power per degree l spectrum."""
        return self._powerperdegree()

    def get_powerperband(self, bandwidth):
        """Return the power per log_{bandwidth} l spectrum"""
        ls = self.get_degrees()
        return self._powerperdegree() * ls * np.log(bandwidth)

    #---- conversions ----
    def get_coeffs(self, normalization='4pi', kind='real'):
        """
        Return complex or real spherical harmonics coefficients 
        in a different normalization.
        """
        # add output kind and output normalization parameters. 
        # get input values from input class.
        # consider renaming this to "convert" and use only for converting between 
        # normalizations (make_complex and make_real already exist)
        coeffs = self._coeffs(kind)
        raise NotImplementedError('Not yet implemented')

    #---- rotation ----
    def rotate(self, alpha, beta, gamma, degrees=True):
        """
        Rotate the spherical harmonics coefficients by
        the Euler angles alpha, beta, gamma.
        """
        if degrees:
            angles = np.radians([alpha, beta, gamma])
        else:
            angles = np.array([alpha, beta, gamma])
        self._rotate(angles)

    #---- expansion ----
    def expand(self, **kwargs):
        """
        Evaluate the coefficients on a spherical grid. 
        Available grid types:
           kind = 'DH1' or 'DH': equisampled lat/lon grid with nlat=nlon
           kind = 'DH2': equidistant lat/lon grid with nlon=2*nlat
           kind = 'GLQ': Gauss Legendre quadrature grid
        """
        if kwargs['kind'] == 'DH1' or kwargs['kind'] == 'DH':
            grid = self._expandDH(sampling=1)
        elif kwargs['kind'] == 'DH2':
            grid = self._expandDH(sampling=2)
        elif kwargs['kind'] == 'GLQ':
            grid = self._expandGLQ()
        else:
            raise NotImplementedError('Grid type {:s} not implemented'.format(kind))
        return grid

    #---- plotting routines ----
    def plot_powerperdegree(self, loglog=True, show=True, fname=None):
        """
        Plot the power per degree spectrum.
        """
        power = self.get_powerperdegree()
        ls = self.get_degrees()

        fig, ax = plt.subplots(1, 1)
        ax.set_xlabel('degree l')
        ax.set_ylabel('power per degree')
        if loglog:
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.grid(True, which='both')
            ax.plot(ls[1:], power[1:], label='power per degree l')
        if show:
            plt.show()
        if fname is not None:
            fig.savefig(fname)

    def plot_powerperband(self, bandwidth=2, show=True, fname=None):
        """
        Plots the power per log_{bandwidth}(degree) spectrum.
        """
        power = self.get_powerperband(bandwidth)
        ls = self.get_degrees()

        fig, ax = plt.subplots(1, 1)
        ax.set_xlabel('degree l')
        ax.set_ylabel('bandpower')
        ax.set_xscale('log', basex=bandwidth)
        ax.set_yscale('log', basey=bandwidth)
        ax.grid(True, which='both')
        ax.plot(ls[1:], power[1:], label='power per degree l')
        fig.tight_layout(pad=0.1)
        if show:
            plt.show()
        if fname is not None:
            fig.savefig(fname)


#================== REAL SPHERICAL HARMONICS ================

class SHRealCoefficients(SHCoeffs):
    """
    Real Spherical Harmonics Coefficient class.
    """
    
    @staticmethod
    def istype(kind):
        return kind == 'real'

    def __init__(self, coeffs, normalization='4pi', csphase=False):
        #---- create mask to filter out m<=l ----
        lmax = coeffs.shape[1] - 1
        mask = np.zeros((2, lmax + 1, lmax + 1), dtype=np.bool)
        mask[0, 0, 0] = True
        for l in np.arange(lmax + 1):
            mask[:, l, :l + 1] = True
        mask[1, :, 0] = False

        self.lmax = lmax
        self.coeffs = np.copy(coeffs)
        self.coeffs[np.invert(mask)] = 0.
        self.kind = 'real'
        self.normalization = normalization
        self.csphase = csphase

    def make_complex(self, convention=1, switchcs=0):
        """Convert the real coefficient class to the complex harmonic coefficient class."""
        complex_coeffs = SHrtoc(self.coeffs, convention=convention, switchcs=switchcs)
        return SHCoeffs.from_array(complex_coeffs, kind='complex')

    def _powerperdegree(self):
        """Return the power per degree l spectrum."""
        return SHPowerSpectrum(self.coeffs)

    def _get_coeffs(self, kind='real', convention=1, swithchcs=0):
        """
        Return complex or real spherical harmonics coefficients 
        in a different normalization.
        """
        if kind=='real': 
            return self.coeffs
             # need to account for switching normalizations
        elif kind=='complex':
            return SHrtoc(self.coeffs, convention=convention, switchcs=switchcs)

    def _rotate(self, angles, dj_matrix=None):
        """
        Rotate the spherical harmonics coefficients by
        the Euler angles alpha, beta, gamma.
        """
        if dj_matrix is None:
            dj_matrix = djpi2(self.lmax + 1)
        self.coeffs = SHRotateRealCoef(self.coeffs, angles, dj_matrix)

    def _expandDH(self, sampling):
        """Evaluate the coefficients on a DH grid.""" 
        data = MakeGridDH(self.coeffs, sampling=sampling)
        grid = SHGrid.from_array(data, 'DH')
        return grid

    def _expandGLQ(self, zeros=None):
        """Evaluate the coefficients on a GLQ grid.""" 
        if zeros is None:
            zeros, weights = SHGLQ(self.lmax)
        data = MakeGridGLQ(self.coeffs, zeros)
        grid = SHGrid.from_array(data, 'GLQ')
        return grid


#=============== COMPLEX SPHERICAL HARMONICS ================

class SHComplexCoefficients(SHCoeffs):
    """
    Complex Spherical Harmonics Coefficients class.
    """

    @staticmethod
    def istype(kind):
        return kind == 'complex'

    def __init__(self, coeffs, normalization='4pi', csphase=False):
        #---- create mask to filter out m<=l ---- 
        # not implemented
        self.coeffs = coeffs
        self.lmax = coeffs.shape[1] - 1
        self.kind = 'complex'
        self.normalization = normalization
        self.csphase = csphase

    def make_real(self, convention=1, switchcs=0):
        """Convert the complex coefficient class to the real harmonic coefficient class."""
        complex_coeffs = SHctor(self.coeffs, convention=convention, switchcs=switchcs)
        return SHCoeffs.from_array(complex_coeffs, kind='real')

    def _powerperdegree(self):
        """Return the power per degree l spectrum."""
        return SHCPowerSpectrum(self.coeffs)

    def _get_coeffs(self, kind='complex', convention=1, swithchcs=0)):
         """
        Return complex or real spherical harmonics coefficients 
        in a different normalization.
        """
        if kind=='real': 
            return SHctor(self.coeffs, convention=convention, switchcs=switchcs)
        elif kind=='complex':
            return self.coeffs
            # need to account for switching normalizations

    def _rotate(self, angles, dj_matrix=None):
        """
        Rotate the spherical harmonics coefficients by
        the Euler angles alpha, beta, gamma.
        """
        if dj_matrix is None:
            dj_matrix = djpi2(self.lmax + 1)
        self.coeffs = SHRotateRealCoef(self.coeffs, angles, dj_matrix)

    def _expandDH(self, sammpling):
        """Evaluate the coefficients on a DH grid."""
        data = MakeGridDHC(self.coeffs, sampling=sampling)
        grid = SHGrid.from_array(data)
        return grid

    def _expandGLQ(self, zeros=None):
        """Evaluate the coefficients on a GLQ grid."""
        if zeros is None:
            zeros, weights = SHGLQ(self.lmax)
        data = MakeGridGLQ(self.coeffs, zeros)
        grid = SHGrid.from_array(data, 'GLQ')
        return grid


#========================================================================
#======      GRID CLASSES      ==========================================
#========================================================================

class SHGrid(object):

    """
    EXPERIMENTAL:
    Spherical Grid Class that can deal with spatial data on the sphere that is
    defined on different grids. Can be constructed from:

    >> SphericalGrid.from_array(...)
    >> SphericalGrid.from_file(...)
    """
    def __init__():
        print('use one of the following methods to initialize the grid:\n\n' +
              '>> SphericalGrid.from_array(...)\n' +
              '>> SphericalGrid.from_file(...)')

    #---- constructors ----
    @classmethod
    def from_array(self, array, kind='DH'):
        for cls in self.__subclasses__():
            if cls.istype(kind):
                return cls(array)

    #---- extract data ----
    def get_lats():
        return self._get_lats()

    def get_lons():
        return self._get_lats()

    #---- plotting routines ----
    def plot_rawdata(self, show=True, fname=None):
        fig,ax = self._plot_rawdata()
        if show:
            plt.show()
        if fname is not None:
            fig.savefig(fname)

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
        self.nlat, self.nlon = array.shape
        if self.nlat == 2 * self.nlon:
            self.sampling = 2
        if self.nlat == self.nlon:
            self.sampling = 1
        else:
            raise ValueError('input array with shape (nlat={:d},nlon={:d})\n' +
                             'it needs nlat=nlon or nlat=2*nlon'.format(self.nlat, self.nlon))
        self.data = array

    def _getlats(self):
        dlat = 360. / self.nlat
        lats = np.linspace(0. + dlat / 2., 360. - dlat / 2., self.nlat)
        return lats

    def _getlons(self):
        dlon = 360. / self.nlon
        lons = np.linspace(0. + dlon / 2., 360. - dlon / 2., self.nlon)
        return lons

    def _expand(self):
        """-> use expand instead of _expand"""
        cilm = SHExpandDH(self.data)
        coeffs = SHCoeffs.from_array(cilm, kind='DH')
        return coeffs

    def _plot_rawdata(self):
        """-> use plot_rawdata instead of _plot_rawdata"""
        fig, ax = plt.subplots(1, 1)
        ax.imshow(self.data, origin='top', extent=(0., 360., -90., 90.))
        ax.set_title('Driscoll Healy Grid')
        ax.set_xlabel('longitude')
        ax.set_ylabel('latitude')
        fig.tight_layout(pad=0.5)
        return fig,ax

#---- implementation of the Gauss-Legendre Grid class ----


class GLQGrid(SHGrid):
    """
    Gauss Legendre Class
    """
    @staticmethod
    def istype(kind):
        return kind == 'GLQ'

    def __init__(self, array, zeros=None):
        """use superclass constructors"""
        #---- check if input is correct ----
        self.nlat, self.nlon = array.shape
        assert self.nlon - 1 == 2 * (self.nlat - 1), 'nlon should equal 2*nlat for GLQ grid'

        #---- store data in class ----
        self.lmax = self.nlat - 1
        self.data = array
        if zeros is None:
            self.zeros, weights = SHGLQ(self.lmax)
        else:
            self.zeros = zeros

    def _getlats(self):
        """-> use getlats instead of _getlats"""
        lats = 90. - np.degrees(self.zeros)
        return lats

    def _getlons(self):
        """-> use getlons instead of _getlons"""
        dlon = 360. / self.nlon
        lons = np.linspace(0. + dlon / 2., 360. - dlon / 2., self.nlon)
        return lons

    def _expand(self):
        """-> use expand instead of _expand"""
        cilm = SHExpandGLQ(self.data)
        coeffs = SHCoeffs.from_array(cilm, kind='GLQ')
        return coeffs

    def _plot_rawdata(self):
        """use plot_rawdata instead of _rawdata"""
        fig, ax = plt.subplots(1, 1)
        ax.imshow(self.data, origin='top')
        ax.set_title('Gauss-Legendre Quadrature Grid')
        ax.set_xlabel('longitude index')
        ax.set_ylabel('latitude  index')
        fig.tight_layout(pad=0.5)
        return fig,ax


#==== SPHERICAL HARMONICS WINDOW FUNCTION CLASS ====
class SHWindow(object):
    """
    EXPERIMENTAL:
    This class contains collections of spherical harmonics windows that
    provide spectral estimates about a specific region
    """
    def __init__(self):
        print("use one of the following constructors: [...]")

    @classmethod
    def from_cap(self, lmax, nwins, theta, clat=0., clon=0., degrees=True):
        """
        constructs a spherical cap window
        """
        if degrees:
            theta = np.radians(theta)
        tapers, eigenvalues, taper_order = SHReturnTapers(theta, lmax)
        return SHSymmetricWindow(tapers,eigenvalues,taper_order,clat=clat,clon=clon)

    @classmethod
    def from_mask(self, lmax, nwins, dh_mask, sampling=1):
        """
        constructs optimal window functions in a masked region (needs dh grid)
        """
        tapers, eigenvalues = SHReturnTapersMap(dh_mask, lmax, sampling=sampling, Ntapers=nwins)
        return SHAsymmetricWindow(tapers,eigenvalues)

    def plot(self,nwins,show=True,fname=None):
        """
        plots the best concentrated spherical harmonics taper functions
        """
        #setup figure and axes ...
        maxcolumns = 5
        ncolumns = min(maxcolumns,nwins)
        nrows    = np.ceil(nwins/ncolumns).astype(int)
        figsize  = ncolumns * 1.2, nrows * 1.2 + 0.5
        fig, axes = plt.subplots(nrows, ncolumns, figsize=figsize)
        for ax in axes[:-1,:].flatten():
            for xlabel_i in ax.get_xticklabels():
                xlabel_i.set_visible(False)
        for ax in axes[:,1:].flatten():
            for ylabel_i in ax.get_yticklabels():
                ylabel_i.set_visible(False)

        #loop through tapers and plot them
        for itaper in range( min(self.nwins, nwins) ):
            evalue = self.eigenvalues[itaper]
            coeffs = self._coeffs(itaper)
            ax = axes.flatten()[itaper]
            grid = MakeGridDH(coeffs)
            ax.imshow(grid)
            ax.set_title('concentration: {:2.2f}'.format(evalue))
        fig.tight_layout(pad=0.5)

        if show: plt.show()
        if fname is not None:
            fig.savefig(fname)


    def get_spectrum(self, shcoeffs, nwins):
        """Returns the regional spherical harmonics spectrum"""
        for itaper in range(nwins):
            tapercoeffs = self._coeffs(itaper)
            modelcoeffs = shcoeffs.get_coeffs(normalization='4pi',kind='real')
            coeffs = SHMultiply(tapercoeffs,modelcoeffs)

    def get_couplingmatrix(self,lmax,nwins):
        """returns the coupling matrix of the first nwins tapers"""
        #store sqrt of taper power in 'tapers' array:
        if nwins>self.nwins: nwins = self.nwins
        tapers = np.zeros( (self.nl, nwins) )
        for itaper in range(nwins):
            tapers[:,itaper] = np.sqrt(SHPowerSpectrum(self._coeffs(itaper)))

        #compute coupling matrix of the first nwins tapers:
        coupling_matrix = SHMTCouplingMatrix(lmax,tapers[:,:nwins])
        return coupling_matrix

    def plot_couplingmatrix(self,lmax,nwins,show=True,fname=None):
        """plots the window's coupling strength"""
        figsize = mpl.rcParams['figure.figsize']
        figsize[0] = figsize[1]
        fig = plt.figure(figsize=figsize)
        ax  = fig.add_subplot(111)
        coupling_matrix = self.get_couplingmatrix(lmax,nwins)
        ax.imshow(coupling_matrix)
        ax.set_xlabel('output power')
        ax.set_ylabel('input power')
        fig.tight_layout(pad=0.1)

        if show: plt.show()
        if fname is not None:
            fig.savefig(fname)

    def info(self):
        """print meta information about the tapers"""
        self._info()

class SHSymmetricWindow(SHWindow):
    """
    This class saves a symmetric spherical window function. It needs to
    save only the m=0 coefficients
    """
    @staticmethod
    def istype(kind):
        return kind == 'Symmetric'

    def __init__(self,tapers, eigenvalues, orders, clat=0., clon =0.):
        self.clat,self.clon = clat,clon    #center of cap window
        self.nl, self.nwins = tapers.shape #nl: number of degrees, nwins: number of windows
        self.lmax           = self.nl - 1  #lmax: maximum degree
        self.tapers         = tapers       #tapers[nl,nwins]: ith window coefs with m=orders[iwin]
        self.eigenvalues    = eigenvalues  #concentration factor of the ith taper
        self.orders         = orders      #order m of the ith taper

    def _coeffs(self,itaper):
        taperm = self.orders[itaper]
        coeffs = np.zeros( (2, self.nl, self.nl) )
        if taperm<0:
            coeffs[1,:,abs(taperm)] = self.tapers[:,itaper]
        else:
            coeffs[0,:,abs(taperm)] = self.tapers[:,itaper]
        return coeffs

    def _info(self):
        print('Cap window with {:d} tapers'.format(self.nwins))

class SHAsymmetricWindow(SHWindow):
    """
    This class saves a asymmetric spherical window function and is much like a
    set of real sherical harmonics. It could maybe be merged at some point...
    """
    @staticmethod
    def istype(kind):
        return kind == 'Asymmetric'

    def __init__(self,tapers,eigenvalues):
        ncoeffs,self.nwins = tapers.shape
        self.nl = np.sqrt(ncoeffs).astype(int)
        self.lmax = self.nl-1
        self.tapers = tapers
        self.eigenvalues = eigenvalues

    def _coeffs(self,itaper):
        return SHVectorToCilm(self.tapers[:,itaper],self.lmax)

    def _info(self):
        print('Asymmetric window with {:d} tapers'.format(self.nwins))
