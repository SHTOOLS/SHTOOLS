# ========================================================================
# ======      GRID CLASSES      ==========================================
# ========================================================================

from __future__ import absolute_import, division, print_function

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from .. import _SHTOOLS as shtools
#from .SHCoeffs import SHCoeffs

class SHGrid(object):
    """
    Grid Class for global gridded data on the sphere. Grids can be
    initialized from:

    >> x = SHGrid.from_array(array)
    >> x = SHGrid.from_file('fname.dat')

    The class instance defines the following class attributes:

    data       : Gridded array of the data.
    nlat, nlon : The number of latitude and longitude bands in the grid.
    lmax       : The maximum spherical harmonic degree that can be resolved
                 by the grid sampling.
    sampling   : For Driscoll and Healy grids, the longitudinal sampling
                 of the grid. Either nlong = nlat or nlong = 2 * nlat.
    kind       : Either 'complex' or 'real' for the data type.
    grid       : Either 'DH' or 'GLQ' for Driscoll and Healy grids or Gauss-
                 Legendre Quadrature grids.
    zeros      : The cos(colatitude) nodes used with Gauss-Legendre
                 Quadrature grids. Default is None.
    weights    : The latitudinal weights used with Gauss-Legendre
                 Quadrature grids. Default is None.

    Each class instance provides the following methods:

    get_lats()     : Return a vector containing the latitudes of each row
                     of the gridded data.
    get_lons()     : Return a vector containing the longitudes of each column
                     of the gridded data.
    expand()       : Expand the grid into spherical harmonics.
    plot_rawdata() : Plot the raw data using a simple cylindrical projection.
    """

    def __init__():
        pass

    # ---- factory methods
    @classmethod
    def from_array(self, array, grid='DH'):
        """
        Initialize the grid of the class instance from an input array.

        Usage
        -----

        x = SHGrid.from_array(array, [grid])

        Parameters
        ----------

        array : numpy array of size (nlat, nlon)
        grid : 'DH' (default) or 'GLQ' for Driscoll and Healy grids or Gauss
                Legendre Quadrature grids, respectively.
        """
        if np.iscomplexobj(array):
            kind = 'complex'
        else:
            kind = 'real'

        if type(grid) != str:
            raise ValueError('grid must be a string. ' +
                             'Input type was {:s}'
                             .format(str(type(grid))))

        if grid.upper() not in set(['DH', 'GLQ']):
            raise ValueError(
                "grid must be 'DH' or 'GLQ'. Input value was {:s}."
                .format(repr(grid))
                )

        for cls in self.__subclasses__():
            if cls.istype(kind) and cls.isgrid(grid):
                return cls(array)

    @classmethod
    def from_file(self, fname, kind='real', grid='DH'):
        """Initialize the grid of the object from a file."""
        raise NotImplementedError('Not implemented yet')

    # ---- Extract grid properties ----
    def get_lats(self):
        """
        Return a vector containing the latitudes (in degrees) of each row
        of the gridded data.

        Usage
        -----

        lats = x.get_lats()

        Returns
        -------

        lats : numpy array of size nlat containing the latitude (in degrees)
               of each row of the gridded data.
        """
        return self._get_lats()

    def get_lons(self):
        """
        Return a vector containing the longitudes (in degrees) of each
        column of the gridded data.

        Usage
        -----

        lons = x.get_lon()

        Returns
        -------

        lons : numpy array of size nlon containing the longitude (in degrees)
               of each column of the gridded data.
        """
        return self._get_lons()

    # ---- Plotting routines ----
    def plot_rawdata(self, show=True, fname=None):
        """
        Plot the raw data using a simple cylindrical projection.

        Usage
        -----

        x.plot_rawdata([show, fname])

        Parameters
        ----------

        show   : If True (default), plot the image to the screen.
        fname  : If present, save the image to the file.
        """
        fig, ax = self._plot_rawdata()
        if show:
            plt.show()
        if fname is not None:
            fig.savefig(fname)

    def expand(self, normalization='4pi', csphase=1, **kwargs):
        """
        Expand the grid into spherical harmonics.

        Usage
        -----

        SHCoeffsInstance = x.expand([normalization, csphase, lmax_calc])

        Parameters
        ----------

        normalization : '4pi' (default), geodesy 4-pi normalized
                      : 'ortho', orthonormalized
                      : 'schmidt', Schmidt semi-normalized)
        csphase       : 1  (default), exlcude the Condon-Shortley phase factor
        lmax_calc     : maximum spherical harmonic degree to return.
                        Default is x.lmax.
        """
        if type(normalization) != str:
            raise ValueError('normalization must be a string. ' +
                             'Input type was {:s}'
                             .format(str(type(normalization))))

        if normalization.lower() not in set(['4pi', 'ortho', 'schmidt']):
            raise ValueError(
                "The normalization must be '4pi', 'ortho' " +
                "or 'schmidt'. Input value was {:s}."
                .format(repr(normalization))
                )

        if csphase != 1 and csphase != -1:
            raise ValueError(
                "csphase must be either 1 or -1. Input value was {:s}."
                .format(repr(csphase))
                )

        return self._expand(normalization=normalization, csphase=csphase,
                            **kwargs)


# ---- Real Driscoll and Healy grid class ----

class DHRealGrid(SHGrid):
    """
    Class for real Driscoll and Healy (1994) grids.
    """
    @staticmethod
    def istype(kind):
        return kind == 'real'

    @staticmethod
    def isgrid(grid):
        return grid == 'DH'

    def __init__(self, array):
        self.nlat, self.nlon = array.shape

        if self.nlat % 2 != 0:
            raise ValueError('Input arrays for DH grids must have an even ' +
                             'number of latitudes: nlat = {:d}'
                             .format(self.nlat)
                             )
        if self.nlon == 2 * self.nlat:
            self.sampling = 2
        elif self.nlat == self.nlon:
            self.sampling = 1
        else:
            raise ValueError('Input array has shape (nlat={:d},nlon={:d})\n' +
                             'but needs nlat=nlon or nlat=2*nlon'
                             .format(self.nlat, self.nlon)
                             )

        self.lmax = int(self.nlat / 2 - 1)
        self.data = array
        self.grid = 'DH'
        self.kind = 'real'

    def _get_lats(self):
        """
        Return a vector containing the latitudes (in degrees) of each row
        of the gridded data.
        """
        lats = np.linspace(90.0, -90.0 + 180.0 / self.nlat, num=self.nlat)
        return lats

    def _get_lons(self):
        """
        Return a vector containing the longitudes (in degrees) of each row
        of the gridded data.
        """
        lons = np.linspace(0.0, 360.0 - 360.0 / self.nlon, num=self.nlon)
        return lons

    def _expand(self, normalization, csphase, **kwargs):
        """Expand the grid into real spherical harmonics."""
        if normalization.lower() == '4pi':
            norm = 1
        elif normalization.lower() == 'schmidt':
            norm = 2
        elif normalization.lower() == 'ortho':
            norm = 4
        else:
            raise ValueError(
                "The normalization must be '4pi', 'ortho' " +
                "or 'schmidt'. Input value was {:s}."
                .format(repr(normalization))
                )

        cilm = shtools.SHExpandDH(self.data, norm=norm, csphase=csphase,
                                  **kwargs)
        coeffs = SHCoeffs.from_array(cilm,
                                     normalization=normalization.lower(),
                                     csphase=csphase)
        return coeffs

    def _plot_rawdata(self):
        """Plot the raw data using a simply cylindrical projection."""
        fig, ax = plt.subplots(1, 1)
        ax.imshow(self.data, origin='top', extent=(0., 360., -90., 90.))
        ax.set_title('Driscoll and Healy Grid')
        ax.set_xlabel('longitude')
        ax.set_ylabel('latitude')
        fig.tight_layout(pad=0.5)
        return fig, ax


# ---- Complex Driscoll and Healy grid class ----

class DHComplexGrid(SHGrid):
    """
    Class for complex Driscoll and Healy (1994) grids.
    """
    @staticmethod
    def istype(kind):
        return kind == 'complex'

    @staticmethod
    def isgrid(grid):
        return grid == 'DH'

    def __init__(self, array):
        self.nlat, self.nlon = array.shape

        if self.nlat % 2 != 0:
            raise ValueError('Input arrays for DH grids must have an even ' +
                             'number of latitudes: nlat = {:d}'
                             .format(self.nlat)
                             )
        if self.nlon == 2 * self.nlat:
            self.sampling = 2
        elif self.nlat == self.nlon:
            self.sampling = 1
        else:
            raise ValueError('Input array has shape (nlat={:d},nlon={:d})\n' +
                             'but needs nlat=nlon or nlat=2*nlon'
                             .format(self.nlat, self.nlon)
                             )

        self.lmax = int(self.nlat / 2 - 1)
        self.data = array
        self.grid = 'DH'
        self.kind = 'complex'

    def _get_lats(self):
        """
        Return a vector containing the latitudes (in degrees) of each row
        of the gridded data.
        """
        lats = np.linspace(90.0, -90.0 + 180.0 / self.nlat, num=self.nlat)
        return lats

    def _get_lons(self):
        """
        Return a vector containing the longitudes (in degrees) of each row
        of the gridded data.
        """
        lons = np.linspace(0., 360.0 - 360.0 / self.nlon, num=self.nlon)
        return lons

    def _expand(self, normalization, csphase, **kwargs):
        """Expand the grid into real spherical harmonics."""
        if normalization.lower() == '4pi':
            norm = 1
        elif normalization.lower() == 'schmidt':
            norm = 2
        elif normalization.lower() == 'ortho':
            norm = 4
        else:
            raise ValueError(
                "The normalization must be '4pi', 'ortho' " +
                "or 'schmidt'. Input value was {:s}."
                .format(repr(normalization))
                )

        cilm = shtools.SHExpandDHC(self.data, norm=norm, csphase=csphase,
                                   **kwargs)
        coeffs = SHCoeffs.from_array(cilm,
                                     normalization=normalization.lower(),
                                     csphase=csphase)
        return coeffs

    def _plot_rawdata(self):
        """Plot the raw data using a simply cylindrical projection."""
        fig, ax = plt.subplots(2, 1)
        ax.flat[0].imshow(self.data.real, origin='top',
                          extent=(0., 360., -90., 90.))
        ax.flat[0].set_title('Driscoll and Healy Grid (real component)')
        ax.flat[0].set_xlabel('longitude')
        ax.flat[0].set_ylabel('latitude')
        ax.flat[1].imshow(self.data.imag, origin='top',
                          extent=(0., 360., -90., 90.))
        ax.flat[1].set_title('Driscoll and Healy Grid (imaginary component)')
        ax.flat[1].set_xlabel('longitude')
        ax.flat[1].set_ylabel('latitude')
        fig.tight_layout(pad=0.5)
        return fig, ax


# ---- Real Gaus Legendre Quadrature grid class ----

class GLQRealGrid(SHGrid):
    """
    Class for real Gauss-Legendre Quadrature grids.
    """
    @staticmethod
    def istype(kind):
        return kind == 'real'

    @staticmethod
    def isgrid(grid):
        return grid == 'GLQ'

    def __init__(self, array, zeros=None, weights=None):
        self.nlat, self.nlon = array.shape
        self.lmax = self.nlat - 1

        if self.nlat != self.lmax + 1 or self.nlon != 2 * self.lmax + 1:
            raise ValueError('Input array has shape (nlat={:d}, nlon={:d})\n' +
                             'but needs (nlat={:d}, {:d})'
                             .format(self.nlat, self.nlon, self.lmax+1,
                                     2*self.lmax+1)
                             )

        if zeros is None or weights is None:
            self.zeros, self.weights = shtools.SHGLQ(self.lmax)
        else:
            self.zeros = zeros
            self.weights = weights

        self.data = array
        self.grid = 'GLQ'
        self.kind = 'real'

    def _get_lats(self):
        """
        Return a vector containing the latitudes (in degrees) of each row
        of the gridded data.
        """
        lats = 90. - np.arccos(self.zeros) * 180. / np.pi
        return lats

    def _get_lons(self):
        """
        Return a vector containing the longitudes (in degrees) of each column
        of the gridded data.
        """
        lons = np.linspace(0.0, 360.0 - 360.0 / self.nlon, num=self.nlon)
        return lons

    def _expand(self, normalization, csphase, **kwargs):
        """Expand the grid into real spherical harmonics."""
        if normalization.lower() == '4pi':
            norm = 1
        elif normalization.lower() == 'schmidt':
            norm = 2
        elif normalization.lower() == 'ortho':
            norm = 4
        else:
            raise ValueError(
                "The normalization must be '4pi', 'ortho' " +
                "or 'schmidt'. Input value was {:s}."
                .format(repr(normalization))
                )

        cilm = shtools.SHExpandGLQ(self.data, self.weights, self.zeros,
                                   norm=norm, csphase=csphase, **kwargs)
        coeffs = SHCoeffs.from_array(cilm,
                                     normalization=normalization.lower(),
                                     csphase=csphase)
        return coeffs

    def _plot_rawdata(self):
        """Plot the raw data using a simply cylindrical projection."""

        fig, ax = plt.subplots(1, 1)
        ax.imshow(self.data, origin='top')
        ax.set_title('Gauss-Legendre Quadrature Grid')
        ax.set_xlabel('longitude index')
        ax.set_ylabel('latitude index')
        fig.tight_layout(pad=0.5)
        return fig, ax


# ---- Complex Gaus Legendre Quadrature grid class ----

class GLQComplexGrid(SHGrid):
    """
    Class for complex Gauss Legendre Quadrature grids.
    """
    @staticmethod
    def istype(kind):
        return kind == 'complex'

    @staticmethod
    def isgrid(grid):
        return grid == 'GLQ'

    def __init__(self, array, zeros=None, weights=None):
        self.nlat, self.nlon = array.shape
        self.lmax = self.nlat - 1

        if self.nlat != self.lmax + 1 or self.nlon != 2 * self.lmax + 1:
            raise ValueError('Input array has shape (nlat={:d}, nlon={:d})\n' +
                             'but needs (nlat={:d}, {:d})'
                             .format(self.nlat, self.nlon, self.lmax+1,
                                     2*self.lmax+1)
                             )

        if zeros is None or weights is None:
            self.zeros, self.weights = shtools.SHGLQ(self.lmax)
        else:
            self.zeros = zeros
            self.weights = weights

        self.data = array
        self.grid = 'GLQ'
        self.kind = 'complex'

    def _get_lats(self):
        """
        Return a vector containing the latitudes (in degrees) of each row
        of the gridded data.
        """
        lats = 90. - np.arccos(self.zeros) * 180. / np.pi
        return lats

    def _get_lons(self):
        """
        Return a vector containing the longitudes (in degrees) of each column
        of the gridded data.
        """
        lons = np.linspace(0., 360. - 360. / self.nlon, num=self.nlon)
        return lons

    def _expand(self, normalization, csphase, **kwargs):
        """Expand the grid into real spherical harmonics."""
        if normalization.lower() == '4pi':
            norm = 1
        elif normalization.lower() == 'schmidt':
            norm = 2
        elif normalization.lower() == 'ortho':
            norm = 4
        else:
            raise ValueError(
                "The normalization must be '4pi', 'ortho' " +
                "or 'schmidt'. Input value was {:s}."
                .format(repr(normalization))
                )

        cilm = shtools.SHExpandGLQC(self.data, self.weights, self.zeros,
                                    norm=norm, csphase=csphase, **kwargs)
        coeffs = SHCoeffs.from_array(cilm,
                                     normalization=normalization.lower(),
                                     csphase=csphase)
        return coeffs

    def _plot_rawdata(self):
        """Plot the raw data using a simply cylindrical projection."""

        fig, ax = plt.subplots(2, 1)
        ax.flat[0].imshow(self.data.real, origin='top')
        ax.flat[0].set_title('Gauss-Legendre Quadrature Grid (real component)')
        ax.flat[0].set_xlabel('longitude index')
        ax.flat[0].set_ylabel('latitude index')
        ax.flat[1].imshow(self.data.imag, origin='top')
        ax.flat[1].set_title('Gauss-Legendre Quadrature Grid ' +
                             '(imaginary component)')
        ax.flat[1].set_xlabel('longitude index')
        ax.flat[1].set_ylabel('latitude index')
        fig.tight_layout(pad=0.5)
        return fig, ax
