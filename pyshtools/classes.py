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
        DHRealGrid
        DHComplexGrid
        GLQRealGrid
        GLQComplexGrid

    SHWindow
        SymmetricWindow
        AsymmetricWindow
"""

from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from . import _SHTOOLS as _shtools


# =============================================================================
# =========    COEFFICIENT CLASSES    =========================================
# =============================================================================

class SHCoeffs(object):
    """
    Spherical Harmonics Coefficient class. Coefficients can be initialized
    using one of the three constructor methods:

    >> x = SHCoeffs.from_array(np.zeros((2, lmax+1, lmax+1)))
    >> x = SHCoeffs.from_random(np.exp(-ls**2))
    >> x = SHCoeffs.from_file('fname.dat')
    """

    def __init__(self):
        pass

    # ---- factory methods:
    @classmethod
    def from_array(self, coeffs, normalization='4pi', csphase=1):
        """
        Initialize the spherical harmonic coefficients of the object
        using an input ndarray dimensioned as (2, lmax+1, lmax+1).
        """
        if np.iscomplexobj(coeffs):
            kind = 'complex'
        else:
            kind = 'real'

        if normalization.lower() not in set(['4pi', 'ortho', 'schmidt']):
            raise ValueError(
                "The input normalization must be '4pi', 'ortho' " +
                "or 'schmidt'. Provided value was '{:s}'"
                .format(str(normalization))
                )

        if csphase != 1 and csphase != -1:
            raise ValueError(
                "csphase must be 1 or -1. Input value was {:d}"
                .format(csphase)
                )

        for cls in self.__subclasses__():
            if cls.istype(kind):
                return cls(coeffs, normalization=normalization.lower(),
                           csphase=csphase)

    @classmethod
    def from_random(self, power, kind='real', normalization='4pi', csphase=1):
        """
        Initialize the spherical harmonic coefficients using Gaussian
        random variables and a given input power spectrum.
        """
        if normalization.lower() not in set(['4pi', 'ortho', 'schmidt']):
            raise ValueError(
                "The input normalization must be '4pi', 'ortho' " +
                "or 'schmidt'. Provided value was '{:s}'"
                .format(str(normalization))
                )

        if csphase != 1 and csphase != -1:
            raise ValueError(
                "csphase must be 1 or -1. Input value was {:d}"
                .format(csphase)
                )

        if kind not in set(['real', 'complex']):
            raise ValueError(
                "kind must be 'real' or 'complex'. " +
                "Input value was '{:s}'.".format(str(kind)))

        nl = len(power)
        l = np.arange(nl)

        if kind == 'real':
            coeffs = np.random.normal(size=(2, nl, nl))
        elif kind == 'complex':
            coeffs = (np.random.normal(size=(2, nl, nl)) +
                      1j * np.random.normal(size=(2, nl, nl)))

        if normalization.lower() == '4pi':
            coeffs *= np.sqrt(
                power / (2.0 * l + 1.0))[np.newaxis, :, np.newaxis]
        elif normalization.lower() == 'ortho':
            coeffs *= np.sqrt(
                4.0 * np.pi * power / (2.0 * l + 1.0)
                )[np.newaxis, :, np.newaxis]
        elif normalization.lower() == 'schmidt':
            coeffs *= np.sqrt(power)[np.newaxis, :, np.newaxis]

        for cls in self.__subclasses__():
            if cls.istype(kind):
                return cls(coeffs, normalization=normalization.lower(),
                           csphase=csphase)

    @classmethod
    def from_file(self, fname, lmax, format='shtools', kind='real',
                  normalization='4pi', csphase=1):
        """
        Initialize the spherical harmonic coefficients by reading the
        coefficients from a specified file.
        """
        if normalization.lower() not in set(['4pi', 'ortho', 'schmidt']):
            raise ValueError(
                "The input normalization must be '4pi', 'ortho' " +
                "or 'schmidt'. Provided value was '{:s}'"
                .format(str(normalization))
                )
        if csphase != 1 and csphase != -1:
            raise ValueError(
                "csphase must be 1 or -1. Input value was {:d}"
                .format(csphase)
                )

        if format == 'shtools':
            if kind == 'real':
                coeffs, lmax = _shtools.SHRead(fname, lmax)
            else:
                raise NotImplementedError(
                    "kind='{:s}' not yet implemented".format(str(kind)))
        else:
            raise NotImplementedError(
                "format='{:s}' not yet implemented".format(str(format)))

        for cls in self.__subclasses__():
            if cls.istype(kind):
                return cls(coeffs, normalization=normalization.lower(),
                           csphase=csphase)

    # ---- Extract data ----
    def get_degrees(self):
        """
        Return an array listing the spherical harmonic degrees
        from 0 to lmax.
        """
        return np.arange(self.lmax + 1)

    def get_powerperdegree(self):
        """Return the power per degree l spectrum."""
        return self._powerperdegree()

    def get_powerperband(self, bandwidth):
        """Return the power per log_{bandwidth} l spectrum"""
        ls = self.get_degrees()
        return self._powerperdegree() * ls * np.log(bandwidth)

    # ---- Return coefficients with a different normalization convention ----
    def get_coeffs(self, normalization='4pi', csphase=1):
        """
        Return spherical harmonics coefficients as an ndarray with
        a different normalization convention.
        """
        if normalization.lower() not in set(['4pi', 'ortho', 'schmidt']):
            raise ValueError(
                "normalization must be '4pi', 'ortho' " +
                "or 'schmidt'. Provided value was '{:s}'"
                .format(str(output_normalization))
                )
        if csphase != 1 and csphase != -1:
            raise ValueError(
                "csphase must be 1 or -1. Input value was {:d}"
                .format(csphase)
                )

        return self._get_coeffs(
            output_normalization=normalization.lower(),
            output_csphase=csphase)

    # ---- Rotate the coordinate system ----
    def rotate(self, alpha, beta, gamma, degrees=True, dj_matrix=None):
        """
        Rotate the coordinate system used to express the spherical
        harmonics coefficients by the Euler angles alpha, beta, gamma,
        and output as a new class instance.
        """
        if degrees:
            angles = np.radians([alpha, beta, gamma])
        else:
            angles = np.array([alpha, beta, gamma])

        self._rotate(angles, dj_matrix)

    # ---- Expand the coefficients onto a grid ----
    def expand(self, grid='DH'):
        """
        Evaluate the coefficients on a spherical grid.
        Available grid types are:
           grid = 'DH' or 'DH1': equisampled lat/lon grid with nlat=nlon
           grid = 'DH2': equidistant lat/lon grid with nlon=2*nlat
           grid = 'GLQ': Gauss-Legendre quadrature grid
        """
        if grid == 'DH' or grid == 'DH1':
            gridout = self._expandDH(sampling=1)
        elif grid == 'DH2':
            gridout = self._expandDH(sampling=2)
        elif grid == 'GLQ':
            gridout = self._expandGLQ(zero=None)
        else:
            raise ValueError(
                "grid must be 'DH', 'DH1', 'DH2', or 'GLQ'. " +
                "Input value was '{:s}'".format(grid))

        return gridout

    # ---- plotting routines ----
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


# ================== REAL SPHERICAL HARMONICS ================

class SHRealCoefficients(SHCoeffs):
    """
    Real Spherical Harmonics Coefficient class.
    """

    @staticmethod
    def istype(kind):
        return kind == 'real'

    def __init__(self, coeffs, normalization='4pi', csphase=1):
        lmax = coeffs.shape[1] - 1
        # ---- create mask to filter out m<=l ----
        mask = np.zeros((2, lmax + 1, lmax + 1), dtype=np.bool)
        mask[0, 0, 0] = True
        for l in np.arange(lmax + 1):
            mask[:, l, :l + 1] = True
        mask[1, :, 0] = False

        self.mask = mask
        self.lmax = lmax
        self.coeffs = np.copy(coeffs)
        self.coeffs[np.invert(mask)] = 0.
        self.kind = 'real'
        self.normalization = normalization.lower()
        self.csphase = csphase

    def make_complex(self, convention=1, switchcs=0):
        """
        Convert the real coefficient class to the complex harmonic
        coefficient class.
        """
        raise NotImplementedError('Not implemented yet!')

        complex_coeffs = _shtools.SHrtoc(self.coeffs, convention=convention,
                                         switchcs=switchcs)
        # NOT DONE. THESE COEFFICIENTS ARE STILL REAL FLOATS!
#        return SHCoeffs.from_array(complex_coeffs, kind='complex')

    def _powerperdegree(self):
        """Return the power per degree l spectrum."""
        if self.normalization == '4pi':
            return _shtools.SHPowerSpectrum(self.coeffs)
        elif self.normalization == 'schmidt':
            power = _shtools.SHPowerSpectrum(self.coeffs)
            l = self.get_degrees()
            power /= (2.0 * l + 1.0)
            return power
        elif self.normalization == 'ortho':
            return _shtools.SHPowerSpectrum(self.coeffs) / (4.0 * np.pi)
        else:
            raise ValueError(
                "Normalization must be '4pi', 'ortho', or 'schmidt'. " +
                "Input value was '{:s}'".format(str(self.normalization)))

    def _get_coeffs(self, output_normalization, output_csphase):
        """
        Return real spherical harmonic coefficients with a
        different normalization convention.
        """
        coeffs = np.copy(self.coeffs)

        if self.normalization == output_normalization:
            pass
        elif (self.normalization == '4pi' and
              output_normalization == 'schmidt'):
            for l in self.get_degrees():
                coeffs[:, l, :l+1] *= np.sqrt(2.0 * l + 1.0)
        elif self.normalization == '4pi' and output_normalization == 'ortho':
            coeffs *= np.sqrt(4.0 * np.pi)
        elif (self.normalization == 'schmidt' and
              output_normalization == '4pi'):
            for l in self.get_degrees():
                coeffs[:, l, :l+1] /= np.sqrt(2.0 * l + 1.0)
        elif (self.normalization == 'schmidt' and
              output_normalization == 'ortho'):
            for l in self.get_degrees():
                coeffs[:, l, :l+1] *= np.sqrt(4.0 * np.pi / (2.0 * l + 1.0))
        elif self.normalization == 'ortho' and output_normalization == '4pi':
            coeffs /= np.sqrt(4.0 * np.pi)
        elif (self.normalization == 'ortho' and
              output_normalization == 'schmidt'):
            for l in self.get_degrees():
                coeffs[:, l, :l+1] *= np.sqrt((2.0 * l + 1.0) / (4.0 * np.pi))

        if output_csphase != self.csphase:
            for m in self.get_degrees():
                if m % 2 == 1:
                    coeffs[:, :, m] = - coeffs[:, :, m]

        return coeffs

    def _rotate(self, angles, dj_matrix):
        """
        Rotate the coordinate system used to express the spherical
        harmonics coefficients by the Euler angles alpha, beta, gamma.
        """
        if dj_matrix is None:
            dj_matrix = _shtools.djpi2(self.lmax + 1)

        # The coefficients need to be 4pi normalized with csphase = 1
        coeffs = _shtools.SHRotateRealCoef(
            self.get_coeffs(normalization='4pi', csphase=1), angles, dj_matrix)

        # Convert 4pi normalized coefficients to the same normalization
        # as the unrotated coefficients.
        if self.normalization != '4pi' or csphase != 1:
            temp = SHCoeffs.from_array(coeffs, kind='real')
            tempcoeffs = tempcoeffs.get_coeffs(
                normalization=self.normalization, csphase=self.csphase)
            return SHCoeffs.from_array(
                tempcoeffs, normalization=self.normalization,
                csphase=self.csphase)
        else:
            return SHCoeffs.from_array(coeffs)

    def _expandDH(self, sampling):
        """
        Evaluate the coefficients on a Driscoll and Healy (1994)
        sampled grid.
        """
        if self.normalization == '4pi':
            norm = 1
        elif self.normalization == 'schmidt':
            norm = 2
        elif self.normalization == 'ortho':
            norm = 4
        else:
            raise ValueError(
                "Normalization must be '4pi', 'ortho', or 'schmidt'")

        data = _shtools.MakeGridDH(self.coeffs, sampling=sampling, norm=norm,
                                   csphase=self.csphase)
        gridout = SHGrid.from_array(data, grid='DH')
        return gridout

    def _expandGLQ(self, zeros):
        """Evaluate the coefficients on a Gauss Legendre quadrature grid."""
        if self.normalization == '4pi':
            norm = 1
        elif self.normalization == 'schmidt':
            norm = 2
        elif self.normalization == 'ortho':
            norm = 4
        else:
            raise ValueError(
                "Normalization must be '4pi', 'ortho', or 'schmidt'")

        if zeros is None:
            zeros, weights = _shtools.SHGLQ(self.lmax)

        data = _shtools.MakeGridGLQ(self.coeffs, zeros, norm=norm,
                                    csphase=self.csphase)
        gridout = SHGrid.from_array(data, grid='GLQ')
        return gridout


# =============== COMPLEX SPHERICAL HARMONICS ================

class SHComplexCoefficients(SHCoeffs):
    """
    Complex Spherical Harmonics Coefficients class.
    """

    @staticmethod
    def istype(kind):
        return kind == 'complex'

    def __init__(self, coeffs, normalization='4pi', csphase=1):
        lmax = coeffs.shape[1] - 1
        # ---- create mask to filter out m<=l ----
        mask = np.zeros((2, lmax + 1, lmax + 1), dtype=np.bool)
        mask[0, 0, 0] = True
        for l in np.arange(lmax + 1):
            mask[:, l, :l + 1] = True
        mask[1, :, 0] = False

        self.mask = mask
        self.lmax = coeffs.shape[1] - 1
        self.coeffs = np.copy(coeffs)
        self.coeffs[np.invert(mask)] = 0.
        self.kind = 'complex'
        self.normalization = normalization.lower()
        self.csphase = csphase

    def make_real(self, convention=1, switchcs=0):
        """
        Convert the complex coefficient class to the real harmonic
        coefficient class.
        """
        raise NotImplementedError('Not implemented yet!')

        # NOT CORRECT. ONLY WORKS IF THE GRID IS REAL!
        # First need to check that the grid is in fact real.
        complex_coeffs = _shtools.SHctor(self.coeffs, convention=convention,
                                         switchcs=switchcs)
        return SHCoeffs.from_array(complex_coeffs, kind='real')

    def _powerperdegree(self):
        """Return the power per degree l spectrum."""
        if self.normalization == '4pi':
            return _shtools.SHPowerSpectrumC(self.coeffs)

        elif self.normalization == 'schmidt':
            power = _shtools.SHPowerSpectrumC(self.coeffs)
            for l in range(self.lmax + 1):
                power[l] /= (2.0 * l + 1.0)
            return power

        elif self.normalization == 'ortho':
            return _shtools.SHCPowerSpectrumC(self.coeffs) / (4.0 * np.pi)

        else:
            raise ValueError(
                "Normalization must be '4pi', 'ortho', or 'schmidt'. " +
                "Input value was '{:s}'".format(str(self.normalization)))

    def _get_coeffs(self, output_normalization, output_csphase):
        """
        Return complext spherical harmonics coefficients with a
        different normalization convention.
        """
        coeffs = np.copy(self.coeffs)

        if self.normalization == output_normalization:
            pass
        elif (self.normalization == '4pi' and
              output_normalization == 'schmidt'):
            for l in self.get_degrees():
                coeffs[:, l, :l+1] *= np.sqrt(2.0 * l + 1.0)
        elif self.normalization == '4pi' and output_normalization == 'ortho':
            coeffs *= np.sqrt(4.0 * np.pi)
        elif (self.normalization == 'schmidt' and
              output_normalization == '4pi'):
            for l in self.get_degrees():
                coeffs[:, l, :l+1] /= np.sqrt(2.0 * l + 1.0)
        elif (self.normalization == 'schmidt' and
              output_normalization == 'ortho'):
            for l in self.get_degrees():
                coeffs[:, l, :l+1] *= np.sqrt(4.0 * np.pi / (2.0 * l + 1.0))
        elif self.normalization == 'ortho' and output_normalization == '4pi':
            coeffs /= np.sqrt(4.0 * np.pi)
        elif (self.normalization == 'ortho' and
              output_normalization == 'schmidt'):
            for l in self.get_degrees():
                coeffs[:, l, :l+1] *= np.sqrt((2.0 * l + 1.0) / (4.0 * np.pi))

        if output_csphase != self.csphase:
            for m in self.get_degrees():
                if m % 2 == 1:
                    coeffs[:, :, m] = - coeffs[:, :, m]

        return coeffs

    def _rotate(self, angles, dj_matrix):
        """
        Rotate the coordinate system used to express the spherical
        harmonics coefficients by the Euler angles alpha, beta, gamma.
        """
        if dj_matrix is None:
            dj_matrix = _shtools.djpi2(self.lmax + 1)

        raise NotImplementedError('Not implemented yet')
        self.coeffs = _shtools.SHRotateRealCoef(self.coeffs, angles, dj_matrix)

    def _expandDH(self, sampling):
        """
        Evaluate the coefficients on a Driscoll and Healy (1994)
        sampled grid.
        """
        if self.normalization == '4pi':
            norm = 1
        elif self.normalization == 'schmidt':
            norm = 2
        elif self.normalization == 'ortho':
            norm = 4
        else:
            raise NotImplementedError(
                "Normalization must be '4pi', 'ortho', or 'schmidt'")

        data = _shtools.MakeGridDHC(self.coeffs, sampling=sampling,
                                    norm=norm, csphase=self.csphase)
        gridout = SHGrid.from_array(data, grid='DH')
        return gridout

    def _expandGLQ(self, zeros):
        """Evaluate the coefficients on a Gauss Legendre quadrature grid."""
        if self.normalization == '4pi':
            norm = 1
        elif self.normalization == 'schmidt':
            norm = 2
        elif normalization == 'ortho':
            norm = 4
        else:
            raise NotImplementedError(
                "Normalization must be '4pi', 'ortho', or 'schmidt'")

        if zeros is None:
            zeros, weights = _shtools.SHGLQ(self.lmax)

        data = _shtools.MakeGridGLQ(self.coeffs, zeros, norm=norm,
                                    csphase=csphase)
        gridout = SHGrid.from_array(data, grid='GLQ')
        return gridout


# ========================================================================
# ======      GRID CLASSES      ==========================================
# ========================================================================

class SHGrid(object):
    """
    Grid Class for global gridded data on the sphere. Grids can be
    initialized from:

    >> x = SHGrid.from_array( numpy.ndarray )
    >> x = SHGrid.from_file( 'fname.dat' )
    """

    def __init__():
        pass

    # ---- factory methods
    @classmethod
    def from_array(self, array, grid='DH'):
        if np.iscomplexobj(array):
            kind = 'complex'
        else:
            kind = 'real'

        for cls in self.__subclasses__():
            if cls.istype(kind) and cls.isgrid(grid):
                return cls(array)

    @classmethod
    def from_file(self, fname, kind='real', grid='DH'):

        raise NotImplementedError('Not implemented yet')

        # need to open and read grid, and specify binary, ascii
        # nlat and nlong
        for cls in self.__subclasses__():
            if cls.istype(kind) and cls.isgrid(grid):
                return cls(grid)

    # ---- Extract grid properties ----
    def get_lats():
        """
        Return a vector containing the latitudes (in degrees) of each row
        of the gridded data.
        """
        return self._get_lats()

    def get_lons():
        """
        Return a vector containing the longitudes (in degrees) of each
        column of the gridded data.
        """
        return self._get_lats()

    # ---- Plotting routines ----
    def plot_rawdata(self, show=True, fname=None):
        """
        Plot the raw data using a simply cylindrical projection.
        """
        fig, ax = self._plot_rawdata()
        if show:
            plt.show()
        if fname is not None:
            fig.savefig(fname)

    def expand(self):
        """Expand the grid into spherical harmonics."""
        return self._expand(normalization='4pi', csphase=1)


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

        if self.nlat == 2 * self.nlon:
            self.sampling = 2
        elif self.nlat == self.nlon:
            self.sampling = 1
        else:
            raise ValueError('Input array has shape (nlat={:d},nlon={:d})\n' +
                             'but needs nlat=nlon or nlat=2*nlon'
                             .format(self.nlat, self.nlon)
                             )

        self.data = array
        self.grid = 'DH'
        self.kind = 'real'

    def _getlats(self):
        """
        Return a vector containing the latitudes (in degrees) of each row
        of the gridded data.
        """
        lats = np.linspace(90., -90.+180./self.nlat, num=self.nlat)
        return lats

    def _getlons(self):
        """
        Return a vector containing the longitudes (in degrees) of each row
        of the gridded data.
        """
        lons = np.linspace(0., 360.-360./self.nlon, num=self.nlon)
        return lons

    def _expand(self, normalization, csphase):
        """Expand the grid into real spherical harmonics."""
        if normalization == '4pi':
            norm = 1
        elif normalization == 'schmidt' or normalization == 'Schmidt':
            norm = 2
        elif normalization == 'ortho':
            norm = 4
        else:
            raise NotImplementedError(
                "Normalization must be '4pi', 'ortho', or 'schmidt'")

        cilm = _shtools.SHExpandDH(self.data, norm=norm, csphase=csphase)
        coeffs = SHCoeffs.from_array(cilm, kind='real',
                                     normalization=normalization,
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

        if self.nlat == 2 * self.nlon:
            self.sampling = 2
        elif self.nlat == self.nlon:
            self.sampling = 1
        else:
            raise ValueError('Input array has shape (nlat={:d},nlon={:d})\n' +
                             'but needs nlat=nlon or nlat=2*nlon'
                             .format(self.nlat, self.nlon)
                             )

        self.data = array
        self.grid = 'DH'
        self.kind = 'complex'

    def _getlats(self):
        """
        Return a vector containing the latitudes (in degrees) of each row
        of the gridded data.
        """
        lats = np.linspace(90., -90.+180./self.nlat, num=self.nlat)
        return lats

    def _getlons(self):
        """
        Return a vector containing the longitudes (in degrees) of each row
        of the gridded data.
        """
        lons = np.linspace(0., 360.-360./self.nlon, num=self.nlon)
        return lons

    def _expand(self, normalization, csphase):
        """Expand the grid into real spherical harmonics."""
        if normalization == '4pi':
            norm = 1
        elif normalization == 'schmidt' or normalization == 'Schmidt':
            norm = 2
        elif normalization == 'ortho':
            norm = 4
        else:
            raise NotImplementedError(
                "Normalization must be '4pi', 'ortho', or 'schmidt'")

        cilm = _shtools.SHExpandDHC(self.data, norm=norm, csphase=csphase)
        coeffs = SHCoeffs.from_array(cilm, kind='complex',
                                     normalization=normalization,
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
    Class for real Gauss Legendre Quadrature grids.
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

        if zeros is None and weights is None:
            self.zeros, weights = _shtools.SHGLQ(self.lmax)
        else:
            self.zeros = zeros
            self.weights.weights

        self.data = array
        self.grid = 'GLQ'
        self.kind = 'real'

    def _getlats(self):
        """
        Return a vector containing the latitudes (in degrees) of each row
        of the gridded data.
        """
        lats = 90. - np.arccos(self.zeros) * 180. / np.pi
        return lats

    def _getlons(self):
        """
        Return a vector containing the longitudes (in degrees) of each column
        of the gridded data.
        """
        lons = np.linspace(0., 360.-360./self.nlon, num=self.nlon)
        return lons

    def _expand(self, normalization, csphase):
        """Expand the grid into real spherical harmonics."""
        if normalization == '4pi':
            norm = 1
        elif normalization == 'schmidt' or normalization == 'Schmidt':
            norm = 2
        elif normalization == 'ortho':
            norm = 4
        else:
            raise NotImplementedError(
                "Normalization must be '4pi', 'ortho', or 'schmidt'")

        cilm = _shtools.SHExpandGLQ(self.data, self.weights, self.zeros,
                                    norm=norm, csphase=csphase)
        coeffs = SHCoeffs.from_array(cilm, kind='real',
                                     normalization=normalization,
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

        if zeros is None and weights is None:
            self.zeros, weights = _shtools.SHGLQ(self.lmax)
        else:
            self.zeros = zeros
            self.weights.weights

        self.data = array
        self.grid = 'GLQ'
        self.kind = 'complex'

    def _getlats(self):
        """
        Return a vector containing the latitudes (in degrees) of each row
        of the gridded data.
        """
        lats = 90. - np.arccos(self.zeros) * 180. / np.pi
        return lats

    def _getlons(self):
        """
        Return a vector containing the longitudes (in degrees) of each column
        of the gridded data.
        """
        lons = np.linspace(0., 360.-360./self.nlon, num=self.nlon)
        return lons

    def _expand(self, normalization, csphase):
        """Expand the grid into real spherical harmonics."""
        if normalization == '4pi':
            norm = 1
        elif normalization == 'schmidt' or normalization == 'Schmidt':
            norm = 2
        elif normalization == 'ortho':
            norm = 4
        else:
            raise NotImplementedError(
                "Normalization must be '4pi', 'ortho', or 'schmidt'")

        cilm = _shtools.SHExpandGLQC(self.data, self.weights, self.zeros,
                                     norm=norm, csphase=csphase)
        coeffs = SHCoeffs.from_array(cilm, kind='real',
                                     normalization=normalization,
                                     csphase=csphase)
        return coeffs

    def _plot_rawdata(self):
        """Plot the raw data using a simply cylindrical projection."""

        fig, ax = plt.subplots(2, 1)
        ax.flat[0].imshow(self.data, origin='top')
        ax.flat[0].set_title('Gauss-Legendre Quadrature Grid (real component')
        ax.flat[0].set_xlabel('longitude index')
        ax.flat[0].set_ylabel('latitude index')
        ax.flat[1].imshow(self.data, origin='top')
        ax.flat[1].set_title('Gauss-Legendre Quadrature Grid ' +
                             '(imaginary component')
        ax.flat[1].set_xlabel('longitude index')
        ax.flat[1].set_ylabel('latitude index')
        fig.tight_layout(pad=0.5)
        return fig, ax


# ==== SPHERICAL HARMONICS WINDOW FUNCTION CLASS ====
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
        return _shtools.SHSymmetricWindow(tapers, eigenvalues, taper_order,
                                          clat=clat, clon=clon)

    @classmethod
    def from_mask(self, lmax, nwins, dh_mask, sampling=1):
        """
        constructs optimal window functions in a masked region (needs dh grid)
        """
        tapers, eigenvalues = _shtools.SHReturnTapersMap(
            dh_mask, lmax, sampling=sampling, Ntapers=nwins)
        return _shtools.SHAsymmetricWindow(tapers, eigenvalues)

    def plot(self, nwins, show=True, fname=None):
        """
        plots the best concentrated spherical harmonics taper functions
        """
        # ---- setup figure and axes
        maxcolumns = 5
        ncolumns = min(maxcolumns, nwins)
        nrows = np.ceil(nwins / ncolumns).astype(int)
        figsize = ncolumns * 1.2, nrows * 1.2 + 0.5
        fig, axes = plt.subplots(nrows, ncolumns, figsize=figsize)
        for ax in axes[:-1, :].flatten():
            for xlabel_i in ax.get_xticklabels():
                xlabel_i.set_visible(False)
        for ax in axes[:, 1:].flatten():
            for ylabel_i in ax.get_yticklabels():
                ylabel_i.set_visible(False)

        # loop through tapers and plot them
        for itaper in range(min(self.nwins, nwins)):
            evalue = self.eigenvalues[itaper]
            coeffs = self._coeffs(itaper)
            ax = axes.flatten()[itaper]
            grid = MakeGridDH(coeffs)
            ax.imshow(grid)
            ax.set_title('concentration: {:2.2f}'.format(evalue))
        fig.tight_layout(pad=0.5)

        if show:
            plt.show()
        if fname is not None:
            fig.savefig(fname)

    def get_spectrum(self, shcoeffs, nwins):
        """Returns the regional spherical harmonics spectrum"""
        for itaper in range(nwins):
            tapercoeffs = self._coeffs(itaper)
            modelcoeffs = shcoeffs.get_coeffs(normalization='4pi', kind='real')
            coeffs = _shtools.SHMultiply(tapercoeffs, modelcoeffs)

    def get_couplingmatrix(self, lmax, nwins):
        """returns the coupling matrix of the first nwins tapers"""
        # store sqrt of taper power in 'tapers' array:
        if nwins > self.nwins:
            nwins = self.nwins
        tapers = np.zeros((self.nl, nwins))
        for itaper in range(nwins):
            tapers[:, itaper] = np.sqrt(_shtools.SHPowerSpectrum(
                self._coeffs(itaper)))

        # compute coupling matrix of the first nwins tapers:
        coupling_matrix = _shtools.SHMTCouplingMatrix(lmax, tapers[:, :nwins])
        return coupling_matrix

    def plot_couplingmatrix(self, lmax, nwins, show=True, fname=None):
        """plots the window's coupling strength"""
        figsize = mpl.rcParams['figure.figsize']
        figsize[0] = figsize[1]
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        coupling_matrix = self.get_couplingmatrix(lmax, nwins)
        ax.imshow(coupling_matrix)
        ax.set_xlabel('output power')
        ax.set_ylabel('input power')
        fig.tight_layout(pad=0.1)

        if show:
            plt.show()
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

    def __init__(self, tapers, eigenvalues, orders, clat=0., clon=0.):
        # center of cap window
        self.clat, self.clon = clat, clon
        # nl: number of degrees, nwins: number of windows
        self.nl, self.nwins = tapers.shape
        # lmax: maximum degree
        self.lmax = self.nl - 1
        # tapers[nl,nwins]: ith window coefs with m=orders[iwin]
        self.tapers = tapers
        # concentration factor of the ith taper
        self.eigenvalues = eigenvalues
        # order m of the ith taper
        self.orders = orders

    def _coeffs(self, itaper):
        taperm = self.orders[itaper]
        coeffs = np.zeros((2, self.nl, self.nl))
        if taperm < 0:
            coeffs[1, :, abs(taperm)] = self.tapers[:, itaper]
        else:
            coeffs[0, :, abs(taperm)] = self.tapers[:, itaper]
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

    def __init__(self, tapers, eigenvalues):
        ncoeffs, self.nwins = tapers.shape
        self.nl = np.sqrt(ncoeffs).astype(int)
        self.lmax = self.nl-1
        self.tapers = tapers
        self.eigenvalues = eigenvalues

    def _coeffs(self, itaper):
        return _shtools.SHVectorToCilm(self.tapers[:, itaper], self.lmax)

    def _info(self):
        print('Asymmetric window with {:d} tapers'.format(self.nwins))
