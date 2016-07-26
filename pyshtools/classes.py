"""
pyshtools defines several classes that facilitate the interactive
examination of geographical gridded data and their associated
spherical harmonic coefficients. Subclasses are used to handle different
internal data types and superclasses are used to implement interface
functions and the documentation.

For more information, see the documentation for the following classes
and subclasses:

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
    Spherical Harmonics Coefficient class.

    The coefficients of this class can be initialized using one of the
    three constructor methods:

    >> x = SHCoeffs.from_array(numpy.zeros((2, lmax+1, lmax+1)))
    >> x = SHCoeffs.from_random(powerspectrum[0:lmax+1])
    >> x = SHCoeffs.from_file('fname.dat')

    The normalization convention of the input coefficents is specified
    by the normalization and csphase parameters, which can take the following
    values:

    normalization : '4pi' (default), geodesy 4-pi normalized
                  : 'ortho', orthonormalized
                  : 'schmidt', Schmidt semi-normalized

    csphase       : 1 (default), exlcude the Condon-Shortley phase factor
                  : -1, include the Condon-Shortley phase factor

    See the documentation for each constructor method for further options.

    The class instance defines the following class attributes:

    lmax          : The maximum spherical harmonic degree of the
                    coefficients.
    coeffs        : The raw coefficients with the specified normalization and
                    phase conventions.
    mask          : A boolean mask that is True for the allowable values of
                    degree l and order m.
    kind          : Either 'complex' or 'real' for the coefficient data type.
    normalization : The normalization of the coefficients, which is '4pi',
                    'ortho', or 'schmidt'.
    csphse        : Defines whether the Condon-Shortley phase is used (1)
                    or not (-1).

    Each class instance also provides the following methods:

    get_degrees()         : Return an array listing the spherical harmonic
                            degrees from 0 to lmax.
    get_powerperdegree()  : Return the power per degree l spectrum.
    get_powerperband()    : Return the power per log_{bandwidth} l spectrum.
    get_coeffs()          : Return spherical harmonics coefficients with
                            a different normalization convention.
    rotate()              : Rotate the coordinate system used to express the
                            spherical harmonics coefficients.
    expand()              : Evaluate the coefficients on a spherical grid.
    plot_powerperdegree() : Plot the power per degree spectrum.
    plot_powerperband()   : Plot the power per log_{bandwidth} l spectrum.
    make_complex()        : Convert real coefficients to complex form.
    make_real()           : Convert complex coefficients to real form.
    """

    def __init__(self):
        pass

    # ---- factory methods:
    @classmethod
    def from_array(self, coeffs, normalization='4pi', csphase=1):
        """
        Initialize the spherical harmonic coefficients of the object from
        an input array.

        Usage
        -----

        x = SHCoeffs.from_array(array, [normalization, csphase])

        Parameters
        ----------

        array         : numpy array of size (2, lmax+1, lmax+1)
        normalization : '4pi' (default), geodesy 4-pi normlized
                      : 'ortho', orthonormalized
                      : 'schmidt', Schmidt semi-normalized
        csphase       : 1 (default), exlcude the Condon-Shortley phase factor
                      : -1, include the Condon-Shortley phase factor
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

        Usage
        -----

        x = SHCoeffs.from_random(power, [kind, normalization, csphase])

        Parameters
        ----------

        power         : numpy array of the power spectrum of size (lmax+1)
        kind          : 'real' (default) or 'complex' output coefficients
        normalization : '4pi' (default), geodesy 4-pi normalized
                      : 'ortho', orthonormalized
                      : 'schmidt', Schmidt semi-normalized
        csphase       : 1 (default) exlcude the Condon-Shortley phase factor
                      : -1, include the Condon-Shortley phase factor
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

        if kind.lower() not in set(['real', 'complex']):
            raise ValueError(
                "kind must be 'real' or 'complex'. " +
                "Input value was '{:s}'.".format(str(kind)))

        nl = len(power)
        l = np.arange(nl)

        if kind.lower() == 'real':
            coeffs = np.random.normal(size=(2, nl, nl))
        elif kind.lower() == 'complex':
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
                  normalization='4pi', csphase=1, **kwargs):
        """
        Initialize the spherical harmonic coefficients by reading the
        coefficients from a specified file.

        Usage
        -----

        x = SHCoeffs.from_file(filename, lmax, [format, kind, normalization,
                                                csphase, skip])

        Parameters
        ----------

        filename      : name of the file
        lmax          : maximum spherical harmonic degree to read from the file
        format        : 'shtools' (default)
        kind          : Output 'real' (default) or 'complex' coefficients
        normalization : '4pi' (default), geodesy 4-pi normalized
                      : 'ortho', orthonormalized
                      : 'schmidt', Schmidt semi-normalized)
        csphase       : 1  (default), exlcude the Condon-Shortley phase factor
                      : -1, include the Condon-Shortley phase factor
        skip          : Number of lines to skip at the beginning of the file

        Description
        -----------

        If format='shtools' spherical harmonic coefficients from an
        ascii-formatted file. The maximum spherical harmonic degree that is
        read is determined by the input value lmax. If the optional value skip
        is specified, parsing of the file will commence after the first skip
        lines.

        The spherical harmonic coefficients in the file are assumed to be
        ordered by increasing degree l and angular order m according to the
        format

        l, m, cilm[0, l, m], cilm[1, l, m]

        For each value of increasing l, all the angular orders are listed
        in inceasing order.
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
                coeffs, lmax = _shtools.SHRead(fname, lmax, **kwargs)
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

        Usage
        -----

        degrees = x.get_degrees()

        Returns
        -------

        degrees : ndarray of size (lmax+1)
        """
        return np.arange(self.lmax + 1)

    def get_powerperdegree(self):
        """
        Return the power per degree l spectrum.

        Usage
        -----

        power = x.get_powerperdegree()

        Returns
        -------

        power : ndarray of size (lmax+1)
        """
        return self._powerperdegree()

    def get_powerperband(self, bandwidth):
        """
        Return the power per log_{bandwidth} l spectrum.

        Usage
        -----

        power = x.get_powerperband()

        Returns
        -------

        power : ndarray of size (lmax+1)
        """
        ls = self.get_degrees()
        return self._powerperdegree() * ls * np.log(bandwidth)

    # ---- Return coefficients with a different normalization convention ----
    def get_coeffs(self, normalization='4pi', csphase=1):
        """
        Return spherical harmonics coefficients with a different
        normalization convention.

        Usage
        -----

        power = x.get_coeffs([normalization, csphase])

        Returns
        -------

        power : ndarray of size (lmax+1)

        Parameters
        ----------

        normalization : '4pi' (default), geodesy 4-pi normalized
                      : 'ortho', orthonormalized
                      : 'schmidt', Schmidt semi-normalized)
        csphase       : 1  (default), exlcude the Condon-Shortley phase factor
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

        Usage
        -----

        SHCoeffsInstance = x.rotate(alpha, beta, gamma, [degrees, dj_matrix])

        Parameters
        ----------

        alpha, beta, gamma : Euler rotation angles in degrees.
        degrees            : True (default) to use degrees, False for radians
        dj_matrix          : The djpi2 rotation matrix (default=None)
        """
        if degrees:
            angles = np.radians([alpha, beta, gamma])
        else:
            angles = np.array([alpha, beta, gamma])

        rot = self._rotate(angles, dj_matrix)
        return rot

    # ---- Expand the coefficients onto a grid ----
    def expand(self, grid='DH', **kwargs):
        """
        Evaluate the coefficients on a spherical grid.

        Usage
        -----

        SHGridInstance = x.expand([grid, lmax, lmax_calc, zeros])

        Parameters
        ----------

        grid      : 'DH' or 'DH1', equisampled lat/lon grid with nlat=nlon
                  : 'DH2', equidistant lat/lon grid with nlon=2*nlat
                  : 'GLQ', Gauss-Legendre quadrature grid
        lmax      : maximum spherical harmonic degree, which determines the
                    grid spacing of the output grid. Default is x.lmax.
        lmax_calc : maximum spherical harmonic degree to use when evaluating
                    the function. Default is x.lmax.
        zeros      : The cos(colatitude) nodes used in Gauss-Legendre
                    Quadrature grids. Default is None
        """
        if grid.upper() == 'DH' or grid.upper() == 'DH1':
            gridout = self._expandDH(sampling=1, **kwargs)
        elif grid.upper() == 'DH2':
            gridout = self._expandDH(sampling=2, **kwargs)
        elif grid.upper() == 'GLQ':
            gridout = self._expandGLQ(zeros=None, **kwargs)
        else:
            raise ValueError(
                "grid must be 'DH', 'DH1', 'DH2', or 'GLQ'. " +
                "Input value was '{:s}'".format(grid))

        return gridout

    # ---- plotting routines ----
    def plot_powerperdegree(self, loglog=True, show=True, fname=None):
        """
        Plot the power per degree spectrum.

        Usage
        -----

        x.plot_powerperdegree([loglog, show, fname])

        Parameters
        ----------

        loglog : if True (default), use log-log axis
        show   : if True (default), plot to the screen.
        fname  : if present, image will be saved to the file
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

        Usage
        -----

        x.plot_powerperdegree([loglog, show, fname])

        Parameters
        ----------

        loglog : if True (default), use log-log axis
        show   : if True (default), plot to the screen.
        fname  : if present, image will be saved to the file
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

    def make_complex(self):
        """
        Convert the real coefficient class to the complex harmonic
        coefficient class with the same normalization and csphase
        conventions.
        """
        rcomplex_coeffs = _shtools.SHrtoc(self.coeffs,
                                          convention=1, switchcs=0)

        # These coefficients are using real floats, and need to be
        # converted to complex form.
        complex_coeffs = np.zeros((2, self.lmax+1, self.lmax+1),
                                  dtype='complex')
        complex_coeffs[0, :, :] = (rcomplex_coeffs[0, :, :] + 1j *
                                   rcomplex_coeffs[1, :, :])
        complex_coeffs[1, :, :] = complex_coeffs[0, :, :].conjugate()
        for m in self.get_degrees():
            if m % 2 == 1:
                complex_coeffs[1, :, m] = - complex_coeffs[1, :, m]

        return SHCoeffs.from_array(complex_coeffs,
                                   normalization=self.normalization,
                                   csphase=self.csphase)

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
        if self.normalization != '4pi' or self.csphase != 1:
            temp = SHCoeffs.from_array(coeffs, kind='real')
            tempcoeffs = temp.get_coeffs(
                normalization=self.normalization, csphase=self.csphase)
            return SHCoeffs.from_array(
                tempcoeffs, normalization=self.normalization,
                csphase=self.csphase)
        else:
            return SHCoeffs.from_array(coeffs)

    def _expandDH(self, sampling, **kwargs):
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
                                   csphase=self.csphase, **kwargs)
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
                                    csphase=self.csphase, **kwargs)
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
        self.lmax = lmax
        self.coeffs = np.copy(coeffs)
        self.coeffs[np.invert(mask)] = 0.
        self.kind = 'complex'
        self.normalization = normalization.lower()
        self.csphase = csphase

    def make_real(self):
        """
        Convert the complex coefficient class to the real harmonic
        coefficient class.
        """
        # First test if the coefficients correspond to a real grid.
        # This is not very elegant. Also, the equality condition for
        # is probably not robust to round off errors.
        for l in self.get_degrees():
            if self.coeffs[0, l, 0] != self.coeffs[0, l, 0].conjugate():
                raise Error('Complex coefficients do not correspond to a ' +
                            'real field. l={:d}, m=0: {:c}'
                            .format(l, self.coeffs[0, l, 0]))
            for m in np.arange(1, l+1):
                if m % 2 == 1:
                    if (self.coeffs[0, l, m] != -
                            self.coeffs[1, l, m].conjugate()):
                        raise Error('Complex coefficients do not correspond ' +
                                    'to a real field. ' +
                                    'l={:d}, m={:d}: {:c}, {:c}'
                                    .format(l, m, self.coeffs[0, l, 0],
                                            self.coeffs[1, l, 0]))
                else:
                    if (self.coeffs[0, l, m] !=
                            self.coeffs[1, l, m].conjugate()):
                        raise Error('Complex coefficients do not correspond ' +
                                    'to a real field. ' +
                                    'l={:d}, m={:d}: {:c}, {:c}'
                                    .format(l, m, self.coeffs[0, l, 0],
                                            self.coeffs[1, l, 0]))

            coeffs_rc = np.zeros((2, self.lmax + 1, self.lmax + 1))
            coeffs_rc[0, :, :] = self.coeffs[0, :, :].real
            coeffs_rc[1, :, :] = self.coeffs[0, :, :].imag
            real_coeffs = _shtools.SHctor(coeffs_rc, convention=1,
                                          switchcs=0)
            return SHCoeffs.from_array(real_coeffs)

    def _powerperdegree(self):
        """Return the power per degree l spectrum."""
        if self.normalization == '4pi':
            return _shtools.SHPowerSpectrumC(self.coeffs)
        elif self.normalization == 'schmidt':
            power = _shtools.SHPowerSpectrumC(self.coeffs)
            l = self.get_degrees()
            power /= (2.0 * l + 1.0)
            return power
        elif self.normalization == 'ortho':
            return _shtools.SHPowerSpectrumC(self.coeffs) / (4.0 * np.pi)
        else:
            raise ValueError(
                "Normalization must be '4pi', 'ortho', or 'schmidt'. " +
                "Input value was '{:s}'".format(str(self.normalization)))

    def _get_coeffs(self, output_normalization, output_csphase):
        """
        Return complex spherical harmonics coefficients with a
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
        # Note that the current method is EXTREMELY inefficient. The complex
        # coefficients are expanded onto real and imaginary grids, each of
        # the two components are rotated separately as real data, they rotated
        # real data are re-expanded on new real and complex grids, they are
        # combined to make a complex grid, and the resultant is expanded
        # in complex spherical harmonics.
        if dj_matrix is None:
            dj_matrix = _shtools.djpi2(self.lmax + 1)

        cgrid = self.expand(grid='DH')
        rgrid, igrid = cgrid.data.real, cgrid.data.imag
        rgridcoeffs = _shtools.SHExpandDH(rgrid, norm=1, sampling=1, csphase=1)
        igridcoeffs = _shtools.SHExpandDH(igrid, norm=1, sampling=1, csphase=1)

        rgridcoeffs_rot = _shtools.SHRotateRealCoef(
            rgridcoeffs, angles, dj_matrix)
        igridcoeffs_rot = _shtools.SHRotateRealCoef(
            igridcoeffs, angles, dj_matrix)

        rgrid_rot = _shtools.MakeGridDH(rgridcoeffs_rot, norm=1,
                                        sampling=1, csphase=1)
        igrid_rot = _shtools.MakeGridDH(igridcoeffs_rot, norm=1,
                                        sampling=1, csphase=1)
        grid_rot = rgrid_rot + 1j * igrid_rot

        if self.normalization == '4pi':
            norm = 1
        elif self.normalization == 'schmidt':
            norm = 2
        elif self.normalization == 'ortho':
            norm = 4
        else:
            raise ValueError(
                "Normalization must be '4pi', 'ortho', or 'schmidt'")

        coeffs_rot = _shtools.SHExpandDHC(grid_rot, norm=norm,
                                          csphase=self.csphase)

        return SHCoeffs.from_array(coeffs_rot,
                                   normalization=self.normalization,
                                   csphase=self.csphase)

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
                                    norm=norm, csphase=self.csphase, **kwargs)
        gridout = SHGrid.from_array(data, grid='DH')
        return gridout

    def _expandGLQ(self, zeros):
        """Evaluate the coefficients on a Gauss-Legendre quadrature grid."""
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
                                    csphase=self.csphase, **kwargs)
        gridout = SHGrid.from_array(data, grid='GLQ')
        return gridout


# ========================================================================
# ======      GRID CLASSES      ==========================================
# ========================================================================

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
    grid       : 'DH' or 'GLQ' for Driscoll and Healy grids or Gauss
                 Legendre quadrature grids.
    zeros      : The cos(colatitude) nodes used with Gauss-Legendre
                 Quadrature grids. Default is None.
    weights    : The latitudinal weights used with Gauss-Legendre
                 Quadrature grids. Default is None.

    Each class instance also provides the following methods:

    get_lats()     : Return a vector containing the latitudes of each row
                     of the gridded data.
    get_lons()     : Return a vector containing the longitudes of each row
                     of the gridded data.
    plot_rawdata() : Plot the raw data using a simply cylindrical projection.
    expand()       : Expand the grid into spherical harmonics.
    """

    def __init__():
        pass

    # ---- factory methods
    @classmethod
    def from_array(self, array, grid='DH'):
        """
        Initialize the grid of the object from an input array.

        Usage
        -----

        x = SHGrid.from_array(array, [grid])

        Parameters
        ----------

        array : numpy array of size (nlat, nlon)
        grid : 'DH' (default) or 'GLQ' for Driscoll and Healy grids or Gauss
                Legendre quadrature grids, respectively.
        """
        if np.iscomplexobj(array):
            kind = 'complex'
        else:
            kind = 'real'

        for cls in self.__subclasses__():
            if cls.istype(kind) and cls.isgrid(grid):
                return cls(array)

    @classmethod
    def from_file(self, fname, kind='real', grid='DH'):
        """Initialize the grid of the object from a file."""
        raise NotImplementedError('Not implemented yet')

    # ---- Extract grid properties ----
    def get_lats():
        """
        Return a vector containing the latitudes (in degrees) of each row
        of the gridded data.

        Usage
        -----

        lats = SHGrid.get_lats()

        Returns
        -------

        lats : numpy array of size nlat containing the latitude (in degrees)
               of each row of the gridded data.
        """
        return self._get_lats()

    def get_lons():
        """
        Return a vector containing the longitudes (in degrees) of each
        column of the gridded data.

        Usage
        -----

        lons = SHGrid.get_lon()

        Returns
        -------

        lons : numpy array of size nlon containing the longitude (in degrees)
               of each column of the gridded data.
        """
        return self._get_lats()

    # ---- Plotting routines ----
    def plot_rawdata(self, show=True, fname=None):
        """
        Plot the raw data using a simply cylindrical projection.

        Usage
        -----

        x.plot_rawdata([show, fname])

        Parameters
        ----------

        show   : if True (default), plot to the screen
        fname  : if present, image will be saved to the file
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

        self.lmax = self.nlat / 2 -1
        self.data = array
        self.grid = 'DH'
        self.kind = 'real'

    def _getlats(self):
        """
        Return a vector containing the latitudes (in degrees) of each row
        of the gridded data.
        """
        lats = np.linspace(90.0, -90.0 + 180.0 / self.nlat, num=self.nlat)
        return lats

    def _getlons(self):
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
            raise NotImplementedError(
                "Normalization must be '4pi', 'ortho', or 'schmidt'")

        cilm = _shtools.SHExpandDH(self.data, norm=norm, csphase=csphase,
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

        self.lmax = self.nlat / 2 -1
        self.data = array
        self.grid = 'DH'
        self.kind = 'complex'

    def _getlats(self):
        """
        Return a vector containing the latitudes (in degrees) of each row
        of the gridded data.
        """
        lats = np.linspace(90.0, -90.0 + 180.0 / self.nlat, num=self.nlat)
        return lats

    def _getlons(self):
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
            raise NotImplementedError(
                "Normalization must be '4pi', 'ortho', or 'schmidt'")

        cilm = _shtools.SHExpandDHC(self.data, norm=norm, csphase=csphase,
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
            self.zeros, self.weights = _shtools.SHGLQ(self.lmax)
        else:
            self.zeros = zeros
            self.weights = weights

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
            raise NotImplementedError(
                "Normalization must be '4pi', 'ortho', or 'schmidt'")

        cilm = _shtools.SHExpandGLQ(self.data, self.weights, self.zeros,
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
            self.zeros, self.weights = _shtools.SHGLQ(self.lmax)
        else:
            self.zeros = zeros
            self.weights = weights

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
            raise NotImplementedError(
                "Normalization must be '4pi', 'ortho', or 'schmidt'")

        cilm = _shtools.SHExpandGLQC(self.data, self.weights, self.zeros,
                                     norm=norm, csphase=csphase, **kwargs)
        coeffs = SHCoeffs.from_array(cilm,
                                     normalization=normalization.lower(),
                                     csphase=csphase)
        return coeffs

    def _plot_rawdata(self):
        """Plot the raw data using a simply cylindrical projection."""

        fig, ax = plt.subplots(2, 1)
        ax.flat[0].imshow(self.data.real, origin='top')
        ax.flat[0].set_title('Gauss-Legendre Quadrature Grid (real component')
        ax.flat[0].set_xlabel('longitude index')
        ax.flat[0].set_ylabel('latitude index')
        ax.flat[1].imshow(self.data.imag, origin='top')
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
