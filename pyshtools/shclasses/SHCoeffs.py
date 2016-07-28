# =============================================================================
# =========    COEFFICIENT CLASSES    =========================================
# =============================================================================

from __future__ import absolute_import, division, print_function

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from .. import _SHTOOLS as shtools
from .SHGrid import SHGrid

class SHCoeffs(object):
    """
    Spherical Harmonics Coefficient class.

    The coefficients of this class can be initialized using one of the
    three constructor methods:

    >> x = SHCoeffs.from_array(numpy.zeros((2, lmax+1, lmax+1)))
    >> x = SHCoeffs.from_random(powerspectrum[0:lmax+1])
    >> x = SHCoeffs.from_file('fname.dat')

    The normalization convention of the input coefficents is specified
    by the normalization and csphase parameters, which take the following
    values:

    normalization : '4pi' (default), geodesy 4-pi normalized.
                  : 'ortho', orthonormalized.
                  : 'schmidt', Schmidt semi-normalized.

    csphase       : 1 (default), exlcude the Condon-Shortley phase factor.
                  : -1, include the Condon-Shortley phase factor.

    See the documentation for each constructor method for further options.

    Once initialized, each class instance defines the following class
    attributes:

    lmax          : The maximum spherical harmonic degree of the coefficients.
    coeffs        : The raw coefficients with the specified normalization and
                    phase conventions.
    normalization : The normalization of the coefficients: '4pi', 'ortho', or
                    'schmidt'.
    csphse        : Defines whether the Condon-Shortley phase is used (1)
                    or not (-1).
    mask          : A boolean mask that is True for the permissible values of
                    degree l and order m.
    kind          : The coefficient data type: either 'complex' or 'real'.

    Each class instance provides the following methods:

    get_degrees()         : Return an array listing the spherical harmonic
                            degrees from 0 to lmax.
    get_powerperdegree()  : Return an array with the power per degree spectrum.
    get_powerperband()    : Return an array with the power per log_{bandwidth}
                            spectrum.
    get_coeffs()          : Return an array of spherical harmonics coefficients
                            with a different normalization convention.
    rotate()              : Rotate the coordinate system used to express the
                            spherical harmonics coefficients and return a new
                            class instance.
    convert()             : Convert the spherical harmonic coefficients to a
                            different normalization and return a new class
                            instance.
    expand()              : Evaluate the coefficients on a spherical grid and
                            return a new SHGrid class instance.
    make_complex()        : Convert a real SHCoeffs class instance to a complex
                            class instance.
    make_real()           : Convert a complex SHCoeffs class instance to a real
                            class instance.
    plot_powerperdegree() : Plot the power per degree spectrum.
    plot_powerperband()   : Plot the power per log_{bandwidth}(degree)
                            spectrum.
    """

    def __init__(self):
        pass

    # ---- factory methods:
    @classmethod
    def from_array(self, coeffs, normalization='4pi', csphase=1):
        """
        Initialize the spherical harmonic coefficients of the class instance
        from an input numpy array.

        Usage
        -----

        x = SHCoeffs.from_array(array, [normalization, csphase])

        Parameters
        ----------

        array         : numpy array of size (2, lmax+1, lmax+1).
        normalization : '4pi' (default), 'ortho' or 'schmidt' for geodesy 4pi
                        normalized, orthonormalized, or Schmidt semi-normalized
                        coefficients, respectively.
        csphase       : 1 (default) if the coefficients exclude the Condon-
                        Shortley phase factor, or -1 if they include it.
        """
        if np.iscomplexobj(coeffs):
            kind = 'complex'
        else:
            kind = 'real'

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

        for cls in self.__subclasses__():
            if cls.istype(kind):
                return cls(coeffs, normalization=normalization.lower(),
                           csphase=csphase)

    @classmethod
    def from_random(self, power, kind='real', normalization='4pi', csphase=1):
        """
        Initialize the spherical harmonic coefficients of the class instance
        using Gaussian random variables with a given power spectrum.

        Usage
        -----

        x = SHCoeffs.from_random(power, [kind, normalization, csphase])

        Parameters
        ----------

        power         : numpy array of the power spectrum of size (lmax+1).
        kind          : 'real' (default) or 'complex' output coefficients.
        normalization : '4pi' (default), 'ortho' or 'schmidt' for geodesy 4pi
                        normalized, orthonormalized, or Schmidt semi-normalized
                        coefficients, respectively.
        csphase       : 1 (default) if the coefficients exclude the Condon-
                        Shortley phase factor, or -1 if they include it.
        """
        if type(normalization) != str:
            raise ValueError('normalization must be a string. ' +
                             'Input type was {:s}'
                             .format(str(type(normalization))))

        if normalization.lower() not in set(['4pi', 'ortho', 'schmidt']):
            raise ValueError(
                "The input normalization must be '4pi', 'ortho' " +
                "or 'schmidt'. Provided value was {:s}"
                .format(repr(normalization))
                )

        if csphase != 1 and csphase != -1:
            raise ValueError(
                "csphase must be 1 or -1. Input value was {:s}"
                .format(repr(csphase))
                )

        if kind.lower() not in set(['real', 'complex']):
            raise ValueError(
                "kind must be 'real' or 'complex'. " +
                "Input value was {:s}.".format(repr(kind)))

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
        Initialize the spherical harmonic coefficients of the class instance
        by reading the coefficients from a specified file.

        Usage
        -----

        x = SHCoeffs.from_file(filename, lmax, [format, kind, normalization,
                                                csphase, skip])

        Parameters
        ----------

        filename      : Name of the file, including path.
        lmax          : Maximum spherical harmonic degree to read from the
                        file.
        format        : 'shtools' (default).
        kind          : Output 'real' (default) or 'complex' coefficients.
        normalization : '4pi' (default), 'ortho' or 'schmidt' for geodesy 4pi
                        normalized, orthonormalized, or Schmidt semi-normalized
                        coefficients, respectively.
        csphase       : 1 (default) if the coefficients exclude the Condon-
                        Shortley phase factor, or -1 if they include it.
        skip          : Number of lines to skip at the beginning of the file.

        Description
        -----------

        If format='shtools' spherical harmonic coefficients are read from an
        ascii-formatted file. The maximum spherical harmonic degree that is
        read is determined by the input value lmax. If the optional value skip
        is specified, parsing of the file will commence after the first skip
        lines.

        Each line of the file must contain

        l, m, cilm[0, l, m], cilm[1, l, m]

        For each value of increasing l, increasing from zero, all the angular
        orders are listed in inceasing order, from 0 to l.

        For more information, see SHRead.

        """
        if type(normalization) != str:
            raise ValueError('normalization must be a string. ' +
                             'Input type was {:s}'
                             .format(str(type(normalization))))

        if normalization.lower() not in set(['4pi', 'ortho', 'schmidt']):
            raise ValueError(
                "The input normalization must be '4pi', 'ortho' " +
                "or 'schmidt'. Provided value was {:s}"
                .format(repr(normalization))
                )
        if csphase != 1 and csphase != -1:
            raise ValueError(
                "csphase must be 1 or -1. Input value was {:s}"
                .format(repr(csphase))
                )

        if format.lower() == 'shtools':
            if kind.lower() == 'real':
                coeffs, lmax = shtools.SHRead(fname, lmax, **kwargs)
            else:
                raise NotImplementedError(
                    "kind={:s} not yet implemented".format(repr(kind)))
        else:
            raise NotImplementedError(
                "format={:s} not yet implemented".format(repr(format)))

        for cls in self.__subclasses__():
            if cls.istype(kind):
                return cls(coeffs, normalization=normalization.lower(),
                           csphase=csphase)

    # ---- Extract data ----
    def get_degrees(self):
        """
        Return a numpy array listing the spherical harmonic degrees
        from 0 to lmax.

        Usage
        -----

        degrees = x.get_degrees()

        Returns
        -------

        degrees : numpy ndarray of size (lmax+1).
        """
        return np.arange(self.lmax + 1)

    def get_powerperdegree(self):
        """
        Return a numpy array with the power per degree l spectrum.

        Usage
        -----

        power = x.get_powerperdegree()

        Returns
        -------

        power : numpy ndarray of size (lmax+1).
        """
        return self._powerperdegree()

    def get_powerperband(self, bandwidth):
        """
        Return a numpy array with the power per log_{bandwidth}(degree)
        spectrum.

        Usage
        -----

        power = x.get_powerperband()

        Returns
        -------

        power : numpy ndarray of size (lmax+1).
        """
        ls = self.get_degrees()
        return self._powerperdegree() * ls * np.log(bandwidth)

    # ---- Return coefficients with a different normalization convention ----
    def get_coeffs(self, normalization='4pi', csphase=1):
        """
        Return spherical harmonics coefficients as a numpy array with a
        different normalization convention.

        Usage
        -----

        coeffs = x.get_coeffs([normalization, csphase])

        Returns
        -------

        coeffs : numpy ndarray of size (lmax+1).

        Parameters
        ----------

        normalization : Normalization of the output coefficients:
                        '4pi' (default), 'ortho' or 'schmidt' for geodesy 4pi
                        normalized, orthonormalized, or Schmidt semi-normalized
                        coefficients, respectively.
        csphase       : Output Condon-Shortley phase convention: 1 (default)
                        to exlcude the phase factor, or -1 to include it.
        """
        if type(normalization) != str:
            raise ValueError('normalization must be a string. ' +
                             'Input type was {:s}'
                             .format(str(type(normalization))))

        if normalization.lower() not in set(['4pi', 'ortho', 'schmidt']):
            raise ValueError(
                "normalization must be '4pi', 'ortho' " +
                "or 'schmidt'. Provided value was {:s}"
                .format(repr(output_normalization))
                )
        if csphase != 1 and csphase != -1:
            raise ValueError(
                "csphase must be 1 or -1. Input value was {:s}"
                .format(repr(csphase))
                )

        return self._get_coeffs(
            output_normalization=normalization.lower(),
            output_csphase=csphase)

    # ---- Rotate the coordinate system ----
    def rotate(self, alpha, beta, gamma, degrees=True, dj_matrix=None):
        """
        Rotate the body or coordinate system used to express the spherical
        harmonic coefficients and output as a new class instance.

        Usage
        -----

        SHCoeffsInstance = x.rotate(alpha, beta, gamma, [degrees, dj_matrix])

        Parameters
        ----------

        alpha, beta, gamma : The three Euler rotation angles in degrees.
        degrees            : True (default) if the Euler angles are in degrees,
                             False if they are in radians.
        dj_matrix          : The djpi2 rotation matrix (default=None), computed
                             by a call to djpi2.

        Description
        -----------
        This method will take the spherical harmonic coefficients of a
        function, rotate the coordinate frame by the three Euler anlges, and
        output the spherical harmonic coefficients of the rotated function.

        The rotation of a coordinate system or body can be viewed in two
        complementary ways involving three successive rotations. Both methods
        have the same initial and final configurations, and the angles listed
        in both schemes are the same.

        Scheme A:

        (I) Rotation about the z axis by alpha.
        (II) Rotation about the new y axis by beta.
        (III) Rotation about the new z axis by gamma.

        Scheme B:

        (I) Rotation about the z axis by gamma.
        (II) Rotation about the initial y axis by beta.
        (III) Rotation about the initial z axis by alpha.

        The rotations can further be viewed either as a rotation of the
        coordinate system or the physical body. For a rotation of the
        coordinate system without rotation of the physical body, use

        (alpha, beta, gamma).

        For a rotation of the physical body without rotation of the coordinate
        system, use

        (-gamma, -beta, -alpha).

        To perform the inverse transform of (alpha, beta, gamma), use

        (-gamma, -beta, -alpha).

        Note that this routine uses the "y convention", where the second
        rotation is with respect to the new y axis. If alpha, beta, and gamma
        were orginally defined in terms of the "x convention", where the second
        rotation was with respect to the new x axis, the Euler angles according
        to the y convention would be

        alpha_y=alpha_x-pi/2, beta_x=beta_y, and gamma_y=gamma_x+pi/2.
        """
        if degrees:
            angles = np.radians([alpha, beta, gamma])
        else:
            angles = np.array([alpha, beta, gamma])

        rot = self._rotate(angles, dj_matrix)
        return rot

    # ---- Convert spherical harmonic coefficients to a different normalization
    def convert(self, normalization='4pi', csphase=1):
        """
        Convert the spherical harmonic coefficients to a different
        normalization convention, and return a new class instance.

        Usage
        -----

        SHCoeffsInstance = x.convert([normalization, csphase])

        Parameters
        ----------

        normalization : Normalization of the output class: '4pi' (default),
                        'ortho' or 'schmidt' for geodesy 4pi normalized,
                        orthonormalized, or Schmidt semi-normalized
                        coefficients, respectively.
        csphase       : Output Condon-Shortley phase convention: 1 (default)
                        to exlcude the phase factor, or -1 to include it.
        """
        if type(normalization) != str:
            raise ValueError('normalization must be a string. ' +
                             'Input type was {:s}'
                             .format(str(type(normalization))))

        if normalization.lower() not in set(['4pi', 'ortho', 'schmidt']):
            raise ValueError(
                "normalization must be '4pi', 'ortho' " +
                "or 'schmidt'. Provided value was {:s}"
                .format(repr(normalization))
                )
        if csphase != 1 and csphase != -1:
            raise ValueError(
                "csphase must be 1 or -1. Input value was {:s}"
                .format(repr(csphase))
                )

        coeffs = self.get_coeffs(normalization=normalization.lower(),
                                 csphase=csphase)
        return SHCoeffs.from_array(coeffs,
                                   normalization=normalization.lower(),
                                   csphase=csphase)

    # ---- Expand the coefficients onto a grid ----
    def expand(self, grid='DH', **kwargs):
        """
        Evaluate the coefficients on a spherical grid.

        Usage
        -----

        SHGridInstance = x.expand([grid, lmax, lmax_calc, zeros])

        Parameters
        ----------

        grid      : 'DH' or 'DH1' for an equisampled lat/lon grid with
                    nlat=nlon, 'DH2' for an equidistant lat/lon grid with
                    nlon=2*nlat, or 'GLQ' for a Gauss-Legendre quadrature grid.
        lmax      : The maximum spherical harmonic degree, which determines the
                    grid spacing of the output grid. Default is x.lmax.
        lmax_calc : The maximum spherical harmonic degree to use when
                    evaluating the function. Default is x.lmax.
        zeros     : The cos(colatitude) nodes used in the Gauss-Legendre
                    Quadrature grids. Default is None.

        Description
        -----------

        For more information concerning the spherical harmonic expansions, and
        the properties of the output grids, see the documentation for
        SHExpandDH, SHExpandDHC, SHExpandGLQ and SHExpandGLQC.
        """
        if type(grid) != str:
            raise ValueError('grid must be a string. ' +
                             'Input type was {:s}'
                             .format(str(type(grid))))

        if grid.upper() == 'DH' or grid.upper() == 'DH1':
            gridout = self._expandDH(sampling=1, **kwargs)
        elif grid.upper() == 'DH2':
            gridout = self._expandDH(sampling=2, **kwargs)
        elif grid.upper() == 'GLQ':
            gridout = self._expandGLQ(zeros=None, **kwargs)
        else:
            raise ValueError(
                "grid must be 'DH', 'DH1', 'DH2', or 'GLQ'. " +
                "Input value was {:s}".format(repr(grid)))

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

        loglog : If True (default), use log-log axis.
        show   : If True (default), plot to the screen.
        fname  : If present, save the image to the file.
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
        Plot the power per log_{bandwidth}(degree) spectrum.

        Usage
        -----

        x.plot_powerperband([loglog, show, fname])

        Parameters
        ----------

        bandwidth : The bandwidth, default = 2.
        loglog    : If True (default), use log-log axis.
        show      : If True (default), plot to the screen.
        fname     : If present, save the image to the file
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

class SHRealCoeffs(SHCoeffs):
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
        self.normalization = normalization
        self.csphase = csphase

    def make_complex(self):
        """
        Convert the real spherical harmonic coefficient class to the complex
        spherical harmonic coefficient class with the same normalization and
        csphase conventions.

        Usage
        -----

        SHComplexCoeffsInstance = x.make_complex()
        """
        rcomplex_coeffs = shtools.SHrtoc(self.coeffs,
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
            return shtools.SHPowerSpectrum(self.coeffs)
        elif self.normalization == 'schmidt':
            power = shtools.SHPowerSpectrum(self.coeffs)
            l = self.get_degrees()
            power /= (2.0 * l + 1.0)
            return power
        elif self.normalization == 'ortho':
            return shtools.SHPowerSpectrum(self.coeffs) / (4.0 * np.pi)
        else:
            raise ValueError(
                "Normalization must be '4pi', 'ortho', or 'schmidt'. " +
                "Input value was {:s}".format(repr(self.normalization)))

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
            dj_matrix = shtools.djpi2(self.lmax + 1)

        # The coefficients need to be 4pi normalized with csphase = 1
        coeffs = shtools.SHRotateRealCoef(
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
                "Normalization must be '4pi', 'ortho', or 'schmidt'. " +
                "Input value was {:s}".format(repr(self.normalization)))

        data = shtools.MakeGridDH(self.coeffs, sampling=sampling, norm=norm,
                                  csphase=self.csphase, **kwargs)
        gridout = SHGrid.from_array(data, grid='DH')
        return gridout

    def _expandGLQ(self, zeros, **kwargs):
        """Evaluate the coefficients on a Gauss Legendre quadrature grid."""
        if self.normalization == '4pi':
            norm = 1
        elif self.normalization == 'schmidt':
            norm = 2
        elif self.normalization == 'ortho':
            norm = 4
        else:
            raise ValueError(
                "Normalization must be '4pi', 'ortho', or 'schmidt'. " +
                "Input value was {:s}".format(repr(self.normalization)))

        if zeros is None:
            zeros, weights = shtools.SHGLQ(self.lmax)

        data = shtools.MakeGridGLQ(self.coeffs, zeros, norm=norm,
                                   csphase=self.csphase, **kwargs)
        gridout = SHGrid.from_array(data, grid='GLQ')
        return gridout


# =============== COMPLEX SPHERICAL HARMONICS ================

class SHComplexCoeffs(SHCoeffs):
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
        self.normalization = normalization
        self.csphase = csphase

    def make_real(self):
        """
        Convert the complex coefficient class to the real harmonic
        coefficient class.

        Usage
        -----

        SHRealCoeffsInstance = x.make_real()
        """
        # First test if the coefficients correspond to a real grid.
        # This is not very elegant. Also, the equality condition
        # is probably not robust to round off errors.
        for l in self.get_degrees():
            if self.coeffs[0, l, 0] != self.coeffs[0, l, 0].conjugate():
                raise RuntimeError('Complex coefficients do not correspond ' +
                                   'to a real field. l = {:d}, m = 0: {:e}'
                                   .format(l, self.coeffs[0, l, 0]))

            for m in np.arange(1, l + 1):
                if m % 2 == 1:
                    if (self.coeffs[0, l, m] != -
                            self.coeffs[1, l, m].conjugate()):
                        raise RuntimeError('Complex coefficients do not ' +
                                           'correspond to a real field. ' +
                                           'l = {:d}, m = {:d}: {:e}, {:e}'
                                           .format(l, m, self.coeffs[0, l, 0],
                                                   self.coeffs[1, l, 0]))
                else:
                    if (self.coeffs[0, l, m] !=
                            self.coeffs[1, l, m].conjugate()):
                        raise RuntimeError('Complex coefficients do not ' +
                                           'correspond to a real field. ' +
                                           'l = {:d}, m = {:d}: {:e}, {:e}'
                                           .format(l, m, self.coeffs[0, l, 0],
                                                   self.coeffs[1, l, 0]))

        coeffs_rc = np.zeros((2, self.lmax + 1, self.lmax + 1))
        coeffs_rc[0, :, :] = self.coeffs[0, :, :].real
        coeffs_rc[1, :, :] = self.coeffs[0, :, :].imag
        real_coeffs = shtools.SHctor(coeffs_rc, convention=1,
                                     switchcs=0)
        return SHCoeffs.from_array(real_coeffs,
                                   normalization=self.normalization,
                                   csphase=self.csphase)

    def _powerperdegree(self):
        """Return the power per degree l spectrum."""
        if self.normalization == '4pi':
            return shtools.SHPowerSpectrumC(self.coeffs)
        elif self.normalization == 'schmidt':
            power = shtools.SHPowerSpectrumC(self.coeffs)
            l = self.get_degrees()
            power /= (2.0 * l + 1.0)
            return power
        elif self.normalization == 'ortho':
            return shtools.SHPowerSpectrumC(self.coeffs) / (4.0 * np.pi)
        else:
            raise ValueError(
                "Normalization must be '4pi', 'ortho', or 'schmidt'. " +
                "Input value was {:s}".format(repr(self.normalization)))

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
            dj_matrix = shtools.djpi2(self.lmax + 1)

        cgrid = self.expand(grid='DH')
        rgrid, igrid = cgrid.data.real, cgrid.data.imag
        rgridcoeffs = shtools.SHExpandDH(rgrid, norm=1, sampling=1, csphase=1)
        igridcoeffs = shtools.SHExpandDH(igrid, norm=1, sampling=1, csphase=1)

        rgridcoeffs_rot = shtools.SHRotateRealCoef(
            rgridcoeffs, angles, dj_matrix)
        igridcoeffs_rot = shtools.SHRotateRealCoef(
            igridcoeffs, angles, dj_matrix)

        rgrid_rot = shtools.MakeGridDH(rgridcoeffs_rot, norm=1,
                                       sampling=1, csphase=1)
        igrid_rot = shtools.MakeGridDH(igridcoeffs_rot, norm=1,
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

        coeffs_rot = shtools.SHExpandDHC(grid_rot, norm=norm,
                                         csphase=self.csphase)

        return SHCoeffs.from_array(coeffs_rot,
                                   normalization=self.normalization,
                                   csphase=self.csphase)

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
                "Normalization must be '4pi', 'ortho', or 'schmidt'. " +
                "Input value was {:s}".format(repr(self.normalization)))

        data = shtools.MakeGridDHC(self.coeffs, sampling=sampling,
                                   norm=norm, csphase=self.csphase, **kwargs)
        gridout = SHGrid.from_array(data, grid='DH')
        return gridout

    def _expandGLQ(self, zeros, **kwargs):
        """Evaluate the coefficients on a Gauss-Legendre quadrature grid."""
        if self.normalization == '4pi':
            norm = 1
        elif self.normalization == 'schmidt':
            norm = 2
        elif self.normalization == 'ortho':
            norm = 4
        else:
            raise ValueError(
                "Normalization must be '4pi', 'ortho', or 'schmidt'. " +
                "Input value was {:s}".format(repr(self.normalization)))

        if zeros is None:
            zeros, weights = shtools.SHGLQ(self.lmax)

        data = shtools.MakeGridGLQC(self.coeffs, zeros, norm=norm,
                                    csphase=self.csphase, **kwargs)
        gridout = SHGrid.from_array(data, grid='GLQ')
        return gridout
