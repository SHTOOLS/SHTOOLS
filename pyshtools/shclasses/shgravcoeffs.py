"""
    Class for spherical harmonic coefficients of the gravitational potential.
"""
from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

import numpy as _np
import matplotlib as _mpl
import matplotlib.pyplot as _plt
import copy as _copy
import warnings as _warnings
from scipy.special import factorial as _factorial

from .shcoeffsgrid import SHCoeffs as _SHCoeffs
from .shcoeffsgrid import SHRealCoeffs as _SHRealCoeffs
from .shcoeffsgrid import DHRealGrid as _DHRealGrid
from .shgravgrid import SHGravGrid as _SHGravGrid
from .shtensor import SHGravTensor as _SHGravTensor
from .shgeoid import SHGeoid as _SHGeoid

from ..constant import G as _G
from ..spectralanalysis import spectrum as _spectrum
from ..shio import convert as _convert
from ..shio import shread as _shread
from ..shtools import CilmPlusRhoHDH as _CilmPlusRhoHDH
from ..shtools import CilmPlusDH as _CilmPlusDH
from ..shtools import MakeGravGridDH as _MakeGravGridDH
from ..shtools import MakeGravGradGridDH as _MakeGravGradGridDH
from ..shtools import MakeGeoidGridDH as _MakeGeoidGridDH


# =============================================================================
# =========    SHGravCoeffs class    =========================================
# =============================================================================

class SHGravCoeffs(object):
    """
    Spherical harmonic coefficients class for the gravitational potential.

    The coefficients of this class can be initialized using one of the four
    constructor methods:

        x = SHGravCoeffs.from_array(array, gm, r0)
        x = SHGravCoeffs.from_random(powerspectrum, gm, r0)
        x = SHGravCoeffs.from_zeros(lmax, gm, r0)
        x = SHGravCoeffs.from_file('fname.dat')
        x = SHGravCoeffs.from_shape(grid, rho, gm)

    The normalization convention of the input coefficents is specified
    by the optional normalization and csphase parameters, which take the
    following values:

    normalization : '4pi' (default), geodesy 4-pi normalized.
                  : 'ortho', orthonormalized.
                  : 'schmidt', Schmidt semi-normalized.
                  : 'unnorm', unnormalized.

    csphase       : 1 (default), exlcude the Condon-Shortley phase factor.
                  : -1, include the Condon-Shortley phase factor.

    See the documentation for each constructor method for further options.

    Once initialized, each class instance defines the following class
    attributes:

    lmax          : The maximum spherical harmonic degree of the coefficients.
    coeffs        : The raw coefficients with the specified normalization and
                    csphase conventions.
    errors        : The uncertainties of the spherical harmonic coefficients.
    gm            : The gravitational constant times the mass times that is
                    associated with the gravitational potential coefficients.
    r0            : The reference radius of the gravitational potential
                    coefficients.
    omega         : The angular rotation rate of the body.
    normalization : The normalization of the coefficients: '4pi', 'ortho',
                    'schmidt', or 'unnorm'.
    csphase       : Defines whether the Condon-Shortley phase is used (1)
                    or not (-1).
    mask          : A boolean mask that is True for the permissible values of
                    degree l and order m.
    kind          : The coefficient data type (only 'real' is permissible).
    header        : A list of values (of type str) from the header line of the
                    input file used to initialize the class (for 'shtools'
                    formatted files only).

    Each class instance provides the following methods:

    degrees()             : Return an array listing the spherical harmonic
                            degrees from 0 to lmax.
    spectrum()            : Return the spectrum of the function as a function
                            of spherical harmonic degree.
    set_omega()           : Set the angular rotation rate of the body.
    set_coeffs()          : Set coefficients in-place to specified values.
    change_ref()          : Return a new class instance referenced to a
                            different gm, or r0.
    rotate()              : Rotate the coordinate system used to express the
                            spherical harmonic coefficients and return a new
                            class instance.
    convert()             : Return a new class instance using a different
                            normalization convention.
    pad()                 : Return a new class instance that is zero padded or
                            truncated to a different lmax.
    expand()              : Calculate the three vector components of the
                            gravity field, the total field, and the
                            gravitational potential, and return an SHGravGrid
                            class instance.
    tensor()              : Calculate the 9 components of the gravity tensor
                            and return an SHGravTensor class instance.
    geoid()               : Calculate the height of the geoid and return an
                            SHGeoid class instance.
    plot_spectrum()       : Plot the spectrum as a function of spherical
                            harmonic degree.
    plot_spectrum2d()     : Plot the 2D spectrum of all spherical harmonic
                            degrees and orders.
    to_array()            : Return an array of spherical harmonic coefficients
                            with a different normalization convention.
    to_file()             : Save the spherical harmonic coefficients as a file.
    copy()                : Return a copy of the class instance.
    info()                : Print a summary of the data stored in the
                            SHGravCoeffs instance.
    """

    def __init__(self):
        """Unused constructor of the super class."""
        print('Initialize the class using one of the class methods:\n'
              '>>> pyshtools.SHGravCoeffs.from_array\n'
              '>>> pyshtools.SHGravCoeffs.from_random\n'
              '>>> pyshtools.SHGravCoeffs.from_zeros\n'
              '>>> pyshtools.SHGravCoeffs.from_file\n'
              '>>> pyshtools.SHGravCoeffs.from_shape\n')

    # ---- Factory methods ----
    @classmethod
    def from_array(self, coeffs, gm, r0, omega=None, errors=None,
                   normalization='4pi', csphase=1, lmax=None, copy=True):
        """
        Initialize the class with spherical harmonic coefficients from an input
        array.

        Usage
        -----
        x = SHGravCoeffs.from_array(array, gm, r0, [omega, errors,
                                                    normalization, csphase,
                                                    lmax, copy])

        Returns
        -------
        x : SHGravCoeffs class instance.

        Parameters
        ----------
        array : ndarray, shape (2, lmaxin+1, lmaxin+1).
            The input spherical harmonic coefficients.
        gm : float
            The gravitational constant times the mass that is associated with
            the gravitational potential coefficients.
        mass : float
            The mass of the planet in kg.
        r0 : float
            The reference radius of the spherical harmonic coefficients.
        omega : float, optional, default = None
            The angular rotation rate of the body.
        errors : ndarray, optional, default = None
            The uncertainties of the spherical harmonic coefficients.
        normalization : str, optional, default = '4pi'
            '4pi', 'ortho', 'schmidt', or 'unnorm' for geodesy 4pi normalized,
            orthonormalized, Schmidt semi-normalized, or unnormalized
            coefficients, respectively.
        csphase : int, optional, default = 1
            Condon-Shortley phase convention: 1 to exclude the phase factor,
            or -1 to include it.
        lmax : int, optional, default = None
            The maximum spherical harmonic degree to include in the returned
            class instance. This must be less than or equal to lmaxin.
        copy : bool, optional, default = True
            If True, make a copy of array when initializing the class instance.
            If False, initialize the class instance with a reference to array.

        Notes
        -----
        If the degree-0 term of the input array is equal to zero, it will be
        set to 1.
        """
        if _np.iscomplexobj(coeffs):
            raise TypeError('The input array must be real.')

        if type(normalization) != str:
            raise ValueError('normalization must be a string. '
                             'Input type was {:s}'
                             .format(str(type(normalization))))

        if normalization.lower() not in ('4pi', 'ortho', 'schmidt', 'unnorm'):
            raise ValueError(
                "The normalization must be '4pi', 'ortho', 'schmidt', "
                "or 'unnorm'. Input value was {:s}."
                .format(repr(normalization))
                )

        if csphase != 1 and csphase != -1:
            raise ValueError(
                "csphase must be either 1 or -1. Input value was {:s}."
                .format(repr(csphase))
                )

        if errors is not None:
            if coeffs.shape != errors.shape:
                raise ValueError(
                    "The shape of coeffs and errors must be the same."
                    "Shape of coeffs = {:s}, shape of errors = {:s}"
                    .format(repr(coeffs.shape), repr(coeffs.errors))
                    )

        lmaxin = coeffs.shape[1] - 1
        if lmax is None:
            lmax = lmaxin
        else:
            if lmax > lmaxin:
                lmax = lmaxin

        if normalization.lower() == 'unnorm' and lmax > 85:
            _warnings.warn("Calculations using unnormalized coefficients "
                           "are stable only for degrees less than or equal "
                           "to 85. lmax for the coefficients will be set to "
                           "85. Input value was {:d}.".format(lmax),
                           category=RuntimeWarning)
            lmax = 85

        if coeffs[0, 0, 0] == 0:
            _warnings.warn('The degree 0 term of the array was not set. This, '
                           'will be set to 1.', category=RuntimeWarning)
            coeffs[0, 0, 0] = 1.0

        if errors is not None:
            clm = SHGravRealCoeffs(coeffs[:, 0:lmax+1, 0:lmax+1], gm=gm, r0=r0,
                                   omega=omega, errors=errors[:, 0:lmax+1,
                                                              0:lmax+1],
                                   normalization=normalization.lower(),
                                   csphase=csphase, copy=copy)
        else:
            clm = SHGravRealCoeffs(coeffs[:, 0:lmax+1, 0:lmax+1], gm=gm, r0=r0,
                                   omega=omega,
                                   normalization=normalization.lower(),
                                   csphase=csphase, copy=copy)
        return clm

    @classmethod
    def from_zeros(self, lmax, gm, r0, omega=None, errors=False,
                   normalization='4pi', csphase=1):
        """
        Initialize the class with spherical harmonic coefficients set to zero
        from degree 1 to lmax, and set the degree 0 term to 1.

        Usage
        -----
        x = SHGravCoeffs.from_zeros(lmax, gm, r0, [omega, errors,
                                                   normalization, csphase])

        Returns
        -------
        x : SHGravCoeffs class instance.

        Parameters
        ----------
        lmax : int
            The maximum spherical harmonic degree l of the coefficients.
        gm : float
            The gravitational constant times the mass that is associated with
            the gravitational potential coefficients.
        r0 : float
            The reference radius of the spherical harmonic coefficients.
        omega : float, optional, default = None
            The angular rotation rate of the body.
        errors : bool, optional, default = False
            If True, initialize the attribute errors with zeros.
        normalization : str, optional, default = '4pi'
            '4pi', 'ortho', 'schmidt', or 'unnorm' for geodesy 4pi normalized,
            orthonormalized, Schmidt semi-normalized, or unnormalized
             coefficients, respectively.
        csphase : int, optional, default = 1
            Condon-Shortley phase convention: 1 to exclude the phase factor,
            or -1 to include it.
        """
        if normalization.lower() not in ('4pi', 'ortho', 'schmidt', 'unnorm'):
            raise ValueError(
                "The normalization must be '4pi', 'ortho', 'schmidt', "
                "or 'unnorm'. Input value was {:s}."
                .format(repr(normalization))
                )

        if csphase != 1 and csphase != -1:
            raise ValueError(
                "csphase must be either 1 or -1. Input value was {:s}."
                .format(repr(csphase))
                )

        if normalization.lower() == 'unnorm' and lmax > 85:
            _warnings.warn("Calculations using unnormalized coefficients "
                           "are stable only for degrees less than or equal "
                           "to 85. lmax for the coefficients will be set to "
                           "85. Input value was {:d}.".format(lmax),
                           category=RuntimeWarning)
            lmax = 85

        coeffs = _np.zeros((2, lmax + 1, lmax + 1))
        coeffs[0, 0, 0] = 1.0

        if errors is False:
            clm = SHGravRealCoeffs(coeffs, gm=gm, r0=r0, omega=omega,
                                   normalization=normalization.lower(),
                                   csphase=csphase)
        else:
            clm = SHGravRealCoeffs(coeffs, gm=gm, r0=r0, omega=omega,
                                   errors=_np.zeros((2, lmax + 1, lmax + 1)),
                                   normalization=normalization.lower(),
                                   csphase=csphase)
        return clm

    @classmethod
    def from_file(self, fname, format='shtools', gm=None, r0=None,
                  omega=None, lmax=None, normalization='4pi', skip=0,
                  header=True, errors=False, csphase=1, r0_index=0, gm_index=1,
                  omega_index=None, header_units='m', **kwargs):
        """
        Initialize the class with spherical harmonic coefficients from a file.

        Usage
        -----
        x = SHGravCoeffs.from_file(filename, [format='shtools', gm, r0, omega,
                                              lmax, normalization, csphase,
                                              skip, header, errors, gm_index,
                                              r0_index, omega_index,
                                              header_units])
        x = SHGravCoeffs.from_file(filename, format='npy', gm, r0,
                                   [omega, normalization, csphase, **kwargs])

        Returns
        -------
        x : SHGravCoeffs class instance.

        Parameters
        ----------
        filename : str
            Name of the file, including path.
        format : str, optional, default = 'shtools'
            'shtools' format or binary numpy 'npy' format.
        lmax : int, optional, default = None
            The maximum spherical harmonic degree to read from 'shtools'
            formatted files.
        header : bool, optional, default = True
            If True, read a list of values from the header line of an 'shtools'
            formatted file.
        errors : bool, optional, default = False
            If True, read errors from the file (for 'shtools' formatted files
            only).
        r0_index : int, optional, default = 0
            For shtools formatted files, if header is True, r0 will be set
            using the value from the header line with this index.
        gm_index : int, optional, default = 1
            For shtools formatted files, if header is True, gm will be set
            using the value from the header line with this index.
        omega_index : int, optional, default = None
            For shtools formatted files, if header is True, omega will be set
            using the value from the header line with this index.
        gm : float, optional, default = None
            The gravitational constant time the mass that is associated with
            the gravitational potential coefficients.
        r0 : float, optional, default = None
            The reference radius of the spherical harmonic coefficients.
        omega : float, optional, default = None
            The angular rotation rate of the body.
        header_units : str, optional, default = 'm'
            The units used for r0 and gm in the header line of an shtools
            formatted file: 'm' or 'km'. If 'km', the values of r0 and gm will
            be converted to meters.
        normalization : str, optional, default = '4pi'
            '4pi', 'ortho', 'schmidt', or 'unnorm' for geodesy 4pi normalized,
            orthonormalized, Schmidt semi-normalized, or unnormalized
            coefficients, respectively.
        csphase : int, optional, default = 1
            Condon-Shortley phase convention: 1 to exclude the phase factor,
            or -1 to include it.
        skip : int, optional, default = 0
            Number of lines to skip at the beginning of the file when format is
            'shtools'.
        **kwargs : keyword argument list, optional for format = 'npy'
            Keyword arguments of numpy.load() when format is 'npy'.

        Description
        -----------
        If format='shtools', spherical harmonic coefficients will be read from
        a text file. The optional parameter `skip` specifies how many lines
        should be skipped before attempting to parse the file, the optional
        parameter `header` specifies whether to read a list of values from a
        header line, and the optional parameter `lmax` specifies the maximum
        degree to read from the file. If a header line is read, r0_index,
        gm_index, and omega_index, are used as the indices to set r0, gm, and
        omega. If header_unit is specified as 'km', the values of r0 and gm
        that are read from the header will be converted to meters.

        For shtools formatted files, all lines that do not start with 2
        integers and that are less than 3 words long will be treated as
        comments and ignored. For this format, each line of the file must
        contain

        l, m, coeffs[0, l, m], coeffs[1, l, m]

        where l and m are the spherical harmonic degree and order,
        respectively. The terms coeffs[1, l, 0] can be neglected as they are
        zero. For more information, see `shio.shread()`. If errors are read,
        each line must contain:

        l, m, coeffs[0, l, m], coeffs[1, l, m], error[0, l, m], error[1, l, m]

        If format='npy', a binary numpy 'npy' file will be read using
        numpy.load().

        If the degree 0 term of the file is zero (or not specified), this will
        be set to 1.
        """
        error = None

        if type(normalization) != str:
            raise ValueError('normalization must be a string. '
                             'Input type was {:s}'
                             .format(str(type(normalization))))

        if normalization.lower() not in ('4pi', 'ortho', 'schmidt', 'unnorm'):
            raise ValueError(
                "The input normalization must be '4pi', 'ortho', 'schmidt', "
                "or 'unnorm'. Provided value was {:s}"
                .format(repr(normalization))
                )

        if csphase != 1 and csphase != -1:
            raise ValueError(
                "csphase must be 1 or -1. Input value was {:s}"
                .format(repr(csphase))
                )

        if header_units.lower() not in ('m', 'km'):
            raise ValueError("header_units can be only 'm', or 'km'. Input "
                             "value is {:s}.".format(repr(header_units)))

        if format == 'shtools':
            if r0_index is not None and r0 is not None:
                raise ValueError('Can not specify both r0_index and r0')
            if gm_index is not None and gm is not None:
                raise ValueError('Can not specify both gm_index and gm')
            if omega_index is not None and omega is not None:
                raise ValueError('Can not specify both omega_index and omega')
            if header is False and (r0 is None or gm is None):
                raise ValueError('If header is False, r0 and gm must be '
                                 'specified.')

        header_list = None
        if format.lower() == 'shtools':
            if header is True:
                if errors is True:
                    coeffs, error, lmaxout, header_list = _shread(
                        fname, lmax=lmax, skip=skip, header=True, error=True)
                else:
                    coeffs, lmaxout, header_list = _shread(
                        fname, lmax=lmax, skip=skip, header=True)
            else:
                if errors is True:
                    coeffs, error, lmaxout = _shread(
                        fname, lmax=lmax, error=True, skip=skip)
                else:
                    coeffs, lmaxout = _shread(fname, lmax=lmax, skip=skip)

        elif format.lower() == 'npy':
            if gm is None or r0 is None:
                raise ValueError('For binary npy files, gm and r0 must be '
                                 'specified.')
            coeffs = _np.load(fname, **kwargs)
            lmaxout = coeffs.shape[1] - 1
        else:
            raise NotImplementedError(
                'format={:s} not implemented'.format(repr(format)))

        if _np.iscomplexobj(coeffs):
            raise TypeError('The input coefficients must be real.')

        if normalization.lower() == 'unnorm' and lmaxout > 85:
            _warnings.warn("Calculations using unnormalized coefficients "
                           "are stable only for degrees less than or equal "
                           "to 85. lmax for the coefficients will be set to "
                           "85. Input value was {:d}.".format(lmaxout),
                           category=RuntimeWarning)
            lmaxout = 85

        if coeffs[0, 0, 0] == 0:
            warnstr = ("The degree 0 term of the file was not set. "
                       "This will be set to 1.")
            _warnings.warn(warnstr, category=RuntimeWarning)
            coeffs[0, 0, 0] = 1.0

        if format.lower() == 'shtools' and header is True:
            if r0_index is not None:
                r0 = float(header_list[r0_index])
            if gm_index is not None:
                gm = float(header_list[gm_index])
            if omega_index is not None:
                omega = float(header_list[omega_index])
            if header_units.lower() == 'km':
                r0 *= 1.e3
                gm *= 1.e9

        clm = SHGravRealCoeffs(coeffs, gm=gm, r0=r0, omega=omega,
                               errors=error,
                               normalization=normalization.lower(),
                               csphase=csphase, header=header_list)
        return clm

    @classmethod
    def from_random(self, power, gm, r0, omega=None, function='geoid',
                    lmax=None, normalization='4pi', csphase=1,
                    exact_power=False):
        """
        Initialize the class of gravitational potential spherical harmonic
        coefficients as random variables with a given spectrum.

        Usage
        -----
        x = SHGravCoeffs.from_random(power, gm, r0, [omega, function, lmax,
                                                     normalization,
                                                     csphase, exact_power])

        Returns
        -------
        x : SHGravCoeffs class instance.

        Parameters
        ----------
        power : ndarray, shape (L+1)
            numpy array of shape (L+1) that specifies the expected power per
            degree l, where L is the maximum spherical harmonic bandwidth.
        gm : float
            The gravitational constant times the mass that is associated with
            the gravitational potential coefficients.
        r0 : float
            The reference radius of the spherical harmonic coefficients.
        omega : float, optional, default = None
            The angular rotation rate of the body.
        function : str, optional, default = 'geoid'
            The type of input power spectrum: 'potential' for the gravitational
            potential, 'geoid' for the geoid, 'radial' for the radial gravity,
            or 'total' for the total gravity field.
        lmax : int, optional, default = len(power) - 1
            The maximum spherical harmonic degree l of the output coefficients.
            The coefficients will be set to zero for degrees greater than L.
        normalization : str, optional, default = '4pi'
            '4pi', 'ortho', 'schmidt', or 'unnorm' for geodesy 4pi normalized,
            orthonormalized, Schmidt semi-normalized, or unnormalized
            coefficients, respectively.
        csphase : int, optional, default = 1
            Condon-Shortley phase convention: 1 to exclude the phase factor,
            or -1 to include it.
        exact_power : bool, optional, default = False
            The total variance of the coefficients is set exactly to the input
            power. The distribution of power at degree l amongst the angular
            orders is random, but the total power is fixed.

        Description
        -----------
        This routine returns a random realization of spherical harmonic
        gravitational potential coefficients obtained from a normal
        distribution. The variance of each coefficient at degree l is equal to
        the total power at degree l divided by the number of coefficients at
        that degree (2l+1). These coefficients are then divided by a prefactor
        that depends upon the function being used to calculate the spectrum:
        gm/r0 for the gravitiational potential, r0 for the geoid, and
        (l+1)*gm/(r**2) for the radial gravity. The power spectrum of the
        random realization can be fixed exactly to the input spectrum by
        setting exact_power to True.

        Note that the degree 0 term is set to 1, and the degree-1 terms are
        set to 0.
        """
        if type(normalization) != str:
            raise ValueError('normalization must be a string. '
                             'Input type was {:s}'
                             .format(str(type(normalization))))

        if function.lower() not in ('potential', 'geoid', 'radial', 'total'):
            raise ValueError(
                "function must be of type 'potential', "
                "'geoid', 'radial', or 'total'. Provided value was {:s}"
                .format(repr(function))
                )

        if normalization.lower() not in ('4pi', 'ortho', 'schmidt', 'unnorm'):
            raise ValueError(
                "The input normalization must be '4pi', 'ortho', 'schmidt', "
                "or 'unnorm'. Provided value was {:s}"
                .format(repr(normalization))
                )

        if csphase != 1 and csphase != -1:
            raise ValueError(
                "csphase must be 1 or -1. Input value was {:s}"
                .format(repr(csphase))
                )

        if lmax is None:
            nl = len(power)
            lmax = nl - 1
        else:
            if lmax <= len(power) - 1:
                nl = lmax + 1
            else:
                nl = len(power)
        degrees = _np.arange(nl)

        if normalization.lower() == 'unnorm' and nl - 1 > 85:
            _warnings.warn("Calculations using unnormalized coefficients "
                           "are stable only for degrees less than or equal "
                           "to 85. lmax for the coefficients will be set to "
                           "85. Input value was {:d}.".format(nl-1),
                           category=RuntimeWarning)
            nl = 85 + 1
            lmax = 85

        # Create coefficients with unit variance, which returns an expected
        # total power per degree of (2l+1) for 4pi normalized harmonics.
        coeffs = _np.empty((2, nl, nl))
        for l in degrees:
            coeffs[:2, l, :l+1] = _np.random.normal(size=(2, l+1))

        if exact_power:
            power_per_l = _spectrum(coeffs, normalization='4pi', unit='per_l')
            coeffs *= _np.sqrt(
                power[0:nl] / power_per_l)[_np.newaxis, :, _np.newaxis]
        else:
            coeffs *= _np.sqrt(
                power[0:nl] / (2 * degrees + 1))[_np.newaxis, :, _np.newaxis]

        if normalization.lower() == '4pi':
            pass
        elif normalization.lower() == 'ortho':
            coeffs = _convert(coeffs, normalization_in='4pi',
                              normalization_out='ortho')
        elif normalization.lower() == 'schmidt':
            coeffs = _convert(coeffs, normalization_in='4pi',
                              normalization_out='schmidt')
        elif normalization.lower() == 'unnorm':
            coeffs = _convert(coeffs, normalization_in='4pi',
                              normalization_out='unnorm')

        if function.lower() == 'potential':
            coeffs /= (gm / r0)
        elif function.lower() == 'geoid':
            coeffs /= r0
        elif function.lower() == 'radial':
            for l in degrees:
                coeffs[:, l, :l+1] /= (gm * (l + 1) / r0**2)
        elif function.lower() == 'total':
            for l in degrees:
                coeffs[:, l, :l+1] /= (gm / r0**2) * _np.sqrt((l + 1) *
                                                              (2 * l + 1))

        if lmax > nl - 1:
            coeffs = _np.pad(coeffs, ((0, 0), (0, lmax - nl + 1),
                             (0, lmax - nl + 1)), 'constant')

        coeffs[0, 0, 0] = 1.0
        coeffs[:, 1, :] = 0.0

        clm = SHGravRealCoeffs(coeffs, gm=gm, r0=r0, omega=omega,
                               normalization=normalization.lower(),
                               csphase=csphase)
        return clm

    @classmethod
    def from_shape(self, shape, rho, gm, nmax=7, lmax=None, lmax_grid=None,
                   lmax_calc=None, omega=None):
        """
        Initialize a class of gravitational potential spherical harmonic
        coefficients by calculuting the gravitational potential associatiated
        with relief along an interface.

        Usage
        -----
        x = SHGravCoeffs.from_shape(shape, rho, gm, [nmax, lmax, lmax_grid,
                                                     lmax_calc, omega])

        Returns
        -------
        x : SHGravCoeffs class instance.

        Parameters
        ----------
        shape : SHGrid or SHCoeffs class instance
            The shape of the interface, either as an SHGrid or SHCoeffs class
            instance. If the input is an SHCoeffs class instance, this will be
            expaned on a grid using the optional parameters lmax_grid and
            lmax_calc.
        rho : int, float, or ndarray, or an SHGrid or SHCoeffs class instance
            The density contrast associated with the interface in kg / m3. If
            the input is a scalar, the density contrast is constant. If
            the input is an SHCoeffs or SHGrid class instance, the density
            contrast will vary laterally.
        gm : float
            The gravitational constant times the mass that is associated with
            the gravitational potential coefficients.
        nmax : integer, optional, default = 7
             The maximum order used in the Taylor-series expansion when
             calculating the potential coefficients.
        lmax : int, optional, shape.lmax
            The maximum spherical harmonic degree of the output spherical
            harmonic coefficients.
        lmax_grid : int, optional, default = lmax
            If shape or rho is of type SHCoeffs, this parameter determines the
            maximum spherical harmonic degree that is resolvable when expanded
            onto a grid.
        lmax_calc : optional, integer, default = lmax
            If shape or rho is of type SHCoeffs, this parameter determines the
            maximum spherical harmonic degree that will be used when expanded
            onto a grid.
        omega : float, optional, default = None
            The angular rotation rate of the body.

        Description
        -----------
        Initialize an SHGravCoeffs class instance by calculating the spherical
        harmonic coefficients of the gravitational potential associated with
        the shape of a density interface. The potential is calculated using the
        finite-amplitude technique of Wieczorek and Phillips (1998) for a
        constant density contrast and Wieczorek (2007) for a density contrast
        that varies laterally. The output coefficients are referenced to the
        mean radius of shape, and the potential is strictly valid only when it
        is evaluated at a radius greater than the maximum radius of shape.

        The input shape (and density contrast rho for variable density) can be
        either an SHGrid or SHCoeffs class instance. The routine makes direct
        use of gridded versions of these quantities, so if the input is of type
        SHCoeffs, it will first be expanded onto a grid. This exansion will be
        performed on a grid that can resolve degrees up to lmax_grid, with only
        the first lmax_calc coefficients being used. The input shape must
        correspond to absolute radii as the degree 0 term determines the
        reference radius of the coefficients.

        As an intermediate step, this routine calculates the spherical harmonic
        coefficients of the interface raised to the nth power, i.e.,
        (shape-r0)**n, where r0 is the mean radius of shape. If the input shape
        is bandlimited to degree L, the resulting function will thus be
        bandlimited to degree L*nmax. This subroutine assumes implicitly that
        the maximum spherical harmonic degree of the input shape (when
        SHCoeffs) or maximum resolvable spherical harmonic degree of shape
        (when SHGrid) is greater or equal to this value. If this is not the
        case, aliasing will occur. In practice, for accurate results, the
        effective bandwidth needs only to be about three times the size of L,
        though this should be verified for each application. The effective
        bandwidth of shape (when SHCoeffs) can be increased by preprocessing
        with the method pad(), or by increaesing the value of lmax_grid (when
        SHGrid).
        """
        mass = gm / _G.value

        if type(shape) is not _SHRealCoeffs and type(shape) is not _DHRealGrid:
            raise ValueError('shape must be of type SHRealCoeffs '
                             'or DHRealGrid. Input type is {:s}'
                             .format(repr(type(shape))))

        if (not issubclass(type(rho), float) and type(rho) is not int
                and type(rho) is not _np.ndarray and
                type(rho) is not _SHRealCoeffs and
                type(rho is not _DHRealGrid)):
            raise ValueError('rho must be of type float, int, ndarray, '
                             'SHRealCoeffs or DHRealGrid. Input type is {:s}'
                             .format(repr(type(rho))))

        if type(shape) is _SHRealCoeffs:
            shape = shape.expand(lmax=lmax_grid, lmax_calc=lmax_calc)

        if type(rho) is _SHRealCoeffs:
            rho = rho.expand(lmax=lmax_grid, lmax_calc=lmax_calc)

        if type(rho) is _DHRealGrid:
            if shape.lmax != rho.lmax:
                raise ValueError('The grids for shape and rho must have the '
                                 'same size. '
                                 'lmax of shape = {:d}, lmax of rho = {:d}'
                                 .format(shape.lmax, rho.lmax))
            cilm, d = _CilmPlusRhoHDH(shape.data, nmax, mass, rho.data,
                                      lmax=lmax)

        else:
            cilm, d = _CilmPlusDH(shape.data, nmax, mass, rho, lmax=lmax)

        clm = SHGravRealCoeffs(cilm, gm=gm, r0=d, omega=omega,
                               normalization='4pi', csphase=1)
        return clm

    @property
    def mass(self):
        """Return the mass of the planet in kg.
        """
        return self.gm / _G.value

    # ---- Define methods that modify internal variables ----
    def set_omega(self, omega):
        """
        Set the angular rotation rate of the class instance.

        Usage
        -----
        x.set_omega(omega)

        Parameters
        ----------
        omega : float
            The angular rotation rate of the body.
        """
        self.omega = omega

    def set_coeffs(self, values, ls, ms):
        """
        Set spherical harmonic coefficients in-place to specified values.

        Usage
        -----
        x.set_coeffs(values, ls, ms)

        Parameters
        ----------
        values : float (list)
            The value(s) of the spherical harmonic coefficient(s).
        ls : int (list)
            The degree(s) of the coefficient(s) that should be set.
        ms : int (list)
            The order(s) of the coefficient(s) that should be set. Positive
            and negative values correspond to the cosine and sine
            components, respectively.

        Examples
        --------
        x.set_coeffs(10., 1, 1)                   # x.coeffs[0, 1, 1] = 10.
        x.set_coeffs(5., 1, -1)                   # x.coeffs[1, 1, 1] = 5.
        x.set_coeffs([1., 2], [1, 2], [0, -2])    # x.coeffs[0, 1, 0] = 1.
                                                  # x.coeffs[1, 2, 2] = 2.
        """
        # Ensure that the type is correct
        values = _np.array(values)
        ls = _np.array(ls)
        ms = _np.array(ms)

        mneg_mask = (ms < 0).astype(_np.int)
        self.coeffs[mneg_mask, ls, _np.abs(ms)] = values

    # ---- IO routines ----
    def to_file(self, filename, format='shtools', header=None, errors=False,
                **kwargs):
        """
        Save spherical harmonic coefficients to a file.

        Usage
        -----
        x.to_file(filename, [format='shtools', header, errors])
        x.to_file(filename, [format='npy', **kwargs])

        Parameters
        ----------
        filename : str
            Name of the output file.
        format : str, optional, default = 'shtools'
            'shtools' or 'npy'. See method from_file() for more information.
        header : str, optional, default = None
            A header string written to an 'shtools'-formatted file directly
            before the spherical harmonic coefficients.
        errors : bool, optional, default = False
            If True, save the errors in the file (for 'shtools' formatted
            files only).
        **kwargs : keyword argument list, optional for format = 'npy'
            Keyword arguments of numpy.save().

        Description
        -----------
        If format='shtools', the coefficients and meta-data will be written to
        an ascii formatted file. The first line is an optional user provided
        header line, and the following line provides the attributes r0, gm,
        omega, and lmax. The spherical harmonic coefficients are then listed,
        with increasing degree and order, with the format

        l, m, coeffs[0, l, m], coeffs[1, l, m]

        where l and m are the spherical harmonic degree and order,
        respectively. If the errors are to be saved, the format of each line
        will be

        l, m, coeffs[0, l, m], coeffs[1, l, m], error[0, l, m], error[1, l, m]

        If format='npy', the spherical harmonic coefficients (but not the
        meta-data nor errors) will be saved to a binary numpy 'npy' file using
        numpy.save().
        """
        if format is 'shtools':
            if errors is True and self.errors is None:
                raise ValueError('Can not save errors when then have not been '
                                 'initialized.')

            if self.omega is None:
                omega = 0.
            else:
                omega = self.omega

            with open(filename, mode='w') as file:
                if header is not None:
                    file.write(header + '\n')
                file.write('{:.16e}, {:.16e}, {:.16e}, {:d}\n'.format(
                    self.r0, self.gm, omega, self.lmax))
                for l in range(self.lmax+1):
                    for m in range(l+1):
                        if errors is True:
                            file.write('{:d}, {:d}, {:.16e}, {:.16e}, '
                                       '{:.16e}, {:.16e}\n'
                                       .format(l, m, self.coeffs[0, l, m],
                                               self.coeffs[1, l, m],
                                               self.errors[0, l, m],
                                               self.errors[1, l, m]))
                        else:
                            file.write('{:d}, {:d}, {:.16e}, {:.16e}\n'
                                       .format(l, m, self.coeffs[0, l, m],
                                               self.coeffs[1, l, m]))
        elif format is 'npy':
            _np.save(filename, self.coeffs, **kwargs)
        else:
            raise NotImplementedError(
                'format={:s} not implemented'.format(repr(format)))

    def to_array(self, normalization=None, csphase=None, lmax=None):
        """
        Return spherical harmonic coefficients (and errors) as a numpy array.

        Usage
        -----
        coeffs, [errors] = x.to_array([normalization, csphase, lmax])

        Returns
        -------
        coeffs : ndarry, shape (2, lmax+1, lmax+1)
            numpy ndarray of the spherical harmonic coefficients.
        errors : ndarry, shape (2, lmax+1, lmax+1)
            numpy ndarray of the errors of the spherical harmonic
            coefficients if they are not None.

        Parameters
        ----------
        normalization : str, optional, default = x.normalization
            Normalization of the output coefficients: '4pi', 'ortho',
            'schmidt', or 'unnorm' for geodesy 4pi normalized, orthonormalized,
            Schmidt semi-normalized, or unnormalized coefficients,
            respectively.
        csphase : int, optional, default = x.csphase
            Condon-Shortley phase convention: 1 to exclude the phase factor,
            or -1 to include it.
        lmax : int, optional, default = x.lmax
            Maximum spherical harmonic degree to output. If lmax is greater
            than x.lmax, the array will be zero padded.

        Description
        -----------
        This method will return an array of the spherical harmonic coefficients
        using a different normalization and Condon-Shortley phase convention,
        and a different maximum spherical harmonic degree. If the maximum
        degree is smaller than the maximum degree of the class instance, the
        coefficients will be truncated. Conversely, if this degree is larger
        than the maximum degree of the class instance, the output array will be
        zero padded. If the errors of the coefficients are set, they will be
        output as a separate array.
        """
        if normalization is None:
            normalization = self.normalization
        if csphase is None:
            csphase = self.csphase
        if lmax is None:
            lmax = self.lmax

        coeffs = _convert(self.coeffs, normalization_in=self.normalization,
                          normalization_out=normalization,
                          csphase_in=self.csphase, csphase_out=csphase,
                          lmax=lmax)

        if self.errors is not None:
            errors = _convert(self.errors, normalization_in=self.normalization,
                              normalization_out=normalization,
                              csphase_in=self.csphase, csphase_out=csphase,
                              lmax=lmax)
            return coeffs, errors
        else:
            return coeffs

    def copy(self):
        """
        Return a deep copy of the class instance.

        Usage
        -----
        copy = x.copy()
        """
        return _copy.deepcopy(self)

    def info(self):
        """
        Print a summary of the data stored in the SHGravCoeffs class instance.

        Usage
        -----
        x.info()
        """
        print(repr(self))

    # -------------------------------------------------------------------------
    #    Mathematical operators
    #
    #    Operations that involve a change of units are not permitted, such as
    #    SHGravCoeffs*SHGravCoeffs, SHGravCoeffs/SHGravCoeffs, and
    #    SHGravCoeffs+SHCoeffs. All operations ignore the errors of the
    #    coefficients.
    # -------------------------------------------------------------------------
    def __add__(self, other):
        """
        Add two similar sets of gravitational potential coefficients:
        self + other.
        """
        if isinstance(other, SHGravCoeffs):
            if (self.gm == other.gm and self.r0 == other.r0 and
                    self.normalization == other.normalization and
                    self.csphase == other.csphase and self.kind == other.kind
                    and self.lmax == other.lmax):
                coeffs = _np.empty([2, self.lmax+1, self.lmax+1],
                                   dtype=self.coeffs.dtype)
                coeffs[self.mask] = (self.coeffs[self.mask] +
                                     other.coeffs[self.mask])
                return SHGravCoeffs.from_array(
                    coeffs, gm=self.gm, r0=self.r0, omega=self.omega,
                    csphase=self.csphase, normalization=self.normalization)
            else:
                raise ValueError('Addition is permitted only when the two '
                                 'SHGravCoeffs instances have the same kind, '
                                 'normalization, csphase, gm, r0, and lmax.')
        else:
            raise TypeError('Addition is permitted only for two SHGravCoeffs '
                            'instances. Type of other is {:s}'
                            .format(repr(type(other))))

    def __radd__(self, other):
        """
        Add two similar sets of gravitational potential coefficients:
        other + self.
        """
        return self.__add__(other)

    def __sub__(self, other):
        """
        Subtract two similar sets of gravitational potential coefficients:
        self - other.
        """
        if isinstance(other, SHGravCoeffs):
            if (self.gm == other.gm and self.r0 == other.r0 and
                    self.normalization == other.normalization and
                    self.csphase == other.csphase and self.kind == other.kind
                    and self.lmax == other.lmax):
                coeffs = _np.empty([2, self.lmax+1, self.lmax+1],
                                   dtype=self.coeffs.dtype)
                coeffs[self.mask] = (self.coeffs[self.mask] -
                                     other.coeffs[self.mask])
                return SHGravCoeffs.from_array(
                    coeffs, gm=self.gm, r0=self.r0, omega=self.omega,
                    csphase=self.csphase, normalization=self.normalization)
            else:
                raise ValueError('Subtraction is permitted only when the two '
                                 'SHGravCoeffs instances have the same kind, '
                                 'normalization, csphase, gm, r0, and lmax.')
        else:
            raise TypeError('Subtraction is permitted only for two '
                            'SHGravCoeffs instances. Type of other is {:s}'
                            .format(repr(type(other))))

    def __rsub__(self, other):
        """
        Subtract two similar sets of gravitational potential coefficients:
        other - self.
        """
        if isinstance(other, SHGravCoeffs):
            if (self.gm == other.gm and self.r0 == other.r0 and
                    self.normalization == other.normalization and
                    self.csphase == other.csphase and self.kind == other.kind
                    and self.lmax == other.lmax):
                coeffs = _np.empty([2, self.lmax+1, self.lmax+1],
                                   dtype=self.coeffs.dtype)
                coeffs[self.mask] = (other.coeffs[self.mask] -
                                     self.coeffs[self.mask])
                return SHGravCoeffs.from_array(
                    coeffs, gm=self.gm, r0=self.r0, omega=self.omega,
                    csphase=self.csphase, normalization=self.normalization)
            else:
                raise ValueError('Subtraction is permitted only when the two '
                                 'SHGravCoeffs instances have the same kind, '
                                 'normalization, csphase, gm, r0, and lmax.')
        else:
            raise TypeError('Subtraction is permitted only for two '
                            'SHGravCoeffs instances. Type of other is {:s}'
                            .format(repr(type(other))))

    def __mul__(self, other):
        """
        Multiply an SHGravCoeffs instance by an SHCoeffs instance or scalar:
        self * other.
        """
        if isinstance(other, _SHCoeffs):
            if (self.normalization == other.normalization and
                    self.csphase == other.csphase and self.kind == other.kind
                    and self.lmax == other.lmax):
                coeffs = _np.empty([2, self.lmax+1, self.lmax+1],
                                   dtype=self.coeffs.dtype)
                coeffs[self.mask] = (self.coeffs[self.mask] *
                                     other.coeffs[self.mask])
                return SHGravCoeffs.from_array(
                    coeffs, gm=self.gm, r0=self.r0, omega=self.omega,
                    csphase=self.csphase, normalization=self.normalization)
            else:
                raise ValueError('The two sets of coefficients must have the '
                                 'same kind, normalization, csphase, and '
                                 'lmax.')
        elif _np.isscalar(other) is True:
            if self.kind == 'real' and _np.iscomplexobj(other):
                raise ValueError('Can not multiply real gravitational '
                                 'potential coefficients by a complex '
                                 'constant.')
            coeffs = _np.empty([2, self.lmax+1, self.lmax+1],
                               dtype=self.coeffs.dtype)
            coeffs[self.mask] = self.coeffs[self.mask] * other
            return SHGravCoeffs.from_array(
                coeffs, gm=self.gm, r0=self.r0, omega=self.omega,
                csphase=self.csphase, normalization=self.normalization)
        else:
            raise TypeError('Multiplication of an SHGravCoeffs instance is '
                            'permitted only with either an SHCoeffs instance '
                            'or a scalar. '
                            'Type of other is {:s}'.format(repr(type(other))))

    def __rmul__(self, other):
        """
        Multiply an SHGravCoeffs instance by an SHCoeffs instance or scalar:
        other * self.
        """
        return self.__mul__(other)

    def __div__(self, other):
        """
        Divide an SHGravCoeffs instance by an SHCoeffs instance or scalar
        when __future__.division is not in effect: self / other.
        """
        if isinstance(other, _SHCoeffs):
            if (self.normalization == other.normalization and
                    self.csphase == other.csphase and self.kind == other.kind
                    and self.lmax == other.lmax):
                coeffs = _np.empty([2, self.lmax+1, self.lmax+1],
                                   dtype=self.coeffs.dtype)
                coeffs[self.mask] = (self.coeffs[self.mask] /
                                     other.coeffs[self.mask])
                return SHGravCoeffs.from_array(
                    coeffs, gm=self.gm, r0=self.r0, omega=self.omega,
                    csphase=self.csphase, normalization=self.normalization)
            else:
                raise ValueError('The two sets of coefficients must have the '
                                 'same kind, normalization, csphase, and '
                                 'lmax.')
        elif _np.isscalar(other) is True:
            if self.kind == 'real' and _np.iscomplexobj(other):
                raise ValueError('Can not divide real gravitational '
                                 'potential coefficients by a complex '
                                 'constant.')
            coeffs = _np.empty([2, self.lmax+1, self.lmax+1],
                               dtype=self.coeffs.dtype)
            coeffs[self.mask] = self.coeffs[self.mask] / other
            return SHGravCoeffs.from_array(
                coeffs, gm=self.gm, r0=self.r0, omega=self.omega,
                csphase=self.csphase, normalization=self.normalization)
        else:
            raise TypeError('Division of an SHGravCoeffs instance is '
                            'permitted only with either an SHCoeffs instance '
                            'or a scalar. '
                            'Type of other is {:s}'.format(repr(type(other))))

    def __truediv__(self, other):
        """
        Divide an SHGravCoeffs instance by an SHCoeffs instance or scalar
        when __future__.division is in effect: self / other.
        """
        if isinstance(other, _SHCoeffs):
            if (self.normalization == other.normalization and
                    self.csphase == other.csphase and self.kind == other.kind
                    and self.lmax == other.lmax):
                coeffs = _np.empty([2, self.lmax+1, self.lmax+1],
                                   dtype=self.coeffs.dtype)
                coeffs[self.mask] = (self.coeffs[self.mask] /
                                     other.coeffs[self.mask])
                return SHGravCoeffs.from_array(
                    coeffs, gm=self.gm, r0=self.r0, omega=self.omega,
                    csphase=self.csphase, normalization=self.normalization)
            else:
                raise ValueError('The two sets of coefficients must have the '
                                 'same kind, normalization, csphase, and '
                                 'lmax.')
        elif _np.isscalar(other) is True:
            if self.kind == 'real' and _np.iscomplexobj(other):
                raise ValueError('Can not divide real gravitational '
                                 'potential coefficients by a complex '
                                 'constant.')
            coeffs = _np.empty([2, self.lmax+1, self.lmax+1],
                               dtype=self.coeffs.dtype)
            coeffs[self.mask] = self.coeffs[self.mask] / other
            return SHGravCoeffs.from_array(
                coeffs, gm=self.gm, r0=self.r0, omega=self.omega,
                csphase=self.csphase, normalization=self.normalization)
        else:
            raise TypeError('Division of an SHGravCoeffs instance is '
                            'permitted only with either an SHCoeffs instance '
                            'or a scalar. '
                            'Type of other is {:s}'.format(repr(type(other))))

    # ---- Extract data ----
    def degrees(self):
        """
        Return a numpy array with the spherical harmonic degrees from 0 to
        lmax.

        Usage
        -----
        degrees = x.degrees()

        Returns
        -------
        degrees : ndarray, shape (lmax+1)
            1-D numpy ndarray listing the spherical harmonic degrees, where
            lmax is the maximum spherical harmonic degree.
        """
        return _np.arange(self.lmax + 1)

    def spectrum(self, function='geoid', lmax=None, unit='per_l', base=10.):
        """
        Return the spectrum as a function of spherical harmonic degree.

        Usage
        -----
        spectrum, [error_spectrum] = x.spectrum([function, lmax, unit, base])

        Returns
        -------
        spectrum : ndarray, shape (lmax+1)
            1-D numpy ndarray of the spectrum, where lmax is the maximum
            spherical harmonic degree.
        error_spectrum : ndarray, shape (lmax+1)
            1-D numpy ndarray of the error_spectrum (if the attribute errors
            is not None).

        Parameters
        ----------
        function : str, optional, default = 'geoid'
            The type of power spectrum to return: 'potential' for the
            gravitational potential in m2/s2, 'geoid' for the geoid in m,
            'radial' for the radial gravity in m/s2, or 'total' for the total
            gravitational field in m/s2.
        lmax : int, optional, default = x.lmax
            Maximum spherical harmonic degree of the spectrum to return.
        unit : str, optional, default = 'per_l'
            If 'per_l', return the total contribution to the spectrum for each
            spherical harmonic degree l. If 'per_lm', return the average
            contribution to the spectrum for each coefficient at spherical
            harmonic degree l. If 'per_dlogl', return the spectrum per log
            interval dlog_a(l).
        base : float, optional, default = 10.
            The logarithm base when calculating the 'per_dlogl' spectrum.

        Description
        -----------
        This method returns the power spectrum of the class instance, where the
        type of function is defined by the function parameter: 'potential' for
        the gravitational potential, 'geoid' for the geoid, 'radial' for
        the radial gravity, or 'total' for the total gravitational field. In
        all cases, the total power of the function is defined as the integral
        of the function squared over all space, divided by the area the
        function spans. If the mean of the function is zero, this is equivalent
        to the variance of the function.

        The output spectrum can be expresed using one of three units. 'per_l'
        returns the contribution to the total spectrum from all angular orders
        at degree l. 'per_lm' returns the average contribution to the total
        spectrum from a single coefficient at degree l, which is equal to the
        'per_l' spectrum divided by (2l+1). 'per_dlogl' returns the
        contribution to the total spectrum from all angular orders over an
        infinitessimal logarithmic degree band. The contrubution in the band
        dlog_a(l) is spectrum(l, 'per_dlogl')*dlog_a(l), where a is the base,
        and where spectrum(l, 'per_dlogl) is equal to
        spectrum(l, 'per_l')*l*log(a).
        """
        if function.lower() not in ('potential', 'geoid', 'radial', 'total'):
            raise ValueError(
                "function must be of type 'potential', 'geoid', 'radial', or "
                "'total'. Provided value was {:s}".format(repr(function))
                )

        s = _spectrum(self.coeffs, normalization=self.normalization,
                      convention='power', unit=unit, base=base, lmax=lmax)

        if self.errors is not None:
            es = _spectrum(self.errors, normalization=self.normalization,
                           convention='power', unit=unit, base=base, lmax=lmax)

        if function.lower() == 'potential':
            s *= (self.gm / self.r0)**2
            if self.errors is not None:
                es *= (self.gm / self.r0)**2
        elif function.lower() == 'geoid':
            s *= self.r0**2
            if self.errors is not None:
                es *= self.r0**2
        elif function.lower() == 'radial':
            degrees = _np.arange(len(s))
            s *= (self.gm * (degrees + 1) / self.r0**2)**2
            if self.errors is not None:
                es *= (self.gm * (degrees + 1) / self.r0**2)**2
        elif function.lower() == 'total':
            degrees = _np.arange(len(s))
            s *= (self.gm / self.r0**2)**2 * (degrees + 1) * (2 * degrees + 1)
            if self.errors is not None:
                es *= (self.gm / self.r0**2)**2 * (degrees + 1) * \
                    (2 * degrees + 1)

        if self.errors is not None:
            return s, es
        else:
            return s

    # ---- Operations that return a new SHGravCoeffs class instance ----
    def rotate(self, alpha, beta, gamma, degrees=True, convention='y',
               body=False, dj_matrix=None):
        """
        Rotate either the coordinate system used to express the spherical
        harmonic coefficients or the physical body, and return a new class
        instance.

        Usage
        -----
        x_rotated = x.rotate(alpha, beta, gamma, [degrees, convention,
                             body, dj_matrix])

        Returns
        -------
        x_rotated : SHGravCoeffs class instance

        Parameters
        ----------
        alpha, beta, gamma : float
            The three Euler rotation angles in degrees.
        degrees : bool, optional, default = True
            True if the Euler angles are in degrees, False if they are in
            radians.
        convention : str, optional, default = 'y'
            The convention used for the rotation of the second angle, which
            can be either 'x' or 'y' for a rotation about the x or y axes,
            respectively.
        body : bool, optional, default = False
            If true, rotate the physical body and not the coordinate system.
        dj_matrix : ndarray, optional, default = None
            The djpi2 rotation matrix computed by a call to djpi2.

        Description
        -----------
        This method will take the spherical harmonic coefficients of a
        function, rotate the coordinate frame by the three Euler anlges, and
        output the spherical harmonic coefficients of the new function. If
        the optional parameter body is set to True, then the physical body will
        be rotated instead of the coordinate system.

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

        Here, the 'y convention' is employed, where the second rotation is with
        respect to the y axis. When using the 'x convention', the second
        rotation is instead with respect to the x axis. The relation between
        the Euler angles in the x and y conventions is given by

        alpha_y=alpha_x-pi/2, beta_y=beta_x, and gamma_y=gamma_x+pi/2.

        To perform the inverse transform associated with the three angles
        (alpha, beta, gamma), one would perform an additional rotation using
        the angles (-gamma, -beta, -alpha).

        The rotations can be viewed either as a rotation of the coordinate
        system or the physical body. To rotate the physical body without
        rotation of the coordinate system, set the optional parameter body to
        True. This rotation is accomplished by performing the inverse rotation
        using the angles (-gamma, -beta, -alpha).
        """
        if type(convention) != str:
            raise ValueError('convention must be a string. '
                             'Input type was {:s}'
                             .format(str(type(convention))))

        if convention.lower() not in ('x', 'y'):
            raise ValueError(
                "convention must be either 'x' or 'y'. "
                "Provided value was {:s}".format(repr(convention))
                )

        if convention is 'y':
            if body is True:
                angles = _np.array([-gamma, -beta, -alpha])
            else:
                angles = _np.array([alpha, beta, gamma])
        elif convention is 'x':
            if body is True:
                angles = _np.array([-gamma - _np.pi/2, -beta,
                                    -alpha + _np.pi/2])
            else:
                angles = _np.array([alpha - _np.pi/2, beta, gamma + _np.pi/2])

        if degrees:
            angles = _np.radians(angles)

        if self.lmax > 1200:
            _warnings.warn("The rotate() method is accurate only to about"
                           " spherical harmonic degree 1200. "
                           "lmax = {:d}".format(self.lmax),
                           category=RuntimeWarning)

        rot = self._rotate(angles, dj_matrix, gm=self.gm, r0=self.r0,
                           omega=self.omega)
        return rot

    def convert(self, normalization=None, csphase=None, lmax=None):
        """
        Return an SHGravCoeffs class instance with a different normalization
        convention.

        Usage
        -----
        clm = x.convert([normalization, csphase, lmax])

        Returns
        -------
        clm : SHGravCoeffs class instance

        Parameters
        ----------
        normalization : str, optional, default = x.normalization
            Normalization of the output class: '4pi', 'ortho', 'schmidt', or
            'unnorm', for geodesy 4pi normalized, orthonormalized, Schmidt
            semi-normalized, or unnormalized coefficients, respectively.
        csphase : int, optional, default = x.csphase
            Condon-Shortley phase convention for the output class: 1 to exclude
            the phase factor, or -1 to include it.
        lmax : int, optional, default = x.lmax
            Maximum spherical harmonic degree to output.

        Description
        -----------
        This method will return a new class instance of the spherical
        harmonic coefficients using a different normalization and
        Condon-Shortley phase convention. A different maximum spherical
        harmonic degree of the output coefficients can be specified, and if
        this maximum degree is smaller than the maximum degree of the original
        class, the coefficients will be truncated. Conversely, if this degree
        is larger than the maximum degree of the original class, the
        coefficients of the new class will be zero padded.
        """
        if normalization is None:
            normalization = self.normalization
        if csphase is None:
            csphase = self.csphase
        if lmax is None:
            lmax = self.lmax

        # check argument consistency
        if type(normalization) != str:
            raise ValueError('normalization must be a string. '
                             'Input type was {:s}'
                             .format(str(type(normalization))))
        if normalization.lower() not in ('4pi', 'ortho', 'schmidt', 'unnorm'):
            raise ValueError(
                "normalization must be '4pi', 'ortho', 'schmidt', or "
                "'unnorm'. Provided value was {:s}"
                .format(repr(normalization)))
        if csphase != 1 and csphase != -1:
            raise ValueError(
                "csphase must be 1 or -1. Input value was {:s}"
                .format(repr(csphase)))

        if self.errors is not None:
            coeffs, errors = self.to_array(normalization=normalization.lower(),
                                           csphase=csphase, lmax=lmax)
            return SHGravCoeffs.from_array(
                coeffs, gm=self.gm, r0=self.r0, omega=self.omega,
                errors=errors, normalization=normalization.lower(),
                csphase=csphase, copy=False)
        else:
            coeffs = self.to_array(normalization=normalization.lower(),
                                   csphase=csphase, lmax=lmax)
            return SHGravCoeffs.from_array(
                coeffs, gm=self.gm, r0=self.r0, omega=self.omega,
                normalization=normalization.lower(), csphase=csphase,
                copy=False)

    def pad(self, lmax):
        """
        Return an SHGravCoeffs class where the coefficients are zero padded or
        truncated to a different lmax.

        Usage
        -----
        clm = x.pad(lmax)

        Returns
        -------
        clm : SHGravCoeffs class instance

        Parameters
        ----------
        lmax : int
            Maximum spherical harmonic degree to output.
        """
        clm = self.copy()

        if lmax <= self.lmax:
            clm.coeffs = clm.coeffs[:, :lmax+1, :lmax+1]
            clm.mask = clm.mask[:, :lmax+1, :lmax+1]
            if self.errors is not None:
                clm.errors = clm.errors[:, :lmax+1, :lmax+1]
        else:
            clm.coeffs = _np.pad(clm.coeffs, ((0, 0), (0, lmax - self.lmax),
                                 (0, lmax - self.lmax)), 'constant')
            if self.errors is not None:
                clm.errors = _np.pad(
                    clm.errors, ((0, 0), (0, lmax - self.lmax),
                                 (0, lmax - self.lmax)), 'constant')
            mask = _np.zeros((2, lmax + 1, lmax + 1), dtype=_np.bool)
            for l in _np.arange(lmax + 1):
                mask[:, l, :l + 1] = True
            mask[1, :, 0] = False
            clm.mask = mask

        clm.lmax = lmax
        return clm

    def change_ref(self, gm=None, r0=None, lmax=None):
        """
        Return a new SHGravCoeffs class instance with a different reference gm
        or r0.

        Usage
        -----
        clm = x.change_ref([gm, r0, lmax])

        Returns
        -------
        clm : SHGravCoeffs class instance.

        Parameters
        ----------
        gm : float, optional, default = self.gm
            The gravitational constant time the mass that is associated with
            the gravitational potential coefficients.
        r0 : float, optional, default = self.r0
            The reference radius of the spherical harmonic coefficients.
        lmax : int, optional, default = self.lmax
            Maximum spherical harmonic degree to output.

        Description
        -----------
        This method returns a new class instance of the gravitational
        potential, but using a difference reference gm or r0. When
        changing the reference radius r0, the spherical harmonic coefficients
        will be upward or downward continued under the assumption that the
        reference radius is exterior to the body.
        """
        if lmax is None:
            lmax = self.lmax

        clm = self.pad(lmax)

        if gm is not None and gm != self.gm:
            clm.coeffs *= self.gm / gm
            clm.gm = gm
            if self.errors is not None:
                clm.errors *= self.gm / gm

        if r0 is not None and r0 != self.r0:
            for l in _np.arange(lmax+1):
                clm.coeffs[:, l, :l+1] *= (self.r0 / r0)**l
                if self.errors is not None:
                    clm.errors[:, l, :l+1] *= (self.r0 / r0)**l
            clm.r0 = r0

        return clm

    # ---- Routines that return different gravity-related class instances ----
    def expand(self, a=None, f=None, lmax=None, lmax_calc=None,
               normal_gravity=True, sampling=2):
        """
        Create 2D cylindrical maps on a flattened and rotating ellipsoid of all
        three components of the gravity field, the gravity disturbance, and the
        gravitational potential, and return as a SHGravGrid class instance.

        Usage
        -----
        grav = x.expand([a, f, lmax, lmax_calc, normal_gravity, sampling])

        Returns
        -------
        grav : SHGravGrid class instance.

        Parameters
        ----------
        a : optional, float, default = self.r0
            The semi-major axis of the flattened ellipsoid on which the field
            is computed.
        f : optional, float, default = 0
            The flattening of the reference ellipsoid: f=(a-b)/a.
        lmax : optional, integer, default = self.lmax
            The maximum spherical harmonic degree, which determines the number
            of samples of the output grids, n=2lmax+2, and the latitudinal
            sampling interval, 90/(lmax+1).
        lmax_calc : optional, integer, default = lmax
            The maximum spherical harmonic degree used in evaluating the
            functions. This must be less than or equal to lmax.
        normal_gravity : optional, bool, default = True
            If True, the normal gravity (the gravitational acceleration on the
            ellipsoid) will be subtracted from the total gravity, yielding the
            "gravity disturbance." This is done using Somigliana's formula
            (after converting geocentric to geodetic coordinates).
        sampling : optional, integer, default = 2
            If 1 the output grids are equally sampled (n by n). If 2 (default),
            the grids are equally spaced in degrees (n by 2n).

        Description
        -----------
        This method will create 2-dimensional cylindrical maps of the three
        components of the gravity field, the total field, and the gravitational
        potential, and return these as an SHGravGrid class instance. Each
        map is stored as an SHGrid class instance using Driscoll and Healy
        grids that are either equally sampled (n by n) or equally spaced
        (n by 2n) in latitude and longitude. All grids use geocentric
        coordinates, the output is in SI units, and the sign of the radial
        components is positive when directed upwards. If the optional angular
        rotation rate omega is specified in the SHGravCoeffs instance, the
        potential and radial gravitational acceleration will be calculated in a
        body-fixed rotating reference frame. If normal_gravity is set to True,
        the normal gravity will be removed from the total field, yielding the
        gravity disturbance.

        The gravitational potential is given by

            V = GM/r Sum_{l=0}^lmax (r0/r)^l Sum_{m=-l}^l C_{lm} Y_{lm},

        and the gravitational acceleration is

            B = Grad V.

        The coefficients are referenced to the radius r0, and the function is
        computed on a flattened ellipsoid with semi-major axis a (i.e., the
        mean equatorial radius) and flattening f. To convert m/s^2 to mGals,
        multiply the gravity grids by 10^5.
        """
        if a is None:
            a = self.r0
        if f is None:
            f = 0.
        if normal_gravity is True:
            ng = 1
        else:
            ng = 0
        if lmax is None:
            lmax = self.lmax
        if lmax_calc is None:
            lmax_calc = lmax

        if self.errors is not None:
            coeffs, errors = self.to_array(normalization='4pi', csphase=1)
        else:
            coeffs = self.to_array(normalization='4pi', csphase=1)

        rad, theta, phi, total, pot = _MakeGravGridDH(
            coeffs, self.gm, self.r0, a=a, f=f, lmax=lmax,
            lmax_calc=lmax_calc, sampling=sampling, omega=self.omega,
            normal_gravity=ng)

        return _SHGravGrid(rad, theta, phi, total, pot, self.gm, a, f,
                           self.omega, normal_gravity, lmax, lmax_calc)

    def tensor(self, a=None, f=None, lmax=None, lmax_calc=None, degree0=False,
               sampling=2):
        """
        Create 2D cylindrical maps on a flattened ellipsoid of the 9
        components of the gravity "gradient" tensor in a local north-oriented
        reference frame, and return an SHGravTensor class instance.

        Usage
        -----
        tensor = x.tensor([a, f, lmax, lmax_calc, sampling])

        Returns
        -------
        tensor : SHGravTensor class instance.

        Parameters
        ----------
        a : optional, float, default = self.r0
            The semi-major axis of the flattened ellipsoid on which the field
            is computed.
        f : optional, float, default = 0
            The flattening of the reference ellipsoid: f=(a-b)/a.
        lmax : optional, integer, default = self.lmax
            The maximum spherical harmonic degree that determines the number of
            samples of the output grids, n=2lmax+2, and the latitudinal
            sampling interval, 90/(lmax+1).
        lmax_calc : optional, integer, default = lmax
            The maximum spherical harmonic degree used in evaluating the
            functions. This must be less than or equal to lmax.
        degree0 : optional, default = False
            If True, include the degree-0 term when calculating the tensor. If
            False, set the degree-0 term to zero.
        sampling : optional, integer, default = 2
            If 1 the output grids are equally sampled (n by n). If 2 (default),
            the grids are equally spaced in degrees (n by 2n).

        Description
        -----------
        This method will create 2-dimensional cylindrical maps for the 9
        components of the gravity 'gradient' tensor and return an SHGravTensor
        class instance. The components are

            (Vxx, Vxy, Vxz)
            (Vyx, Vyy, Vyz)
            (Vzx, Vzy, Vzz)

        where the reference frame is north-oriented, where x points north, y
        points west, and z points upward (all tangent or perpendicular to a
        sphere of radius r, where r is the local radius of the flattened
        ellipsoid). The gravitational potential is defined as

            V = GM/r Sum_{l=0}^lmax (r0/r)^l Sum_{m=-l}^l C_{lm} Y_{lm},

        where r0 is the reference radius of the spherical harmonic coefficients
        Clm, and the gravitational acceleration is

            B = Grad V.

        The components of the gravity tensor are calculated according to eq. 1
        in Petrovskaya and Vershkov (2006), which is based on eq. 3.28 in Reed
        (1973) (noting that Reed's equations are in terms of latitude and that
        the y axis points east):

            Vzz = Vrr
            Vxx = 1/r Vr + 1/r^2 Vtt
            Vyy = 1/r Vr + 1/r^2 /tan(t) Vt + 1/r^2 /sin(t)^2 Vpp
            Vxy = 1/r^2 /sin(t) Vtp - cos(t)/sin(t)^2 /r^2 Vp
            Vxz = 1/r^2 Vt - 1/r Vrt
            Vyz = 1/r^2 /sin(t) Vp - 1/r /sin(t) Vrp

        where r, t, p stand for radius, theta, and phi, respectively, and
        subscripts on V denote partial derivatives. The output grids are in
        units of Eotvos (10^-9 s^-2).

        References
        ----------
        Reed, G.B., Application of kinematical geodesy for determining
        the short wave length components of the gravity field by satellite
        gradiometry, Ohio State University, Dept. of Geod. Sciences, Rep. No.
        201, Columbus, Ohio, 1973.

        Petrovskaya, M.S. and A.N. Vershkov, Non-singular expressions for the
        gravity gradients in the local north-oriented and orbital reference
        frames, J. Geod., 80, 117-127, 2006.
        """
        if a is None:
            a = self.r0
        if f is None:
            f = 0.
        if lmax is None:
            lmax = self.lmax
        if lmax_calc is None:
            lmax_calc = lmax

        if self.errors is not None:
            coeffs, errors = self.to_array(normalization='4pi', csphase=1)
        else:
            coeffs = self.to_array(normalization='4pi', csphase=1)

        if degree0 is False:
            coeffs[0, 0, 0] = 0.

        vxx, vyy, vzz, vxy, vxz, vyz = _MakeGravGradGridDH(
            coeffs, self.gm, self.r0, a=a, f=f, lmax=lmax,
            lmax_calc=lmax_calc, sampling=sampling)

        return _SHGravTensor(1.e9*vxx, 1.e9*vyy, 1.e9*vzz, 1.e9*vxy, 1.e9*vxz,
                             1.e9*vyz, self.gm, a, f, lmax, lmax_calc)

    def geoid(self, potref, a=None, f=None, r=None, omega=None, order=2,
              lmax=None, lmax_calc=None, grid='DH2'):
        """
        Create a global map of the height of the geoid and return an SHGeoid
        class instance.

        Usage
        -----
        geoid = x.geoid(potref, [a, f, r, omega, order, lmax, lmax_calc, grid])

        Returns
        -------
        geoid : SHGeoid class instance.

        Parameters
        ----------
        potref : float
            The value of the potential on the chosen geoid, in m2 / s2.
        a : optional, float, default = self.r0
            The semi-major axis of the flattened ellipsoid on which the field
            is computed.
        f : optional, float, default = 0
            The flattening of the reference ellipsoid: f=(a-b)/a.
        r : optional, float, default = self.r0
            The radius of the reference sphere that the Taylor expansion of the
            potential is calculated on.
        order : optional, integer, default = 2
            The order of the Taylor series expansion of the potential about the
            reference radius r. This can be either 1, 2, or 3.
        omega : optional, float, default = self.omega
            The angular rotation rate of the planet.
        lmax : optional, integer, default = self.lmax
            The maximum spherical harmonic degree that determines the number
            of samples of the output grid, n=2lmax+2, and the latitudinal
            sampling interval, 90/(lmax+1).
        lmax_calc : optional, integer, default = lmax
            The maximum spherical harmonic degree used in evaluating the
            functions. This must be less than or equal to lmax.
        grid : str, optional, default = 'DH2'
            'DH' or 'DH1' for an equisampled lat/lon grid with nlat=nlon, or
            'DH2' for an equidistant lat/lon grid with nlon=2*nlat.

        Description
        -----------
        This method will create a global map of the geoid height, accurate to
        either first, second, or third order, using the method described in
        Wieczorek (2007; equation 19-20). The algorithm expands the potential
        in a Taylor series on a spherical interface of radius r, and computes
        the height above this interface to the potential potref exactly from
        the linear, quadratic, or cubic equation at each grid point. If the
        optional parameters a and f are specified, the geoid height will be
        referenced to a flattened ellipsoid with semi-major axis a and
        flattening f. The pseudo-rotational potential is explicitly accounted
        for by using the angular rotation rate omega of the planet in the
        SHGravCoeffs class instance. If omega is explicitly specified for this
        method, it will override the value present in the class instance.

        Reference
        ----------
        Wieczorek, M. A. Gravity and topography of the terrestrial planets,
        Treatise on Geophysics, 10, 165-206, 2007.
        """
        if a is None:
            a = self.r0
        if f is None:
            f = 0.
        if r is None:
            r = self.r0
        if lmax is None:
            lmax = self.lmax
        if lmax_calc is None:
            lmax_calc = lmax

        if grid.upper() in ('DH', 'DH1'):
            sampling = 1
        elif grid.upper() == 'DH2':
            sampling = 2
        else:
            raise ValueError(
                    "grid must be 'DH', 'DH1', or 'DH2'. "
                    "Input value was {:s}".format(repr(grid)))

        if self.errors is not None:
            coeffs, errors = self.to_array(normalization='4pi', csphase=1)
        else:
            coeffs = self.to_array(normalization='4pi', csphase=1)

        if omega is None:
            omega = self.omega

        geoid = _MakeGeoidGridDH(coeffs, self.r0, self.gm, potref, lmax=lmax,
                                 omega=omega, r=r, order=order,
                                 lmax_calc=lmax_calc, a=a, f=f,
                                 sampling=sampling)

        return _SHGeoid(geoid, self.gm, potref, a, f, omega, r, order,
                        lmax, lmax_calc)

    # ---- Plotting routines ----
    def plot_spectrum(self, function='geoid', unit='per_l', base=10.,
                      lmax=None, xscale='lin', yscale='log', grid=True,
                      legend=None, axes_labelsize=None, tick_labelsize=None,
                      show=True, ax=None, fname=None, **kwargs):
        """
        Plot the spectrum as a function of spherical harmonic degree.

        Usage
        -----
        x.plot_spectrum([function, unit, base, lmax, xscale, yscale, grid,
                         legend, axes_labelsize, tick_labelsize, show, ax,
                         fname, **kwargs])

        Parameters
        ----------
        function : str, optional, default = 'geoid'
            The type of power spectrum to calculate: 'potential' for the
            gravitational potential, 'geoid' for the geoid, 'radial' for
            the radial gravity, or 'total' for the total gravitational field.
        unit : str, optional, default = 'per_l'
            If 'per_l', plot the total contribution to the spectrum for each
            spherical harmonic degree l. If 'per_lm', plot the average
            contribution to the spectrum for each coefficient at spherical
            harmonic degree l. If 'per_dlogl', plot the spectrum per log
            interval dlog_a(l).
        base : float, optional, default = 10.
            The logarithm base when calculating the 'per_dlogl' spectrum, and
            the base to use for logarithmic axes.
        lmax : int, optional, default = self.lmax
            The maximum spherical harmonic degree to plot.
        xscale : str, optional, default = 'lin'
            Scale of the x axis: 'lin' for linear or 'log' for logarithmic.
        yscale : str, optional, default = 'log'
            Scale of the y axis: 'lin' for linear or 'log' for logarithmic.
        grid : bool, optional, default = True
            If True, plot grid lines.
        legend : str, optional, default = None
            Text to use for the legend.
        axes_labelsize : int, optional, default = None
            The font size for the x and y axes labels.
        tick_labelsize : int, optional, default = None
            The font size for the x and y tick labels.
        show : bool, optional, default = True
            If True, plot to the screen.
        ax : matplotlib axes object, optional, default = None
            A single matplotlib axes object where the plot will appear.
        fname : str, optional, default = None
            If present, and if axes is not specified, save the image to the
            specified file.
        **kwargs : keyword arguments, optional
            Keyword arguments for pyplot.plot().

        Description
        -----------
        This method plots the power (and error) spectrum of the class instance,
        where the type of spectrum is defined by the parameter function:
        'potential' for the gravitational potential, 'geoid' for the geoid,
        'radial' for the radial gravity, or 'total' for the total gravitational
        field. The power for the degree 0 and 1 terms are not plotted. In all
        cases, the total power of the function is defined as the integral of
        the function squared over all space, divided by the area the function
        spans. If the mean of the function is zero, this is equivalent to the
        variance of the function.

        The output spectrum can be expresed using one of three units. 'per_l'
        returns the contribution to the total spectrum from all angular orders
        at degree l. 'per_lm' returns the average contribution to the total
        spectrum from a single coefficient at degree l, which is equal to the
        'per_l' spectrum divided by (2l+1). 'per_dlogl' returns the
        contribution to the total spectrum from all angular orders over an
        infinitessimal logarithmic degree band. The contrubution in the band
        dlog_a(l) is spectrum(l, 'per_dlogl')*dlog_a(l), where a is the base,
        and where spectrum(l, 'per_dlogl) is equal to
        spectrum(l, 'per_l')*l*log(a).
        """
        if lmax is None:
            lmax = self.lmax

        if self.errors is not None:
            spectrum, error_spectrum = self.spectrum(function=function,
                                                     unit=unit, base=base,
                                                     lmax=lmax)
            if function == 'radial' or function == 'total':
                spectrum *= 1.e10
                error_spectrum *= 1.e10
        else:
            spectrum = self.spectrum(function=function, unit=unit, base=base,
                                     lmax=lmax)
            if function == 'radial' or function == 'total':
                spectrum *= 1.e10

        ls = _np.arange(lmax + 1)

        if ax is None:
            fig, axes = _plt.subplots(1, 1)
        else:
            axes = ax

        if axes_labelsize is None:
            axes_labelsize = _mpl.rcParams['axes.labelsize']
        if tick_labelsize is None:
            tick_labelsize = _mpl.rcParams['xtick.labelsize']

        axes.set_xlabel('Spherical harmonic degree', fontsize=axes_labelsize)

        if function == 'geoid':
            axes.set_ylabel('Power, m$^2$', fontsize=axes_labelsize)
        elif function == 'potential':
            axes.set_ylabel('Power, m$^4$ s$^{-4}$', fontsize=axes_labelsize)
        elif function == 'radial':
            axes.set_ylabel('Power, mGal$^2$', fontsize=axes_labelsize)
        elif function == 'total':
            axes.set_ylabel('Power, mGal$^2$', fontsize=axes_labelsize)

        if legend is None:
            if (unit == 'per_l'):
                legend = 'power per degree'
            elif (unit == 'per_lm'):
                legend = 'power per coefficient'
            elif (unit == 'per_dlogl'):
                legend = 'power per log bandwidth'

        if xscale == 'log':
            axes.set_xscale('log', basex=base)
        if yscale == 'log':
            axes.set_yscale('log', basey=base)

        if self.errors is not None:
            axes.plot(ls[2:lmax + 1], spectrum[2:lmax + 1], label=legend,
                      **kwargs)
            axes.plot(ls[2:lmax + 1], error_spectrum[2:lmax + 1],
                      label='error', **kwargs)
        else:
            axes.plot(ls[2:lmax + 1], spectrum[2: lmax + 1], label=legend,
                      **kwargs)

        if xscale == 'lin':
            axes.set(xlim=(ls[0], ls[lmax]))

        axes.grid(grid, which='major')
        axes.minorticks_on()
        axes.tick_params(labelsize=tick_labelsize)
        axes.legend()

        if ax is None:
            fig.tight_layout(pad=0.5)
            if show:
                fig.show()
            if fname is not None:
                fig.savefig(fname)
            return fig, axes

    def plot_spectrum2d(self, function='geoid', xscale='lin', yscale='lin',
                        grid=True, axes_labelsize=None, tick_labelsize=None,
                        vscale='log', vrange=None, vmin=None, vmax=None,
                        lmax=None, errors=False, show=True, ax=None,
                        fname=None):
        """
        Plot the spectrum as a function of spherical harmonic degree and order.

        Usage
        -----
        x.plot_spectrum2d([function, xscale, yscale, grid, axes_labelsize,
                           tick_labelsize, vscale, vrange, vmin, vmax, lmax,
                           errors, show, ax, fname])

        Parameters
        ----------
        function : str, optional, default = 'geoid'
            The type of power spectrum to calculate: 'potential' for the
            gravitational potential, 'geoid' for the geoid, 'radial' for
            the radial gravity, or 'total' for the total gravitational field.
        xscale : str, optional, default = 'lin'
            Scale of the l axis: 'lin' for linear or 'log' for logarithmic.
        yscale : str, optional, default = 'lin'
            Scale of the m axis: 'lin' for linear or 'log' for logarithmic.
        grid : bool, optional, default = True
            If True, plot grid lines.
        axes_labelsize : int, optional, default = None
            The font size for the x and y axes labels.
        tick_labelsize : int, optional, default = None
            The font size for the x and y tick labels.
        vscale : str, optional, default = 'log'
            Scale of the color axis: 'lin' for linear or 'log' for logarithmic.
        vrange : (float, float), optional, default = None
            Colormap range (min, max), relative to the maximum value. If None,
            scale the image to the maximum and minimum values.
        vmin : float, optional, default=None
            The minmum range of the colormap. If None, the minimum value of the
            spectrum will be used.
        vmax : float, optional, default=None
            The maximum range of the colormap. If None, the maximum value of
            the spectrum will be used.
        lmax : int, optional, default = self.lmax
            The maximum spherical harmonic degree to plot.
        errors : bool, optional, default = False
            If True, plot the spectrum of the errors.
        show : bool, optional, default = True
            If True, plot to the screen.
        ax : matplotlib axes object, optional, default = None
            A single matplotlib axes object where the plot will appear.
        fname : str, optional, default = None
            If present, and if axes is not specified, save the image to the
            specified file.

        Description
        -----------
        This method plots the power of the class instance for each spherical
        harmonic degree and order, where the type of spectrum is defined by
        the parameter function: 'potential' for the gravitational potential,
        'geoid' for the geoid, 'radial' for the radial gravity, or 'total' for
        the total gravitational field. In all cases, the total power of the
        function is defined as the integral of the function squared over all
        space, divided by the area the function spans. If the mean of the
        function is zero, this is equivalent to the variance of the function.
        """
        if axes_labelsize is None:
            axes_labelsize = _mpl.rcParams['axes.labelsize']
        if tick_labelsize is None:
            tick_labelsize = _mpl.rcParams['xtick.labelsize']

        if lmax is None:
            lmax = self.lmax
        degrees = _np.arange(lmax + 1)

        # Create the matrix of the spectrum for each coefficient
        if errors is True:
            if self.errors is None:
                raise ValueError('Can not plot the error spectrum when the '
                                 'errors are not set.')
            coeffs = self.errors
        else:
            coeffs = self.coeffs

        spectrum = _np.empty((lmax + 1, 2 * lmax + 1))
        mpositive = _np.abs(coeffs[0, :lmax + 1, :lmax + 1])**2
        mpositive[0, 0] = 0.
        mpositive[~self.mask[0, :lmax + 1, :lmax + 1]] = _np.nan
        mnegative = _np.abs(coeffs[1, :lmax + 1, :lmax + 1])**2
        mnegative[~self.mask[1, :lmax + 1, :lmax + 1]] = _np.nan

        spectrum[:, :lmax] = _np.fliplr(mnegative)[:, :lmax]
        spectrum[:, lmax:] = mpositive

        if self.normalization == '4pi':
            pass
        elif self.normalization == 'schmidt':
            for l in degrees:
                spectrum[l, :] /= (2. * l + 1.)
        elif self.normalization == 'ortho':
            for l in degrees:
                spectrum[l, :] /= (4. * _np.pi)
        elif self.normalization == 'unnorm':
            for l in degrees:
                ms = _np.arange(l+1)
                conv = _factorial(l+ms) / (2. * l + 1.) / _factorial(l-ms)
                if self.kind == 'real':
                    conv[1:l + 1] = conv[1:l + 1] / 2.
                spectrum[l, lmax-l:lmax] *= conv[::-1][0:l]
                spectrum[l, lmax:lmax+l+1] *= conv[0:l+1]
        else:
            raise ValueError(
                "normalization must be '4pi', 'ortho', 'schmidt', " +
                "or 'unnorm'. Input value was {:s}"
                .format(repr(self.normalization)))

        if function == 'geoid':
            spectrum *= self.r0**2
        elif function == 'potential':
            spectrum *= (self.gm / self.r0)**2
        elif function == 'radial':
            for l in degrees:
                spectrum[l, :] *= 1.e10 * (self.gm * (l + 1) /
                                           self.r0**2)**2
        elif function == 'total':
            for l in degrees:
                spectrum[l, :] *= 1.e10 * (self.gm / self.r0**2)**2 * \
                    (l + 1) * (2 * l + 1)

        spectrum_masked = _np.ma.masked_invalid(spectrum)

        # need to add one extra value to each in order for pcolormesh
        # to plot the last row and column.
        ls = _np.arange(lmax+2).astype(_np.float)
        ms = _np.arange(-lmax, lmax + 2, dtype=_np.float)
        lgrid, mgrid = _np.meshgrid(ls, ms, indexing='ij')
        lgrid -= 0.5
        mgrid -= 0.5

        if ax is None:
            fig, axes = _plt.subplots()
        else:
            axes = ax

        if vrange is not None:
            vmin = _np.nanmax(spectrum) * vrange[0]
            vmax = _np.nanmax(spectrum) * vrange[1]
        else:
            if vmin is None:
                _temp = spectrum
                _temp[_temp == 0] = _np.NaN
                vmin = _np.nanmin(_temp)
            if vmax is None:
                vmax = _np.nanmax(spectrum)

        if vscale.lower() == 'log':
            norm = _mpl.colors.LogNorm(vmin, vmax, clip=True)
            # Clipping is required to avoid an invalid value error
        elif vscale.lower() == 'lin':
            norm = _plt.Normalize(vmin, vmax)
        else:
            raise ValueError(
                "vscale must be 'lin' or 'log'. " +
                "Input value was {:s}".format(repr(vscale)))

        if (xscale == 'lin'):
            cmesh = axes.pcolormesh(lgrid, mgrid, spectrum_masked,
                                    norm=norm, cmap='viridis')
            axes.set(xlim=(-0.5, lmax + 0.5))
        elif (xscale == 'log'):
            cmesh = axes.pcolormesh(lgrid[1:], mgrid[1:], spectrum_masked[1:],
                                    norm=norm, cmap='viridis')
            axes.set(xscale='log', xlim=(1., lmax + 0.5))
        else:
            raise ValueError(
                "xscale must be 'lin' or 'log'. " +
                "Input value was {:s}".format(repr(xscale)))

        if (yscale == 'lin'):
            axes.set(ylim=(-lmax - 0.5, lmax + 0.5))
        elif (yscale == 'log'):
            axes.set(yscale='symlog', ylim=(-lmax - 0.5, lmax + 0.5))
        else:
            raise ValueError(
                "yscale must be 'lin' or 'log'. " +
                "Input value was {:s}".format(repr(yscale)))

        cb = _plt.colorbar(cmesh, ax=ax)

        if function == 'geoid':
            cb.set_label('Power, m$^2$', fontsize=axes_labelsize)
        elif function == 'potential':
            cb.set_label('Power, m$^4$ s$^{-4}$', fontsize=axes_labelsize)
        elif function == 'radial':
            cb.set_label('Power, mGal$^2$', fontsize=axes_labelsize)
        elif function == 'total':
            cb.set_label('Power, mGal$^2$', fontsize=axes_labelsize)

        cb.ax.tick_params(labelsize=tick_labelsize)
        axes.set_xlabel('Spherical harmonic degree', fontsize=axes_labelsize)
        axes.set_ylabel('Spherical harmonic order', fontsize=axes_labelsize)
        axes.minorticks_on()
        axes.grid(grid, which='major')

        if ax is None:
            fig.tight_layout(pad=0.5)
            if show:
                fig.show()
            if fname is not None:
                fig.savefig(fname)
            return fig, axes


class SHGravRealCoeffs(SHGravCoeffs):
    """
    Real spherical harmonic coefficient class for the gravitational potential.
    """

    def __init__(self, coeffs, gm=None, r0=None, omega=None, errors=None,
                 normalization='4pi', csphase=1, copy=True, header=None):
        """Initialize real gravitational potential coefficients class."""
        lmax = coeffs.shape[1] - 1
        # ---- create mask to filter out m<=l ----
        mask = _np.zeros((2, lmax + 1, lmax + 1), dtype=_np.bool)
        mask[0, 0, 0] = True
        for l in _np.arange(lmax + 1):
            mask[:, l, :l + 1] = True
        mask[1, :, 0] = False
        self.mask = mask
        self.lmax = lmax
        self.kind = 'real'
        self.normalization = normalization
        self.csphase = csphase
        self.header = header
        self.gm = gm
        self.r0 = r0
        self.omega = omega

        if copy:
            self.coeffs = _np.copy(coeffs)
            self.coeffs[~mask] = 0.
        else:
            self.coeffs = coeffs

        if errors is not None:
            if copy:
                self.errors = _np.copy(errors)
                self.errors[~mask] = 0.
            else:
                self.errors = errors
        else:
            self.errors = None

    def __repr__(self):
        if self.errors is not None:
            err_set = True
        else:
            err_set = False

        return ('kind = {:s}\n'
                'normalization = {:s}\n'
                'csphase = {:d}\n'
                'lmax = {:d}\n'
                'GM (m3 / s2) = {:s}\n'
                'r0 (m) = {:s}\n'
                'Omega (rad / s) = {:s}\n'
                'errors are set: {:s}\n'
                'header = {:s}'
                .format(repr(self.kind), repr(self.normalization),
                        self.csphase, self.lmax, repr(self.gm), repr(self.r0),
                        repr(self.omega), repr(err_set), repr(self.header)))

    def _rotate(self, angles, dj_matrix, gm=None, r0=None, omega=None):
        """Rotate the coefficients by the Euler angles alpha, beta, gamma."""
        if dj_matrix is None:
            dj_matrix = _shtools.djpi2(self.lmax + 1)

        # The coefficients need to be 4pi normalized with csphase = 1
        coeffs = _shtools.SHRotateRealCoef(
            self.to_array(normalization='4pi', csphase=1), angles, dj_matrix)

        # Convert 4pi normalized coefficients to the same normalization
        # as the unrotated coefficients.
        if self.normalization != '4pi' or self.csphase != 1:
            temp = _convert(coeffs, normalization_in='4pi', csphase=1,
                            normalization_out=self.normalization,
                            csphase_out=self.csphase)
            return SHGravCoeffs.from_array(
                temp, normalization=self.normalization,
                csphase=self.csphase, copy=False, gm=gm, r0=r0, omega=omega)
        else:
            return SHGravCoeffs.from_array(coeffs, gm=gm, r0=r0, omega=omega,
                                           copy=False)
