"""
    Spherical Harmonic Coefficients classes
"""
import numpy as _np
import matplotlib as _mpl
import matplotlib.pyplot as _plt
import copy as _copy
import gzip as _gzip
import shutil as _shutil
import warnings as _warnings
from scipy.special import factorial as _factorial
import xarray as _xr

from .. import shtools as _shtools
from ..spectralanalysis import spectrum as _spectrum
from ..spectralanalysis import cross_spectrum as _cross_spectrum
from ..shio import convert as _convert
from ..shio import shread as _shread
from ..shio import shwrite as _shwrite
from ..shio import read_dov as _read_dov
from ..shio import write_dov as _write_dov
from ..shio import read_bshc as _read_bshc
from ..shio import write_bshc as _write_bshc


class SHCoeffs(object):
    """
    Spherical Harmonic Coefficients class.

    The coefficients of this class can be initialized using one of the four
    constructor methods:

        x = SHCoeffs.from_array(array)
        x = SHCoeffs.from_random(powerspectrum)
        x = SHCoeffs.from_zeros(lmax)
        x = SHCoeffs.from_file('fname.dat')
        x = SHCoeffs.from_netcdf('ncname.nc')
        x = SHCoeffs.from_cap(theta, lmax)

    The normalization convention of the input coefficents is specified
    by the normalization and csphase parameters, which take the following
    values:

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
    error_kind    : An arbitrary string describing the kind of errors, such as
                    'unknown', 'unspecified', 'calibrated', 'formal' or None.
    normalization : The normalization of the coefficients: '4pi', 'ortho',
                    'schmidt', or 'unnorm'.
    csphase       : Defines whether the Condon-Shortley phase is used (1)
                    or not (-1).
    mask          : A boolean mask that is True for the permissible values of
                    degree l and order m.
    kind          : The coefficient data type: either 'complex' or 'real'.
    units         : The units of the spherical harmonic coefficients.
    header        : A list of values (of type str) from the header line of the
                    input file used to initialize the class (for 'shtools'
                    and 'dov' formatted files).
    header2       : A list of values (of type str) from the second header line
                    of the input file used to initialize the class (for
                    'shtools' and 'dov' formatted files only).

    Each class instance provides the following methods:

    degrees()             : Return an array listing the spherical harmonic
                            degrees from 0 to lmax.
    spectrum()            : Return the spectrum of the function as a function
                            of spherical harmonic degree.
    cross_spectrum()      : Return the cross-spectrum of two functions as a
                            function of spherical harmonic degree.
    admittance()          : Return the admittance with another function.
    correlation()         : Return the spectral correlation with another
                            function.
    admitcorr()           : Return the admittance and spectral correlation with
                            another function.
    volume()              : Calculate the volume of the body.
    centroid()            : Compute the centroid of the body.
    set_coeffs()          : Set coefficients in-place to specified values.
    rotate()              : Rotate the coordinate system used to express the
                            spherical harmonic coefficients and return a new
                            class instance.
    convert()             : Return a new class instance using a different
                            normalization convention.
    pad()                 : Return a new class instance that is zero padded or
                            truncated to a different lmax.
    expand()              : Evaluate the coefficients either on a spherical
                            grid and return an SHGrid class instance, or for
                            a list of latitude and longitude coordinates.
    plot_spectrum()       : Plot the spectrum as a function of spherical
                            harmonic degree.
    plot_cross_spectrum() : Plot the cross-spectrum of two functions.
    plot_admittance()     : Plot the admittance with another function.
    plot_correlation()    : Plot the spectral correlation with another
                            function.
    plot_admitcorr()      : Plot the admittance and spectral correlation with
                            another function.
    plot_spectrum2d()     : Plot the 2D spectrum of all spherical harmonic
                            degrees and orders.
    plot_cross_spectrum2d() : Plot the 2D cross-spectrum of all spherical
                              harmonic degrees and orders.
    to_array()            : Return an array of spherical harmonic coefficients
                            with a different normalization convention.
    to_file()             : Save raw spherical harmonic coefficients as a file.
    to_netcdf()           : Save raw spherical harmonic coefficients as a
                            netcdf file.
    copy()                : Return a copy of the class instance.
    info()                : Print a summary of the data stored in the SHCoeffs
                            instance.
    """

    def __init__(self):
        """Unused constructor of the super class."""
        print('Initialize the class using one of the class methods:\n'
              '>>> pyshtools.SHCoeffs.from_array\n'
              '>>> pyshtools.SHCoeffs.from_random\n'
              '>>> pyshtools.SHCoeffs.from_zeros\n'
              '>>> pyshtools.SHCoeffs.from_file\n'
              '>>> pyshtools.SHCoeffs.from_netcdf\n'
              '>>> pyshtools.SHCoeffs.from_cap\n'
              )

    # ---- Factory methods ----
    @classmethod
    def from_array(self, coeffs, errors=None, error_kind=None,
                   normalization='4pi', csphase=1, lmax=None, units=None,
                   copy=True):
        """
        Initialize the class with spherical harmonic coefficients from an input
        array.

        Usage
        -----
        x = SHCoeffs.from_array(array, [errors, error_kind, normalization,
                                        csphase, lmax, units, copy])

        Returns
        -------
        x : SHCoeffs class instance.

        Parameters
        ----------
        array : ndarray, shape (2, lmaxin+1, lmaxin+1).
            The input spherical harmonic coefficients.
        errors : ndarray, optional, default = None
            The uncertainties of the spherical harmonic coefficients.
        error_kind : str, optional, default = None
            An arbitrary string describing the kind of errors, such as None,
            'unspecified', 'calibrated' or 'formal'.
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
        units : str, optional, default = None
            The units of the spherical harmonic coefficients.
        copy : bool, optional, default = True
            If True, make a copy of array when initializing the class instance.
            If False, initialize the class instance with a reference to array.
        """
        if _np.iscomplexobj(coeffs):
            kind = 'complex'
        else:
            kind = 'real'

        if type(normalization) != str:
            raise ValueError('normalization must be a string. ' +
                             'Input type is {:s}.'
                             .format(str(type(normalization))))

        if normalization.lower() not in ('4pi', 'ortho', 'schmidt', 'unnorm'):
            raise ValueError(
                "The normalization must be '4pi', 'ortho', 'schmidt', " +
                "or 'unnorm'. Input value is {:s}."
                .format(repr(normalization))
                )

        if csphase != 1 and csphase != -1:
            raise ValueError(
                "csphase must be either 1 or -1. Input value is {:s}."
                .format(repr(csphase))
                )

        if errors is not None:
            if coeffs.shape != errors.shape:
                raise ValueError(
                    "The shape of coeffs and errors must be the same."
                    "Shape of coeffs = {:s}, shape of errors = {:s}."
                    .format(repr(coeffs.shape), repr(coeffs.errors))
                    )
            if error_kind is None:
                error_kind = 'unspecified'

        lmaxin = coeffs.shape[1] - 1
        if lmax is None:
            lmax = lmaxin
        else:
            if lmax > lmaxin:
                lmax = lmaxin

        if normalization.lower() == 'unnorm' and lmax > 85:
            _warnings.warn("Calculations using unnormalized coefficients " +
                           "are stable only for degrees less than or equal " +
                           "to 85. lmax for the coefficients will be set to " +
                           "85. Input value is {:d}.".format(lmax),
                           category=RuntimeWarning)
            lmax = 85

        for cls in self.__subclasses__():
            if cls.istype(kind):
                if errors is not None:
                    return cls(coeffs[:, 0:lmax+1, 0:lmax+1],
                               errors=errors[:, 0:lmax+1, 0:lmax+1],
                               error_kind=error_kind,
                               normalization=normalization.lower(),
                               csphase=csphase, units=units, copy=copy)
                else:
                    return cls(coeffs[:, 0:lmax+1, 0:lmax+1],
                               normalization=normalization.lower(),
                               csphase=csphase, units=units, copy=copy)

    @classmethod
    def from_zeros(self, lmax, errors=None, error_kind=None, kind='real',
                   normalization='4pi', csphase=1, units=None):
        """
        Initialize class with spherical harmonic coefficients set to zero from
        degree 0 to lmax.

        Usage
        -----
        x = SHCoeffs.from_zeros(lmax, [errors, error_kind, normalization,
                                       csphase, kind, units])

        Returns
        -------
        x : SHCoeffs class instance.

        Parameters
        ----------
        lmax : int
            The highest spherical harmonic degree l of the coefficients.
        errors : bool, optional, default = None
            If True, initialize the attribute errors with zeros.
        error_kind : str, optional, default = None
            An arbitrary string describing the kind of errors, such as None,
            'unspecified', 'calibrated' or 'formal'.
        normalization : str, optional, default = '4pi'
            '4pi', 'ortho', 'schmidt', or 'unnorm' for geodesy 4pi normalized,
            orthonormalized, Schmidt semi-normalized, or unnormalized
             coefficients, respectively.
        csphase : int, optional, default = 1
            Condon-Shortley phase convention: 1 to exclude the phase factor,
            or -1 to include it.
        kind : str, optional, default = 'real'
            'real' or 'complex' spherical harmonic coefficients.
        units : str, optional, default = None
            The units of the spherical harmonic coefficients.
        """
        error_coeffs = None
        if kind.lower() not in ('real', 'complex'):
            raise ValueError(
                "Kind must be 'real' or 'complex'. Input value is {:s}."
                .format(repr(kind))
                )

        if normalization.lower() not in ('4pi', 'ortho', 'schmidt', 'unnorm'):
            raise ValueError(
                "The normalization must be '4pi', 'ortho', 'schmidt', " +
                "or 'unnorm'. Input value is {:s}."
                .format(repr(normalization))
                )

        if csphase != 1 and csphase != -1:
            raise ValueError(
                "csphase must be either 1 or -1. Input value is {:s}."
                .format(repr(csphase))
                )

        if normalization.lower() == 'unnorm' and lmax > 85:
            _warnings.warn("Calculations using unnormalized coefficients " +
                           "are stable only for degrees less than or equal " +
                           "to 85. lmax for the coefficients will be set to " +
                           "85. Input value is {:d}.".format(lmax),
                           category=RuntimeWarning)
            lmax = 85

        if kind.lower() == 'real':
            coeffs = _np.zeros((2, lmax + 1, lmax + 1))
            if errors:
                error_coeffs = _np.zeros((2, lmax + 1, lmax + 1))
        else:
            coeffs = _np.zeros((2, lmax + 1, lmax + 1), dtype=complex)
            if errors:
                error_coeffs = _np.zeros((2, lmax + 1, lmax + 1),
                                         dtype=complex)
        if errors is True and error_kind is None:
            error_kind = 'unspecified'

        for cls in self.__subclasses__():
            if cls.istype(kind):
                return cls(coeffs, errors=error_coeffs, error_kind=error_kind,
                           normalization=normalization.lower(),
                           csphase=csphase, units=units)

    @classmethod
    def from_file(self, fname, lmax=None, format='shtools', kind='real',
                  errors=None, error_kind=None, normalization='4pi', skip=0,
                  header=False, header2=False, csphase=1, units=None,
                  **kwargs):
        """
        Initialize the class with spherical harmonic coefficients from a file.

        Usage
        -----
        x = SHCoeffs.from_file(filename, [format='shtools' or 'dov', lmax,
                               errors, error_kind, normalization, csphase,
                               skip, header, header2, units])
        x = SHCoeffs.from_file(filename, format='bshc', [lmax, normalization,
                               csphase, units])
        x = SHCoeffs.from_file(filename, format='npy', [lmax, normalization,
                               csphase, units, **kwargs])

        Returns
        -------
        x : SHCoeffs class instance.

        Parameters
        ----------
        filename : str
            File name or URL containing the spherical harmonic coefficients.
            filename will be treated as a URL if it starts with 'http://',
            'https://', or 'ftp://'. For 'shtools' and 'bshc' formatted files,
            if filename ends with '.gz' or '.zip', the file will be
            uncompressed before parsing.
        format : str, optional, default = 'shtools'
            'shtools' for generic text files, 'dov' for [degree, order, value]
            text files, 'bshc' for binary spherical harmonic coefficient
            files, or 'npy' for binary numpy files.
        lmax : int, optional, default = None
            The maximum spherical harmonic degree to read from the file. The
            default is to read the entire file.
        errors : bool, optional, default = None
            If True, read errors from the file (for 'shtools' and 'dov'
            formatted files only).
        error_kind : str, optional, default = None
            For 'shtools' and 'dov' formatted files: An arbitrary string
            describing the kind of errors, such as None, 'unspecified',
            'calibrated' or 'formal'.
        normalization : str, optional, default = '4pi'
            '4pi', 'ortho', 'schmidt', or 'unnorm' for geodesy 4pi normalized,
            orthonormalized, Schmidt semi-normalized, or unnormalized
            coefficients, respectively.
        csphase : int, optional, default = 1
            Condon-Shortley phase convention: 1 to exclude the phase factor,
            or -1 to include it.
        skip : int, optional, default = 0
            Number of lines to skip at the beginning of the file for 'shtools'
            formatted files.
        header : bool, optional, default = False
            If True, read a list of values from the header line of an 'shtools'
            or 'dov' formatted file.
        header2 : bool, optional, default = False
            If True, read a list of values from a second header line of an
            'shtools' or 'dov' formatted file.
        units : str, optional, default = None
            The units of the spherical harmonic coefficients.
        **kwargs : keyword argument list, optional for format = 'npy'
            Keyword arguments of numpy.load() when format is 'npy'.

        Notes
        -----
        Supported file formats:
            'shtools' (see pyshtools.shio.shread)
            'dov' (see pyshtools.shio.read_dov)
            'bshc' (see pyshtools.shio.read_bshc)
            'npy' (see numpy.load)

        For 'shtools', 'dov' or 'bshc' formatted files, if filename starts with
        'http://', 'https://', or 'ftp://', the file will be treated as a URL.
        In this case, the file will be downloaded in its entirety before it is
        parsed. If the filename ends with '.gz' or '.zip', the file will be
        automatically uncompressed before parsing. For zip files, archives with
        only a single file are supported. Note that reading '.gz' and '.zip'
        files will be extremely slow if lmax is not specified.

        For 'shtools' and 'dov' formatted files, the optional parameter `skip`
        specifies how many lines should be skipped before attempting to parse
        the file, the optional parameter `header` specifies whether to read a
        list of values from a header line, and the optional parameter `lmax`
        specifies the maximum degree to read from the file.
        """
        error_coeffs = None
        header_list = None
        header2_list = None

        if type(normalization) != str:
            raise ValueError('normalization must be a string. '
                             'Input type is {:s}.'
                             .format(str(type(normalization))))

        if normalization.lower() not in ('4pi', 'ortho', 'schmidt', 'unnorm'):
            raise ValueError(
                "The input normalization must be '4pi', 'ortho', 'schmidt', "
                "or 'unnorm'. Provided value is {:s}."
                .format(repr(normalization))
                )

        if csphase != 1 and csphase != -1:
            raise ValueError(
                "csphase must be 1 or -1. Input value is {:s}."
                .format(repr(csphase))
                )

        if format.lower() == 'shtools' or format.lower() == 'dov':
            if format.lower() == 'shtools':
                read_func = _shread
            else:
                read_func = _read_dov

            if header is True:
                if errors:
                    if header2:
                        coeffs, error_coeffs, lmaxout, header_list, \
                            header2_list = read_func(fname, lmax=lmax,
                                                     skip=skip, header=True,
                                                     header2=True, error=True)
                    else:
                        coeffs, error_coeffs, lmaxout, header_list = read_func(
                            fname, lmax=lmax, skip=skip, header=True,
                            error=True)
                else:
                    if header2:
                        coeffs, lmaxout, header_list, header2_list = read_func(
                            fname, lmax=lmax, skip=skip, header=True,
                            header2=True)
                    else:
                        coeffs, lmaxout, header_list = read_func(
                            fname, lmax=lmax, skip=skip, header=True)
            else:
                if errors:
                    coeffs, error_coeffs, lmaxout = read_func(fname, lmax=lmax,
                                                              skip=skip,
                                                              error=True)
                else:
                    coeffs, lmaxout = read_func(fname, lmax=lmax, skip=skip)

            if errors is True and error_kind is None:
                error_kind = 'unspecified'

        elif format.lower() == 'bshc':
            coeffs, lmaxout = _read_bshc(fname, lmax=lmax)

        elif format.lower() == 'npy':
            coeffs = _np.load(fname, **kwargs)
            lmaxout = coeffs.shape[1] - 1
            if lmax is not None:
                if lmax < lmaxout:
                    coeffs = coeffs[:, :lmax+1, :lmax+1]
                    lmaxout = lmax

        else:
            raise NotImplementedError(
                'format={:s} not implemented.'.format(repr(format)))

        if normalization.lower() == 'unnorm' and lmaxout > 85:
            _warnings.warn("Calculations using unnormalized coefficients "
                           "are stable only for degrees less than or equal "
                           "to 85. lmax for the coefficients will be set to "
                           "85. Input value is {:d}.".format(lmaxout),
                           category=RuntimeWarning)
            lmaxout = 85
            coeffs = coeffs[:, :lmaxout+1, :lmaxout+1]

        if _np.iscomplexobj(coeffs):
            kind = 'complex'
        else:
            kind = 'real'

        for cls in self.__subclasses__():
            if cls.istype(kind):
                return cls(coeffs, errors=error_coeffs, error_kind=error_kind,
                           normalization=normalization.lower(),
                           csphase=csphase, header=header_list,
                           header2=header2_list, units=units)

    @classmethod
    def from_random(self, power, lmax=None, kind='real', normalization='4pi',
                    csphase=1, units=None, exact_power=False, seed=None):
        """
        Initialize the class with spherical harmonic coefficients as random
        variables with a given spectrum.

        Usage
        -----
        x = SHCoeffs.from_random(power, [lmax, kind, normalization, csphase,
                                         units, exact_power, seed])

        Returns
        -------
        x : SHCoeffs class instance.

        Parameters
        ----------
        power : ndarray, shape (L+1)
            numpy array of shape (L+1) that specifies the expected power per
            degree l of the random coefficients, where L is the maximum
            spherical harmonic bandwidth.
        lmax : int, optional, default = len(power) - 1
            The maximum spherical harmonic degree l of the output coefficients.
            The coefficients will be set to zero for degrees greater than L.
        kind : str, optional, default = 'real'
            'real' or 'complex' spherical harmonic coefficients.
        normalization : str, optional, default = '4pi'
            '4pi', 'ortho', 'schmidt', or 'unnorm' for geodesy 4pi normalized,
            orthonormalized, Schmidt semi-normalized, or unnormalized
            coefficients, respectively.
        csphase : int, optional, default = 1
            Condon-Shortley phase convention: 1 to exclude the phase factor,
            or -1 to include it.
        units : str, optional, default = None
            The units of the spherical harmonic coefficients.
        exact_power : bool, optional, default = False
            The total variance of the coefficients is set exactly to the input
            power. The distribution of power at degree l amongst the angular
            orders is random, but the total power is fixed.
        seed : int, optional, default = None
            Set the seed for the numpy random number generator.

        Notes
        -----
        This routine returns a random realization of spherical harmonic
        coefficients obtained from a normal distribution. The variance of
        each coefficient at degree l is equal to the total power at degree
        l divided by the number of coefficients at that degree. The power
        spectrum of the random realization can be fixed exactly to the input
        spectrum by setting exact_power to True.
        """
        # check if all arguments are correct
        if type(normalization) != str:
            raise ValueError('normalization must be a string. ' +
                             'Input type is {:s}.'
                             .format(str(type(normalization))))

        if normalization.lower() not in ('4pi', 'ortho', 'schmidt', 'unnorm'):
            raise ValueError(
                "The input normalization must be '4pi', 'ortho', 'schmidt', " +
                "or 'unnorm'. Provided value is {:s}."
                .format(repr(normalization))
                )

        if csphase != 1 and csphase != -1:
            raise ValueError(
                "csphase must be 1 or -1. Input value is {:s}."
                .format(repr(csphase))
                )

        if kind.lower() not in ('real', 'complex'):
            raise ValueError(
                "kind must be 'real' or 'complex'. " +
                "Input value is {:s}.".format(repr(kind)))

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
            _warnings.warn("Calculations using unnormalized coefficients " +
                           "are stable only for degrees less than or equal " +
                           "to 85. lmax for the coefficients will be set to " +
                           "85. Input value is {:d}.".format(nl-1),
                           category=RuntimeWarning)
            nl = 85 + 1
            lmax = 85

        # Create coefficients with unit variance, which returns an expected
        # total power per degree of (2l+1) for 4pi normalized harmonics.
        if seed is not None:
            _np.random.seed(seed=seed)
        if kind.lower() == 'real':
            coeffs = _np.zeros((2, nl, nl))
            for l in degrees:
                coeffs[:2, l, :l+1] = _np.random.normal(size=(2, l+1))
        elif kind.lower() == 'complex':
            # - need to divide by sqrt 2 as there are two terms for each coeff.
            coeffs = _np.zeros((2, nl, nl), dtype=complex)
            for l in degrees:
                coeffs[:2, l, :l+1] = (_np.random.normal(size=(2, l+1)) +
                                       1j * _np.random.normal(size=(2, l+1))
                                       ) / _np.sqrt(2.)

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

        if lmax > nl - 1:
            coeffs = _np.pad(coeffs, ((0, 0), (0, lmax - nl + 1),
                             (0, lmax - nl + 1)), 'constant')

        for cls in self.__subclasses__():
            if cls.istype(kind):
                return cls(coeffs, errors=None,
                           normalization=normalization.lower(),
                           csphase=csphase, units=units)

    @classmethod
    def from_netcdf(self, filename, lmax=None, normalization='4pi', csphase=1,
                    units=None):
        """
        Initialize the class with spherical harmonic coefficients from a
        netcdf file.

        Usage
        -----
        x = SHCoeffs.from_netcdf(filename, [lmax, normalization, csphase,
                                            units])

        Returns
        -------
        x : SHCoeffs class instance.

        Parameters
        ----------
        filename : str
            Name of the file, including path.
        lmax : int, optional, default = None
            The maximum spherical harmonic degree to read.
        normalization : str, optional, default = '4pi'
            Spherical harmonic normalization if not specified in the netcdf
            file: '4pi', 'ortho', 'schmidt', or 'unnorm' for geodesy 4pi
            normalized, orthonormalized, Schmidt semi-normalized, or
            unnormalized coefficients, respectively.
        csphase : int, optional, default = 1
            Condon-Shortley phase convention if not specified in the netcdf
            file: 1 to exclude the phase factor, or -1 to include it.
        units : str, optional, default = None
            The units of the spherical harmonic coefficients.

        Description
        -----------
        The format of the netcdf file has to be exactly as the format that is
        used in SHCoeffs.to_netcdf().
        """
        ds = _xr.open_dataset(filename)

        try:
            normalization = ds.coeffs.normalization
        except:
            pass

        if type(normalization) != str:
            raise ValueError('normalization must be a string. '
                             'Input type was {:s}'
                             .format(str(type(normalization))))
        if normalization.lower() not in ('4pi', 'ortho', 'schmidt', 'unnorm'):
            raise ValueError(
                "The input normalization must be '4pi', 'ortho', "
                "'schmidt', or 'unnorm'. Provided value was {:s}"
                .format(repr(normalization))
                )

        try:
            csphase = ds.coeffs.csphase
        except:
            pass

        if csphase != 1 and csphase != -1:
            raise ValueError(
                "csphase must be 1 or -1. Input value was {:s}"
                .format(repr(csphase))
                )

        try:
            units = ds.coeffs.units
        except:
            pass

        lmaxout = ds.dims['degree'] - 1
        c = _np.tril(ds.coeffs.data)
        s = _np.triu(ds.coeffs.data, k=1)
        s = _np.vstack([s[-1], s[:-1]])
        s = _np.transpose(s)
        if isinstance(lmax, int):
            c, s = c[:lmax+1, :lmax+1], s[:lmax+1, :lmax+1]
            lmaxout = lmax

        if normalization.lower() == 'unnorm' and lmaxout > 85:
            _warnings.warn("Calculations using unnormalized coefficients " +
                           "are stable only for degrees less than or equal " +
                           "to 85. lmax for the coefficients will be set to " +
                           "85. Input value was {:d}.".format(lmaxout),
                           category=RuntimeWarning)
            lmaxout = 85
            c, s = c[:lmaxout+1, :lmaxout+1], s[:lmaxout+1, :lmaxout+1]
        coeffs = _np.array([c, s])

        try:
            cerrors = _np.tril(ds.errors.data)
            serrors = _np.triu(ds.errors.data, k=1)
            serrors = _np.vstack([serrors[-1], serrors[:-1]])
            serrors = _np.transpose(serrors)
            cerrors = cerrors[:lmaxout+1, :lmaxout+1]
            serrors = serrors[:lmaxout+1, :lmaxout+1]
            errors = _np.array([cerrors, serrors])
            error_kind = ds.errors.error_kind
        except:
            errors = None
            error_kind = None

        if _np.iscomplexobj(coeffs):
            kind = 'complex'
        else:
            kind = 'real'

        for cls in self.__subclasses__():
            if cls.istype(kind):
                return cls(coeffs, errors=errors, error_kind=error_kind,
                           normalization=normalization.lower(),
                           csphase=csphase, units=units)

    @classmethod
    def from_cap(self, theta, lmax, clat=None, clon=None, normalization='4pi',
                 csphase=1, kind='real', units=None, degrees=True, copy=True):
        """
        Initialize the class with spherical harmonic coefficients of a
        spherical cap centered at the north pole.

        Usage
        -----
        x = SHCoeffs.from_cap(theta, lmax, [clat, clon, normalization, csphase,
                                            kind, units, degrees, copy])

        Returns
        -------
        x : SHCoeffs class instance.

        Parameters
        ----------
        theta : float
            The angular radius of the spherical cap, default in degrees.
        lmax : int
            The maximum spherical harmonic degree of the coefficients.
        clat, clon : float, optional, default = None
            Latitude and longitude of the center of the rotated spherical cap
            (default in degrees).
        normalization : str, optional, default = '4pi'
            '4pi', 'ortho', 'schmidt', or 'unnorm' for geodesy 4pi normalized,
            orthonormalized, Schmidt semi-normalized, or unnormalized
            coefficients, respectively.
        csphase : int, optional, default = 1
            Condon-Shortley phase convention: 1 to exclude the phase factor,
            or -1 to include it.
        kind : str, optional, default = 'real'
            'real' or 'complex' spherical harmonic coefficients.
        units : str, optional, default = None
            The units of the spherical harmonic coefficients.
        degrees : bool, optional = True
            If True, theta, clat, and clon are in degrees.
        copy : bool, optional, default = True
            If True, make a copy of array when initializing the class instance.
            If False, initialize the class instance with a reference to array.

        Notes
        -----
        The spherical harmonic coefficients are normalized such that the
        average value of the function is equal to 1. To rotate the cap to a
        specified latitude and longitude, specify the optional parameters clat
        and clon.
        """
        if type(normalization) != str:
            raise ValueError('normalization must be a string. ' +
                             'Input type is {:s}.'
                             .format(str(type(normalization))))

        if normalization.lower() not in ('4pi', 'ortho', 'schmidt', 'unnorm'):
            raise ValueError(
                "The normalization must be '4pi', 'ortho', 'schmidt', " +
                "or 'unnorm'. Input value is {:s}."
                .format(repr(normalization))
                )

        if csphase != 1 and csphase != -1:
            raise ValueError(
                "csphase must be either 1 or -1. Input value is {:s}."
                .format(repr(csphase))
                )

        if kind.lower() not in ('real', 'complex'):
            raise ValueError(
                "kind must be 'real' or 'complex'. " +
                "Input value is {:s}.".format(repr(kind)))

        if (clat is None and clon is not None) or \
                (clat is not None and clon is None):
            raise ValueError('clat and clon must both be input. ' +
                             'clat = {:s}, clon = {:s}.'
                             .format(repr(clat), repr(clon)))

        if degrees is True:
            theta = _np.deg2rad(theta)

        cl = _shtools.SphericalCapCoef(theta, lmax)
        coeffs = _np.zeros((2, lmax+1, lmax+1))
        coeffs[0, 0:lmax+1, 0] = cl[0:lmax+1]
        coeffs = _convert(coeffs, normalization_in='4pi',
                          normalization_out=normalization,
                          csphase_in=1, csphase_out=csphase
                          )

        if kind == 'complex':
            coeffs = _shtools.SHrtoc(coeffs)

        for cls in self.__subclasses__():
            if cls.istype(kind):
                temp = cls(coeffs[:, 0:lmax+1, 0:lmax+1],
                           normalization=normalization.lower(),
                           csphase=csphase, units=units, copy=copy)

        if clat is not None and clon is not None:
            if degrees is True:
                temp = temp.rotate(0., -90 + clat, -clon, degrees=True)
            else:
                temp = temp.rotate(0., -_np.pi/2. + clat, -clon,
                                   degrees=False)

        return temp

    # ---- Define methods that modify internal variables ----
    def set_coeffs(self, values, ls, ms):
        """
        Set spherical harmonic coefficients in-place to specified values.

        Usage
        -----
        x.set_coeffs(values, ls, ms)

        Parameters
        ----------
        values : float or complex (list)
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

    # ---- IO Routines
    def to_file(self, filename, format='shtools', header=None, header2=None,
                errors=True, lmax=None, **kwargs):
        """
        Save raw spherical harmonic coefficients to a file.

        Usage
        -----
        x.to_file(filename, [format='shtools', header, header2, errors, lmax])
        x.to_file(filename, format='dov', [header, header2, errors, lmax])
        x.to_file(filename, format='bshc', [lmax])
        x.to_file(filename, format='npy', [**kwargs])

        Parameters
        ----------
        filename : str
            Name of the output file. If the filename ends with '.gz', the file
            will be compressed using gzip.
        format : str, optional, default = 'shtools'
            'shtools', 'dov', 'bshc' or 'npy'.
        header : str, optional, default = None
            A header string written to an 'shtools' or 'dov' formatted file
            directly before the spherical harmonic coefficients.
        header2 : str, optional, default = None
            A second header string written to an 'shtools' or 'dov' formatted
            file directly before the spherical harmonic coefficients.
        errors : bool, optional, default = False
            If True, save the errors in the file (for 'shtools' formatted
            files only).
        lmax : int, optional, default = self.lmax
            The maximum spherical harmonic degree to write to the file.
        **kwargs : keyword argument list, optional for format = 'npy'
            Keyword arguments of numpy.save().

        Notes
        -----
        Supported file formats:
            'shtools' (see pyshtools.shio.shwrite)
            'dov' (see pyshtools.shio.write_dov)
            'bshc' (see pyshtools.shio.write_bshc)
            'npy' (see numpy.save)

        If the filename end with '.gz', the file will be compressed using gzip.

        'shtools': The coefficients will be written to an ascii formatted file.
        The first line of the file is an optional user provided header line,
        and the spherical harmonic coefficients (and optionally the errors)
        are then listed, with increasing degree and order, with the format

        l, m, coeffs[0, l, m], coeffs[1, l, m], error[0, l, m], error[1, l, m]

        where l and m are the spherical harmonic degree and order,
        respectively.

        'dov': This format is nearly the same as 'shtools', with the exception
        that each line contains a single coefficient (and optionally an error)
        for each degree and order:

        l, m, coeffs[0, l, m], error[0, l, m]
        l, -m, coeffs[1, l, m], error[1, l, m]

        'bshc': The coefficients will be written to a binary file composed
        solely of 8-byte floats. The file starts with the minimum and maximum
        degree, and is followed by the cosine coefficients and then sine
        coefficients (with all orders being listed, one degree at a time).

        'npy': The spherical harmonic coefficients will be saved to a binary
        numpy 'npy' file.
        """
        if errors is True and self.errors is None:
            errors = False

        if filename[-3:] == '.gz':
            filebase = filename[:-3]
        else:
            filebase = filename

        if format.lower() == 'shtools':
            if errors:
                _shwrite(filebase, self.coeffs, errors=self.errors,
                         header=header, header2=header2, lmax=lmax)
            else:
                _shwrite(filebase, self.coeffs, errors=None,
                         header=header, header2=header2, lmax=lmax)

        elif format.lower() == 'dov':
            if errors:
                _write_dov(filebase, self.coeffs, errors=self.errors,
                           header=header, header2=header2, lmax=lmax)
            else:
                _write_dov(filebase, self.coeffs, errors=None,
                           header=header, header2=header2, lmax=lmax)

        elif format.lower() == 'bshc':
            _write_bshc(filebase, self.coeffs, lmax=lmax)

        elif format.lower() == 'npy':
            _np.save(filebase, self.coeffs, **kwargs)

        else:
            raise NotImplementedError(
                'format={:s} not implemented.'.format(repr(format)))

        if filename[-3:] == '.gz':
            with open(filebase, 'rb') as f_in:
                with _gzip.open(filename, 'wb') as f_out:
                    _shutil.copyfileobj(f_in, f_out)

    def to_netcdf(self, filename, title='', description='', lmax=None):
        """
        Return the coefficient data as a netcdf formatted file or object.

        Usage
        -----
        x.to_netcdf(filename, [title, description, lmax])

        Parameters
        ----------
        filename : str
            Name of the output file.
        title : str, optional, default = ''
            Title of the dataset
        description : str, optional, default = ''
            Description of the data.
        lmax : int, optional, default = self.lmax
            The maximum spherical harmonic degree to output.
        """
        if lmax is None:
            lmax = self.lmax

        ds = _xr.Dataset()
        ds.coords['degree'] = ('degree', _np.arange(lmax+1))
        ds.coords['order'] = ('order', _np.arange(lmax+1))
        # c coeffs as lower triangular matrix
        c = self.coeffs[0, :lmax+1, :lmax+1]
        # s coeffs as upper triangular matrix
        s = _np.transpose(self.coeffs[1, :lmax+1, :lmax+1])
        s = _np.vstack([s[1:], s[0]])
        ds['coeffs'] = (('degree', 'order'), c + s)
        ds['coeffs'].attrs['title'] = title
        ds['coeffs'].attrs['description'] = description
        ds['coeffs'].attrs['normalization'] = self.normalization
        ds['coeffs'].attrs['csphase'] = self.csphase
        if self.units is not None:
            ds['coeffs'].attrs['units'] = self.units

        if self.errors is not None:
            cerrors = self.errors[0, :lmax+1, :lmax+1]
            serrors = _np.transpose(self.errors[1, :lmax+1, :lmax+1])
            serrors = _np.vstack([serrors[1:], serrors[0]])
            ds['errors'] = (('degree', 'order'), cerrors + serrors)
            ds['errors'].attrs['normalization'] = self.normalization
            ds['errors'].attrs['csphase'] = self.csphase
            if self.units is not None:
                ds['errors'].attrs['units'] = self.units
            if self.error_kind is not None:
                ds['errors'].attrs['error_kind'] = self.error_kind

        ds.to_netcdf(filename)

    def to_array(self, normalization=None, csphase=None, lmax=None,
                 errors=True):
        """
        Return spherical harmonic coefficients (and errors) as a numpy array.

        Usage
        -----
        coeffs, [errors] = x.to_array([normalization, csphase, lmax, errors])

        Returns
        -------
        coeffs : ndarry, shape (2, lmax+1, lmax+1)
            numpy ndarray of the spherical harmonic coefficients.
        errors : ndarry, shape (2, lmax+1, lmax+1)
            numpy ndarray of the errors of the spherical harmonic
            coefficients if they are not None and errors is True.

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
        errors : bool, optional, default = True
            If True, return separate arrays of the coefficients and errors. If
            False, return only the coefficients.

        Notes
        -----
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

        if self.errors is not None and errors:
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
        Print a summary of the data stored in the SHCoeffs instance.

        Usage
        -----
        x.info()
        """
        print(repr(self))

    # -------------------------------------------------------------------------
    #    Mathematical operators
    #
    #    All operations ignore the errors of the coefficients.
    #    All operations ignore the units of the coefficients, with the
    #    exception of multiplying and dividing by a scalar.
    # -------------------------------------------------------------------------
    def __add__(self, other):
        """
        Add two similar sets of coefficients or coefficients and a scalar:
        self + other. For the addition of a scalar, only the degree 0
        term is modified.
        """
        if isinstance(other, SHCoeffs):
            if (self.normalization == other.normalization and self.csphase ==
                    other.csphase and self.kind == other.kind and
                    self.lmax == other.lmax):
                coeffs = _np.empty([2, self.lmax+1, self.lmax+1],
                                   dtype=self.coeffs.dtype)
                coeffs[self.mask] = (self.coeffs[self.mask] +
                                     other.coeffs[self.mask])
                return SHCoeffs.from_array(coeffs, csphase=self.csphase,
                                           normalization=self.normalization)
            else:
                raise ValueError('The two sets of coefficients must have the '
                                 'same kind, normalization, csphase, and '
                                 'lmax.')
        elif _np.isscalar(other) is True:
            if self.kind == 'real' and _np.iscomplexobj(other):
                raise ValueError('Can not add a complex constant to real '
                                 'coefficients.')
            coeffs = self.coeffs.copy()
            coeffs[0, 0, 0] += other
            return SHCoeffs.from_array(coeffs, csphase=self.csphase,
                                       normalization=self.normalization)
        else:
            raise NotImplementedError('Mathematical operator not implemented '
                                      'for these operands.')

    def __radd__(self, other):
        """
        Add two similar sets of coefficients or coefficients and a scalar:
        other + self. For the addition of a scalar, only the degree 0
        term is modified.
        """
        return self.__add__(other)

    def __sub__(self, other):
        """
        Subtract two similar sets of coefficients or coefficients and a scalar:
        self - other. For the subtraction of a scalar, only the degree 0
        term is modified.
        """
        if isinstance(other, SHCoeffs):
            if (self.normalization == other.normalization and self.csphase ==
                    other.csphase and self.kind == other.kind and
                    self.lmax == other.lmax):
                coeffs = _np.empty([2, self.lmax+1, self.lmax+1],
                                   dtype=self.coeffs.dtype)
                coeffs[self.mask] = (self.coeffs[self.mask] -
                                     other.coeffs[self.mask])
                return SHCoeffs.from_array(coeffs, csphase=self.csphase,
                                           normalization=self.normalization)
            else:
                raise ValueError('The two sets of coefficients must have the '
                                 'same kind, normalization, csphase, and '
                                 'lmax.')
        elif _np.isscalar(other) is True:
            if self.kind == 'real' and _np.iscomplexobj(other):
                raise ValueError('Can not subtract a complex constant from '
                                 'real coefficients.')
            coeffs = self.coeffs.copy()
            coeffs[0, 0, 0] -= other
            return SHCoeffs.from_array(coeffs, csphase=self.csphase,
                                       normalization=self.normalization)
        else:
            raise NotImplementedError('Mathematical operator not implemented '
                                      'for these operands.')

    def __rsub__(self, other):
        """
        Subtract two similar sets of coefficients or coefficients and a scalar:
        other - self. For the subtraction from a scalar, self is multiplied by
        -1 and then other is added to the degree 0 coefficient.
        """
        if isinstance(other, SHCoeffs):
            if (self.normalization == other.normalization and self.csphase ==
                    other.csphase and self.kind == other.kind and
                    self.lmax == other.lmax):
                coeffs = _np.empty([2, self.lmax+1, self.lmax+1],
                                   dtype=self.coeffs.dtype)
                coeffs[self.mask] = (other.coeffs[self.mask] -
                                     self.coeffs[self.mask])
                return SHCoeffs.from_array(coeffs, csphase=self.csphase,
                                           normalization=self.normalization)
            else:
                raise ValueError('The two sets of coefficients must have the '
                                 'same kind, normalization, csphase, and '
                                 'lmax.')
        elif _np.isscalar(other) is True:
            if self.kind == 'real' and _np.iscomplexobj(other):
                raise ValueError('Can not subtract a complex constant from '
                                 'real coefficients.')
            coeffs = - self.coeffs.copy()
            coeffs[0, 0, 0] += other
            return SHCoeffs.from_array(coeffs, csphase=self.csphase,
                                       normalization=self.normalization)
        else:
            raise NotImplementedError('Mathematical operator not implemented '
                                      'for these operands.')

    def __mul__(self, other):
        """
        Multiply two similar sets of coefficients or coefficients and a scalar:
        self * other.
        """
        if isinstance(other, SHCoeffs):
            if (self.normalization == other.normalization and self.csphase ==
                    other.csphase and self.kind == other.kind and
                    self.lmax == other.lmax):
                coeffs = _np.empty([2, self.lmax+1, self.lmax+1],
                                   dtype=self.coeffs.dtype)
                coeffs[self.mask] = (self.coeffs[self.mask] *
                                     other.coeffs[self.mask])
                return SHCoeffs.from_array(coeffs, csphase=self.csphase,
                                           normalization=self.normalization)
            else:
                raise ValueError('The two sets of coefficients must have the '
                                 'same kind, normalization, csphase and '
                                 'lmax.')
        elif _np.isscalar(other) is True:
            coeffs = _np.empty([2, self.lmax+1, self.lmax+1],
                               dtype=self.coeffs.dtype)
            if self.kind == 'real' and _np.iscomplexobj(other):
                raise ValueError('Can not multiply real coefficients by '
                                 'a complex constant.')
            coeffs[self.mask] = self.coeffs[self.mask] * other
            return SHCoeffs.from_array(coeffs, csphase=self.csphase,
                                       normalization=self.normalization,
                                       units=self.units)
        else:
            raise NotImplementedError('Mathematical operator not implemented '
                                      'for these operands.')

    def __rmul__(self, other):
        """
        Multiply two similar sets of coefficients or coefficients and a scalar:
        other * self.
        """
        return self.__mul__(other)

    def __truediv__(self, other):
        """
        Divide two similar sets of coefficients or coefficients and a scalar:
        self / other.
        """
        if isinstance(other, SHCoeffs):
            if (self.normalization == other.normalization and self.csphase ==
                    other.csphase and self.kind == other.kind and
                    self.lmax == other.lmax):
                coeffs = _np.empty([2, self.lmax+1, self.lmax+1],
                                   dtype=self.coeffs.dtype)
                coeffs[self.mask] = (self.coeffs[self.mask] /
                                     other.coeffs[self.mask])
                return SHCoeffs.from_array(coeffs, csphase=self.csphase,
                                           normalization=self.normalization)
            else:
                raise ValueError('The two sets of coefficients must have the '
                                 'same kind, normalization, csphase and '
                                 'lmax.')
        elif _np.isscalar(other) is True:
            coeffs = _np.empty([2, self.lmax+1, self.lmax+1],
                               dtype=self.coeffs.dtype)
            if self.kind == 'real' and _np.iscomplexobj(other):
                raise ValueError('Can not multiply real coefficients by '
                                 'a complex constant.')
            coeffs[self.mask] = self.coeffs[self.mask] / other
            return SHCoeffs.from_array(coeffs, csphase=self.csphase,
                                       normalization=self.normalization,
                                       units=self.units)
        else:
            raise NotImplementedError('Mathematical operator not implemented '
                                      'for these operands.')

    def __pow__(self, other):
        """
        Raise the spherical harmonic coefficients to a scalar power:
        pow(self, other).
        """
        if _np.isscalar(other) is True:
            return SHCoeffs.from_array(pow(self.coeffs, other),
                                       csphase=self.csphase,
                                       normalization=self.normalization)
        else:
            raise NotImplementedError('Mathematical operator not implemented '
                                      'for these operands.')

    def __repr__(self):
        return ('kind = {:s}\n'
                'normalization = {:s}\n'
                'csphase = {:d}\n'
                'lmax = {:d}\n'
                'error_kind = {:s}\n'
                'header = {:s}\n'
                'header2 = {:s}\n'
                'units = {:s}'
                .format(
                    repr(self.kind), repr(self.normalization), self.csphase,
                    self.lmax, repr(self.error_kind), repr(self.header),
                    repr(self.header2), repr(self.units)))

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

    def spectrum(self, lmax=None, convention='power', unit='per_l', base=10.):
        """
        Return the spectrum as a function of spherical harmonic degree.

        Usage
        -----
        spectrum, [error_spectrum] = x.spectrum([lmax, convention, unit, base])

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
        lmax : int, optional, default = x.lmax
            Maximum spherical harmonic degree of the spectrum to output.
        convention : str, optional, default = 'power'
            The type of spectrum to return: 'power' for power spectrum,
            'energy' for energy spectrum, and 'l2norm' for the l2 norm
            spectrum.
        unit : str, optional, default = 'per_l'
            If 'per_l', return the total contribution to the spectrum for each
            spherical harmonic degree l. If 'per_lm', return the average
            contribution to the spectrum for each coefficient at spherical
            harmonic degree l. If 'per_dlogl', return the spectrum per log
            interval dlog_a(l).
        base : float, optional, default = 10.
            The logarithm base when calculating the 'per_dlogl' spectrum.

        Notes
        -----
        This method returns either the power spectrum, energy spectrum, or
        l2-norm spectrum. Total power is defined as the integral of the
        function squared over all space, divided by the area the function
        spans. If the mean of the function is zero, this is equivalent to the
        variance of the function. The total energy is the integral of the
        function squared over all space and is 4pi times the total power. For
        normalized coefficients ('4pi', 'ortho', or 'schmidt'), the l2-norm is
        the sum of the magnitude of the coefficients squared.

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
        s = _spectrum(self.coeffs, normalization=self.normalization,
                      convention=convention, unit=unit, base=base,
                      lmax=lmax)

        if self.errors is not None:
            es = _spectrum(self.errors, normalization=self.normalization,
                           convention=convention, unit=unit, base=base,
                           lmax=lmax)
            return s, es
        else:
            return s

    def cross_spectrum(self, clm, lmax=None, convention='power', unit='per_l',
                       base=10.):
        """
        Return the cross-spectrum of two functions.

        Usage
        -----
        cross_spectrum = x.cross_spectrum(clm, [lmax, convention, unit, base])

        Returns
        -------
        cross_spectrum : ndarray, shape (lmax+1)
            1-D numpy ndarray of the cross-spectrum, where lmax is the maximum
            spherical harmonic degree.

        Parameters
        ----------
        clm : SHCoeffs class instance.
            The second function used in computing the cross-spectrum.
        lmax : int, optional, default = x.lmax
            Maximum spherical harmonic degree of the spectrum to output.
        convention : str, optional, default = 'power'
            The type of spectrum to return: 'power' for power spectrum,
            'energy' for energy spectrum, and 'l2norm' for the l2 norm
            spectrum.
        unit : str, optional, default = 'per_l'
            If 'per_l', return the total contribution to the spectrum for each
            spherical harmonic degree l. If 'per_lm', return the average
            contribution to the spectrum for each coefficient at spherical
            harmonic degree l. If 'per_dlogl', return the spectrum per log
            interval dlog_a(l).
        base : float, optional, default = 10.
            The logarithm base when calculating the 'per_dlogl' spectrum.

        Notes
        -----
        This method returns either the cross-power spectrum, cross-energy
        spectrum, or l2-cross-norm spectrum. Total cross-power is defined as
        the integral of the function times the conjugate of clm over all space,
        divided by the area the functions span. If the means of the functions
        are zero, this is equivalent to the covariance of the two functions.
        The total cross-energy is the integral of this function times the
        conjugate of clm over all space and is 4pi times the total power. The
        l2-cross-norm is the sum of this function times the conjugate of clm
        over all angular orders as a function of spherical harmonic degree.

        The output spectrum can be expresed using one of three units. 'per_l'
        returns the contribution to the total spectrum from all angular orders
        at degree l. 'per_lm' returns the average contribution to the total
        spectrum from a single coefficient at degree l, and is equal to the
        'per_l' spectrum divided by (2l+1). 'per_dlogl' returns the
        contribution to the total spectrum from all angular orders over an
        infinitessimal logarithmic degree band. The contrubution in the band
        dlog_a(l) is spectrum(l, 'per_dlogl')*dlog_a(l), where a is the base,
        and where spectrum(l, 'per_dlogl) is equal to
        spectrum(l, 'per_l')*l*log(a).
        """
        if not isinstance(clm, SHCoeffs):
            raise ValueError('clm must be an SHCoeffs class instance. Input '
                             'type is {:s}.'.format(repr(type(clm))))

        if lmax is None:
            lmax = min(self.lmax, clm.lmax)

        return _cross_spectrum(self.coeffs,
                               clm.to_array(normalization=self.normalization,
                                            csphase=self.csphase,
                                            lmax=lmax),
                               normalization=self.normalization,
                               convention=convention, unit=unit, base=base,
                               lmax=lmax)

    def admittance(self, hlm, errors=True, lmax=None):
        """
        Return the admittance of two functions, Sgh(l) / Shh(l).

        Usage
        -----
        admittance = g.admittance(hlm, [errors, lmax])

        Returns
        -------
        admittance : ndarray, shape (lmax+1) or (2, lmax+1)
            1-D array of the admittance (errors=False) or 2-D array of the
            admittance and its uncertainty (errors=True), where lmax is the
            maximum spherical harmonic degree.

        Parameters
        ----------
        hlm : SHCoeffs class instance.
            The function h used in computing the admittance Sgh / Shh.
        errors : bool, optional, default = True
            Return the uncertainty of the admittance.
        lmax : int, optional, default = g.lmax
            Maximum spherical harmonic degree of the spectrum to output.

        Notes
        -----
        If two functions g and h are related by the equation

            glm = Z(l) hlm + nlm

        where nlm is a zero-mean random variable, the admittance Z(l) can be
        estimated using

            Z(l) = Sgh(l) / Shh(l),

        where Sgh, Shh and Sgg are the cross-power and power spectra of the
        functions g (self) and h (input).
        """
        if not isinstance(hlm, SHCoeffs):
            raise ValueError('hlm must be an SHCoeffs class instance. Input '
                             'type is {:s}.'.format(repr(type(hlm))))

        if lmax is None:
            lmax = min(self.lmax, hlm.lmax)

        sgh = _cross_spectrum(self.coeffs,
                              hlm.to_array(normalization=self.normalization,
                                           csphase=self.csphase, lmax=lmax,
                                           errors=False),
                              normalization=self.normalization,
                              lmax=lmax)
        shh = _spectrum(hlm.coeffs, normalization=hlm.normalization, lmax=lmax)

        with _np.errstate(invalid='ignore', divide='ignore'):
            admit = sgh / shh
            if errors:
                sgg = _spectrum(self.coeffs, normalization=self.normalization,
                                lmax=lmax)
                sigma = (sgg / shh) * (1. - sgh**2 / sgg / shh) / \
                    _np.arange(lmax+1) / 2.
                admit = _np.column_stack((admit, _np.sqrt(sigma)))
        return admit

    def correlation(self, hlm, lmax=None):
        """
        Return the spectral correlation with another function.

        Usage
        -----
        correlation = g.correlation(hlm, [lmax])

        Returns
        -------
        correlation : ndarray, shape (lmax+1)
            1-D numpy ndarray of the spectral correlation, where lmax is the
            maximum spherical harmonic degree.

        Parameters
        ----------
        hlm : SHCoeffs class instance.
            The function h used in computing the spectral correlation.
        lmax : int, optional, default = g.lmax
            Maximum spherical harmonic degree of the spectrum to output.

        Notes
        -----
        The spectral correlation is defined as

            gamma(l) = Sgh(l) / sqrt( Sgg(l) Shh(l) )

        where Sgh, Shh and Sgg are the cross-power and power spectra of the
        functions g (self) and h (input).
        """
        from .shgravcoeffs import SHGravCoeffs as _SHGravCoeffs
        from .shmagcoeffs import SHMagCoeffs as _SHMagCoeffs
        if not isinstance(hlm, (SHCoeffs, _SHMagCoeffs, _SHGravCoeffs)):
            raise ValueError('hlm must be an SHCoeffs, SHMagCoeffs or '
                             'SHGravCoeffs class instance. Input type is {:s}.'
                             .format(repr(type(hlm))))

        if lmax is None:
            lmax = min(self.lmax, hlm.lmax)

        sgg = _spectrum(self.coeffs, normalization=self.normalization,
                        lmax=lmax)
        shh = _spectrum(hlm.coeffs, normalization=hlm.normalization, lmax=lmax)
        sgh = _cross_spectrum(self.coeffs,
                              hlm.to_array(normalization=self.normalization,
                                           csphase=self.csphase, lmax=lmax,
                                           errors=False),
                              normalization=self.normalization,
                              lmax=lmax)

        with _np.errstate(invalid='ignore', divide='ignore'):
            return sgh / _np.sqrt(sgg * shh)

    def admitcorr(self, hlm, errors=True, lmax=None):
        """
        Return the admittance and spectral correlation of two functions.

        Usage
        -----
        admittance, correlation = g.admitcorr(hlm, [errors, lmax])

        Returns
        -------
        admittance : ndarray, shape (lmax+1) or (2, lmax+1)
            1-D array of the admittance (errors=False) or 2-D array of the
            admittance and its uncertainty (errors=True), where lmax is the
            maximum spherical harmonic degree.
        correlation : ndarray, shape (lmax+1)
            1-D numpy ndarray of the spectral correlation, where lmax is the
            maximum spherical harmonic degree.

        Parameters
        ----------
        hlm : SHCoeffs class instance.
            The function h used in computing the admittance and spectral
            correlation.
        errors : bool, optional, default = True
            Return the uncertainty of the admittance.
        lmax : int, optional, default = g.lmax
            Maximum spherical harmonic degree of the spectrum to output.

        Notes
        -----
        If two functions g and h are related by the equation

            glm = Z(l) hlm + nlm

        where nlm is a zero-mean random variable, the admittance and spectral
        correlation gamma(l) can be estimated using

            Z(l) = Sgh(l) / Shh(l)
            gamma(l) = Sgh(l) / sqrt( Sgg(l) Shh(l) )

        where Sgh, Shh and Sgg are the cross-power and power spectra of the
        functions g (self) and h (input).
        """
        if not isinstance(hlm, SHCoeffs):
            raise ValueError('hlm must be an SHCoeffs class instance. Input '
                             'type is {:s}.'.format(repr(type(hlm))))

        if lmax is None:
            lmax = min(self.lmax, hlm.lmax)

        sgg = _spectrum(self.coeffs, normalization=self.normalization,
                        lmax=lmax)
        shh = _spectrum(hlm.coeffs, normalization=hlm.normalization, lmax=lmax)
        sgh = _cross_spectrum(self.coeffs,
                              hlm.to_array(normalization=self.normalization,
                                           csphase=self.csphase, lmax=lmax,
                                           errors=False),
                              normalization=self.normalization,
                              lmax=lmax)

        with _np.errstate(invalid='ignore', divide='ignore'):
            admit = sgh / shh
            corr = sgh / _np.sqrt(sgg * shh)
            if errors:
                sigma = (sgg / shh) * (1. - corr**2) / _np.arange(lmax+1) / 2.
                admit = _np.column_stack((admit, _np.sqrt(sigma)))
            return admit, corr

    def volume(self, lmax=None):
        """
        If the function is the real shape of an object, calculate the volume
        of the body.

        Usage
        -----
        volume = x.volume([lmax])

        Returns
        -------
        volume : float
            The volume of the object.

        Parameters
        ----------
        lmax : int, optional, default = x.lmax
            The maximum spherical harmonic degree to use when calculating the
            volume.

        Notes
        -----
        If the function is the real shape of an object, this method will
        calculate the volume of the body exactly by integration. This routine
        raises the function to the nth power, with n from 1 to 3, and
        calculates the spherical harmonic degree and order 0 term. To avoid
        aliases, the function is first expand on a grid that can resolve
        spherical harmonic degrees up to 3*lmax.
        """
        if self.coeffs[0, 0, 0] == 0:
            raise ValueError('The volume of the object can not be calculated '
                             'when the degree and order 0 term is equal to '
                             'zero.')

        if self.kind == 'complex':
            raise ValueError('The volume of the object can not be calculated '
                             'for complex functions.')

        if lmax is None:
            lmax = self.lmax

        r0 = self.coeffs[0, 0, 0]
        grid = self.expand(lmax=min(3*lmax, 2800)) - r0
        h200 = (grid**2).expand(lmax_calc=0).coeffs[0, 0, 0]
        h300 = (grid**3).expand(lmax_calc=0).coeffs[0, 0, 0]

        volume = 4 * _np.pi / 3 * (h300 + 3 * r0 * h200 + r0**3)
        return volume

    def centroid(self):
        """
        Compute the centroid of the body in Cartesian coordinates.

        Usage
        -----
        centroid = x.centroid()

        Returns
        -------
        [x, y, z] : ndarray
            The centroid of the object in meters.

        Notes
        -----
        The centroid is computed as the center of mass of a homogeneous body.
        The units of the input function must be in meters.
        """
        from .shgravcoeffs import SHGravCoeffs as _SHGravCoeffs
        from ..constant import G as _G

        density = 1.
        gm = density * _G.value * self.volume()
        potential = _SHGravCoeffs.from_shape(self, density, gm)
        return potential.center_of_mass

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
        x_rotated : SHCoeffs class instance

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

        Notes
        -----
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
            raise ValueError('convention must be a string. Input type is {:s}.'
                             .format(str(type(convention))))

        if convention.lower() not in ('x', 'y'):
            raise ValueError(
                "convention must be either 'x' or 'y'. " +
                "Provided value is {:s}.".format(repr(convention))
                )

        if convention == 'y':
            if body is True:
                angles = _np.array([-gamma, -beta, -alpha])
            else:
                angles = _np.array([alpha, beta, gamma])
        elif convention == 'x':
            if body is True:
                angles = _np.array([-gamma - _np.pi/2, -beta,
                                    -alpha + _np.pi/2])
            else:
                angles = _np.array([alpha - _np.pi/2, beta, gamma + _np.pi/2])

        if degrees:
            angles = _np.radians(angles)

        if self.lmax > 1200:
            _warnings.warn("The rotate() method is accurate only to about" +
                           " spherical harmonic degree 1200. " +
                           "lmax = {:d}".format(self.lmax),
                           category=RuntimeWarning)

        rot = self._rotate(angles, dj_matrix)
        return rot

    def convert(self, normalization=None, csphase=None, lmax=None, kind=None,
                check=True):
        """
        Return a SHCoeffs class instance with a different normalization
        convention.

        Usage
        -----
        clm = x.convert([normalization, csphase, lmax, kind, check])

        Returns
        -------
        clm : SHCoeffs class instance

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
        kind : str, optional, default = clm.kind
            'real' or 'complex' spherical harmonic coefficients for the output
            class.
        check : bool, optional, default = True
            When converting complex coefficients to real coefficients, if True,
            check if function is entirely real.

        Notes
        -----
        This method will return a new class instance of the spherical
        harmonic coefficients using a different normalization and
        Condon-Shortley phase convention. The coefficients can be converted
        between real and complex form, and a different maximum spherical
        harmonic degree of the output coefficients can be specified. If this
        maximum degree is smaller than the maximum degree of the original
        class, the coefficients will be truncated. Conversely, if this degree
        is larger than the maximum degree of the original class, the
        coefficients of the new class will be zero padded.
        """
        error_coeffs = None
        if normalization is None:
            normalization = self.normalization
        if csphase is None:
            csphase = self.csphase
        if lmax is None:
            lmax = self.lmax
        if kind is None:
            kind = self.kind

        # check argument consistency
        if type(normalization) != str:
            raise ValueError('normalization must be a string. ' +
                             'Input type is {:s}.'
                             .format(str(type(normalization))))
        if normalization.lower() not in ('4pi', 'ortho', 'schmidt', 'unnorm'):
            raise ValueError(
                "normalization must be '4pi', 'ortho', 'schmidt', or " +
                "'unnorm'. Provided value is {:s}."
                .format(repr(normalization)))
        if csphase != 1 and csphase != -1:
            raise ValueError(
                "csphase must be 1 or -1. Input value is {:s}."
                .format(repr(csphase)))

        if (kind != self.kind):
            if (kind == 'complex'):
                temp = self._make_complex()
            else:
                temp = self._make_real(check=check)
            if self.errors is not None:
                coeffs, error_coeffs = temp.to_array(
                    normalization=normalization.lower(), csphase=csphase,
                    lmax=lmax)
            else:
                coeffs = temp.to_array(normalization=normalization.lower(),
                                       csphase=csphase, lmax=lmax)
        else:
            if self.errors is not None:
                coeffs, error_coeffs = temp.to_array(
                    normalization=normalization.lower(), csphase=csphase,
                    lmax=lmax)
            else:
                coeffs = self.to_array(normalization=normalization.lower(),
                                       csphase=csphase, lmax=lmax)

        return SHCoeffs.from_array(coeffs, errors=error_coeffs,
                                   error_kind=self.error_kind,
                                   normalization=normalization.lower(),
                                   csphase=csphase, units=self.units,
                                   copy=False)

    def pad(self, lmax, copy=True):
        """
        Return a SHCoeffs class where the coefficients are zero padded or
        truncated to a different lmax.

        Usage
        -----
        clm = x.pad(lmax, [copy])

        Returns
        -------
        clm : SHCoeffs class instance

        Parameters
        ----------
        lmax : int
            Maximum spherical harmonic degree to output.
        copy : bool, optional, default = True
            If True, make a copy of x when initializing the class instance.
            If False, modify x itself.
        """
        if copy:
            clm = self.copy()
        else:
            clm = self

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

    # ---- Expand the coefficients onto a grid ----
    def expand(self, grid='DH2', lat=None, colat=None, lon=None, degrees=True,
               zeros=None, lmax=None, lmax_calc=None, extend=True):
        """
        Evaluate the spherical harmonic coefficients either on a global grid
        or for a list of coordinates.

        Usage
        -----
        f = x.expand([grid, lmax, lmax_calc, zeros])
        g = x.expand(lat=lat, lon=lon, [lmax_calc, degrees])
        g = x.expand(colat=colat, lon=lon, [lmax_calc, degrees])

        Returns
        -------
        f : SHGrid class instance
        g : float, ndarray, or list

        Parameters
        ----------
        lat : int, float, ndarray, or list, optional, default = None
            Latitude coordinates where the function is to be evaluated.
        colat : int, float, ndarray, or list, optional, default = None
            Colatitude coordinates where the function is to be evaluated.
        lon : int, float, ndarray, or list, optional, default = None
            Longitude coordinates where the function is to be evaluated.
        degrees : bool, optional, default = True
            True if lat, colat and lon are in degrees, False if in radians.
        grid : str, optional, default = 'DH2'
            'DH' or 'DH1' for an equisampled lat/lon grid with nlat=nlon,
            'DH2' for an equidistant lat/lon grid with nlon=2*nlat, or 'GLQ'
            for a Gauss-Legendre quadrature grid.
        lmax : int, optional, default = x.lmax
            The maximum spherical harmonic degree, which determines the grid
            spacing of the output grid.
        lmax_calc : int, optional, default = x.lmax
            The maximum spherical harmonic degree to use when evaluating the
            function.
        extend : bool, optional, default = True
            If True, compute the longitudinal band for 360 E (DH and GLQ grids)
            and the latitudinal band for 90 S (DH grids only).
        zeros : ndarray, optional, default = None
            The cos(colatitude) nodes used in the Gauss-Legendre Quadrature
            grids.

        Notes
        -----
        This method either (1) evaluates the spherical harmonic coefficients on
        a global grid and returns an SHGrid class instance, or (2) evaluates
        the spherical harmonic coefficients for a list of (co)latitude and
        longitude coordinates. For the first case, the grid type is defined
        by the optional parameter grid, which can be 'DH', 'DH2' or 'GLQ'.For
        the second case, the optional parameters lon and either colat or lat
        must be provided.
        """
        if lat is not None and colat is not None:
            raise ValueError('lat and colat can not both be specified.')

        if lat is not None and lon is not None:
            if lmax_calc is None:
                lmax_calc = self.lmax

            values = self._expand_coord(lat=lat, lon=lon, degrees=degrees,
                                        lmax_calc=lmax_calc)
            return values

        if colat is not None and lon is not None:
            if lmax_calc is None:
                lmax_calc = self.lmax

            if type(colat) is list:
                lat = list(map(lambda x: 90 - x, colat))
            else:
                lat = 90 - colat

            values = self._expand_coord(lat=lat, lon=lon, degrees=degrees,
                                        lmax_calc=lmax_calc)
            return values

        else:
            if lmax is None:
                lmax = self.lmax
            if lmax_calc is None:
                lmax_calc = lmax

            if type(grid) != str:
                raise ValueError('grid must be a string. Input type is {:s}.'
                                 .format(str(type(grid))))

            if grid.upper() in ('DH', 'DH1'):
                gridout = self._expandDH(sampling=1, lmax=lmax,
                                         lmax_calc=lmax_calc, extend=extend)
            elif grid.upper() == 'DH2':
                gridout = self._expandDH(sampling=2, lmax=lmax,
                                         lmax_calc=lmax_calc, extend=extend)
            elif grid.upper() == 'GLQ':
                gridout = self._expandGLQ(zeros=zeros, lmax=lmax,
                                          lmax_calc=lmax_calc, extend=extend)
            else:
                raise ValueError(
                    "grid must be 'DH', 'DH1', 'DH2', or 'GLQ'. " +
                    "Input value is {:s}.".format(repr(grid)))

            return gridout

    # ---- Plotting routines ----
    def plot_spectrum(self, convention='power', unit='per_l', base=10.,
                      lmax=None, xscale='lin', yscale='log', grid=True,
                      legend=None, legend_error='error', legend_loc='best',
                      axes_labelsize=None, tick_labelsize=None, show=True,
                      ax=None, fname=None, **kwargs):
        """
        Plot the spectrum as a function of spherical harmonic degree.

        Usage
        -----
        x.plot_spectrum([convention, unit, base, lmax, xscale, yscale, grid,
                         axes_labelsize, tick_labelsize, legend, legend_error,
                         legend_loc, show, ax, fname, **kwargs])

        Parameters
        ----------
        convention : str, optional, default = 'power'
            The type of spectrum to plot: 'power' for power spectrum,
            'energy' for energy spectrum, and 'l2norm' for the l2 norm
            spectrum.
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
        legend_error : str, optional, default = 'error'
            Text to use for the legend of the error spectrum.
        legend_loc : str, optional, default = 'best'
            Location of the legend, such as 'upper right' or 'lower center'
            (see pyplot.legend for all options).
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

        Notes
        -----
        This method plots either the power spectrum, energy spectrum, or
        l2-norm spectrum. If the error coefficients are specified, the error
        spectrum will also be plotted. Total power is defined as the integral
        of the function squared over all space, divided by the area the
        function spans. If the mean of the function is zero, this is equivalent
        to the variance of the function. The total energy is the integral of
        the function squared over all space and is 4pi times the total power.
        For normalized coefficients ('4pi', 'ortho', or 'schmidt'), the l2-norm
        is the sum of the magnitude of the coefficients squared.

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
            spectrum, error_spectrum = self.spectrum(convention=convention,
                                                     unit=unit, base=base,
                                                     lmax=lmax)
        else:
            spectrum = self.spectrum(convention=convention, unit=unit,
                                     base=base, lmax=lmax)

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
        if convention == 'Energy':
            axes.set_ylabel('Energy', fontsize=axes_labelsize)
            if legend is None:
                if (unit == 'per_l'):
                    legend = 'Energy per degree'
                elif (unit == 'per_lm'):
                    legend = 'Energy per coefficient'
                elif (unit == 'per_dlogl'):
                    legend = 'Energy per log bandwidth'
        elif convention == 'l2norm':
            axes.set_ylabel('l2 norm', fontsize=axes_labelsize)
            if legend is None:
                if (unit == 'per_l'):
                    legend = 'l2 norm per degree'
                elif (unit == 'per_lm'):
                    legend = 'l2 norm per coefficient'
                elif (unit == 'per_dlogl'):
                    legend = 'l2 norm per log bandwidth'
        else:
            axes.set_ylabel('Power', fontsize=axes_labelsize)
            if legend is None:
                if (unit == 'per_l'):
                    legend = 'Power per degree'
                elif (unit == 'per_lm'):
                    legend = 'Power per coefficient'
                elif (unit == 'per_dlogl'):
                    legend = 'Power per log bandwidth'

        if xscale == 'log':
            axes.set_xscale('log', base=base)
        if yscale == 'log':
            axes.set_yscale('log', base=base)

        if xscale == 'log':
            axes.plot(ls[1:lmax+1], spectrum[1:lmax+1], label=legend, **kwargs)
            if self.errors is not None:
                axes.plot(ls[1:lmax+1], error_spectrum[1:lmax+1],
                          label='error', **kwargs)
        else:
            axes.plot(ls[:lmax+1], spectrum[:lmax+1], label=legend, **kwargs)
            if self.errors is not None:
                axes.plot(ls[:lmax+1], error_spectrum[:lmax+1],
                          label=legend_error, **kwargs)
            if ax is None:
                axes.set(xlim=(ls[0], ls[lmax]))
            else:
                axes.set(xlim=(ls[0], max(ls[lmax], ax.get_xbound()[1])))

        axes.grid(grid, which='major')
        axes.minorticks_on()
        axes.tick_params(labelsize=tick_labelsize)
        axes.legend(loc=legend_loc)

        if ax is None:
            fig.tight_layout(pad=0.5)
            if show:
                fig.show()
            if fname is not None:
                fig.savefig(fname)
            return fig, axes

    def plot_cross_spectrum(self, clm, convention='power', unit='per_l',
                            base=10., lmax=None, xscale='lin', yscale='log',
                            grid=True, legend=None, legend_loc='best',
                            axes_labelsize=None, tick_labelsize=None,
                            show=True, ax=None, fname=None, **kwargs):
        """
        Plot the cross-spectrum of two functions.

        Usage
        -----
        x.plot_cross_spectrum(clm, [convention, unit, base, lmax, xscale,
                                    yscale, grid, axes_labelsize,
                                    tick_labelsize, legend, legend_loc, show,
                                    ax, fname, **kwargs])

        Parameters
        ----------
        clm : SHCoeffs class instance.
            The second function used in computing the cross-spectrum.
        convention : str, optional, default = 'power'
            The type of spectrum to plot: 'power' for power spectrum,
            'energy' for energy spectrum, and 'l2norm' for the l2 norm
            spectrum.
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
        legend_loc : str, optional, default = 'best'
            Location of the legend, such as 'upper right' or 'lower center'
            (see pyplot.legend for all options).
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

        Notes
        -----
        This method plots either the cross-power spectrum, cross-energy
        spectrum, or l2-cross-norm spectrum. Total cross-power is defined as
        the integral of the function times the conjugate of clm over all space,
        divided by the area the functions span. If the means of the functions
        are zero, this is equivalent to the covariance of the two functions.
        The total cross-energy is the integral of this function times the
        conjugate of clm over all space and is 4pi times the total power. The
        l2-cross-norm is the sum of this function times the conjugate of clm
        over all angular orders as a function of spherical harmonic degree.

        The output spectrum can be expresed using one of three units. 'per_l'
        returns the contribution to the total spectrum from all angular orders
        at degree l. 'per_lm' returns the average contribution to the total
        spectrum from a single coefficient at degree l, and is equal to the
        'per_l' spectrum divided by (2l+1). 'per_dlogl' returns the
        contribution to the total spectrum from all angular orders over an
        infinitessimal logarithmic degree band. The contrubution in the band
        dlog_a(l) is spectrum(l, 'per_dlogl')*dlog_a(l), where a is the base,
        and where spectrum(l, 'per_dlogl) is equal to
        spectrum(l, 'per_l')*l*log(a). If the input fields are complex, the
        absolute value of the cross-spectrum will be plotted.
        """
        if not isinstance(clm, SHCoeffs):
            raise ValueError('clm must be an SHCoeffs class instance. Input '
                             'type is {:s}.'.format(repr(type(clm))))

        if lmax is None:
            lmax = min(self.lmax, clm.lmax)

        spectrum = self.cross_spectrum(clm, convention=convention, unit=unit,
                                       base=base, lmax=lmax)
        spectrum = abs(spectrum)

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
        if convention == 'Energy':
            axes.set_ylabel('Energy', fontsize=axes_labelsize)
            if legend is None:
                if (unit == 'per_l'):
                    legend = 'Energy per degree'
                elif (unit == 'per_lm'):
                    legend = 'Energy per coefficient'
                elif (unit == 'per_dlogl'):
                    legend = 'Energy per log bandwidth'
        elif convention == 'l2norm':
            axes.set_ylabel('l2 norm', fontsize=axes_labelsize)
            if legend is None:
                if (unit == 'per_l'):
                    legend = 'l2 norm per degree'
                elif (unit == 'per_lm'):
                    legend = 'l2 norm per coefficient'
                elif (unit == 'per_dlogl'):
                    legend = 'l2 norm per log bandwidth'
        else:
            axes.set_ylabel('Power', fontsize=axes_labelsize)
            if legend is None:
                if (unit == 'per_l'):
                    legend = 'Power per degree'
                elif (unit == 'per_lm'):
                    legend = 'Power per coefficient'
                elif (unit == 'per_dlogl'):
                    legend = 'Power per log bandwidth'

        if xscale == 'log':
            axes.set_xscale('log', base=base)
        if yscale == 'log':
            axes.set_yscale('log', base=base)

        if xscale == 'log':
            axes.plot(ls[1:lmax+1], spectrum[1:lmax+1], label=legend, **kwargs)
        else:
            axes.plot(ls[:lmax+1], spectrum[:lmax+1], label=legend, **kwargs)
            if ax is None:
                axes.set(xlim=(ls[0], ls[lmax]))
            else:
                axes.set(xlim=(ls[0], max(ls[lmax], ax.get_xbound()[1])))

        axes.grid(grid, which='major')
        axes.minorticks_on()
        axes.tick_params(labelsize=tick_labelsize)
        axes.legend(loc=legend_loc)

        if ax is None:
            fig.tight_layout(pad=0.5)
            if show:
                fig.show()
            if fname is not None:
                fig.savefig(fname)
            return fig, axes

    def plot_spectrum2d(self, convention='power', xscale='lin', yscale='lin',
                        grid=True, axes_labelsize=None, tick_labelsize=None,
                        vscale='log', vrange=None, vmin=None, vmax=None,
                        lmax=None, errors=False, show=True, ax=None,
                        fname=None):
        """
        Plot the spectrum as a function of spherical harmonic degree and order.

        Usage
        -----
        x.plot_spectrum2d([convention, xscale, yscale, grid, axes_labelsize,
                           tick_labelsize, vscale, vrange, vmin, vmax, lmax,
                           errors, show, ax, fname])

        Parameters
        ----------
        convention : str, optional, default = 'power'
            The type of spectrum to plot: 'power' for power spectrum,
            'energy' for energy spectrum, and 'l2norm' for the l2 norm
            spectrum.
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
            Colormap range (min, max) relative to the maximum value. If None,
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

        Notes
        -----
        This method plots either the power, energy, or l2-norm for each
        spherical harmonic degree and order of the function. Total power is
        defined as the integral of the function squared over all space,
        divided by the area the function spans. If the mean of the function is
        zero, this is equivalent to the variance of the function. The total
        energy is the integral of the function squared over all space and is
        4pi times the total power. For normalized coefficients ('4pi',
        'ortho', or 'schmidt'), the l2-norm is the sum of the magnitude of the
        coefficients squared.
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
        mpositive = coeffs[0, :lmax + 1, :lmax + 1] * \
            coeffs[0, :lmax + 1, :lmax + 1].conj()
        mpositive[~self.mask[0, :lmax + 1, :lmax + 1]] = _np.nan
        mnegative = coeffs[1, :lmax + 1, :lmax + 1] * \
            coeffs[1, :lmax + 1, :lmax + 1].conj()
        mnegative[~self.mask[1, :lmax + 1, :lmax + 1]] = _np.nan

        spectrum[:, :lmax] = _np.fliplr(mnegative)[:, :lmax]
        spectrum[:, lmax:] = mpositive

        if (convention.lower() == 'l2norm'):
            if self.normalization == 'unnorm':
                raise ValueError("convention can not be set to 'l2norm' " +
                                 "when using unnormalized harmonics.")
            else:
                pass
        elif convention.lower() in ('power', 'energy'):
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
                    "or 'unnorm'. Input value is {:s}."
                    .format(repr(self.normalization)))
        else:
            raise ValueError(
                "convention must be 'power', 'energy', or 'l2norm'. " +
                "Input value is {:s}.".format(repr(convention)))

        if convention == 'energy':
            spectrum *= 4.0 * _np.pi

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
                "Input value is {:s}.".format(repr(vscale)))

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
                "Input value is {:s}.".format(repr(xscale)))

        if (yscale == 'lin'):
            axes.set(ylim=(-lmax - 0.5, lmax + 0.5))
        elif (yscale == 'log'):
            axes.set(yscale='symlog', ylim=(-lmax - 0.5, lmax + 0.5))
        else:
            raise ValueError(
                "yscale must be 'lin' or 'log'. " +
                "Input value is {:s}.".format(repr(yscale)))

        cb = _plt.colorbar(cmesh, ax=ax)

        if (convention == 'energy'):
            cb.set_label('Energy per coefficient', fontsize=axes_labelsize)
        elif (convention == 'power'):
            cb.set_label('Power per coefficient', fontsize=axes_labelsize)
        else:
            cb.set_label('Magnitude-squared coefficient',
                         fontsize=axes_labelsize)

        cb.ax.tick_params(labelsize=tick_labelsize)
        axes.set_xlabel('Spherical harmonic degree', fontsize=axes_labelsize)
        axes.set_ylabel('Spherical harmonic order', fontsize=axes_labelsize)
        axes.tick_params(labelsize=tick_labelsize)
        axes.minorticks_on()
        axes.grid(grid, which='major')

        if ax is None:
            fig.tight_layout(pad=0.5)
            if show:
                fig.show()
            if fname is not None:
                fig.savefig(fname)
            return fig, axes

    def plot_cross_spectrum2d(self, clm, convention='power', xscale='lin',
                              yscale='lin', grid=True, axes_labelsize=None,
                              tick_labelsize=None, vscale='log', vrange=None,
                              vmin=None, vmax=None, lmax=None, show=True,
                              ax=None, fname=None):
        """
        Plot the cross-spectrum of two functions as a function of spherical
        harmonic degree and order.

        Usage
        -----
        x.plot_cross_spectrum2d(clm, [convention, xscale, yscale, grid,
                                      axes_labelsize, tick_labelsize, vscale,
                                      vrange, vmin, vmax, lmax, show, ax,
                                      fname])

        Parameters
        ----------
        clm : SHCoeffs class instance.
            The second function used in computing the cross-spectrum.
        convention : str, optional, default = 'power'
            The type of spectrum to plot: 'power' for power spectrum,
            'energy' for energy spectrum, and 'l2norm' for the l2 norm
            spectrum.
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
            Colormap range (min, max) relative to the maximum value. If None,
            scale the image to the maximum and minimum values.
        vmin : float, optional, default=None
            The minmum range of the colormap. If None, the minimum value of the
            spectrum will be used.
        vmax : float, optional, default=None
            The maximum range of the colormap. If None, the maximum value of
            the spectrum will be used.
        lmax : int, optional, default = self.lmax
            The maximum spherical harmonic degree to plot.
        show : bool, optional, default = True
            If True, plot to the screen.
        ax : matplotlib axes object, optional, default = None
            A single matplotlib axes object where the plot will appear.
        fname : str, optional, default = None
            If present, and if axes is not specified, save the image to the
            specified file.

        Notes
        -----
        This method plots either the power, energy, or l2-norm for each
        spherical harmonic degree and order of the function. Total power is
        defined as the integral of the function squared over all space,
        divided by the area the function spans. If the mean of the function is
        zero, this is equivalent to the variance of the function. The total
        energy is the integral of the function squared over all space and is
        4pi times the total power. For normalized coefficients ('4pi',
        'ortho', or 'schmidt'), the l2-norm is the sum of the magnitude of the
        coefficients squared. If the input fields are complex, the absolute
        value of the cross-spectrum will be plotted.
        """
        if not isinstance(clm, SHCoeffs):
            raise ValueError('clm must be an SHCoeffs class instance. Input '
                             'type is {:s}.'.format(repr(type(clm))))

        if axes_labelsize is None:
            axes_labelsize = _mpl.rcParams['axes.labelsize']
        if tick_labelsize is None:
            tick_labelsize = _mpl.rcParams['xtick.labelsize']

        if lmax is None:
            lmax = min(self.lmax, clm.lmax)

        degrees = _np.arange(lmax + 1)

        coeffs = clm.to_array(normalization=self.normalization,
                              csphase=self.csphase,
                              lmax=self.lmax)

        # Create the matrix of the spectrum for each coefficient
        spectrum = _np.empty((lmax + 1, 2 * lmax + 1))
        mpositive = _np.abs(self.coeffs[0, :lmax + 1, :lmax + 1] *
                            coeffs[0, :lmax + 1, :lmax + 1].conj())
        mpositive[~self.mask[0, :lmax + 1, :lmax + 1]] = _np.nan
        mnegative = _np.abs(self.coeffs[1, :lmax + 1, :lmax + 1] *
                            coeffs[1, :lmax + 1, :lmax + 1].conj())
        mnegative[~self.mask[1, :lmax + 1, :lmax + 1]] = _np.nan

        spectrum[:, :lmax] = _np.fliplr(mnegative)[:, :lmax]
        spectrum[:, lmax:] = mpositive

        if (convention.lower() == 'l2norm'):
            if self.normalization == 'unnorm':
                raise ValueError("convention can not be set to 'l2norm' " +
                                 "when using unnormalized harmonics.")
            else:
                pass
        elif convention.lower() in ('power', 'energy'):
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
                    "or 'unnorm'. Input value is {:s}."
                    .format(repr(self.normalization)))
        else:
            raise ValueError(
                "convention must be 'power', 'energy', or 'l2norm'. " +
                "Input value is {:s}.".format(repr(convention)))

        if convention == 'energy':
            spectrum *= 4.0 * _np.pi

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
                "Input value is {:s}.".format(repr(vscale)))

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
                "Input value is {:s}.".format(repr(xscale)))

        if (yscale == 'lin'):
            axes.set(ylim=(-lmax - 0.5, lmax + 0.5))
        elif (yscale == 'log'):
            axes.set(yscale='symlog', ylim=(-lmax - 0.5, lmax + 0.5))
        else:
            raise ValueError(
                "yscale must be 'lin' or 'log'. " +
                "Input value is {:s}.".format(repr(yscale)))

        cb = _plt.colorbar(cmesh, ax=ax)

        if (convention == 'energy'):
            cb.set_label('Energy per coefficient', fontsize=axes_labelsize)
        elif (convention == 'power'):
            cb.set_label('Power per coefficient', fontsize=axes_labelsize)
        else:
            cb.set_label('Magnitude-squared coefficient',
                         fontsize=axes_labelsize)

        cb.ax.tick_params(labelsize=tick_labelsize)
        axes.set_xlabel('Spherical harmonic degree', fontsize=axes_labelsize)
        axes.set_ylabel('Spherical harmonic order', fontsize=axes_labelsize)
        axes.tick_params(labelsize=tick_labelsize)
        axes.minorticks_on()
        axes.grid(grid, which='major')

        if ax is None:
            fig.tight_layout(pad=0.5)
            if show:
                fig.show()
            if fname is not None:
                fig.savefig(fname)
            return fig, axes

    def plot_admitcorr(self, hlm, errors=True, style='separate', lmax=None,
                       grid=True, legend=None, legend_loc='best',
                       axes_labelsize=None, tick_labelsize=None,
                       elinewidth=0.75, show=True, ax=None, ax2=None,
                       fname=None, **kwargs):
        """
        Plot the admittance and/or correlation with another function.

        Usage
        -----
        x.plot_admitcorr(hlm, [errors, style, lmax, grid, legend, legend_loc,
                               axes_labelsize, tick_labelsize, elinewidth,
                               show, ax, ax2, fname, **kwargs])

        Parameters
        ----------
        hlm : SHCoeffs class instance.
            The second function used in computing the admittance and
            correlation.
        errors : bool, optional, default = True
            Plot the uncertainty of the admittance.
        style : str, optional, default = 'separate'
            Style of the plot. 'separate' to plot the admittance and
            correlation in separate plots, 'combined' to plot the admittance
            and correlation in a single plot, 'admit' to plot only the
            admittance, or 'corr' to plot only the correlation.
        lmax : int, optional, default = self.lmax
            The maximum spherical harmonic degree to plot.
        grid : bool, optional, default = True
            If True, plot grid lines. grid is set to False when style is
            'combined'.
        legend : str, optional, default = None
            Text to use for the legend. If style is 'combined' or 'separate',
            provide a list of two strings for the admittance and correlation,
            respectively.
        legend_loc : str, optional, default = 'best'
            Location of the legend, such as 'upper right' or 'lower center'
            (see pyplot.legend for all options). If style is 'separate',
            provide a list of two strings for the admittance and correlation,
            respectively.
        axes_labelsize : int, optional, default = None
            The font size for the x and y axes labels.
        tick_labelsize : int, optional, default = None
            The font size for the x and y tick labels.
        elinewidth : float, optional, default = 0.75
            Line width of the error bars when errors is True.
        show : bool, optional, default = True
            If True, plot to the screen.
        ax : matplotlib axes object, optional, default = None
            A single matplotlib axes object where the plot will appear.
        ax2 : matplotlib axes object, optional, default = None
            A single matplotlib axes object where the second plot will appear
            when style is 'separate'.
        fname : str, optional, default = None
            If present, and if axes is not specified, save the image to the
            specified file.
        **kwargs : keyword arguments, optional
            Keyword arguments for pyplot.plot() and pyplot.errorbar().

        Notes
        -----
        If two functions g and h are related by the equation

            glm = Z(l) hlm + nlm

        where nlm is a zero-mean random variable, the admittance and spectral
        correlation gamma(l) can be estimated using

            Z(l) = Sgh(l) / Shh(l)
            gamma(l) = Sgh(l) / sqrt( Sgg(l) Shh(l) )

        where Sgh, Shh and Sgg are the cross-power and power spectra of the
        functions g (self) and h (input).
        """
        if lmax is None:
            lmax = min(self.lmax, hlm.lmax)

        if style in ('combined', 'separate'):
            admit, corr = self.admitcorr(hlm, errors=errors, lmax=lmax)
        elif style == 'corr':
            corr = self.correlation(hlm, lmax=lmax)
        elif style == 'admit':
            admit = self.admittance(hlm, errors=errors, lmax=lmax)
        else:
            raise ValueError("style must be 'combined', 'separate', 'admit' "
                             "or 'corr'. Input value is {:s}"
                             .format(repr(style)))

        ls = _np.arange(lmax + 1)

        if style == 'separate':
            if ax is None:
                scale = 0.4
                figsize = (_mpl.rcParams['figure.figsize'][0],
                           _mpl.rcParams['figure.figsize'][0]*scale)
                fig, (axes, axes2) = _plt.subplots(1, 2, figsize=figsize)
            else:
                axes = ax
                axes2 = ax2
        elif style == 'combined':
            if ax is None:
                fig, axes = _plt.subplots(1, 1)
                axes2 = axes.twinx()
            else:
                axes = ax
                axes2 = axes.twinx()
        else:
            if ax is None:
                fig, axes = _plt.subplots(1, 1)
            else:
                axes = ax

        if style in ('separate', 'combined'):
            admitax = axes
            corrax = axes2
        elif style == 'admit':
            admitax = axes
        elif style == 'corr':
            corrax = axes

        if legend is None:
            legend = [None, None]
        elif style == 'admit':
            legend = [legend, None]
            legend_loc = [legend_loc, None]
        elif style == 'corr':
            legend = [None, legend]
            legend_loc = [None, legend_loc]
        elif style == 'combined':
            legend_loc = [legend_loc, legend_loc]
        else:
            if type(legend_loc) is str:
                legend_loc = [legend_loc, legend_loc]

        if axes_labelsize is None:
            axes_labelsize = _mpl.rcParams['axes.labelsize']
        if tick_labelsize is None:
            tick_labelsize = _mpl.rcParams['xtick.labelsize']

        if style in ('admit', 'separate', 'combined'):
            if errors:
                admitax.errorbar(ls, admit[:, 0], yerr=admit[:, 1],
                                 label=legend[0], elinewidth=elinewidth,
                                 **kwargs)
            else:
                admitax.plot(ls, admit, label=legend[0], **kwargs)
            if ax is None:
                admitax.set(xlim=(0, lmax))
            else:
                admitax.set(xlim=(0, max(lmax, ax.get_xbound()[1])))
            admitax.set_xlabel('Spherical harmonic degree',
                               fontsize=axes_labelsize)
            admitax.set_ylabel('Admittance', fontsize=axes_labelsize)
            admitax.minorticks_on()
            admitax.tick_params(labelsize=tick_labelsize)
            if legend[0] is not None:
                if style != 'combined':
                    admitax.legend(loc=legend_loc[0])
            if style != 'combined':
                admitax.grid(grid, which='major')

        if style in ('corr', 'separate', 'combined'):
            if style == 'combined':
                # plot with next color
                next(corrax._get_lines.prop_cycler)['color']
            corrax.plot(ls, corr, label=legend[1], **kwargs)
            if ax is None:
                corrax.set(xlim=(0, lmax))
                corrax.set(ylim=(-1, 1))
            else:
                corrax.set(xlim=(0, max(lmax, ax.get_xbound()[1])))
            corrax.set_xlabel('Spherical harmonic degree',
                              fontsize=axes_labelsize)
            corrax.set_ylabel('Correlation', fontsize=axes_labelsize)
            corrax.minorticks_on()
            corrax.tick_params(labelsize=tick_labelsize)
            if legend[1] is not None:
                if style == 'combined':
                    lines, labels = admitax.get_legend_handles_labels()
                    lines2, labels2 = corrax.get_legend_handles_labels()
                    corrax.legend(lines + lines2, labels + labels2,
                                  loc=legend_loc[1])
                else:
                    corrax.legend(loc=legend_loc[1])
            if style != 'combined':
                corrax.grid(grid, which='major')

        if ax is None:
            fig.tight_layout(pad=0.5)
            if show:
                fig.show()
            if fname is not None:
                fig.savefig(fname)
            if style in ('separate', 'combined'):
                return fig, (axes, axes2)
            else:
                return fig, axes

    def plot_admittance(self, hlm, errors=True, lmax=None, grid=True,
                        legend=None, legend_loc='best', axes_labelsize=None,
                        tick_labelsize=None, elinewidth=0.75, show=True,
                        ax=None, fname=None, **kwargs):
        """
        Plot the admittance with another function.

        Usage
        -----
        x.plot_admittance(hlm, [errors, lmax, grid, legend, legend_loc,
                                axes_labelsize, tick_labelsize, elinewidth,
                                show, ax, fname, **kwargs])

        Parameters
        ----------
        hlm : SHCoeffs class instance.
            The second function used in computing the admittance.
        errors : bool, optional, default = True
            Plot the uncertainty of the admittance.
        lmax : int, optional, default = self.lmax
            The maximum spherical harmonic degree to plot.
        grid : bool, optional, default = True
            If True, plot grid lines.
        legend : str, optional, default = None
            Text to use for the legend.
        legend_loc : str, optional, default = 'best'
            Location of the legend, such as 'upper right' or 'lower center'
            (see pyplot.legend for all options).
        axes_labelsize : int, optional, default = None
            The font size for the x and y axes labels.
        tick_labelsize : int, optional, default = None
            The font size for the x and y tick labels.
        elinewidth : float, optional, default = 0.75
            Line width of the error bars when errors is True.
        show : bool, optional, default = True
            If True, plot to the screen.
        ax : matplotlib axes object, optional, default = None
            A single matplotlib axes object where the plot will appear.
        fname : str, optional, default = None
            If present, and if axes is not specified, save the image to the
            specified file.
        **kwargs : keyword arguments, optional
            Keyword arguments for pyplot.plot() and pyplot.errorbar().

        Notes
        -----
        If two functions g and h are related by the equation

            glm = Z(l) hlm + nlm

        where nlm is a zero-mean random variable, the admittance can be
        estimated using

            Z(l) = Sgh(l) / Shh(l)

        where Sgh and Shh are the cross-power and power spectra of the
        functions g (self) and h (input).
        """
        return self.plot_admitcorr(hlm, errors=errors, style='admit',
                                   lmax=lmax, grid=grid, legend=legend,
                                   legend_loc=legend_loc,
                                   axes_labelsize=axes_labelsize,
                                   tick_labelsize=tick_labelsize,
                                   elinewidth=elinewidth, show=True,
                                   fname=fname, ax=ax, **kwargs)

    def plot_correlation(self, hlm, lmax=None, grid=True, legend=None,
                         legend_loc='best', axes_labelsize=None,
                         tick_labelsize=None, elinewidth=0.75, show=True,
                         ax=None, fname=None, **kwargs):
        """
        Plot the correlation with another function.

        Usage
        -----
        x.plot_correlation(hlm, [lmax, grid, legend, legend_loc,
                                 axes_labelsize, tick_labelsize, elinewidth,
                                 show, ax, fname, **kwargs])

        Parameters
        ----------
        hlm : SHCoeffs class instance.
            The second function used in computing the correlation.
        lmax : int, optional, default = self.lmax
            The maximum spherical harmonic degree to plot.
        grid : bool, optional, default = True
            If True, plot grid lines.
        legend : str, optional, default = None
            Text to use for the legend.
        legend_loc : str, optional, default = 'best'
            Location of the legend, such as 'upper right' or 'lower center'
            (see pyplot.legend for all options).
        axes_labelsize : int, optional, default = None
            The font size for the x and y axes labels.
        tick_labelsize : int, optional, default = None
            The font size for the x and y tick labels.
        elinewidth : float, optional, default = 0.75
            Line width of the error bars when errors is True.
        show : bool, optional, default = True
            If True, plot to the screen.
        ax : matplotlib axes object, optional, default = None
            A single matplotlib axes object where the plot will appear.
        fname : str, optional, default = None
            If present, and if axes is not specified, save the image to the
            specified file.
        **kwargs : keyword arguments, optional
            Keyword arguments for pyplot.plot() and pyplot.errorbar().

        Notes
        -----
        The spectral correlation is defined as

            gamma(l) = Sgh(l) / sqrt( Sgg(l) Shh(l) )

        where Sgh, Shh and Sgg are the cross-power and power spectra of the
        functions g (self) and h (input).
        """
        return self.plot_admitcorr(hlm, style='corr', lmax=lmax, grid=grid,
                                   legend=legend, legend_loc=legend_loc,
                                   axes_labelsize=axes_labelsize,
                                   tick_labelsize=tick_labelsize,
                                   show=True, fname=fname, ax=ax, **kwargs)


# ================== REAL SPHERICAL HARMONICS ================

class SHRealCoeffs(SHCoeffs):
    """Real Spherical Harmonics Coefficient class."""

    @staticmethod
    def istype(kind):
        """Test if class is Real or Complex."""
        return kind == 'real'

    def __init__(self, coeffs, errors=None, error_kind=None,
                 normalization='4pi', csphase=1, units=None, copy=True,
                 header=None, header2=None):
        """Initialize Real SH Coefficients."""
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
        self.header2 = header2
        self.units = units
        self.error_kind = error_kind

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

    def _make_complex(self):
        """Convert the real SHCoeffs class to the complex class."""
        rcomplex_coeffs = _shtools.SHrtoc(self.coeffs,
                                          convention=1, switchcs=0)

        # These coefficients are using real floats, and need to be
        # converted to complex form.
        complex_coeffs = _np.zeros((2, self.lmax+1, self.lmax+1),
                                   dtype='complex')
        complex_coeffs[0, :, :] = (rcomplex_coeffs[0, :, :] + 1j *
                                   rcomplex_coeffs[1, :, :])
        complex_coeffs[1, :, :] = complex_coeffs[0, :, :].conjugate()
        for m in self.degrees():
            if m % 2 == 1:
                complex_coeffs[1, :, m] = - complex_coeffs[1, :, m]

        # complex_coeffs is initialized in this function and can be
        # passed as reference
        return SHCoeffs.from_array(complex_coeffs,
                                   normalization=self.normalization,
                                   csphase=self.csphase, units=self.units,
                                   copy=False)

    def _rotate(self, angles, dj_matrix):
        """Rotate the coefficients by the Euler angles alpha, beta, gamma."""
        if dj_matrix is None:
            dj_matrix = _shtools.djpi2(self.lmax + 1)

        # The coefficients need to be 4pi normalized with csphase = 1
        coeffs = _shtools.SHRotateRealCoef(
            self.to_array(normalization='4pi', csphase=1), angles, dj_matrix)

        # Convert 4pi normalized coefficients to the same normalization
        # as the unrotated coefficients.
        if self.normalization != '4pi' or self.csphase != 1:
            temp = _convert(coeffs, normalization_in='4pi', csphase_in=1,
                            normalization_out=self.normalization,
                            csphase_out=self.csphase)
            return SHCoeffs.from_array(
                temp, normalization=self.normalization,
                csphase=self.csphase, units=self.units, copy=False)
        else:
            return SHCoeffs.from_array(coeffs, units=self.units, copy=False)

    def _expandDH(self, sampling, lmax, lmax_calc, extend):
        """Evaluate the coefficients on a Driscoll and Healy (1994) grid."""
        from .shgrid import SHGrid
        if self.normalization == '4pi':
            norm = 1
        elif self.normalization == 'schmidt':
            norm = 2
        elif self.normalization == 'unnorm':
            norm = 3
        elif self.normalization == 'ortho':
            norm = 4
        else:
            raise ValueError(
                "Normalization must be '4pi', 'ortho', 'schmidt', or " +
                "'unnorm'. Input value is {:s}."
                .format(repr(self.normalization)))

        data = _shtools.MakeGridDH(self.coeffs, sampling=sampling, norm=norm,
                                   csphase=self.csphase, lmax=lmax,
                                   lmax_calc=lmax_calc, extend=extend)
        gridout = SHGrid.from_array(data, grid='DH', units=self.units,
                                    copy=False)
        return gridout

    def _expandGLQ(self, zeros, lmax, lmax_calc, extend):
        """Evaluate the coefficients on a Gauss Legendre quadrature grid."""
        from .shgrid import SHGrid
        if self.normalization == '4pi':
            norm = 1
        elif self.normalization == 'schmidt':
            norm = 2
        elif self.normalization == 'unnorm':
            norm = 3
        elif self.normalization == 'ortho':
            norm = 4
        else:
            raise ValueError(
                "Normalization must be '4pi', 'ortho', 'schmidt', or " +
                "'unnorm'. Input value is {:s}."
                .format(repr(self.normalization)))

        if zeros is None:
            zeros, weights = _shtools.SHGLQ(self.lmax)

        data = _shtools.MakeGridGLQ(self.coeffs, zeros, norm=norm,
                                    csphase=self.csphase, lmax=lmax,
                                    lmax_calc=lmax_calc, extend=extend)
        gridout = SHGrid.from_array(data, grid='GLQ', units=self.units,
                                    copy=False)
        return gridout

    def _expand_coord(self, lat, lon, lmax_calc, degrees):
        """Evaluate the function at the coordinates lat and lon."""
        if self.normalization == '4pi':
            norm = 1
        elif self.normalization == 'schmidt':
            norm = 2
        elif self.normalization == 'unnorm':
            norm = 3
        elif self.normalization == 'ortho':
            norm = 4
        else:
            raise ValueError(
                "Normalization must be '4pi', 'ortho', 'schmidt', or " +
                "'unnorm'. Input value is {:s}."
                .format(repr(self.normalization)))

        if degrees is True:
            latin = lat
            lonin = lon
        else:
            latin = _np.rad2deg(lat)
            lonin = _np.rad2deg(lon)

        if type(lat) is not type(lon):
            raise ValueError('lat and lon must be of the same type. ' +
                             'Input types are {:s} and {:s}.'
                             .format(repr(type(lat)), repr(type(lon))))

        if type(lat) is int or type(lat) is float or type(lat) is _np.float_:
            return _shtools.MakeGridPoint(self.coeffs, lat=latin, lon=lonin,
                                          lmax=lmax_calc, norm=norm,
                                          csphase=self.csphase)
        elif type(lat) is _np.ndarray:
            values = _np.empty_like(lat, dtype=float)
            for v, latitude, longitude in _np.nditer([values, latin, lonin],
                                                     op_flags=['readwrite']):
                v[...] = _shtools.MakeGridPoint(self.coeffs, lat=latitude,
                                                lon=longitude,
                                                lmax=lmax_calc, norm=norm,
                                                csphase=self.csphase)
            return values
        elif type(lat) is list:
            values = []
            for latitude, longitude in zip(latin, lonin):
                values.append(
                    _shtools.MakeGridPoint(self.coeffs, lat=latitude,
                                           lon=longitude,
                                           lmax=lmax_calc, norm=norm,
                                           csphase=self.csphase))
            return values
        else:
            raise ValueError('lat and lon must be either an int, float, '
                             'ndarray, or list. Input types are {:s} and {:s}.'
                             .format(repr(type(lat)), repr(type(lon))))


# =============== COMPLEX SPHERICAL HARMONICS ================

class SHComplexCoeffs(SHCoeffs):
    """Complex Spherical Harmonics Coefficients class."""

    @staticmethod
    def istype(kind):
        """Check if class has kind 'real' or 'complex'."""
        return kind == 'complex'

    def __init__(self, coeffs, errors=None, error_kind=None,
                 normalization='4pi', csphase=1, units=None, copy=True,
                 header=None, header2=None):
        """Initialize Complex coefficients."""
        lmax = coeffs.shape[1] - 1
        # ---- create mask to filter out m<=l ----
        mask = _np.zeros((2, lmax + 1, lmax + 1), dtype=_np.bool)
        mask[0, 0, 0] = True
        for l in _np.arange(lmax + 1):
            mask[:, l, :l + 1] = True
        mask[1, :, 0] = False

        self.mask = mask
        self.lmax = lmax
        self.kind = 'complex'
        self.normalization = normalization
        self.csphase = csphase
        self.header = header
        self.header2 = header2
        self.units = units
        self.error_kind = error_kind

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

    def _make_real(self, check=True):
        """Convert the complex SHCoeffs class to the real class."""
        # Test if the coefficients correspond to a real grid.
        # This is not very elegant, and the equality condition
        # is probably not robust to round off errors.
        if check:
            for l in self.degrees():
                if self.coeffs[0, l, 0] != self.coeffs[0, l, 0].conjugate():
                    raise RuntimeError('Complex coefficients do not '
                                       'correspond to a real field. '
                                       'l = {:d}, m = 0: {:e}'
                                       .format(l, self.coeffs[0, l, 0]))
                for m in _np.arange(1, l + 1):
                    if m % 2 == 1:
                        if (self.coeffs[0, l, m] != -
                                self.coeffs[1, l, m].conjugate()):
                            raise RuntimeError('Complex coefficients do not '
                                               'correspond to a real field. '
                                               'l = {:d}, m = {:d}: {:e}, {:e}'
                                               .format(
                                                   l, m, self.coeffs[0, l, 0],
                                                   self.coeffs[1, l, 0]))
                    else:
                        if (self.coeffs[0, l, m] !=
                                self.coeffs[1, l, m].conjugate()):
                            raise RuntimeError('Complex coefficients do not '
                                               'correspond to a real field. '
                                               'l = {:d}, m = {:d}: {:e}, {:e}'
                                               .format(
                                                   l, m, self.coeffs[0, l, 0],
                                                   self.coeffs[1, l, 0]))

        coeffs_rc = _np.zeros((2, self.lmax + 1, self.lmax + 1))
        coeffs_rc[0, :, :] = self.coeffs[0, :, :].real
        coeffs_rc[1, :, :] = self.coeffs[0, :, :].imag
        real_coeffs = _shtools.SHctor(coeffs_rc, convention=1,
                                      switchcs=0)
        return SHCoeffs.from_array(real_coeffs,
                                   normalization=self.normalization,
                                   csphase=self.csphase, units=self.units)

    def _rotate(self, angles, dj_matrix):
        """Rotate the coefficients by the Euler angles alpha, beta, gamma."""
        # Note that the current method is EXTREMELY inefficient. The complex
        # coefficients are expanded onto real and imaginary grids, each of
        # the two components are rotated separately as real data, the rotated
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
        elif self.normalization == 'unnorm':
            norm = 3
        elif self.normalization == 'ortho':
            norm = 4
        else:
            raise ValueError(
                "Normalization must be '4pi', 'ortho', 'schmidt', or " +
                "'unnorm'. Input value is {:s}."
                .format(repr(self.normalization)))

        coeffs_rot = _shtools.SHExpandDHC(grid_rot, norm=norm,
                                          csphase=self.csphase)

        return SHCoeffs.from_array(coeffs_rot,
                                   normalization=self.normalization,
                                   csphase=self.csphase, units=self.units,
                                   copy=False)

    def _expandDH(self, sampling, lmax, lmax_calc, extend):
        """Evaluate the coefficients on a Driscoll and Healy (1994) grid."""
        from .shgrid import SHGrid
        if self.normalization == '4pi':
            norm = 1
        elif self.normalization == 'schmidt':
            norm = 2
        elif self.normalization == 'unnorm':
            norm = 3
        elif self.normalization == 'ortho':
            norm = 4
        else:
            raise ValueError(
                "Normalization must be '4pi', 'ortho', 'schmidt', or " +
                "'unnorm'. Input value is {:s}."
                .format(repr(self.normalization)))

        data = _shtools.MakeGridDHC(self.coeffs, sampling=sampling,
                                    norm=norm, csphase=self.csphase, lmax=lmax,
                                    lmax_calc=lmax_calc, extend=extend)
        gridout = SHGrid.from_array(data, grid='DH', units=self.units,
                                    copy=False)
        return gridout

    def _expandGLQ(self, zeros, lmax, lmax_calc, extend):
        """Evaluate the coefficients on a Gauss-Legendre quadrature grid."""
        from .shgrid import SHGrid
        if self.normalization == '4pi':
            norm = 1
        elif self.normalization == 'schmidt':
            norm = 2
        elif self.normalization == 'unnorm':
            norm = 3
        elif self.normalization == 'ortho':
            norm = 4
        else:
            raise ValueError(
                "Normalization must be '4pi', 'ortho', 'schmidt', or " +
                "'unnorm'. Input value is {:s}."
                .format(repr(self.normalization)))

        if zeros is None:
            zeros, weights = _shtools.SHGLQ(self.lmax)

        data = _shtools.MakeGridGLQC(self.coeffs, zeros, norm=norm,
                                     csphase=self.csphase, lmax=lmax,
                                     lmax_calc=lmax_calc, extend=extend)
        gridout = SHGrid.from_array(data, grid='GLQ', units=self.units,
                                    copy=False)
        return gridout

    def _expand_coord(self, lat, lon, lmax_calc, degrees):
        """Evaluate the function at the coordinates lat and lon."""
        if self.normalization == '4pi':
            norm = 1
        elif self.normalization == 'schmidt':
            norm = 2
        elif self.normalization == 'unnorm':
            norm = 3
        elif self.normalization == 'ortho':
            norm = 4
        else:
            raise ValueError(
                "Normalization must be '4pi', 'ortho', 'schmidt', or " +
                "'unnorm'. Input value is {:s}."
                .format(repr(self.normalization)))

        if degrees is True:
            latin = lat
            lonin = lon
        else:
            latin = _np.rad2deg(lat)
            lonin = _np.rad2deg(lon)

        if type(lat) is not type(lon):
            raise ValueError('lat and lon must be of the same type. ' +
                             'Input types are {:s} and {:s}.'
                             .format(repr(type(lat)), repr(type(lon))))

        if type(lat) is int or type(lat) is float or type(lat) is _np.float_:
            return _shtools.MakeGridPointC(self.coeffs, lat=latin, lon=lonin,
                                           lmax=lmax_calc, norm=norm,
                                           csphase=self.csphase)
        elif type(lat) is _np.ndarray:
            values = _np.empty_like(lat, dtype=_np.complex)
            for v, latitude, longitude in _np.nditer([values, latin, lonin],
                                                     op_flags=['readwrite']):
                v[...] = _shtools.MakeGridPointC(self.coeffs, lat=latitude,
                                                 lon=longitude,
                                                 lmax=lmax_calc, norm=norm,
                                                 csphase=self.csphase)
            return values
        elif type(lat) is list:
            values = []
            for latitude, longitude in zip(latin, lonin):
                values.append(
                    _shtools.MakeGridPointC(self.coeffs, lat=latitude,
                                            lon=longitude,
                                            lmax=lmax_calc, norm=norm,
                                            csphase=self.csphase))
            return values
        else:
            raise ValueError('lat and lon must be either an int, float, ' +
                             'ndarray, or list. ' +
                             'Input types are {:s} and {:s}.'
                             .format(repr(type(lat)), repr(type(lon))))
