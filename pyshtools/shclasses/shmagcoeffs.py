"""
    Class for spherical harmonic coefficients of the magnetic potential.
"""
import numpy as _np
import matplotlib as _mpl
import matplotlib.pyplot as _plt
from mpl_toolkits.axes_grid1 import make_axes_locatable as _make_axes_locatable
import copy as _copy
import warnings as _warnings
import xarray as _xr
from scipy.special import factorial as _factorial
import gzip as _gzip
import shutil as _shutil

from .shcoeffs import SHCoeffs as _SHCoeffs
from .shmaggrid import SHMagGrid as _SHMagGrid
from .shtensor import SHMagTensor as _SHMagTensor

from ..spectralanalysis import spectrum as _spectrum
from ..spectralanalysis import cross_spectrum as _cross_spectrum
from ..shio import convert as _convert
from ..shio import shread as _shread
from ..shio import shwrite as _shwrite
from ..shio import read_dov as _read_dov
from ..shio import write_dov as _write_dov
from ..shio import read_bshc as _read_bshc
from ..shio import write_bshc as _write_bshc
from ..shio import read_igrf as _read_igrf
from ..backends.shtools import MakeMagGridDH as _MakeMagGridDH
from ..backends.shtools import MakeMagGradGridDH as _MakeMagGradGridDH
from ..backends.shtools import djpi2 as _djpi2
from ..backends.shtools import MakeMagGridPoint as _MakeMagGridPoint
from ..backends import backend_module
from ..backends import preferred_backend


class SHMagCoeffs(object):
    """
    Spherical harmonic coefficients class for the magnetic potential.

    The coefficients of this class (in units of nT) can be initialized using
    one of the four constructor methods:

        x = SHMagCoeffs.from_array(array, r0)
        x = SHMagCoeffs.from_random(powerspectrum, r0)
        x = SHMagCoeffs.from_zeros(lmax, r0)
        x = SHMagCoeffs.from_file('fname.dat')
        x = SHMagCoeffs.from_netcdf('ncname.nc')

    The normalization convention of the input coefficents is specified
    by the optional normalization and csphase parameters, which take the
    following values:

    normalization : '4pi', geodesy 4-pi normalized.
                  : 'ortho', orthonormalized.
                  : 'schmidt' (default), Schmidt semi-normalized.
                  : 'unnorm', unnormalized.

    csphase       : 1 (default), exlcude the Condon-Shortley phase factor.
                  : -1, include the Condon-Shortley phase factor.

    See the documentation for each constructor method for further options.

    Once initialized, each class instance defines the following class
    attributes:

    lmax          : The maximum spherical harmonic degree of the coefficients.
    coeffs        : The raw coefficients with the specified normalization and
                    csphase conventions. This is a three-dimensional array
                    formatted as coeffs[i, degree, order], where i=0
                    corresponds to positive orders and i=1 to negative orders.
    errors        : The uncertainties of the spherical harmonic coefficients.
    error_kind    : An arbitrary string describing the kind of errors, such as
                    'unknown', 'unspecified', 'calibrated', 'formal' or None.
    r0            : The reference radius of the magnetic potential
                    coefficients.
    normalization : The normalization of the coefficients: '4pi', 'ortho',
                    'schmidt', or 'unnorm'.
    csphase       : Defines whether the Condon-Shortley phase is used (1)
                    or not (-1).
    mask          : A boolean mask that is True for the permissible values of
                    degree l and order m.
    kind          : The coefficient data type (only 'real' is permissible).
    name          : The name of the dataset.
    units         : The units of the spherical harmonic coefficients.
    year          : The year of the time-variable spherical harmonic
                    coefficients.
    header        : A list of values (of type str) from the header line of the
                    input file used to initialize the class (for 'shtools'
                    and 'dov' formatted files only).
    header2       : A list of values (of type str) from the second header line
                    of the input file used to initialize the class (for
                    'shtools' and 'dov' formatted files only).

    Each class instance provides the following methods:

    degrees()             : Return an array listing the spherical harmonic
                            degrees from 0 to lmax.
    spectrum()            : Return the spectrum of the function as a function
                            of spherical harmonic degree.
    correlation()         : Return the spectral correlation with another
                            function.
    set_coeffs()          : Set coefficients in-place to specified values.
    change_ref()          : Return a new class instance referenced to a
                            different reference radius.
    change_units()        : Return a new class instance with the internal
                            coefficients converted to either nT or T.
    rotate()              : Rotate the coordinate system used to express the
                            spherical harmonic coefficients and return a new
                            class instance.
    convert()             : Return a new class instance using a different
                            normalization convention.
    pad()                 : Return a new class instance that is zero padded or
                            truncated to a different lmax.
    expand()              : Calculate the three vector components of the
                            magnetic field, the total field, and the magnetic
                            potential, and return an SHMagGrid class instance.
    tensor()              : Calculate the 9 components of the magnetic field
                            tensor and return an SHMagTensor class instance.
    plot_spectrum()       : Plot the spectrum as a function of spherical
                            harmonic degree.
    plot_spectrum2d()     : Plot the 2D spectrum of all spherical harmonic
                            degrees and orders.
    to_array()            : Return an array of spherical harmonic coefficients
                            with a different normalization convention.
    to_file()             : Save the spherical harmonic coefficients to a file.
    to_netcdf()           : Save raw spherical harmonic coefficients as a
                            netcdf file.
    copy()                : Return a copy of the class instance.
    info()                : Print a summary of the data stored in the
                            SHMagCoeffs instance.
    """

    def __init__(self):
        """Unused constructor of the super class."""
        print('Initialize the class using one of the class methods:\n'
              '>>> pyshtools.SHMagCoeffs.from_array\n'
              '>>> pyshtools.SHMagCoeffs.from_random\n'
              '>>> pyshtools.SHMagCoeffs.from_zeros\n'
              '>>> pyshtools.SHCoeffs.from_file\n'
              '>>> pyshtools.SHMagCoeffs.from_netcdf\n')

    # ---- Factory methods ----
    @classmethod
    def from_array(self, coeffs, r0, errors=None, error_kind=None,
                   normalization='schmidt', csphase=1, lmax=None, name=None,
                   units='nT', year=None, copy=True):
        """
        Initialize the class with spherical harmonic coefficients from an input
        array.

        Usage
        -----
        x = SHMagCoeffs.from_array(array, r0, [errors, error_kind,
                                               normalization, csphase, lmax,
                                               name, units, year, copy])

        Returns
        -------
        x : SHMagCoeffs class instance.

        Parameters
        ----------
        array : ndarray, shape (2, lmaxin+1, lmaxin+1).
            The input spherical harmonic coefficients. This is a three-
            dimensional array formatted as coeffs[i, degree, order], where i=0
            corresponds to positive orders and i=1 to negative orders.
        r0 : float
            The reference radius of the spherical harmonic coefficients.
        errors : ndarray, optional, default = None
            The uncertainties of the spherical harmonic coefficients.
        error_kind : str, optional, default = None
            An arbitrary string describing the kind of errors, such as None,
            'unspecified', 'calibrated' or 'formal'.
        normalization : str, optional, default = 'schmidt'
            '4pi', 'ortho', 'schmidt', or 'unnorm' for geodesy 4pi normalized,
            orthonormalized, Schmidt semi-normalized, or unnormalized
            coefficients, respectively.
        csphase : int, optional, default = 1
            Condon-Shortley phase convention: 1 to exclude the phase factor,
            or -1 to include it.
        lmax : int, optional, default = None
            The maximum spherical harmonic degree to include in the returned
            class instance. This must be less than or equal to lmaxin.
        name : str, optional, default = None
            The name of the dataset.
        units : str, optional, default = 'nT'
            The units of the spherical harmonic coefficients, which can be
            either 'T' or 'nT'.
        year : float, default = None.
            The year of the time-variable spherical harmonic coefficients.
        copy : bool, optional, default = True
            If True, make a copy of array when initializing the class instance.
            If False, initialize the class instance with a reference to array.
        """
        if _np.iscomplexobj(coeffs):
            raise TypeError('The input array must be real.')

        if type(normalization) != str:
            raise ValueError('normalization must be a string. '
                             'Input type is {:s}.'
                             .format(str(type(normalization))))

        if normalization.lower() not in ('4pi', 'ortho', 'schmidt', 'unnorm'):
            raise ValueError(
                "The normalization must be '4pi', 'ortho', 'schmidt', "
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

        if units.lower() not in ('nt', 't'):
            raise ValueError("units can be only 'T' or 'nT'. Input "
                             "value is {:s}.".format(repr(units)))

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
                           "85. Input value is {:d}.".format(lmax),
                           category=RuntimeWarning)
            lmax = 85

        if errors is not None:
            clm = SHMagRealCoeffs(coeffs[:, 0:lmax+1, 0:lmax+1], r0=r0,
                                  errors=errors[:, 0:lmax+1, 0:lmax+1],
                                  error_kind=error_kind,
                                  normalization=normalization.lower(),
                                  csphase=csphase, name=name, units=units,
                                  year=year, copy=copy)
        else:
            clm = SHMagRealCoeffs(coeffs[:, 0:lmax+1, 0:lmax+1], r0=r0,
                                  normalization=normalization.lower(),
                                  csphase=csphase, name=name, units=units,
                                  year=year, copy=copy)
        return clm

    @classmethod
    def from_zeros(self, lmax, r0, errors=None, error_kind=None,
                   normalization='schmidt', csphase=1, name=None, units='nT',
                   year=None):
        """
        Initialize the class with spherical harmonic coefficients set to zero.

        Usage
        -----
        x = SHMagCoeffs.from_zeros(lmax, r0, [errors, error_kind,
                                              normalization, csphase,
                                              name, units, year])

        Returns
        -------
        x : SHMagCoeffs class instance.

        Parameters
        ----------
        lmax : int
            The maximum spherical harmonic degree l of the coefficients.
        r0 : float
            The reference radius of the spherical harmonic coefficients.
        errors : bool, optional, default = None
            If True, initialize the attribute errors with zeros.
        error_kind : str, optional, default = None
            An arbitrary string describing the kind of errors, such as None,
            'unspecified', 'calibrated' or 'formal'.
        normalization : str, optional, default = 'schmidt'
            '4pi', 'ortho', 'schmidt', or 'unnorm' for geodesy 4pi normalized,
            orthonormalized, Schmidt semi-normalized, or unnormalized
             coefficients, respectively.
        csphase : int, optional, default = 1
            Condon-Shortley phase convention: 1 to exclude the phase factor,
            or -1 to include it.
        name : str, optional, default = None
            The name of the dataset.
        units : str, optional, default = 'nT'
            The units of the spherical harmonic coefficients, which can be
            either 'T' or 'nT'.
        year : float, default = None.
            The year of the time-variable spherical harmonic coefficients.
        """
        if normalization.lower() not in ('4pi', 'ortho', 'schmidt', 'unnorm'):
            raise ValueError(
                "The normalization must be '4pi', 'ortho', 'schmidt', "
                "or 'unnorm'. Input value is {:s}."
                .format(repr(normalization))
                )

        if csphase != 1 and csphase != -1:
            raise ValueError(
                "csphase must be either 1 or -1. Input value is {:s}."
                .format(repr(csphase))
                )

        if units.lower() not in ('nt', 't'):
            raise ValueError("units can be only 'T' or 'nT'. Input "
                             "value is {:s}.".format(repr(units)))

        if normalization.lower() == 'unnorm' and lmax > 85:
            _warnings.warn("Calculations using unnormalized coefficients "
                           "are stable only for degrees less than or equal "
                           "to 85. lmax for the coefficients will be set to "
                           "85. Input value is {:d}.".format(lmax),
                           category=RuntimeWarning)
            lmax = 85

        coeffs = _np.zeros((2, lmax + 1, lmax + 1))
        if errors is True:
            error_coeffs = _np.zeros((2, lmax + 1, lmax + 1))
            if error_kind is None:
                error_kind = 'unspecified'
        else:
            error_coeffs = None

        clm = SHMagRealCoeffs(coeffs, r0=r0, errors=error_coeffs,
                              error_kind=error_kind,
                              normalization=normalization.lower(),
                              csphase=csphase, name=name, units=units,
                              year=year)
        return clm

    @classmethod
    def from_file(self, fname, format='shtools', r0=None, lmax=None,
                  normalization='schmidt', skip=0, header=True, header2=False,
                  errors=None, error_kind=None, csphase=1, r0_index=0,
                  header_units='m', file_units='nT', name=None, units='nT',
                  year=None, encoding=None, **kwargs):
        """
        Initialize the class with spherical harmonic coefficients from a file.

        Usage
        -----
        x = SHMagCoeffs.from_file(filename, [format='shtools' or 'dov', r0,
                                  lmax, normalization, csphase, skip, header,
                                  header2, errors, error_kind, r0_index,
                                  header_units, file_units, name, units, year,
                                  encoding])
        x = SHMagCoeffs.from_file(filename, format='igrf', r0, year, [lmax,
                                  normalization, csphase, file_units, name,
                                  units, encoding])
        x = SHMagCoeffs.from_file(filename, format='bshc', r0, [lmax,
                                  normalization, csphase, file_units, name,
                                  units, year])
        x = SHMagCoeffs.from_file(filename, format='npy', r0, [lmax,
                                  normalization, csphase, file_units, name,
                                  units, year, **kwargs])

        Returns
        -------
        x : SHMagCoeffs class instance.

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
            text files, 'igrf' for International Geomagnetic Reference Field
            files, 'bshc' for binary spherical harmonic coefficient files, or
            'npy' for binary numpy files.
        lmax : int, optional, default = None
            The maximum spherical harmonic degree to read from the file. The
            default is to read the entire file.
        header : bool, optional, default = True
            If True, read a list of values from the header line of an 'shtools'
            or 'dov' formatted file. The last header line will contain the
            value for r0.
        header2 : bool, optional, default = False
            If True, read a list of values from a second header line of an
            'shtools' or 'dov' formatted file. The last header line will
            contain the value for r0.
        errors : bool, optional, default = None
            If True, read errors from the file (for 'shtools' and 'dov'
            formatted files only).
        error_kind : str, optional, default = None
            For 'shtools' and 'dov' formatted files: An arbitrary string
            describing the kind of errors, such as None, 'unspecified',
            'calibrated' or 'formal'.
        r0_index : int, optional, default = 0
            For 'shtools' and 'dov' formatted files, r0 will be set using the
            value from the last header line with this index.
        r0 : float, optional, default = None
            The reference radius of the spherical harmonic coefficients.
        header_units : str, optional, default = 'm'
            The units of r0 in the header line of an 'shtools' formatted file:
            'm' or 'km'. If 'km', the value of r0 will be converted to meters.
        file_units : str, optional, default = 'nT'
            The units of the coefficients read from the file: 'nT' or 'T'. If
            these differ from those of the argument units, the spherical
            harmonic coefficients will be converted to those of units.
        name : str, optional, default = None
            The name of the dataset.
        units : str, optional, default = 'nT'
            The units of the spherical harmonic coefficients, which can be
            either 'T' or 'nT'.
        normalization : str, optional, default = 'schmidt'
            '4pi', 'ortho', 'schmidt', or 'unnorm' for geodesy 4pi normalized,
            orthonormalized, Schmidt semi-normalized, or unnormalized
            coefficients, respectively.
        csphase : int, optional, default = 1
            Condon-Shortley phase convention: 1 to exclude the phase factor,
            or -1 to include it.
        skip : int, optional, default = 0
            Number of lines to skip at the beginning of the file for 'shtools'
            formatted files.
        year : float, default = None.
            The year to compute the coefficients for 'igrf' formatted files.
        encoding : str, optional, default = None
            Encoding of the input file when format is 'shtools', 'dov' or
            'igrf'. The default is to use the system default.
        **kwargs : keyword argument list, optional for format = 'npy'
            Keyword arguments of numpy.load() when format is 'npy'.

        Notes
        -----
        Supported file formats:
            'shtools' (see pyshtools.shio.shread)
            'dov' (see pyshtools.shio.read_dov)
            'igrf' (see pyshtools.shio.read_igrf)
            'bshc' (see pyshtools.shio.read_bshc)
            'npy' (see numpy.load)

        The units of the coefficients read from the file can be either 'nT' or
        'T'. If these differ from those of the argument units, the spherical
        harmonic coefficients will be converted.

        For 'shtools', 'dob', 'igrf' and 'bshc' formatted files, if filename
        starts with 'http://', 'https://', or 'ftp://', the file will be
        treated as a URL. In this case, the file will be downloaded in its
        entirety before it is parsed. If the filename ends with '.gz' or
        '.zip', the file will be automatically uncompressed before parsing.
        For zip files, archives with only a single file are supported. Note
        that reading '.gz' and '.zip' files for 'shtools' formatted files will
        be extremely slow if lmax is not specified.

        If format='shtools' or 'dov', spherical harmonic coefficients will be
        read from a text file. The optional parameter `skip` specifies how many
        lines should be skipped before attempting to parse the file, the
        optional parameters `header` and `header2` specify whether to read a
        list of values from one or two header line, and the optional parameter
        `lmax` specifies the maximum degree to read from the file. If headers
        are read, r0_index is used as the indice to set r0 from the last header
        line. If header_unit is specified as 'km', the value of r0 read from
        the header will be converted to meters.
        """
        error_coeffs = None
        header_list = None
        if not header:
            r0_index = None
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

        if format == 'shtools' or format == 'dov':
            if r0_index is not None and r0 is not None:
                raise ValueError('Can not specify both r0_index and r0.')
            if header is False and r0 is None:
                raise ValueError('If header is False, r0 must be specified.')
            if header_units.lower() not in ('m', 'km'):
                raise ValueError("header_units can be only 'm' or 'km'. Input "
                                 "value is {:s}.".format(repr(header_units)))

        if file_units.lower() not in ('nt', 't'):
            raise ValueError("file_units can be only 'T' or 'nT'. Input "
                             "value is {:s}.".format(repr(file_units)))

        if units.lower() not in ('nt', 't'):
            raise ValueError("units can be only 'T' or 'nT'. Input "
                             "value is {:s}.".format(repr(units)))

        if format.lower() == 'shtools' or format.lower() == 'dov':
            if format.lower() == 'shtools':
                read_func = _shread
            else:
                read_func = _read_dov

            if header is True:
                if errors is True:
                    if header2:
                        coeffs, error_coeffs, lmaxout, header_list, \
                            header2_list = read_func(fname, lmax=lmax,
                                                     skip=skip, header=True,
                                                     header2=True, error=True,
                                                     encoding=encoding)
                    else:
                        coeffs, error_coeffs, lmaxout, header_list = read_func(
                            fname, lmax=lmax, skip=skip, header=True,
                            error=True, encoding=encoding)
                else:
                    if header2:
                        coeffs, lmaxout, header_list, header2_list = read_func(
                            fname, lmax=lmax, skip=skip, header=True,
                            header2=True, encoding=encoding)
                    else:
                        coeffs, lmaxout, header_list = read_func(
                            fname, lmax=lmax, skip=skip, header=True,
                            encoding=encoding)

                if r0_index is not None:
                    if header2:
                        r0 = float(header2_list[r0_index])
                    else:
                        r0 = float(header_list[r0_index])
                    if header_units.lower() == 'km':
                        r0 *= 1.e3

            else:
                if errors is True:
                    coeffs, error_coeffs, lmaxout = read_func(
                        fname, lmax=lmax, error=True, skip=skip,
                        encoding=encoding)
                else:
                    coeffs, lmaxout = read_func(fname, lmax=lmax, skip=skip,
                                                encoding=encoding)

            if errors is True and error_kind is None:
                error_kind = 'unspecified'

        elif format.lower() == 'bshc':
            if r0 is None:
                raise ValueError('For binary bshc files, r0 must be '
                                 'specified.')
            coeffs, lmaxout = _read_bshc(fname, lmax=lmax)

        elif format.lower() == 'igrf':
            if r0 is None:
                raise ValueError('For igrf files, r0 must be specified.')
            if year is None:
                raise ValueError('For igrf files, year must be specified.')
            coeffs = _read_igrf(fname, year=year, encoding=encoding)
            lmaxout = coeffs.shape[1] - 1
            if lmax is not None:
                if lmax < lmaxout:
                    coeffs = coeffs[:, :lmax+1, :lmax+1]
                    lmaxout = lmax

        elif format.lower() == 'npy':
            if r0 is None:
                raise ValueError('For binary npy files, r0 must be specified.')
            coeffs = _np.load(fname, **kwargs)
            lmaxout = coeffs.shape[1] - 1
            if lmax is not None:
                if lmax < lmaxout:
                    coeffs = coeffs[:, :lmax+1, :lmax+1]
                    lmaxout = lmax

        else:
            raise NotImplementedError(
                'format={:s} not implemented'.format(repr(format)))

        if _np.iscomplexobj(coeffs):
            raise TypeError('The input coefficients must be real.')

        if normalization.lower() == 'unnorm' and lmaxout > 85:
            _warnings.warn("Calculations using unnormalized coefficients "
                           "are stable only for degrees less than or equal "
                           "to 85. lmax for the coefficients will be set to "
                           "85. Input value is {:d}.".format(lmaxout),
                           category=RuntimeWarning)
            lmaxout = 85
            coeffs = coeffs[:, :lmaxout+1, :lmaxout+1]

        if file_units.lower() == units.lower():
            pass
        elif file_units.lower() == 't' and units.lower() == 'nt':
            coeffs *= 1.e9
            if errors is True:
                error_coeffs *= 1.e9
        elif file_units.lower() == 'nt' and units.lower() == 't':
            coeffs /= 1.e9
            if errors is True:
                error_coeffs /= 1.e9

        clm = SHMagRealCoeffs(coeffs, r0=r0, errors=error_coeffs,
                              error_kind=error_kind,
                              normalization=normalization.lower(),
                              csphase=csphase, header=header_list,
                              header2=header2_list, name=name, units=units,
                              year=year)
        return clm

    @classmethod
    def from_random(self, power, r0, function='total', lmax=None,
                    normalization='schmidt', csphase=1, name=None, units='nT',
                    year=None, exact_power=False, power_unit='per_l'):
        """
        Initialize the class of magnetic potential spherical harmonic
        coefficients as random variables with a given spectrum.

        Usage
        -----
        x = SHMagCoeffs.from_random(power, r0, [function, lmax, normalization,
                                                csphase, name, units, year,
                                                exact_power, power_unit])

        Returns
        -------
        x : SHMagCoeffs class instance.

        Parameters
        ----------
        power : ndarray, shape (L+1)
            numpy array of shape (L+1) that specifies the expected power
            spectrum of the output function, where L is the maximum spherical
            harmonic bandwidth. By default, the input power spectrum represents
            the power of all angular orders of the geoid as a function of
            spherical harmonic degree (see function and power_unit).
        r0 : float
            The reference radius of the spherical harmonic coefficients.
        function : str, optional, default = 'total'
            The type of input power spectrum: 'potential' for the magnetic
            potential in nT m, 'radial' for the radial magnetic field in nT,
            or 'total' for the total magnetic field (Lowes-Mauersberger) in nT.
        lmax : int, optional, default = len(power) - 1
            The maximum spherical harmonic degree l of the output coefficients.
            The coefficients will be set to zero for degrees greater than L.
        normalization : str, optional, default = 'schmidt'
            '4pi', 'ortho', 'schmidt', or 'unnorm' for geodesy 4pi normalized,
            orthonormalized, Schmidt semi-normalized, or unnormalized
            coefficients, respectively.
        csphase : int, optional, default = 1
            Condon-Shortley phase convention: 1 to exclude the phase factor,
            or -1 to include it.
        name : str, optional, default = None
            The name of the dataset.
        units : str, optional, default = 'nT'
            The units of the spherical harmonic coefficients, which can be
            either 'T' or 'nT'.
        year : float, default = None.
            The year of the time-variable spherical harmonic coefficients.
        exact_power : bool, optional, default = False
            If True, the spherical harmonic coefficients of the random
            realization will be rescaled such that the power spectrum is
            exactly equal to the input spectrum.
        power_unit : str, optional, default = 'per_l'
            If 'per_l', the input power spectrum represents the total power of
            all angular orders as a function of spherical harmonic degree. If
            'per_lm', the input power spectrum represents the power per
            coefficient (which is assumed isotropic and varies only as a
            function of spherical harmonic degree).

        Notes
        -----
        This routine returns a random realization of spherical harmonic
        magnetic potential coefficients obtained from a normal distribution.
        The variance of each coefficient is determined by the input power
        spectrum and the type of spectrum (as specified by function and
        power_unit). If power_unit is 'per_l' (default), the variance of each
        coefficient at spherical harmonic degree l is equal to the input total
        power at degree l divided by the number of coefficients at that degree.
        If power_unit is 'per_lm', the variance of each coefficient at degree l
        is equal to the input power at that degree. The power of the input
        function can be either for the potential (default), radial magnetic
        field, or total (Lowes-Mauersberger) magnetic field. The power spectrum
        of the random realization can be fixed exactly to the input spectrum by
        setting exact_power to True.
        """
        if type(normalization) != str:
            raise ValueError('normalization must be a string. '
                             'Input type is {:s}.'
                             .format(str(type(normalization))))

        if function.lower() not in ('potential', 'radial', 'total'):
            raise ValueError(
                "function must be of type 'potential', "
                "'radial', or 'total'. Provided value is {:s}."
                .format(repr(function))
                )

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

        if units.lower() not in ('nt', 't'):
            raise ValueError("units can be only 'T' or 'nT'. Input "
                             "value is {:s}.".format(repr(units)))

        if power_unit.lower() not in ('per_l', 'per_lm'):
            raise ValueError("power_unit must be 'per_l' or 'per_lm'. " +
                             "Input value was {:s}".format(repr(power_unit)))

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
                           "85. Input value is {:d}.".format(nl-1),
                           category=RuntimeWarning)
            nl = 85 + 1
            lmax = 85

        # Create coefficients with unit variance, which returns an expected
        # total power per degree of (2l+1) for 4pi normalized harmonics.
        coeffs = _np.empty((2, nl, nl))
        for l in degrees:
            coeffs[:2, l, :l+1] = _np.random.normal(size=(2, l+1))

        if exact_power:
            power_realization = _spectrum(coeffs, normalization='4pi',
                                          unit=power_unit)
            coeffs *= _np.sqrt(
                power[0:nl] / power_realization)[_np.newaxis, :, _np.newaxis]
        else:
            if power_unit == 'per_l':
                coeffs *= \
                    _np.sqrt(power[0:nl] / (2 * degrees + 1))[_np.newaxis, :,
                                                              _np.newaxis]
            elif power_unit == 'per_lm':
                coeffs *= _np.sqrt(power[0:nl])[_np.newaxis, :, _np.newaxis]

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
            coeffs /= r0
        elif function.lower() == 'radial':
            for l in degrees:
                coeffs[:, l, :l+1] /= (l + 1)
        elif function.lower() == 'total':
            for l in degrees:
                coeffs[:, l, :l+1] /= _np.sqrt((l + 1) * (2 * l + 1))

        if lmax > nl - 1:
            coeffs = _np.pad(coeffs, ((0, 0), (0, lmax - nl + 1),
                             (0, lmax - nl + 1)), 'constant')

        coeffs[0, 0, 0] = 0.0

        clm = SHMagRealCoeffs(coeffs, r0=r0, errors=None,
                              normalization=normalization.lower(),
                              csphase=csphase, name=name, units=units,
                              year=year)
        return clm

    @classmethod
    def from_netcdf(self, filename, lmax=None, normalization='schmidt',
                    csphase=1, name=None, units='nT', year=None):
        """
        Initialize the class with spherical harmonic coefficients from a
        netcdf file.

        Usage
        -----
        x = SHMagCoeffs.from_netcdf(filename, [lmax, normalization, csphase,
                                               name, units, year])

        Returns
        -------
        x : SHMagCoeffs class instance.

        Parameters
        ----------
        filename : str
            Name of the file, including path.
        lmax : int, optional, default = None
            The maximum spherical harmonic degree to read.
        normalization : str, optional, default = 'schmidt'
            Spherical harmonic normalization if not specified in the netcdf
            file: '4pi', 'ortho', 'schmidt', or 'unnorm' for geodesy 4pi
            normalized, orthonormalized, Schmidt semi-normalized, or
            unnormalized coefficients, respectively.
        csphase : int, optional, default = 1
            Condon-Shortley phase convention if not specified in the netcdf
            file: 1 to exclude the phase factor, or -1 to include it.
        name : str, optional, default = None
            The name of the dataset.
        units : str, optional, default = 'nT'
            The units of the spherical harmonic coefficients, which can be
            either 'T' or 'nT'.
        year : float, default = None.
            The year of the time-variable spherical harmonic coefficients.

        Description
        -----------
        The format of the netcdf file has to be exactly as the format that is
        used in SHMagCoeffs.to_netcdf().
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

        if units.lower() not in ('nt', 't'):
            raise ValueError("units can be only 'T' or 'nT'. Input "
                             "value is {:s}.".format(repr(units)))

        try:
            r0 = ds.coeffs.r0
        except:
            raise ValueError("coeffs.r0 must be specified in the netcdf file.")

        try:
            year = ds.coeffs.year
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
            raise ValueError('Magnetic potential coefficients must be '
                             'real. Input coefficients are complex.')

        clm = SHMagRealCoeffs(coeffs, r0=r0, errors=errors,
                              error_kind=error_kind,
                              normalization=normalization.lower(),
                              csphase=csphase, name=name, units=units,
                              year=year)
        return clm

    # ---- Define methods that modify internal variables ----
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

        mneg_mask = (ms < 0).astype(_np.int_)
        self.coeffs[mneg_mask, ls, _np.abs(ms)] = values

    # ---- IO routines ----
    def to_file(self, filename, format='shtools', header=None, errors=True,
                lmax=None, encoding=None, **kwargs):
        """
        Save spherical harmonic coefficients to a file.

        Usage
        -----
        x.to_file(filename, [format='shtools', header, errors, encoding])
        x.to_file(filename, format='dov', [header, errors, encoding])
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
            A header string written to an 'shtools' or 'dov'-formatted file
            directly before the metadata and spherical harmonic coefficients.
        errors : bool, optional, default = True
            If True, save the errors in the file (for 'shtools' and 'dov'
            formatted files only).
        lmax : int, optional, default = self.lmax
            The maximum spherical harmonic degree to write to the file.
        encoding : str, optional, default = None
            Encoding of the output file when format is 'shtools' or 'dov'. The
            default is to use the system default.
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

        'shtools': The coefficients and meta-data will be written to an ascii
        formatted file. The first line is an optional user provided header
        line, and the following line provides the attributes r0, and lmax. The
        spherical harmonic coefficients (and optionally the errors) are then
        listed, with increasing degree and order, with the format

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
        coefficients (with all orders being listed, one degree at a time). This
        format does noe support additional metadata or coefficient errors.

        'npy': The spherical harmonic coefficients (but not the meta-data nor
        errors) will be saved to a binary numpy 'npy' file.
        """
        if lmax is None:
            lmax = self.lmax

        if errors is True and self.errors is None:
            errors = False

        if filename[-3:] == '.gz':
            filebase = filename[:-3]
        else:
            filebase = filename

        if format.lower() == 'shtools' or format.lower() == 'dov':
            if format.lower() == 'shtools':
                write_func = _shwrite
            else:
                write_func = _write_dov

            if errors is True and self.errors is None:
                raise ValueError('Can not save errors when then have not been '
                                 'initialized.')

            header_str = '{:.16e}, {:d}'.format(self.r0, lmax)
            if header is None:
                header = header_str
                header2 = None
            else:
                header2 = header_str

            if errors:
                write_func(filebase, self.coeffs, errors=self.errors,
                           header=header, header2=header2, lmax=lmax,
                           encoding=encoding)
            else:
                write_func(filebase, self.coeffs, errors=None,
                           header=header, header2=header2, lmax=lmax,
                           encoding=encoding)

        elif format.lower() == 'bshc':
            _write_bshc(filebase, self.coeffs, lmax=lmax)

        elif format == 'npy':
            _np.save(filename, self.coeffs, **kwargs)
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
        ds['coeffs'].attrs['r0'] = self.r0
        if self.units is not None:
            ds['coeffs'].attrs['units'] = self.units
        if self.year is not None:
            ds['coeffs'].attrs['year'] = self.year

        if self.errors is not None:
            cerrors = self.errors[0, :lmax+1, :lmax+1]
            serrors = _np.transpose(self.errors[1, :lmax+1, :lmax+1])
            serrors = _np.vstack([serrors[1:], serrors[0]])
            ds['errors'] = (('degree', 'order'), cerrors + serrors)
            ds['errors'].attrs['normalization'] = self.normalization
            ds['errors'].attrs['csphase'] = self.csphase
            ds['errors'].attrs['r0'] = self.r0
            if self.units is not None:
                ds['errors'].attrs['units'] = self.units
            if self.year is not None:
                ds['errors'].attrs['year'] = self.year
            if self.error_kind is not None:
                ds['errors'].attrs['error_kind'] = self.error_kind

        ds.to_netcdf(filename)

    def to_array(self, normalization=None, csphase=None, lmax=None,
                 errors=False):
        """
        Return spherical harmonic coefficients (and errors) as a numpy array.

        Usage
        -----
        coeffs, [errors] = x.to_array([normalization, csphase, lmax, errors])

        Returns
        -------
        coeffs : ndarry, shape (2, lmax+1, lmax+1)
            numpy ndarray of the spherical harmonic coefficients. This is a
            three-dimensional array formatted as coeffs[i, degree, order],
            where i=0 corresponds to positive orders and i=1 to negative
            orders.
        errors : ndarry, shape (2, lmax+1, lmax+1)
            numpy ndarray of the errors of the spherical harmonic coefficients
            if they are not None and errors is True.

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
        errors : bool, optional, default = False
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
        zero padded. If the errors of the coefficients are set, and the
        optional parameter errors is set to True, the errors will be output as
        a separate array.
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
        Print a summary of the data stored in the SHMagCoeffs class instance.

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
        Add two similar sets of magnetic potential coefficients:
        self + other.
        """
        if isinstance(other, SHMagCoeffs):
            if (self.r0 == other.r0 and
                    self.normalization == other.normalization and
                    self.csphase == other.csphase and self.kind == other.kind
                    and self.lmax == other.lmax):
                coeffs = _np.empty([2, self.lmax+1, self.lmax+1],
                                   dtype=self.coeffs.dtype)
                coeffs[self.mask] = (self.coeffs[self.mask] +
                                     other.coeffs[self.mask])
                return SHMagCoeffs.from_array(
                    coeffs, r0=self.r0, csphase=self.csphase,
                    normalization=self.normalization)
            else:
                raise ValueError('Addition is permitted only when the two '
                                 'SHMagCoeffs instances have the same kind, '
                                 'normalization, csphase, r0 and lmax.')
        else:
            raise TypeError('Addition is permitted only for two SHMagCoeffs '
                            'instances. Type of other is {:s}.'
                            .format(repr(type(other))))

    def __radd__(self, other):
        """
        Add two similar sets of magnetic potential coefficients:
        other + self.
        """
        return self.__add__(other)

    def __sub__(self, other):
        """
        Subtract two similar sets of magnetic potential coefficients:
        self - other.
        """
        if isinstance(other, SHMagCoeffs):
            if (self.r0 == other.r0 and
                    self.normalization == other.normalization and
                    self.csphase == other.csphase and self.kind == other.kind
                    and self.lmax == other.lmax):
                coeffs = _np.empty([2, self.lmax+1, self.lmax+1],
                                   dtype=self.coeffs.dtype)
                coeffs[self.mask] = (self.coeffs[self.mask] -
                                     other.coeffs[self.mask])
                return SHMagCoeffs.from_array(
                    coeffs, r0=self.r0, csphase=self.csphase,
                    normalization=self.normalization)
            else:
                raise ValueError('Subtraction is permitted only when the two '
                                 'SHMagCoeffs instances have the same kind, '
                                 'normalization, csphase, r0 and lmax.')
        else:
            raise TypeError('Subtraction is permitted only for two '
                            'SHMagCoeffs instances. Type of other is {:s}.'
                            .format(repr(type(other))))

    def __rsub__(self, other):
        """
        Subtract two similar sets of magnetic potential coefficients:
        other - self.
        """
        if isinstance(other, SHMagCoeffs):
            if (self.r0 == other.r0 and
                    self.normalization == other.normalization and
                    self.csphase == other.csphase and self.kind == other.kind
                    and self.lmax == other.lmax):
                coeffs = _np.empty([2, self.lmax+1, self.lmax+1],
                                   dtype=self.coeffs.dtype)
                coeffs[self.mask] = (other.coeffs[self.mask] -
                                     self.coeffs[self.mask])
                return SHMagCoeffs.from_array(
                    coeffs, r0=self.r0, csphase=self.csphase,
                    normalization=self.normalization)
            else:
                raise ValueError('Subtraction is permitted only when the two '
                                 'SHMagCoeffs instances have the same kind, '
                                 'normalization, csphase, r0 and lmax.')
        else:
            raise TypeError('Subtraction is permitted only for two '
                            'SHMagCoeffs instances. Type of other is {:s}.'
                            .format(repr(type(other))))

    def __mul__(self, other):
        """
        Multiply an SHMagCoeffs instance by an SHCoeffs instance or scalar:
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
                return SHMagCoeffs.from_array(
                    coeffs, r0=self.r0, csphase=self.csphase,
                    normalization=self.normalization)
            else:
                raise ValueError('The two sets of coefficients must have the '
                                 'same kind, normalization, csphase, and '
                                 'lmax.')
        elif _np.isscalar(other) is True:
            if self.kind == 'real' and _np.iscomplexobj(other):
                raise ValueError('Can not multiply real magnetic '
                                 'potential coefficients by a complex '
                                 'constant.')
            coeffs = _np.empty([2, self.lmax+1, self.lmax+1],
                               dtype=self.coeffs.dtype)
            coeffs[self.mask] = self.coeffs[self.mask] * other
            return SHMagCoeffs.from_array(
                coeffs, r0=self.r0, csphase=self.csphase,
                normalization=self.normalization, units=self.units)
        else:
            raise TypeError('Multiplication of an SHMagCoeffs instance is '
                            'permitted only with either an SHCoeffs instance '
                            'or a scalar. '
                            'Type of other is {:s}.'.format(repr(type(other))))

    def __rmul__(self, other):
        """
        Multiply an SHMagCoeffs instance by an SHCoeffs instance or scalar:
        other * self.
        """
        return self.__mul__(other)

    def __truediv__(self, other):
        """
        Divide an SHMagCoeffs instance by an SHCoeffs instance or scalar:
        self / other.
        """
        if isinstance(other, _SHCoeffs):
            if (self.normalization == other.normalization and
                    self.csphase == other.csphase and self.kind == other.kind
                    and self.lmax == other.lmax):
                coeffs = _np.empty([2, self.lmax+1, self.lmax+1],
                                   dtype=self.coeffs.dtype)
                coeffs[self.mask] = (self.coeffs[self.mask] /
                                     other.coeffs[self.mask])
                return SHMagCoeffs.from_array(
                    coeffs, r0=self.r0, csphase=self.csphase,
                    normalization=self.normalization)
            else:
                raise ValueError('The two sets of coefficients must have the '
                                 'same kind, normalization, csphase, and '
                                 'lmax.')
        elif _np.isscalar(other) is True:
            if self.kind == 'real' and _np.iscomplexobj(other):
                raise ValueError('Can not divide real magnetic '
                                 'potential coefficients by a complex '
                                 'constant.')
            coeffs = _np.empty([2, self.lmax+1, self.lmax+1],
                               dtype=self.coeffs.dtype)
            coeffs[self.mask] = self.coeffs[self.mask] / other
            return SHMagCoeffs.from_array(
                coeffs, r0=self.r0, csphase=self.csphase,
                normalization=self.normalization, units=self.units)
        else:
            raise TypeError('Division of an SHMagCoeffs instance is '
                            'permitted only with either an SHCoeffs instance '
                            'or a scalar. '
                            'Type of other is {:s}.'.format(repr(type(other))))

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

    def spectrum(self, function='total', lmax=None, unit='per_l', base=10.):
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
        function : str, optional, default = 'total'
            The type of power spectrum to return: 'potential' for the
            magnetic potential in nT m, 'radial' for the radial magnetic field
            in nT, or 'total' for the total magnetic field (Lowes-Mauersberger)
            in nT.
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

        Notes
        -----
        This method returns the power spectrum of the class instance, where the
        type of function is defined by the function parameter: 'potential' for
        the magnetic potential, 'radial' for the radial magnetic field, or
        'total' for the total magnetic field. The case 'total' corresponds to
        the classic Lowes-Mauersberger spectrum. In all cases, total power of
        the function is defined as the integral of the function squared over
        all space, divided by the area the function spans. If the mean of the
        function is zero, this is equivalent to the variance of the function.

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
        if function.lower() not in ('potential', 'radial', 'total'):
            raise ValueError(
                "function must be of type 'potential', 'radial', or 'total'. "
                "Provided value is {:s}."
                .format(repr(function))
                )

        s = _spectrum(self.coeffs, normalization=self.normalization,
                      convention='power', unit=unit, base=base, lmax=lmax)

        if self.errors is not None:
            es = _spectrum(self.errors, normalization=self.normalization,
                           convention='power', unit=unit, base=base, lmax=lmax)

        if function.lower() == 'potential':
            s *= self.r0**2
            if self.errors is not None:
                es *= self.r0**2
        elif function.lower() == 'radial':
            degrees = _np.arange(len(s))
            s *= (degrees + 1)**2
            if self.errors is not None:
                es *= (degrees + 1)**2
        elif function.lower() == 'total':
            degrees = _np.arange(len(s))
            s *= (degrees + 1) * (2 * degrees + 1)
            if self.errors is not None:
                es *= (degrees + 1) * (2 * degrees + 1)

        if self.errors is not None:
            return s, es
        else:
            return s

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
        hlm : SHCoeffs, SHMagCoeffs or SHGravCoeffs class instance.
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
        if not isinstance(hlm, (_SHCoeffs, SHMagCoeffs, _SHGravCoeffs)):
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

    # ---- Operations that return a new SHMagCoeffs class instance ----
    def rotate(self, alpha, beta, gamma, degrees=True, convention='y',
               body=False, dj_matrix=None, backend=None, nthreads=None):
        """
        Rotate either the coordinate system used to express the spherical
        harmonic coefficients or the physical body, and return a new class
        instance.

        Usage
        -----
        x_rotated = x.rotate(alpha, beta, gamma, [degrees, convention,
                             body, dj_matrix, backend, nthreads])

        Returns
        -------
        x_rotated : SHMagCoeffs class instance

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
            The djpi2 rotation matrix computed by a call to djpi2 (not used if
            the backend is 'ducc').
        backend : str, optional, default = preferred_backend()
            Name of the preferred backend, either 'shtools' or 'ducc'.
        nthreads : int, optional, default = 1
            Number of threads to use for the 'ducc' backend. Setting this
            parameter to 0 will use as many threads as there are hardware
            threads on the system.

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
            raise ValueError('convention must be a string. '
                             'Input type is {:s}.'
                             .format(str(type(convention))))

        if convention.lower() not in ('x', 'y'):
            raise ValueError(
                "convention must be either 'x' or 'y'. "
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
        if backend is None:
            backend = preferred_backend()

        if self.lmax > 1200:
            _warnings.warn("The rotate() method is accurate only to about"
                           " spherical harmonic degree 1200. "
                           "lmax = {:d}.".format(self.lmax),
                           category=RuntimeWarning)

        rot = self._rotate(angles, dj_matrix, r0=self.r0, backend=backend,
                           nthreads=nthreads)
        return rot

    def convert(self, normalization=None, csphase=None, lmax=None):
        """
        Return an SHMagCoeffs class instance with a different normalization
        convention.

        Usage
        -----
        clm = x.convert([normalization, csphase, lmax])

        Returns
        -------
        clm : SHMagCoeffs class instance

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

        Notes
        -----
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
                             'Input type is {:s}.'
                             .format(str(type(normalization))))
        if normalization.lower() not in ('4pi', 'ortho', 'schmidt', 'unnorm'):
            raise ValueError(
                "normalization must be '4pi', 'ortho', 'schmidt', or "
                "'unnorm'. Provided value is {:s}."
                .format(repr(normalization)))
        if csphase != 1 and csphase != -1:
            raise ValueError(
                "csphase must be 1 or -1. Input value is {:s}."
                .format(repr(csphase)))

        if self.errors is not None:
            coeffs, errors = self.to_array(normalization=normalization.lower(),
                                           csphase=csphase, lmax=lmax,
                                           errors=True)
            return SHMagCoeffs.from_array(
                coeffs, r0=self.r0, errors=errors, error_kind=self.error_kind,
                normalization=normalization.lower(),
                csphase=csphase, units=self.units, year=self.year, copy=False)
        else:
            coeffs = self.to_array(normalization=normalization.lower(),
                                   csphase=csphase, lmax=lmax)
            return SHMagCoeffs.from_array(
                coeffs, r0=self.r0, normalization=normalization.lower(),
                csphase=csphase, units=self.units, year=self.year, copy=False)

    def pad(self, lmax, copy=True):
        """
        Return an SHMagCoeffs class where the coefficients are zero padded or
        truncated to a different lmax.

        Usage
        -----
        clm = x.pad(lmax)

        Returns
        -------
        clm : SHMagCoeffs class instance

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
            mask = _np.zeros((2, lmax + 1, lmax + 1), dtype=bool)
            for l in _np.arange(lmax + 1):
                mask[:, l, :l + 1] = True
            mask[1, :, 0] = False
            clm.mask = mask

        clm.lmax = lmax
        return clm

    def change_ref(self, r0=None, lmax=None):
        """
        Return a new SHMagCoeffs class instance with a different reference r0.

        Usage
        -----
        clm = x.change_ref([r0, lmax])

        Returns
        -------
        clm : SHMagCoeffs class instance.

        Parameters
        ----------
        r0 : float, optional, default = self.r0
            The reference radius of the spherical harmonic coefficients.
        lmax : int, optional, default = self.lmax
            Maximum spherical harmonic degree to output.

        Notes
        -----
        This method returns a new class instance of the magnetic potential,
        but using a difference reference r0. When changing the reference
        radius r0, the spherical harmonic coefficients will be upward or
        downward continued under the assumption that the reference radius is
        exterior to the body.
        """
        if lmax is None:
            lmax = self.lmax

        clm = self.pad(lmax)

        if r0 is not None and r0 != self.r0:
            for l in _np.arange(lmax+1):
                clm.coeffs[:, l, :l+1] *= (self.r0 / r0)**(l+2)
                if self.errors is not None:
                    clm.errors[:, l, :l+1] *= (self.r0 / r0)**(l+2)
            clm.r0 = r0

        return clm

    def change_units(self, units=None, lmax=None):
        """
        Return a new SHMagCoeffs class instance with the internal coefficients
        converted to either nT or T.

        Usage
        -----
        clm = x.change_units([units, lmax])

        Returns
        -------
        clm : SHMagCoeffs class instance.

        Parameters
        ----------
        units : str, optional, default = self.units
            The new units of the spherical harmonic coefficients, which can be
            either 'T' or 'nT'.
        lmax : int, optional, default = self.lmax
            Maximum spherical harmonic degree to output.
        """
        if lmax is None:
            lmax = self.lmax
        if units is None:
            units = self.units

        clm = self.pad(lmax)

        if units != self.units:
            if units.lower() == 't' and self.units.lower() == 'nt':
                factor = 1.e-9
            elif units.lower() == 'nt' and self.units.lower() == 't':
                factor = 1.e9
            else:
                raise ValueError("units must be either 'nT' or 'T'. "
                                 "Input value is {:s}".format(repr(units)))

            clm.coeffs *= factor
            clm.units = units
            if clm.errors is not None:
                clm.errors *= factor

        return clm

    # ---- Routines that return different magnetic-related class instances ----
    def expand(self, a=None, f=None, lat=None, colat=None, lon=None,
               degrees=True, lmax=None, lmax_calc=None, sampling=2,
               extend=True):
        """
        Create 2D cylindrical maps on a flattened ellipsoid of all three
        components of the magnetic field, the total magnetic intensity,
        and the magnetic potential. Alternatively, compute the magnetic field
        vector at specified coordinates.

        Usage
        -----
        grids = x.expand([a, f, lmax, lmax_calc, sampling, extend])
        g = x.expand(lat, lon, [a, f, lmax, lmax_calc, degrees])
        g = x.expand(colat, lon, [a, f, lmax, lmax_calc, degrees])

        Returns
        -------
        grids : SHMagGrid class instance.
        g     : (r, theta, phi) components of the magnetic field at the
                specified points.

        Parameters
        ----------
        a : optional, float, default = self.r0
            The semi-major axis of the flattened ellipsoid on which the field
            is computed.
        f : optional, float, default = 0
            The flattening of the reference ellipsoid: f=(a-b)/a.
        lat : int, float, ndarray, or list, optional, default = None
            Latitude coordinates where the magnetic field is to be evaluated.
        colat : int, float, ndarray, or list, optional, default = None
            Colatitude coordinates where the magnetic field is to be evaluated.
        lon : int, float, ndarray, or list, optional, default = None
            Longitude coordinates where the magnetic field is to be evaluated.
        degrees : bool, optional, default = True
            True if lat, colat and lon are in degrees, False if in radians.
        lmax : optional, integer, default = self.lmax
            The maximum spherical harmonic degree, which determines the number
            of samples of the output grids, n=2lmax+2, and the latitudinal
            sampling interval, 90/(lmax+1).
        lmax_calc : optional, integer, default = lmax
            The maximum spherical harmonic degree used in evaluating the
            functions. This must be less than or equal to lmax.
        sampling : optional, integer, default = 2
            If 1 the output grids are equally sampled (n by n). If 2 (default),
            the grids are equally spaced in degrees.
        extend : bool, optional, default = True
            If True, compute the longitudinal band for 360 E and the
            latitudinal band for 90 S.

        Notes
        -----
        This method will create 2-dimensional cylindrical maps of the three
        components of the magnetic field, the total field, and the magnetic
        potential, and return these as an SHMagGrid class instance. Each
        map is stored as an SHGrid class instance using Driscoll and Healy
        grids that are either equally sampled (n by n) or equally spaced
        in degreess latitude and longitude. All grids use geocentric
        coordinates, and the units are either in nT (for the magnetic field),
        or nT m (for the potential). If latitude and longitude coordinates
        are specified, this method will instead return only the magnetic field
        vector.

        The magnetic potential is given by

        V = r0 Sum_{l=1}^lmax (r0/r)^{l+1} Sum_{m=-l}^l g_{lm} Y_{lm}

        and the magnetic field is

        B = - Grad V.

        The coefficients are referenced to a radius r0, and the function is
        computed on a flattened ellipsoid with semi-major axis a (i.e., the
        mean equatorial radius) and flattening f.
        """
        if lat is not None and colat is not None:
            raise ValueError('lat and colat can not both be specified.')

        if a is None:
            a = self.r0
        if f is None:
            f = 0.

        if (lat is not None or colat is not None) and lon is not None:
            if lmax_calc is None:
                lmax_calc = self.lmax

            if colat is not None:
                if degrees:
                    temp = 90.
                else:
                    temp = _np.pi/2.

                if type(colat) is list:
                    lat = list(map(lambda x: temp - x, colat))
                else:
                    lat = temp - colat

            values = self._expand_coord(a=a, f=f, lat=lat, lon=lon,
                                        degrees=degrees, lmax_calc=lmax_calc)
            return values

        else:
            if lmax is None:
                lmax = self.lmax
            if lmax_calc is None:
                lmax_calc = lmax

            coeffs = self.to_array(normalization='schmidt', csphase=1,
                                   errors=False)

            rad, theta, phi, total, pot = _MakeMagGridDH(
                coeffs, self.r0, a=a, f=f, lmax=lmax, lmax_calc=lmax_calc,
                sampling=sampling, extend=extend)

        return _SHMagGrid(rad, theta, phi, total, pot, a, f, lmax, lmax_calc,
                          units=self.units, pot_units=self.units+' m',
                          year=self.year)

    def tensor(self, a=None, f=None, lmax=None, lmax_calc=None, sampling=2,
               extend=True):
        """
        Create 2D cylindrical maps on a flattened ellipsoid of the 9
        components of the magnetic field tensor in a local north-oriented
        reference frame, and return an SHMagTensor class instance.

        Usage
        -----
        tensor = x.tensor([a, f, lmax, lmax_calc, sampling, extend])

        Returns
        -------
        tensor : SHMagTensor class instance.

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
        sampling : optional, integer, default = 2
            If 1 the output grids are equally sampled (n by n). If 2 (default),
            the grids are equally spaced in degrees.
        extend : bool, optional, default = True
            If True, compute the longitudinal band for 360 E and the
            latitudinal band for 90 S.

        Notes
        -----
        This method will create 2-dimensional cylindrical maps for the 9
        components of the magnetic field tensor and return an SHMagTensor
        class instance. The components are

            (Vxx, Vxy, Vxz)
            (Vyx, Vyy, Vyz)
            (Vzx, Vzy, Vzz)

        where the reference frame is north-oriented, where x points north, y
        points west, and z points upward (all tangent or perpendicular to a
        sphere of radius r, where r is the local radius of the flattened
        ellipsoid). The magnetic potential is defined as

            V = r0 Sum_{l=0}^lmax (r0/r)^(l+1) Sum_{m=-l}^l C_{lm} Y_{lm},

        where r0 is the reference radius of the spherical harmonic coefficients
        Clm, and the vector magnetic field is

            B = - Grad V.

        The components of the tensor are calculated according to eq. 1 in
        Petrovskaya and Vershkov (2006), which is based on eq. 3.28 in Reed
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
        units of nT / m.

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

        coeffs = self.to_array(normalization='schmidt', csphase=1,
                               errors=False)

        if self.units.lower() == 'nt':
            units = 'nT/m'
        else:
            units = 'T/m'
        vxx, vyy, vzz, vxy, vxz, vyz = _MakeMagGradGridDH(
            coeffs, self.r0, a=a, f=f, lmax=lmax, lmax_calc=lmax_calc,
            sampling=sampling, extend=extend)

        return _SHMagTensor(vxx, vyy, vzz, vxy, vxz, vyz, a, f, lmax,
                            lmax_calc, units=units, year=self.year)

    # ---- Plotting routines ----
    def plot_spectrum(self, function='total', unit='per_l', base=10.,
                      lmax=None, xscale='lin', yscale='log', grid=True,
                      legend=None, legend_error='error', legend_loc='best',
                      axes_labelsize=None, tick_labelsize=None, ax=None,
                      show=True, fname=None, **kwargs):
        """
        Plot the spectrum as a function of spherical harmonic degree.

        Usage
        -----
        x.plot_spectrum([function, unit, base, lmax, xscale, yscale, grid,
                         legend, legend_loc, axes_labelsize, tick_labelsize,
                         ax, show, fname, **kwargs])

        Parameters
        ----------
        function : str, optional, default = 'total'
            The type of power spectrum to calculate: 'potential' for the
            magnetic potential, 'radial' for the radial magnetic field, or
            'total' for the total magnetic field (Lowes-Mauersberger).
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
        ax : matplotlib axes object, optional, default = None
            A single matplotlib axes object where the plot will appear.
        show : bool, optional, default = True
            If True, plot to the screen.
        fname : str, optional, default = None
            If present, and if ax is not specified, save the image to the
            specified file.
        **kwargs : keyword arguments, optional
            Keyword arguments for pyplot.plot().

        Notes
        -----
        This method plots the power (and error) spectrum of the class instance,
        where the type of function is defined by the function parameter:
        'potential' for the magnetic potential, 'radial' for the radial
        magnetic field, or 'total' for the total magnetic field. The case
        'total' corresponds to the classic Lowes-Mauersberger spectrum. In all
        cases, total power of the function is defined as the integral of the
        function squared over all space, divided by the area the function
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
        else:
            spectrum = self.spectrum(function=function, unit=unit, base=base,
                                     lmax=lmax)

        ls = _np.arange(lmax + 1)

        if ax is None:
            fig, axes = _plt.subplots(1, 1)
        else:
            axes = ax

        if axes_labelsize is None:
            axes_labelsize = _mpl.rcParams['axes.labelsize']
            if type(axes_labelsize) == str:
                axes_labelsize = _mpl.font_manager \
                                 .FontProperties(size=axes_labelsize) \
                                 .get_size_in_points()
        if tick_labelsize is None:
            tick_labelsize = _mpl.rcParams['xtick.labelsize']
            if type(tick_labelsize) == str:
                tick_labelsize = _mpl.font_manager \
                                 .FontProperties(size=tick_labelsize) \
                                 .get_size_in_points()

        axes.set_xlabel('Spherical harmonic degree', fontsize=axes_labelsize)

        if function == 'potential':
            axes.set_ylabel('Power, ' + self.units + '$^2$ m$^2$',
                            fontsize=axes_labelsize)
        elif function == 'radial':
            axes.set_ylabel('Power, ' + self.units + '$^2$',
                            fontsize=axes_labelsize)
        elif function == 'total':
            axes.set_ylabel('Power, ' + self.units + '$^2$',
                            fontsize=axes_labelsize)

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

        if self.errors is not None:
            axes.plot(ls[1:lmax + 1], spectrum[1:lmax + 1], label=legend,
                      **kwargs)
            axes.plot(ls[1:lmax + 1], error_spectrum[1:lmax + 1],
                      label=legend_error, **kwargs)
        else:
            axes.plot(ls[1:lmax + 1], spectrum[1: lmax + 1], label=legend,
                      **kwargs)

        if xscale == 'lin':
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

    def plot_spectrum2d(self, function='total', ticks='WSen',
                        tick_interval=[None, None],
                        minor_tick_interval=[None, None],
                        degree_label='Spherical harmonic degree',
                        order_label='Spherical harmonic order', title=None,
                        colorbar='right', origin='top', cmap='viridis',
                        cmap_limits=None, cmap_rlimits=None,
                        cmap_reverse=False, cmap_scale='log',
                        cb_triangles='neither', cb_label=None, cb_offset=None,
                        cb_width=None, lmax=None, errors=False, xscale='lin',
                        yscale='lin', grid=False, titlesize=None,
                        axes_labelsize=None, tick_labelsize=None, ax=None,
                        show=True, fname=None):
        """
        Plot the spectrum as a function of spherical harmonic degree and order.

        Usage
        -----
        x.plot_spectrum2d([function, ticks, degree_label, order_label, title,
                           colorbar, origin, cmap, cmap_limits, cmap_rlimits,
                           cmap_reverse, cmap_scale, cb_triangles, cb_label,
                           cb_offset, cb_width, lmax, errors, xscale, yscale,
                           grid, titlesize, axes_labelsize, tick_labelsize, ax,
                           show, fname])

        Parameters
        ----------
        function : str, optional, default = 'geoid'
            The type of power spectrum to calculate: 'potential' for the
            magnetic potential, 'radial' for the radial magnetic field, or
            'total' for the total magnetic field (Lowes-Mauersberger).
        ticks : str, optional, default = 'WSen'
            Specify which axes should have ticks drawn and annotated. Capital
            letters plot the ticks and annotations, whereas small letters plot
            only the ticks. 'W', 'S', 'E', and 'N' denote the west, south, east
            and north boundaries of the plot, respectively. Alternatively, use
            'L', 'B', 'R', and 'T' for left, bottom, right, and top.
        tick_interval : list or tuple, optional, default = [None, None]
            Intervals to use when plotting the degree and order ticks,
            respectively (used only when xscale and yscale are 'lin'). If set
            to None, ticks will be generated automatically.
        minor_tick_interval : list or tuple, optional, default = [None, None]
            Intervals to use when plotting the minor degree and order ticks,
            respectively (used only when xscale and yscale are 'lin'). If set
            to None, minor ticks will be generated automatically.
        degree_label : str, optional, default = 'Spherical harmonic degree'
            Label for the spherical harmonic degree axis.
        order_label : str, optional, default = 'Spherical harmonic order'
            Label for the spherical harmonic order axis.
        title : str or list, optional, default = None
            The title of the plot.
        colorbar : str, optional, default = 'right'
            Plot a colorbar along the 'top', 'right', 'bottom', or 'left' axis.
        origin : str, optional, default = 'top'
            Location where the degree 0 coefficient is plotted. Either 'left',
            'right', 'top', or 'bottom'.
        cmap : str, optional, default = 'viridis'
            The color map to use when plotting the data and colorbar.
        cmap_limits : list, optional, default = [self.min(), self.max()]
            Set the lower and upper limits of the data used by the colormap,
            and optionally an interval for each color band. If interval is
            specified, the number of discrete colors will be
            (cmap_limits[1]-cmap_limits[0])/cmap_limits[2] for linear scales
            and log10(cmap_limits[1]/cmap_limits[0])*cmap_limits[2] for
            logarithmic scales.
        cmap_rlimits : list, optional, default = None
           Same as cmap_limits, except the provided upper and lower values are
           relative with respect to the maximum value of the data.
        cmap_reverse : bool, optional, default = False
            Set to True to reverse the sense of the color progression in the
            color table.
        cmap_scale : str, optional, default = 'log'
            Scale of the color axis: 'lin' for linear or 'log' for logarithmic.
        cb_triangles : str, optional, default = 'neither'
            Add triangles to the edges of the colorbar for minimum and maximum
            values. Can be 'neither', 'both', 'min', or 'max'.
        cb_label : str, optional, default = None
            Text label for the colorbar.
        cb_offset : float or int, optional, default = None
            Offset of the colorbar from the map edge in points. If None,
            the offset will be calculated automatically.
        cb_width : float, optional, default = None
            Width of the colorbar in percent with respect to the width of the
            respective image axis. Defaults are 2.5 and 5 for vertical and
            horizontal colorbars, respectively.
        lmax : int, optional, default = self.lmax
            The maximum spherical harmonic degree to plot.
        errors : bool, optional, default = False
            If True, plot the spectrum of the errors.
        xscale : str, optional, default = 'lin'
            Scale of the l axis: 'lin' for linear or 'log' for logarithmic.
        yscale : str, optional, default = 'lin'
            Scale of the m axis: 'lin' for linear or 'log' for logarithmic.
        grid : bool, optional, default = False
            If True, plot grid lines.
        axes_labelsize : int, optional, default = None
            The font size for the x and y axes labels.
        tick_labelsize : int, optional, default = None
            The font size for the x and y tick labels.
        titlesize : int, optional, default = None
            The font size of the title.
        ax : matplotlib axes object, optional, default = None
            A single matplotlib axes object where the plot will appear.
        show : bool, optional, default = True
            If True, plot to the screen.
        fname : str, optional, default = None
            If present, and if ax is not specified, save the image to the
            specified file.

        Notes
        -----
        This method plots the power of the class instance for each spherical
        harmonic degree and order, where the type of spectrum is defined by
        the parameter function: 'potential' for the magnetic potential,
        'radial' for the radial magnetic field, or 'total' for the total
        magnetic field. The case 'total' corresponds to the classic
        Lowes-Mauersberger spectrum. In all cases, the total power of the
        function is defined as the integral of the function squared over all
        space, divided by the area the function spans. If the mean of the
        function is zero, this is equivalent to the variance of the function.
        """
        if tick_interval is None:
            tick_interval = [None, None]
        if minor_tick_interval is None:
            minor_tick_interval = [None, None]

        if axes_labelsize is None:
            axes_labelsize = _mpl.rcParams['axes.labelsize']
            if type(axes_labelsize) == str:
                axes_labelsize = _mpl.font_manager \
                                 .FontProperties(size=axes_labelsize) \
                                 .get_size_in_points()
        if tick_labelsize is None:
            tick_labelsize = _mpl.rcParams['xtick.labelsize']
            if type(tick_labelsize) == str:
                tick_labelsize = _mpl.font_manager \
                                 .FontProperties(size=tick_labelsize) \
                                 .get_size_in_points()
        if titlesize is None:
            titlesize = _mpl.rcParams['axes.titlesize']
            if type(titlesize) == str:
                titlesize = _mpl.font_manager \
                                 .FontProperties(size=titlesize) \
                                 .get_size_in_points()

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
                "or 'unnorm'. Input value is {:s}."
                .format(repr(self.normalization)))

        if function == 'potential':
            spectrum *= self.r0**2
        elif function == 'radial':
            for l in degrees:
                spectrum[l, :] *= (l + 1)**2
        elif function == 'total':
            for l in degrees:
                spectrum[l, :] *= (l + 1) * (2 * l + 1)

        if origin in ('top', 'bottom'):
            spectrum = _np.rot90(spectrum, axes=(1, 0))

        spectrum_masked = _np.ma.masked_invalid(spectrum)

        # need to add one extra value to each in order for pcolormesh
        # to plot the last row and column.
        ls = _np.arange(lmax+2).astype(_np.float64)
        ms = _np.arange(-lmax, lmax + 2, dtype=_np.float64)
        if origin in ('left', 'right'):
            xgrid, ygrid = _np.meshgrid(ls, ms, indexing='ij')
        elif origin in ('top', 'bottom'):
            xgrid, ygrid = _np.meshgrid(ms, ls[::-1], indexing='ij')
        else:
            raise ValueError(
                "origin must be 'left', 'right', 'top', or 'bottom'. "
                "Input value is {:s}.".format(repr(origin)))
        xgrid -= 0.5
        ygrid -= 0.5

        if ax is None:
            if colorbar is not None:
                if colorbar in set(['top', 'bottom']):
                    scale = 1.2
                else:
                    scale = 0.9
            else:
                scale = 1.025
            figsize = (_mpl.rcParams['figure.figsize'][0],
                       _mpl.rcParams['figure.figsize'][0] * scale)
            fig = _plt.figure(figsize=figsize)
            axes = fig.add_subplot(111)
        else:
            axes = ax

        # make colormap
        if cmap_limits is None and cmap_rlimits is None:
            if cmap_scale.lower() == 'log':
                _temp = spectrum
                _temp[_temp == 0] = _np.NaN
                vmin = _np.nanmin(_temp)
            else:
                vmin = _np.nanmin(spectrum)
            vmax = _np.nanmax(spectrum)
            cmap_limits = [vmin, vmax]
        elif cmap_rlimits is not None:
            vmin = _np.nanmax(spectrum) * cmap_rlimits[0]
            vmax = _np.nanmax(spectrum) * cmap_rlimits[1]
            cmap_limits = [vmin, vmax]
            if len(cmap_rlimits) == 3:
                cmap_limits.append(cmap_rlimits[2])
        if len(cmap_limits) == 3:
            if cmap_scale.lower() == 'log':
                num = int(_np.log10(cmap_limits[1]/cmap_limits[0])
                          * cmap_limits[2])
            else:
                num = int((cmap_limits[1] - cmap_limits[0]) / cmap_limits[2])
            if isinstance(cmap, _mpl.colors.Colormap):
                cmap_scaled = cmap._resample(num)
            else:
                cmap_scaled = _mpl.cm.get_cmap(cmap, num)
        else:
            cmap_scaled = _mpl.cm.get_cmap(cmap)
        if cmap_reverse:
            cmap_scaled = cmap_scaled.reversed()

        if cmap_scale.lower() == 'log':
            norm = _mpl.colors.LogNorm(cmap_limits[0], cmap_limits[1],
                                       clip=True)
            # Clipping is required to avoid an invalid value error
        elif cmap_scale.lower() == 'lin':
            norm = _plt.Normalize(cmap_limits[0], cmap_limits[1])
        else:
            raise ValueError(
                "cmap_scale must be 'lin' or 'log'. " +
                "Input value is {:s}.".format(repr(cmap_scale)))

        # determine which ticks to plot
        if 'W' in ticks or 'L' in ticks:
            left, labelleft = True, True
        elif 'w' in ticks or 'l' in ticks:
            left, labelleft = True, False
        else:
            left, labelleft = False, False
        if 'S' in ticks or 'B' in ticks:
            bottom, labelbottom = True, True
        elif 's' in ticks or 'b' in ticks:
            bottom, labelbottom = True, False
        else:
            bottom, labelbottom = False, False
        if 'E' in ticks or 'R' in ticks:
            right, labelright = True, True
        elif 'e' in ticks or 'r' in ticks:
            right, labelright = True, False
        else:
            right, labelright = False, False
        if 'N' in ticks or 'T' in ticks:
            top, labeltop = True, True
        elif 'n' in ticks or 't' in ticks:
            top, labeltop = True, False
        else:
            top, labeltop = False, False

        # Set tick intervals (used only for linear axis)
        if tick_interval[0] is not None:
            degree_ticks = _np.linspace(
                0, lmax, num=lmax//tick_interval[0]+1, endpoint=True)
        if tick_interval[1] is not None:
            order_ticks = _np.linspace(
                -lmax, lmax, num=2*lmax//tick_interval[1]+1, endpoint=True)
        if minor_tick_interval[0] is not None:
            degree_minor_ticks = _np.linspace(
                0, lmax, num=lmax//minor_tick_interval[0]+1, endpoint=True)
        if minor_tick_interval[1] is not None:
            order_minor_ticks = _np.linspace(
                -lmax, lmax, num=2*lmax//minor_tick_interval[1]+1,
                endpoint=True)

        if (xscale == 'lin'):
            cmesh = axes.pcolormesh(xgrid, ygrid, spectrum_masked,
                                    norm=norm, cmap=cmap_scaled)
            if origin in ('left', 'right'):
                axes.set(xlim=(-0.5, lmax + 0.5))
                if tick_interval[0] is not None:
                    axes.set(xticks=degree_ticks)
                if minor_tick_interval[0] is not None:
                    axes.set_xticks(degree_minor_ticks, minor=True)
            else:
                axes.set(xlim=(-lmax - 0.5, lmax + 0.5))
                if tick_interval[1] is not None:
                    axes.set(xticks=order_ticks)
                if minor_tick_interval[1] is not None:
                    axes.set_xticks(order_minor_ticks, minor=True)
        elif (xscale == 'log'):
            cmesh = axes.pcolormesh(xgrid[1:], ygrid[1:], spectrum_masked[1:],
                                    norm=norm, cmap=cmap_scaled)
            if origin in ('left', 'right'):
                axes.set(xscale='log', xlim=(1., lmax + 0.5))
            else:
                axes.set(xscale='symlog', xlim=(-lmax - 0.5, lmax + 0.5))
        else:
            raise ValueError(
                "xscale must be 'lin' or 'log'. "
                "Input value is {:s}.".format(repr(xscale)))

        if (yscale == 'lin'):
            if origin in ('left', 'right'):
                axes.set(ylim=(-lmax - 0.5, lmax + 0.5))
                if tick_interval[1] is not None:
                    axes.set(yticks=order_ticks)
                if minor_tick_interval[1] is not None:
                    axes.set_yticks(order_minor_ticks, minor=True)
            else:
                axes.set(ylim=(-0.5, lmax + 0.5))
                if tick_interval[0] is not None:
                    axes.set(yticks=degree_ticks)
                if minor_tick_interval[0] is not None:
                    axes.set_yticks(degree_minor_ticks, minor=True)
        elif (yscale == 'log'):
            if origin in ('left', 'right'):
                axes.set(yscale='symlog', ylim=(-lmax - 0.5, lmax + 0.5))
            else:
                axes.set(yscale='log', ylim=(1., lmax + 0.5))
        else:
            raise ValueError(
                "yscale must be 'lin' or 'log'. "
                "Input value is {:s}.".format(repr(yscale)))

        axes.set_aspect('auto')
        if origin in ('left', 'right'):
            axes.set_xlabel(degree_label, fontsize=axes_labelsize)
            axes.set_ylabel(order_label, fontsize=axes_labelsize)
        else:
            axes.set_xlabel(order_label, fontsize=axes_labelsize)
            axes.set_ylabel(degree_label, fontsize=axes_labelsize)
        if labeltop:
            axes.xaxis.set_label_position('top')
        if labelright:
            axes.yaxis.set_label_position('right')
        axes.tick_params(bottom=bottom, top=top, right=right, left=left,
                         labelbottom=labelbottom, labeltop=labeltop,
                         labelleft=labelleft, labelright=labelright,
                         which='both')
        axes.tick_params(labelsize=tick_labelsize)
        axes.minorticks_on()
        axes.grid(grid, which='major')
        if title is not None:
            axes.set_title(title, fontsize=titlesize)
        if origin == 'right':
            axes.invert_xaxis()
        if origin == 'top':
            axes.invert_yaxis()

        # plot colorbar
        if colorbar is not None:
            if cb_label is None:
                if function == 'potential':
                    cb_label = 'Power, ' + self.units + '$^2$ m$^2$'
                elif function == 'radial':
                    cb_label = 'Power, ' + self.units + '$^2$'
                elif function == 'total':
                    cb_label = 'Power, ' + self.units + '$^2$'

            if cb_offset is None:
                offset = 1.3 * _mpl.rcParams['font.size']
                if (colorbar == 'left' and left) or \
                        (colorbar == 'right' and right) or \
                        (colorbar == 'bottom' and bottom) or \
                        (colorbar == 'top' and top):
                    offset += _mpl.rcParams['xtick.major.size']
                if (colorbar == 'left' and labelleft) or \
                        (colorbar == 'right' and labelright) or \
                        (colorbar == 'bottom' and labelbottom) or \
                        (colorbar == 'top' and labeltop):
                    offset += _mpl.rcParams['xtick.major.pad']
                    offset += tick_labelsize
                if origin in ('left', 'right') and colorbar == 'left' and \
                        order_label != '' and order_label is not None \
                        and labelleft:
                    offset += 1.9 * axes_labelsize
                if origin in ('left', 'right') and colorbar == 'right' \
                        and order_label != '' and order_label is not None \
                        and labelright:
                    offset += 1.9 * axes_labelsize
                if origin in ('bottom', 'top') and colorbar == 'left' \
                        and degree_label != '' \
                        and degree_label is not None and labelleft:
                    offset += 1.9 * axes_labelsize
                if origin in ('bottom', 'top') and colorbar == 'right' \
                        and degree_label != '' \
                        and degree_label is not None and labelright:
                    offset += 1.9 * axes_labelsize
                if origin in ('left', 'right') and colorbar == 'bottom' \
                        and degree_label != '' \
                        and degree_label is not None and labelbottom:
                    offset += axes_labelsize
                if origin in ('left', 'right') and colorbar == 'top' \
                        and degree_label != '' \
                        and degree_label is not None and labeltop:
                    offset += axes_labelsize
                if origin in ('bottom', 'top') and colorbar == 'bottom' \
                        and order_label != '' \
                        and order_label is not None and labelbottom:
                    offset += axes_labelsize
                if origin in ('bottom', 'top') and colorbar == 'top' \
                        and order_label != '' \
                        and order_label is not None and labeltop:
                    offset += axes_labelsize
            else:
                offset = cb_offset

            offset /= 72.  # convert to inches
            divider = _make_axes_locatable(axes)
            if colorbar in set(['left', 'right']):
                orientation = 'vertical'
                extendfrac = 0.025
                if cb_width is None:
                    size = '5%'
                else:
                    size = '{:f}%'.format(cb_width)
            else:
                orientation = 'horizontal'
                extendfrac = 0.025
                if cb_width is None:
                    size = '5%'
                else:
                    size = '{:f}%'.format(cb_width)
            cax = divider.append_axes(colorbar, size=size, pad=offset)
            cbar = _plt.colorbar(cmesh, cax=cax, orientation=orientation,
                                 extend=cb_triangles, extendfrac=extendfrac)
            if colorbar == 'left':
                cbar.ax.yaxis.set_ticks_position('left')
                cbar.ax.yaxis.set_label_position('left')
            if colorbar == 'top':
                cbar.ax.xaxis.set_ticks_position('top')
                cbar.ax.xaxis.set_label_position('top')
            cbar.set_label(cb_label, fontsize=axes_labelsize)
            cbar.ax.tick_params(labelsize=tick_labelsize)

        if ax is None:
            fig.tight_layout(pad=0.5)
            if show:
                fig.show()
            if fname is not None:
                fig.savefig(fname)
            return fig, axes

    def plot_correlation(self, hlm, lmax=None, grid=True, legend=None,
                         legend_loc='best', axes_labelsize=None,
                         tick_labelsize=None, ax=None, show=True, fname=None,
                         **kwargs):
        """
        Plot the correlation with another function.

        Usage
        -----
        x.plot_correlation(hlm, [lmax, grid, legend, legend_loc,
                                 axes_labelsize, tick_labelsize,
                                 ax, show, fname, **kwargs])

        Parameters
        ----------
        hlm : SHCoeffs class instance.
            The second function used in computing the spectral correlation.
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
        ax : matplotlib axes object, optional, default = None
            A single matplotlib axes object where the plot will appear.
        show : bool, optional, default = True
            If True, plot to the screen.
        fname : str, optional, default = None
            If present, and if ax is not specified, save the image to the
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
        if lmax is None:
            lmax = min(self.lmax, hlm.lmax)

        corr = self.correlation(hlm, lmax=lmax)

        ls = _np.arange(lmax + 1)

        if ax is None:
            fig, axes = _plt.subplots(1, 1)
        else:
            axes = ax

        if axes_labelsize is None:
            axes_labelsize = _mpl.rcParams['axes.labelsize']
            if type(axes_labelsize) == str:
                axes_labelsize = _mpl.font_manager \
                                 .FontProperties(size=axes_labelsize) \
                                 .get_size_in_points()
        if tick_labelsize is None:
            tick_labelsize = _mpl.rcParams['xtick.labelsize']
            if type(tick_labelsize) == str:
                tick_labelsize = _mpl.font_manager \
                                 .FontProperties(size=tick_labelsize) \
                                 .get_size_in_points()

        axes.plot(ls, corr, label=legend, **kwargs)
        if ax is None:
            axes.set(xlim=(0, lmax))
            axes.set(ylim=(-1, 1))
        else:
            axes.set(xlim=(0, max(lmax, ax.get_xbound()[1])))

        axes.set_xlabel('Spherical harmonic degree',
                        fontsize=axes_labelsize)
        axes.set_ylabel('Correlation', fontsize=axes_labelsize)
        axes.minorticks_on()
        axes.tick_params(labelsize=tick_labelsize)
        if legend is not None:
            axes.legend(loc=legend_loc)
        axes.grid(grid, which='major')

        if ax is None:
            fig.tight_layout(pad=0.5)
            if show:
                fig.show()
            if fname is not None:
                fig.savefig(fname)
            return fig, axes


class SHMagRealCoeffs(SHMagCoeffs):
    """
    Real spherical harmonic coefficient class for the magnetic potential.
    """

    def __init__(self, coeffs, r0=None, errors=None, error_kind=None,
                 normalization='schmidt', csphase=1, copy=True, header=None,
                 header2=None, name=None, units='nT', year=None):
        """Initialize real magnetic potential coefficients class."""
        lmax = coeffs.shape[1] - 1
        # ---- create mask to filter out m<=l ----
        mask = _np.zeros((2, lmax + 1, lmax + 1), dtype=bool)
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
        self.r0 = r0
        self.name = name
        self.units = units
        self.year = year
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

    def __repr__(self):
        return ('kind = {:s}\n'
                'normalization = {:s}\n'
                'csphase = {:d}\n'
                'lmax = {:d}\n'
                'r0 (m) = {:s}\n'
                'error_kind = {:s}\n'
                'header = {:s}\n'
                'header2 = {:s}\n'
                'name = {:s}\n'
                'units = {:s}\n'
                'year = {:s}'
                .format(repr(self.kind), repr(self.normalization),
                        self.csphase, self.lmax, repr(self.r0),
                        repr(self.error_kind), repr(self.header),
                        repr(self.header2), repr(self.name), repr(self.units),
                        repr(self.year)))

    def _rotate(self, angles, dj_matrix, r0=None, backend=None, nthreads=None):
        """Rotate the coefficients by the Euler angles alpha, beta, gamma."""
        if self.lmax > 1200 and backend.lower() == "shtools":
            _warnings.warn("The rotate() method is accurate only to about" +
                           " spherical harmonic degree 1200 when using the" +
                           " shtools backend. " +
                           "lmax = {:d}".format(self.lmax),
                           category=RuntimeWarning)
        if backend == "shtools" and dj_matrix is None:
            dj_matrix = _djpi2(self.lmax + 1)

        # The coefficients need to be 4pi normalized with csphase = 1
        coeffs = backend_module(
            backend=backend, nthreads=nthreads).SHRotateRealCoef(
                self.to_array(normalization='4pi', csphase=1, errors=False),
                angles, dj_matrix)

        # Convert 4pi normalized coefficients to the same normalization
        # as the unrotated coefficients.
        if self.normalization != '4pi' or self.csphase != 1:
            temp = _convert(coeffs, normalization_in='4pi', csphase_in=1,
                            normalization_out=self.normalization,
                            csphase_out=self.csphase)
            return SHMagCoeffs.from_array(temp, r0=r0, errors=self.errors,
                                          normalization=self.normalization,
                                          csphase=self.csphase,
                                          units=self.units, year=self.year,
                                          copy=False)
        else:
            return SHMagCoeffs.from_array(coeffs, r0=r0, errors=self.errors,
                                          units=self.units, year=self.year,
                                          copy=False)

    def _expand_coord(self, a, f, lat, lon, degrees, lmax_calc):
        """Evaluate the magnetic field at the coordinates lat and lon."""
        coeffs = self.to_array(normalization='schmidt', csphase=1,
                               errors=False)

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
            if f == 0.:
                r = a
            else:
                r = _np.cos(_np.deg2rad(latin))**2 + \
                    _np.sin(_np.deg2rad(latin))**2 / (1.0 - f)**2
                r = a * _np.sqrt(1. / r)

            return _MakeMagGridPoint(coeffs, a=self.r0, r=r, lat=latin,
                                     lon=lonin, lmax=lmax_calc)
        elif type(lat) is _np.ndarray:
            values = _np.empty((len(lat), 3), dtype=_np.float64)
            for i, (latitude, longitude) in enumerate(zip(latin, lonin)):
                if f == 0.:
                    r = a
                else:
                    r = _np.cos(_np.deg2rad(latitude))**2 + \
                        _np.sin(_np.deg2rad(latitude))**2 / (1.0 - f)**2
                    r = a * _np.sqrt(1. / r)

                values[i, :] = _MakeMagGridPoint(coeffs, a=self.r0, r=r,
                                                 lat=latitude, lon=longitude,
                                                 lmax=lmax_calc)
            return values
        elif type(lat) is list:
            values = []
            for latitude, longitude in zip(latin, lonin):
                if f == 0.:
                    r = a
                else:
                    r = _np.cos(_np.deg2rad(latitude))**2 + \
                        _np.sin(_np.deg2rad(latitude))**2 / (1.0 - f)**2
                    r = a * _np.sqrt(1. / r)
                values.append(
                    _MakeMagGridPoint(coeffs, a=self.r0, r=r, lat=latitude,
                                      lon=longitude, lmax=lmax_calc))
            return values
        else:
            raise ValueError('lat and lon must be either an int, float, '
                             'ndarray, or list. Input types are {:s} and {:s}.'
                             .format(repr(type(lat)), repr(type(lon))))
