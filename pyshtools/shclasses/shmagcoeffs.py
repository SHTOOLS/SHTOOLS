"""
    Class for spherical harmonic coefficients of the magnetic potential.
"""
import numpy as _np
import matplotlib as _mpl
import matplotlib.pyplot as _plt
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
from ..shtools import MakeMagGridDH as _MakeMagGridDH
from ..shtools import MakeMagGradGridDH as _MakeMagGradGridDH
from ..shtools import djpi2 as _djpi2
from ..shtools import SHRotateRealCoef as _SHRotateRealCoef


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
                    csphase conventions.
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
                   normalization='schmidt', csphase=1, lmax=None, units='nT',
                   year=None, copy=True):
        """
        Initialize the class with spherical harmonic coefficients from an input
        array.

        Usage
        -----
        x = SHMagCoeffs.from_array(array, r0, [errors, error_kind,
                                               normalization, csphase, lmax,
                                               units, year, copy])

        Returns
        -------
        x : SHMagCoeffs class instance.

        Parameters
        ----------
        array : ndarray, shape (2, lmaxin+1, lmaxin+1).
            The input spherical harmonic coefficients.
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
        units : str, optional, default = 'nT'
            The units of the spherical harmonic coefficients.
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
                                  csphase=csphase, units=units, year=year,
                                  copy=copy)
        else:
            clm = SHMagRealCoeffs(coeffs[:, 0:lmax+1, 0:lmax+1], r0=r0,
                                  normalization=normalization.lower(),
                                  csphase=csphase, units=units, year=year,
                                  copy=copy)
        return clm

    @classmethod
    def from_zeros(self, lmax, r0, errors=None, error_kind=None,
                   normalization='schmidt', csphase=1, units='nT', year=None):
        """
        Initialize the class with spherical harmonic coefficients set to zero.

        Usage
        -----
        x = SHMagCoeffs.from_zeros(lmax, r0, [errors, error_kind,
                                              normalization, csphase,
                                              units, year])

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
        units : str, optional, default = 'nT'
            The units of the spherical harmonic coefficients.
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
                              csphase=csphase, units=units, year=year)
        return clm

    @classmethod
    def from_file(self, fname, format='shtools', r0=None, lmax=None,
                  normalization='schmidt', skip=0, header=True, header2=False,
                  errors=None, error_kind=None, csphase=1, r0_index=0,
                  header_units='m', file_units='nT', units='nT', year=None,
                  **kwargs):
        """
        Initialize the class with spherical harmonic coefficients from a file.

        Usage
        -----
        x = SHMagCoeffs.from_file(filename, [format='shtools' or 'dov', r0,
                                  lmax, normalization, csphase, skip, header,
                                  header2, errors, error_kind, r0_index,
                                  header_units, file_units, units, year])
        x = SHMagCoeffs.from_file(filename, format='igrf', r0, year, [lmax,
                                  normalization, csphase, file_units, units])
        x = SHMagCoeffs.from_file(filename, format='bshc', r0, [lmax,
                                  normalization, csphase, file_units, units,
                                  year])
        x = SHMagCoeffs.from_file(filename, format='npy', r0, [lmax,
                                  normalization, csphase, file_units, units,
                                  year, **kwargs])

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
        units : str, optional, default = 'nT'
            The units of the spherical harmonic coefficients.
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
                        fname, lmax=lmax, error=True, skip=skip)
                else:
                    coeffs, lmaxout = read_func(fname, lmax=lmax, skip=skip)

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
            coeffs = _read_igrf(fname, year=year)
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
                              header2=header2_list, units=units, year=year)
        return clm

    @classmethod
    def from_random(self, power, r0, function='total', lmax=None,
                    normalization='schmidt', csphase=1, units='nT', year=None,
                    exact_power=False):
        """
        Initialize the class of magnetic potential spherical harmonic
        coefficients as random variables with a given spectrum.

        Usage
        -----
        x = SHMagCoeffs.from_random(power, r0, [function, lmax, normalization,
                                                csphase, units, year,
                                                exact_power])

        Returns
        -------
        x : SHMagCoeffs class instance.

        Parameters
        ----------
        power : ndarray, shape (L+1)
            numpy array of shape (L+1) that specifies the expected power per
            degree l, where L is the maximum spherical harmonic bandwidth.
        r0 : float
            The reference radius of the spherical harmonic coefficients.
        function : str, optional, default = 'potential'
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
        units : str, optional, default = 'nT'
            The units of the spherical harmonic coefficients.
        year : float, default = None.
            The year of the time-variable spherical harmonic coefficients.
        exact_power : bool, optional, default = False
            The total variance of the coefficients is set exactly to the input
            power. The distribution of power at degree l amongst the angular
            orders is random, but the total power is fixed.

        Notes
        -----
        This routine returns a random realization of spherical harmonic
        magnetic potential coefficients obtained from a normal distribution.
        The variance of each coefficient at degree l is equal to
        the total power at degree l divided by the number of coefficients at
        that degree (2l+1). These coefficients are then divided by a prefactor
        that depends upon the function being used to calculate the spectrum:
        r0 for the magnetic potential, (l+1) for the radial magnetic field,
        or sqrt((l+1)*(2l+1)). The power spectrum of the random realization can
        be fixed exactly to the input spectrum by setting exact_power to True.
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
                              csphase=csphase, units=units, year=year)
        return clm

    @classmethod
    def from_netcdf(self, filename, lmax=None, normalization='schmidt',
                    csphase=1, units='nT', year=None):
        """
        Initialize the class with spherical harmonic coefficients from a
        netcdf file.

        Usage
        -----
        x = SHMagCoeffs.from_netcdf(filename, [lmax, normalization, csphase,
                                               units, year])

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
        units : str, optional, default = 'nT'
            The units of the spherical harmonic coefficients.
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
                              csphase=csphase, units=units, year=year)
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

        mneg_mask = (ms < 0).astype(_np.int)
        self.coeffs[mneg_mask, ls, _np.abs(ms)] = values

    # ---- IO routines ----
    def to_file(self, filename, format='shtools', header=None, errors=True,
                lmax=None, **kwargs):
        """
        Save spherical harmonic coefficients to a file.

        Usage
        -----
        x.to_file(filename, [format='shtools', header, errors])
        x.to_file(filename, format='dov', [header, errors])
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
                           header=header, header2=header2, lmax=lmax)
            else:
                write_func(filebase, self.coeffs, errors=None,
                           header=header, header2=header2, lmax=lmax)

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

        if self.lmax > 1200:
            _warnings.warn("The rotate() method is accurate only to about"
                           " spherical harmonic degree 1200. "
                           "lmax = {:d}.".format(self.lmax),
                           category=RuntimeWarning)

        rot = self._rotate(angles, dj_matrix, r0=self.r0)
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
                                           csphase=csphase, lmax=lmax)
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
            mask = _np.zeros((2, lmax + 1, lmax + 1), dtype=_np.bool)
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

    # ---- Routines that return different magnetic-related class instances ----
    def expand(self, a=None, f=None, lmax=None, lmax_calc=None, sampling=2,
               extend=True):
        """
        Create 2D cylindrical maps on a flattened and rotating ellipsoid of all
        three components of the magnetic field, the total magnetic intensity,
        and the magnetic potential, and return as a SHMagGrid class instance.

        Usage
        -----
        mag = x.expand([a, f, lmax, lmax_calc, sampling, extend])

        Returns
        -------
        mag : SHMagGrid class instance.

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
        or nT m (for the potential),

        The magnetic potential is given by

        V = r0 Sum_{l=1}^lmax (r0/r)^{l+1} Sum_{m=-l}^l g_{lm} Y_{lm}

        and the magnetic field is

        B = - Grad V.

        The coefficients are referenced to a radius r0, and the function is
        computed on a flattened ellipsoid with semi-major axis a (i.e., the
        mean equatorial radius) and flattening f.
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
            coeffs, errors = self.to_array(normalization='schmidt', csphase=1)
        else:
            coeffs = self.to_array(normalization='schmidt', csphase=1)

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

        if self.errors is not None:
            coeffs, errors = self.to_array(normalization='schmidt', csphase=1)
        else:
            coeffs = self.to_array(normalization='schmidt', csphase=1)

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
                      axes_labelsize=None, tick_labelsize=None, show=True,
                      ax=None, fname=None, **kwargs):
        """
        Plot the spectrum as a function of spherical harmonic degree.

        Usage
        -----
        x.plot_spectrum([function, unit, base, lmax, xscale, yscale, grid,
                         legend, legend_loc, axes_labelsize, tick_labelsize,
                         show, ax, fname, **kwargs])

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
        if tick_labelsize is None:
            tick_labelsize = _mpl.rcParams['xtick.labelsize']

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

    def plot_spectrum2d(self, function='total', xscale='lin', yscale='lin',
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
            magnetic potential, 'radial' for the radial magnetic field, or
            'total' for the total magnetic field (Lowes-Mauersberger).
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

        if function == 'potential':
            cb.set_label('Power, ' + self.units + '$^2$ m$^2$',
                         fontsize=axes_labelsize)
        elif function == 'radial':
            cb.set_label('Power, ' + self.units + '$^2$',
                         fontsize=axes_labelsize)
        elif function == 'total':
            cb.set_label('Power, ' + self.units + '$^2$',
                         fontsize=axes_labelsize)

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

    def plot_correlation(self, hlm, lmax=None, grid=True, legend=None,
                         legend_loc='best', axes_labelsize=None,
                         tick_labelsize=None, show=True, ax=None, fname=None,
                         **kwargs):
        """
        Plot the correlation with another function.

        Usage
        -----
        x.plot_correlation(hlm, [lmax, grid, legend, legend_loc,
                                 axes_labelsize, tick_labelsize,
                                 show, ax, fname, **kwargs])

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
        if tick_labelsize is None:
            tick_labelsize = _mpl.rcParams['xtick.labelsize']

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
                 header2=None, units='nT', year=None):
        """Initialize real magnetic potential coefficients class."""
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
        self.r0 = r0
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
                'units = {:s}\n'
                'year = {:s}'
                .format(repr(self.kind), repr(self.normalization),
                        self.csphase, self.lmax, repr(self.r0),
                        repr(self.error_kind), repr(self.header),
                        repr(self.header2), repr(self.units), repr(self.year)))

    def _rotate(self, angles, dj_matrix, r0=None):
        """Rotate the coefficients by the Euler angles alpha, beta, gamma."""
        if dj_matrix is None:
            dj_matrix = _djpi2(self.lmax + 1)

        # The coefficients need to be 4pi normalized with csphase = 1
        coeffs = _SHRotateRealCoef(
            self.to_array(normalization='4pi', csphase=1), angles, dj_matrix)

        # Convert 4pi normalized coefficients to the same normalization
        # as the unrotated coefficients.
        if self.normalization != '4pi' or self.csphase != 1:
            temp = _convert(coeffs, normalization_in='4pi', csphase_in=1,
                            normalization_out=self.normalization,
                            csphase_out=self.csphase)
            return SHMagCoeffs.from_array(temp, r0=r0,
                                          normalization=self.normalization,
                                          csphase=self.csphase,
                                          units=self.units, year=self.year,
                                          copy=False)
        else:
            return SHMagCoeffs.from_array(coeffs, r0=r0, units=self.units,
                                          year=self.year, copy=False)
