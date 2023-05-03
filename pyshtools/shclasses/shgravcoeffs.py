"""
    Class for spherical harmonic coefficients of the gravitational potential.
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
from .shcoeffs import SHRealCoeffs as _SHRealCoeffs
from .shgrid import DHRealGrid as _DHRealGrid
from .shgravgrid import SHGravGrid as _SHGravGrid
from .shtensor import SHGravTensor as _SHGravTensor
from .shgeoid import SHGeoid as _SHGeoid

from ..constants import G as _G
from ..spectralanalysis import spectrum as _spectrum
from ..spectralanalysis import cross_spectrum as _cross_spectrum
from ..shio import convert as _convert
from ..shio import shread as _shread
from ..shio import shwrite as _shwrite
from ..shio import read_dov as _read_dov
from ..shio import write_dov as _write_dov
from ..shio import read_bshc as _read_bshc
from ..shio import write_bshc as _write_bshc
from ..shio import read_icgem_gfc as _read_icgem_gfc
from ..shio import write_icgem_gfc as _write_icgem_gfc
from ..backends.shtools import CilmPlusRhoHDH as _CilmPlusRhoHDH
from ..backends.shtools import CilmPlusDH as _CilmPlusDH
from ..backends.shtools import MakeGravGridDH as _MakeGravGridDH
from ..backends.shtools import MakeGravGradGridDH as _MakeGravGradGridDH
from ..backends.shtools import MakeGeoidGridDH as _MakeGeoidGridDH
from ..backends.shtools import djpi2 as _djpi2
from ..backends.shtools import MakeGravGridPoint as _MakeGravGridPoint
from ..backends import backend_module
from ..backends import preferred_backend


class SHGravCoeffs(object):
    """
    Spherical harmonic coefficients class for the gravitational potential.

    The coefficients of this class can be initialized using one of the four
    constructor methods:

        x = SHGravCoeffs.from_array(array, gm, r0)
        x = SHGravCoeffs.from_random(powerspectrum, gm, r0)
        x = SHGravCoeffs.from_zeros(lmax, gm, r0)
        x = SHGravCoeffs.from_file('fname.dat')
        x = SHGravCoeffs.from_netcdf('ncname.nc')
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
                    csphase conventions. This is a three-dimensional array
                    formatted as coeffs[i, degree, order], where i=0
                    corresponds to positive orders and i=1 to negative orders.
    errors        : The uncertainties of the spherical harmonic coefficients.
    error_kind    : An arbitrary string describing the kind of errors, such as
                    'unknown', 'unspecified', 'calibrated', 'formal' or None.
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
    name          : The name of the dataset.
    epoch         : The epoch time of the spherical harmonic coefficients.
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
    admittance()          : Return the admittance with an input topography
                            function.
    correlation()         : Return the spectral correlation with another
                            function.
    admitcorr()           : Return the admittance and spectral correlation with
                            an input topography function.
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
    mass                  : Return the mass of the planet.
    center_of_mass        : Return coordinates of the center of mass of the
                            planet.
    inertia_tensor()      : Return an array of the inertia tensor.
    tensor()              : Calculate the 9 components of the gravity tensor
                            and return an SHGravTensor class instance.
    geoid()               : Calculate the height of the geoid and return an
                            SHGeoid class instance.
    plot_spectrum()       : Plot the spectrum as a function of spherical
                            harmonic degree.
    plot_spectrum2d()     : Plot the 2D spectrum of all spherical harmonic
                            degrees and orders.
    plot_correlation()    : Plot the spectral correlation with another
                            function.
    plot_admittance()     : Plot the admittance with an input topography
                            function.
    plot_admitcorr()      : Plot the admittance and spectral correlation with
                            an input topography function.
    to_array()            : Return an array of spherical harmonic coefficients
                            with a different normalization convention.
    to_file()             : Save the spherical harmonic coefficients as a file.
    to_netcdf()           : Save raw spherical harmonic coefficients as a
                            netcdf file.
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
              '>>> pyshtools.SHGravCoeffs.from_netcdf\n'
              '>>> pyshtools.SHGravCoeffs.from_shape\n')

    # ---- Factory methods ----
    @classmethod
    def from_array(self, coeffs, gm, r0, omega=None, errors=None,
                   error_kind=None, normalization='4pi', csphase=1, lmax=None,
                   set_degree0=True, name=None, epoch=None, copy=True):
        """
        Initialize the class with spherical harmonic coefficients from an input
        array.

        Usage
        -----
        x = SHGravCoeffs.from_array(array, gm, r0, [omega, errors, error_kind,
                                                    normalization, csphase,
                                                    lmax, set_degree0, name,
                                                    epoch, copy])

        Returns
        -------
        x : SHGravCoeffs class instance.

        Parameters
        ----------
        array : ndarray, shape (2, lmaxin+1, lmaxin+1).
            The input spherical harmonic coefficients. This is a three-
            dimensional array formatted as coeffs[i, degree, order], where i=0
            corresponds to positive orders and i=1 to negative orders.
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
        set_degree0 : bool, optional, default = True
            If the degree-0 coefficient is zero, set this to 1.
        name : str, optional, default = None
            The name of the dataset.
        epoch : str or float, optional, default = None
            The epoch time of the spherical harmonic coefficients as given by
            the format YYYYMMDD.DD.
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

        if coeffs[0, 0, 0] == 0 and set_degree0:
            coeffs[0, 0, 0] = 1.0

        if errors is not None:
            clm = SHGravRealCoeffs(coeffs[:, 0:lmax+1, 0:lmax+1], gm=gm, r0=r0,
                                   omega=omega, errors=errors[:, 0:lmax+1,
                                                              0:lmax+1],
                                   error_kind=error_kind,
                                   normalization=normalization.lower(),
                                   csphase=csphase, name=name, epoch=epoch,
                                   copy=copy)
        else:
            clm = SHGravRealCoeffs(coeffs[:, 0:lmax+1, 0:lmax+1], gm=gm, r0=r0,
                                   omega=omega,
                                   normalization=normalization.lower(),
                                   csphase=csphase, name=name, epoch=epoch,
                                   copy=copy)
        return clm

    @classmethod
    def from_zeros(self, lmax, gm, r0, omega=None, errors=None,
                   error_kind=None, normalization='4pi', csphase=1,
                   name=None, epoch=None):
        """
        Initialize the class with spherical harmonic coefficients set to zero
        from degree 1 to lmax, and set the degree 0 term to 1.

        Usage
        -----
        x = SHGravCoeffs.from_zeros(lmax, gm, r0, [omega, errors, error_kind,
                                                   normalization, csphase,
                                                   name, epoch])

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
        name : str, optional, default = None
            The name of the dataset.
        epoch : str or float, optional, default = None
            The epoch time of the spherical harmonic coefficients as given by
            the format YYYYMMDD.DD.
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

        if normalization.lower() == 'unnorm' and lmax > 85:
            _warnings.warn("Calculations using unnormalized coefficients "
                           "are stable only for degrees less than or equal "
                           "to 85. lmax for the coefficients will be set to "
                           "85. Input value is {:d}.".format(lmax),
                           category=RuntimeWarning)
            lmax = 85

        coeffs = _np.zeros((2, lmax + 1, lmax + 1))
        coeffs[0, 0, 0] = 1.0
        if errors is True:
            error_coeffs = _np.zeros((2, lmax + 1, lmax + 1))
            if error_kind is None:
                error_kind = 'unspecified'
        else:
            error_coeffs = None

        clm = SHGravRealCoeffs(coeffs, gm=gm, r0=r0, omega=omega,
                               errors=error_coeffs, error_kind=error_kind,
                               normalization=normalization.lower(),
                               csphase=csphase, name=name, epoch=epoch)
        return clm

    @classmethod
    def from_file(self, fname, format='shtools', gm=None, r0=None,
                  omega=None, lmax=None, normalization='4pi', skip=0,
                  header=True, header2=False, errors=None, error_kind=None,
                  csphase=1, r0_index=0, gm_index=1, omega_index=None,
                  header_units='m', set_degree0=True, name=None, epoch=None,
                  encoding=None, quiet=False, **kwargs):
        """
        Initialize the class with spherical harmonic coefficients from a file.

        Usage
        -----
        x = SHGravCoeffs.from_file(filename, [format='shtools' or 'dov', gm,
                                   r0, omega, lmax, normalization, csphase,
                                   skip, header, header2, errors, error_kind,
                                   gm_index, r0_index, omega_index,
                                   header_units, set_degree0, name, encoding])
        x = SHGravCoeffs.from_file(filename, format='icgem', [lmax, omega,
                                   normalization, csphase, errors, set_degree0,
                                   name, name, epoch, encoding, quiet])
        x = SHGravCoeffs.from_file(filename, format='bshc', gm, r0, [lmax,
                                   omega, normalization, csphase, set_degree0,
                                   name])
        x = SHGravCoeffs.from_file(filename, format='npy', gm, r0, [lmax,
                                   omega, normalization, csphase, set_degree0,
                                   name, **kwargs])

        Returns
        -------
        x : SHGravCoeffs class instance.

        Parameters
        ----------
        filename : str
            File name or URL containing the spherical harmonic coefficients.
            filename will be treated as a URL if it starts with 'http://',
            'https://', or 'ftp://'. For 'shtools', 'icgem' and 'bshc'
            formatted files, if filename ends with '.gz' or '.zip' (or if the
            path contains '/zip/'), the file will be uncompressed before
            parsing.
        format : str, optional, default = 'shtools'
            'shtools' for generic text files, 'dov' for [degree, order, value]
            text files, 'icgem' for ICGEM GFC formatted files, 'bshc' for
            binary spherical harmonic coefficient files, or 'npy' for binary
            numpy files.
        lmax : int, optional, default = None
            The maximum spherical harmonic degree to read from the file. The
            default is to read the entire file.
        header : bool, optional, default = True
            If True, read a list of values from the header line of an 'shtools'
            or 'dov' formatted file. If two header lines are present, the
            second contains values for r0, gm, and omega.
        header2 : bool, optional, default = False
            If True, read a list of values from a second header line of an
            'shtools' or 'dov' formatted file. If two header lines are present,
            the second contains values for r0, gm, and omega.
        errors : bool or str, optional, default = None
            For 'shtools' or 'dov' formatted files: if True, read and return
            the spherical harmonic coefficients of the errors. For 'icgem'
            formatted files, specify the type of error to return: 'calibrated'
            or 'formal'.
        error_kind : str, optional, default = None
            For 'shtools' and 'dov' formatted files: An arbitrary string
            describing the kind of errors, such as None, 'unspecified',
            'calibrated' or 'formal'.
        r0_index : int, optional, default = 0
            For 'shtools' and 'dov' formatted files, r0 will be set using the
            value from the last header line with this index.
        gm_index : int, optional, default = 1
            For 'shtools' and 'dov' formatted files, gm will be set using the
            value from the last header line with this index.
        omega_index : int, optional, default = None
            For 'shtools' and 'dov' formatted files, omega will be set using
            the value from the last header line with this index.
        gm : float, optional, default = None
            The gravitational constant time the mass that is associated with
            the gravitational potential coefficients.
        r0 : float, optional, default = None
            The reference radius of the spherical harmonic coefficients.
        omega : float, optional, default = None
            The angular rotation rate of the body.
        header_units : str, optional, default = 'm'
            The units used for r0 and gm in the header line of an 'shtools'
            formatted file: 'm' or 'km'. If 'km', the values of r0 and gm will
            be converted to meters.
        set_degree0 : bool, optional, default = True
            If the degree-0 coefficient is zero, set this to 1.
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
        name : str, optional, default = None
            The name of the dataset.
        epoch : str or float, optional, default = None
            The epoch time of the spherical harmonic coefficients as given by
            the format YYYYMMDD.DD. If format is 'icgem' and epoch is None,
            the reference epoch t0 of the model will be used. Epoch is required
            for 'icgem' v2.0 formatted files.
        encoding : str, optional, default = None
            Encoding of the input file when format is 'shtools', 'dov' or
            'icgem'. The default is to use the system default.
        quiet : bool, default = False
            If True, suppress warnings about undefined keywords when reading
            ICGEM formatted files.
        **kwargs : keyword argument list, optional for format = 'npy'
            Keyword arguments of numpy.load() when format is 'npy'.

        Notes
        -----
        Supported file formats:
            'shtools' (see pyshtools.shio.shread)
            'dov' (see pyshtools.shio.read_dov)
            'icgem' (see pyshtools.shio.read_icgem_gfc)
            'bshc' (see pyshtools.shio.read_bshc)
            'npy' (see numpy.load)

        If the degree 0 term of the file is zero (or not specified), this will
        by default be set to 1.

        For 'shtools', 'dov', 'icgem' and 'bshc' formatted files, if filename
        starts with 'http://', 'https://', or 'ftp://', the file will be
        treated as a URL. In this case, the file will be downloaded in its
        entirety before it is parsed. If the filename ends with '.gz' or '.zip'
        (or if the path contains '/zip/'), the file will be automatically
        uncompressed before parsing. For zip files, archives with only a single
        file are supported. Note that reading '.gz' and '.zip' files in
        'shtools' format will be extremely slow if lmax is not specified.

        For 'shtools' and 'dov' formatted files, the optional parameter `skip`
        specifies how many lines should be skipped before attempting to parse
        the file, the optional parameter `header` and `header2` specifies
        whether to read a list of values from one or two header lines, and the
        optional parameter `lmax` specifies the maximum degree to read from the
        file. If header lines are read, r0_index, gm_index, and omega_index,
        are used as the indices to set r0, gm, and omega from the last header
        line. If header_unit is specified as 'km', the values of r0 and gm that
        are read from the header will be converted to meters.
        """
        error_coeffs = None
        header_list = None
        header2_list = None

        if not header:
            r0_index = None
            gm_index = None

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
            if header_units.lower() not in ('m', 'km'):
                raise ValueError("header_units can be only 'm', or 'km'. "
                                 "Input value is {:s}."
                                 .format(repr(header_units)))
            if r0_index is not None and r0 is not None:
                raise ValueError('Can not specify both r0_index and r0.')
            if gm_index is not None and gm is not None:
                raise ValueError('Can not specify both gm_index and gm.')
            if omega_index is not None and omega is not None:
                raise ValueError('Can not specify both omega_index and omega,')
            if header is False and (r0 is None or gm is None):
                raise ValueError('If header is False, r0 and gm must be '
                                 'specified.')

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
                                                     header2=True,
                                                     error=True,
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
                if gm_index is not None:
                    if header2:
                        gm = float(header2_list[gm_index])
                    else:
                        gm = float(header_list[gm_index])
                if omega_index is not None:
                    if header2:
                        omega = float(header2_list[omega_index])
                    else:
                        omega = float(header_list[omega_index])
                if header_units.lower() == 'km':
                    r0 *= 1.e3
                    gm *= 1.e9

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
            if gm is None or r0 is None:
                raise ValueError('For binary bshc files, gm and r0 must be '
                                 'specified.')
            coeffs, lmaxout = _read_bshc(fname, lmax=lmax)

        elif format.lower() == 'icgem':
            valid_err = ('unknown', 'calibrated', 'formal')
            if errors is False or errors is None:
                coeffs, gm, r0 = _read_icgem_gfc(filename=fname,
                                                 errors=None, lmax=lmax,
                                                 epoch=epoch,
                                                 encoding=encoding,
                                                 quiet=quiet)
            elif errors in valid_err:
                coeffs, gm, r0, error_coeffs = _read_icgem_gfc(
                    filename=fname, errors=errors, lmax=lmax, epoch=epoch,
                    encoding=encoding, quiet=quiet)
                error_kind = errors
            else:
                raise ValueError('errors must be among: {}. '
                                 'Input value is {:s}.'
                                 .format(valid_err, repr(errors)))
            lmaxout = coeffs.shape[1] - 1

        elif format.lower() == 'npy':
            if gm is None or r0 is None:
                raise ValueError('For binary npy files, gm and r0 must be '
                                 'specified.')
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

        if coeffs[0, 0, 0] == 0 and set_degree0:
            coeffs[0, 0, 0] = 1.0

        clm = SHGravRealCoeffs(coeffs, gm=gm, r0=r0, omega=omega,
                               errors=error_coeffs, error_kind=error_kind,
                               normalization=normalization.lower(),
                               csphase=csphase, header=header_list,
                               header2=header2_list, name=name, epoch=epoch)
        return clm

    @classmethod
    def from_random(self, power, gm, r0, omega=None, function='geoid',
                    lmax=None, normalization='4pi', csphase=1,
                    exact_power=False, power_unit='per_l', name=None,
                    epoch=None):
        """
        Initialize the class of gravitational potential spherical harmonic
        coefficients as random variables with a given spectrum.

        Usage
        -----
        x = SHGravCoeffs.from_random(power, gm, r0, [omega, function, lmax,
                                                     normalization,
                                                     csphase, exact_power,
                                                     power_unit, name, epoch])

        Returns
        -------
        x : SHGravCoeffs class instance.

        Parameters
        ----------
        power : ndarray, shape (L+1)
            numpy array of shape (L+1) that specifies the expected power
            spectrum of the output function, where L is the maximum spherical
            harmonic bandwidth. By default, the input power spectrum represents
            the power of all angular orders of the geoid as a function of
            spherical harmonic degree (see function and power_unit).
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
            If True, the spherical harmonic coefficients of the random
            realization will be rescaled such that the power spectrum is
            exactly equal to the input spectrum.
        power_unit : str, optional, default = 'per_l'
            If 'per_l', the input power spectrum represents the total power of
            all angular orders as a function of spherical harmonic degree. If
            'per_lm', the input power spectrum represents the power per
            coefficient (which is assumed isotropic and varies only as a
            function of spherical harmonic degree).
        name : str, optional, default = None
            The name of the dataset.
        epoch : str or float, optional, default = None
            The epoch time of the spherical harmonic coefficients as given by
            the format YYYYMMDD.DD.

        Notes
        -----
        This routine returns a random realization of spherical harmonic
        gravitational potential coefficients obtained from a normal
        distribution. The variance of each coefficient is determined by the
        input power spectrum and the type of spectrum (as specified by
        function and power_unit). If power_unit is 'per_l' (default), the
        variance of each coefficient at spherical harmonic degree l is equal to
        the input total power at degree l divided by the number of coefficients
        at that degree. If power_unit is 'per_lm', the variance of each
        coefficient at degree l is equal to the input power at that degree.
        The power of the input function can be either for the geoid (default),
        potential, radial gravity, or total gravity field. The power spectrum
        of the random realization can be fixed exactly to the input spectrum by
        setting exact_power to True.

        Note that the degree 0 term is set to 1, and the degree-1 terms are
        set to 0.
        """
        if type(normalization) != str:
            raise ValueError('normalization must be a string. '
                             'Input type is {:s}.'
                             .format(str(type(normalization))))

        if function.lower() not in ('potential', 'geoid', 'radial', 'total'):
            raise ValueError(
                "function must be of type 'potential', "
                "'geoid', 'radial', or 'total'. Provided value is {:s}."
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
                               errors=None,
                               normalization=normalization.lower(),
                               csphase=csphase, name=name, epoch=epoch)
        return clm

    @classmethod
    def from_netcdf(self, filename, lmax=None, normalization='4pi', csphase=1,
                    name=None, epoch=None):
        """
        Initialize the class with spherical harmonic coefficients from a
        netcdf file.

        Usage
        -----
        x = SHGravCoeffs.from_netcdf(filename, [lmax, normalization, csphase,
                                                name, epoch])

        Returns
        -------
        x : SHGravCoeffs class instance.

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
        name : str, optional, default = None
            The name of the dataset.
        epoch : str or float, optional, default = None
            The epoch time of the spherical harmonic coefficients as given by
            the format YYYYMMDD.DD.

        Description
        -----------
        The format of the netcdf file has to be exactly as the format that is
        used in SHGravCoeffs.to_netcdf().
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
            gm = ds.coeffs.GM
        except:
            raise ValueError("coeffs.GM must be specified in the netcdf file.")
        try:
            r0 = ds.coeffs.r0
        except:
            raise ValueError("coeffs.r0 must be specified in the netcdf file.")
        try:
            omega = ds.coeffs.omega
        except:
            omega = None
        try:
            epoch = ds.coeffs.epoch
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
            raise ValueError('Gravitational potential coefficients must be '
                             'real. Input coefficients are complex.')

        clm = SHGravRealCoeffs(coeffs, gm=gm, r0=r0, omega=omega,
                               errors=errors, error_kind=error_kind,
                               normalization=normalization.lower(),
                               csphase=csphase, name=name, epoch=epoch)
        return clm

    @classmethod
    def from_shape(self, shape, rho, gm, nmax=7, lmax=None, lmax_grid=None,
                   lmax_calc=None, omega=None, name=None, epoch=None,
                   backend=None, nthreads=None):
        """
        Initialize a class of gravitational potential spherical harmonic
        coefficients by calculuting the gravitational potential associatiated
        with relief along an interface.

        Usage
        -----
        x = SHGravCoeffs.from_shape(shape, rho, gm, [nmax, lmax, lmax_grid,
                                                     lmax_calc, omega, name,
                                                     epoch, backend, nthreads])

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
        name : str, optional, default = None
            The name of the dataset.
        epoch : str or float, optional, default = None
            The epoch time of the spherical harmonic coefficients as given by
            the format YYYYMMDD.DD.
        backend : str, optional, default = preferred_backend()
            Name of the preferred backend, either 'shtools' or 'ducc'.
        nthreads : int, optional, default = 1
            Number of threads to use for the 'ducc' backend. Setting this
            parameter to 0 will use as many threads as there are hardware
            threads on the system.

        Notes
        -----
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

        if backend is None:
            backend = preferred_backend()

        if type(shape) is not _SHRealCoeffs and type(shape) is not _DHRealGrid:
            raise ValueError('shape must be of type SHRealCoeffs '
                             'or DHRealGrid. Input type is {:s}.'
                             .format(repr(type(shape))))

        if (not issubclass(type(rho), float) and type(rho) is not int
                and type(rho) is not _np.ndarray and
                type(rho) is not _SHRealCoeffs and
                type(rho is not _DHRealGrid)):
            raise ValueError('rho must be of type float, int, ndarray, '
                             'SHRealCoeffs or DHRealGrid. Input type is {:s}.'
                             .format(repr(type(rho))))

        if type(shape) is _SHRealCoeffs:
            shape = shape.expand(lmax=lmax_grid, lmax_calc=lmax_calc,
                                 backend=backend, nthreads=nthreads)

        if type(rho) is _SHRealCoeffs:
            rho = rho.expand(lmax=lmax_grid, lmax_calc=lmax_calc,
                             backend=backend, nthreads=nthreads)

        if type(rho) is _DHRealGrid:
            if shape.lmax != rho.lmax:
                raise ValueError('The grids for shape and rho must have the '
                                 'same size. '
                                 'lmax of shape = {:d}, lmax of rho = {:d}.'
                                 .format(shape.lmax, rho.lmax))
            cilm, d = _CilmPlusRhoHDH(shape.data[:shape.nlat-shape.extend,
                                                 :shape.nlon-shape.extend],
                                      nmax, mass,
                                      rho.data[:rho.nlat-rho.extend,
                                               :rho.nlon-rho.extend],
                                      lmax=lmax)

        else:
            cilm, d = _CilmPlusDH(shape.data[:shape.nlat-shape.extend,
                                             :shape.nlon-shape.extend],
                                  nmax, mass, rho, lmax=lmax)

        clm = SHGravRealCoeffs(cilm, gm=gm, r0=d, omega=omega,
                               normalization='4pi', csphase=1, name=name,
                               epoch=epoch)
        return clm

    @property
    def mass(self):
        """Return the mass of the planet in kg.
        """
        return self.gm / _G.value

    @property
    def center_of_mass(self):
        """
        Return the Cartesian coordinates of the center of mass of the planet
        in meters.

        Returns
        -------
        [x, y, z] : numpy ndarray
            Cartesian coordinates of the center of mass in meters.
        """
        coeffs = self.convert(normalization='unnorm', csphase=1, lmax=1).coeffs

        x_cm = coeffs[0, 1, 1] * self.r0
        y_cm = coeffs[1, 1, 1] * self.r0
        z_cm = coeffs[0, 1, 0] * self.r0

        return _np.array([x_cm, y_cm, z_cm])

    def inertia_tensor(self, dynamical_flattening):
        """Return the inertia tensor of the planet in kg * m**2.

        Parameters
        ----------
        dynamical_flattening : float
            Dynamical flattening (or precession constant) of the planet,
            defined as [C-(A+B)/2]/C.

        Returns
        -------
        tensor : ndarray, shape (3, 3)
            Inertia tensor of the planet.

        Notes
        -----
        The moment of inertia tensor is given by 9 components

            (Ixx, Ixy, Ixz)
            (Iyx, Iyy, Iyz)
            (Izx, Izy, Izz)

        The diagonal elements Ixx, Iyy, Izz are the axial moments of inertia,
        and the off-diagonal elements

            Ixy = Iyx, Ixz = Izx, Iyz = Izy

        are the products of inertia.

        References
        ----------
        Heiskanen, W.A. and Moritz, H., 1967. Physical geodesy. San Francisco,
        WH Freeman, 1967.

        Chen, W., Li, J.C., Ray, J., Shen, W.B. and Huang, C.L.,
        Consistent estimates of the dynamic figure parameters of the Earth.
        J. Geod., 89(2), 179-188, 2015.
        """

        coeffs = self.convert(normalization='unnorm', csphase=1, lmax=2).coeffs

        mr02 = self.mass * self.r0**2

        # Products of inertia
        yz = -mr02 * coeffs[1, 2, 1]
        xz = -mr02 * coeffs[0, 2, 1]
        xy = -2 * mr02 * coeffs[1, 2, 2]

        # Axial moments of inertia
        xx = mr02 * ((1 - 1 / dynamical_flattening) * coeffs[0, 2, 0] -
                     2 * coeffs[0, 2, 2])
        yy = mr02 * ((1 - 1 / dynamical_flattening) * coeffs[0, 2, 0] +
                     2 * coeffs[0, 2, 2])
        zz = -mr02 * coeffs[0, 2, 0] / dynamical_flattening

        tensor = _np.array([
            [xx, xy, xz],
            [xy, yy, yz],
            [xz, yz, zz]])

        return tensor

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

        mneg_mask = (ms < 0).astype(_np.int_)
        self.coeffs[mneg_mask, ls, _np.abs(ms)] = values

    # ---- IO routines ----
    def to_file(self, filename, format='shtools', header=None, errors=True,
                lmax=None, modelname=None, tide_system='unknown',
                encoding=None, **kwargs):
        """
        Save spherical harmonic coefficients to a file.

        Usage
        -----
        x.to_file(filename, [format='shtools', header, errors, lmax, encoding])
        x.to_file(filename, format='dov', [header, errors, lmax, encoding])
        x.to_file(filename, format='bshc', [lmax])
        x.to_file(filename, format='icgem', [header, errors, lmax, modelname,
                                             tide_system, encoding])
        x.to_file(filename, format='npy', [**kwargs])

        Parameters
        ----------
        filename : str
            Name of the output file. If the filename ends with '.gz', the file
            will be compressed using gzip.
        format : str, optional, default = 'shtools'
            'shtools', 'dov', 'bshc', 'icgem' or 'npy'.
        header : str, optional, default = None
            A header string written to an 'shtools' or 'dov'-formatted file
            directly before the metadata and spherical harmonic coefficients.
        errors : bool, optional, default = True
            If True, save the errors in the file (for 'shtools', 'dov', and
            'icgem' formatted files only).
        lmax : int, optional, default = self.lmax
            The maximum spherical harmonic degree to write to the file.
        modelname : str, optional, default = None
            The name of the model for 'icgem' formatted files.
        tide_system : str, optional, default = 'unknown'
            The tide system for 'icgem' formatted files: 'zero_tide',
            'tide_free', or 'unknown'.
        encoding : str, optional, default = None
            Encoding of the output file when format is 'shtools', 'dov' or
            'icgem'. The default is to use the system default.
        **kwargs : keyword argument list, optional for format = 'npy'
            Keyword arguments of numpy.save().

        Notes
        -----
        Supported file formats:
            'shtools' (see pyshtools.shio.shwrite)
            'dov' (see pyshtools.shio.write_dov)
            'bshc' (see pyshtools.shio.write_bshc)
            'icgem' (see pyshtools.shio.write_icgem_gfc)
            'npy' (see numpy.save)

        If the filename end with '.gz', the file will be compressed using gzip.

        'shtools': The coefficients and meta-data will be written to an ascii
        formatted file. The first line is an optional user provided header
        line, and the following line provides the attributes r0, gm,
        omega, and lmax. The spherical harmonic coefficients (and optionally
        the errors) are then listed, with increasing degree and order, with the
        format

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

        'icgem': The coefficients will be written to a text file using the
        gfc format of the International Centre for Global Earth Models.

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
            if self.omega is None:
                omega = 0.
            else:
                omega = self.omega

            header_str = '{:.16e}, {:.16e}, {:.16e}, {:d}'.format(
                    self.r0, self.gm, omega, lmax)
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

        elif format.lower() == 'icgem':
            if errors:
                _write_icgem_gfc(filebase, self.coeffs, errors=self.errors,
                                 header=header, lmax=lmax, modelname=modelname,
                                 gm=self.gm, r0=self.r0,
                                 error_kind=self.error_kind,
                                 tide_system=tide_system,
                                 normalization=self.normalization,
                                 encoding=encoding)
            else:
                _write_icgem_gfc(filebase, self.coeffs, header=header,
                                 lmax=lmax, modelname=modelname,
                                 gm=self.gm, r0=self.r0,
                                 tide_system=tide_system,
                                 normalization=self.normalization,
                                 encoding=encoding)

        elif format.lower() == 'npy':
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
        ds['coeffs'].attrs['GM'] = self.gm
        ds['coeffs'].attrs['r0'] = self.r0
        if self.omega is not None:
            ds['coeffs'].attrs['omega'] = self.omega
        if self.epoch is not None:
            ds['coeffs'].attrs['epoch'] = self.epoch

        if self.errors is not None:
            cerrors = self.errors[0, :lmax+1, :lmax+1]
            serrors = _np.transpose(self.errors[1, :lmax+1, :lmax+1])
            serrors = _np.vstack([serrors[1:], serrors[0]])
            ds['errors'] = (('degree', 'order'), cerrors + serrors)
            ds['errors'].attrs['normalization'] = self.normalization
            ds['errors'].attrs['csphase'] = self.csphase
            ds['errors'].attrs['GM'] = self.gm
            ds['errors'].attrs['r0'] = self.r0
            if self.omega is not None:
                ds['errors'].attrs['omega'] = self.omega
            if self.epoch is not None:
                ds['errors'].attrs['epoch'] = self.epoch
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
        Print a summary of the data stored in the SHGravCoeffs class instance.

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
                            'instances. Type of other is {:s}.'
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
                            'SHGravCoeffs instances. Type of other is {:s}.'
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
                            'SHGravCoeffs instances. Type of other is {:s}.'
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
                            'Type of other is {:s}.'.format(repr(type(other))))

    def __rmul__(self, other):
        """
        Multiply an SHGravCoeffs instance by an SHCoeffs instance or scalar:
        other * self.
        """
        return self.__mul__(other)

    def __truediv__(self, other):
        """
        Divide an SHGravCoeffs instance by an SHCoeffs instance or scalar:
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

        Notes
        -----
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
                "'total'. Provided value is {:s}.".format(repr(function))
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

    def admittance(self, hlm, errors=True, function='radial', lmax=None):
        """
        Return the admittance for an input topography function.

        Usage
        -----
        admittance = g.admittance(hlm, [errors, function, lmax])

        Returns
        -------
        admittance : ndarray, shape (lmax+1) or (2, lmax+1)
            1-D array of the admittance (errors=False) or 2-D array of the
            admittance and its uncertainty (errors=True), where lmax is the
            maximum spherical harmonic degree.

        Parameters
        ----------
        hlm : SHCoeffs class instance.
            The topography function h used in computing the admittance
            Sgh / Shh. hlm is assumed to have units of meters.
        errors : bool, optional, default = True
            Return the uncertainty of the admittance.
        function : str, optional, default = 'radial'
            The type of admittance to return: 'geoid' for using the geoid, in
            units of m/km, or 'radial' for using the radial gravity in units
            of mGal/km.
        lmax : int, optional, default = g.lmax
            Maximum spherical harmonic degree of the spectrum to output.

        Notes
        -----
        If gravity g and topography h are related by the equation

            glm = Z(l) hlm + nlm

        where nlm is a zero-mean random variable, the admittance Z(l) can be
        estimated using

            Z(l) = Sgh(l) / Shh(l),

        where Sgh, Shh and Sgg are the cross-power and power spectra of the
        functions g (self) and h (input).
        """
        if not isinstance(hlm, _SHCoeffs):
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

        if function == 'geoid':
            admit *= 1000. * self.r0
        else:
            degrees = _np.arange(lmax+1)
            if errors:
                admit *= 1.e8 * self.gm * (degrees[:, None] + 1) / self.r0**2
            else:
                admit *= 1.e8 * self.gm * (degrees + 1) / self.r0**2
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
        from .shmagcoeffs import SHMagCoeffs as _SHMagCoeffs
        if not isinstance(hlm, (_SHCoeffs, _SHMagCoeffs, SHGravCoeffs)):
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

    def admitcorr(self, hlm, errors=True, function='radial', lmax=None):
        """
        Return the admittance and correlation for an input topography function.

        Usage
        -----
        admittance, correlation = g.admitcorr(hlm, [errors, function, lmax])

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
            The topography function h used in computing the admittance
            Sgh / Shh. hlm is assumed to have units of meters.
        errors : bool, optional, default = True
            Return the uncertainty of the admittance.
        function : str, optional, default = 'radial'
            The type of admittance to return: 'geoid' for using the geoid, in
            units of m/km, or 'radial' for using the radial gravity in units
            of mGal/km.
        lmax : int, optional, default = g.lmax
            Maximum spherical harmonic degree of the spectrum to output.

        Notes
        -----
        If gravity g and topography h are related by the equation

            glm = Z(l) hlm + nlm

        where nlm is a zero-mean random variable, the admittance Z(l) and
        spectral correlation can be estimated using

            Z(l) = Sgh(l) / Shh(l)
            gamma(l) = Sgh(l) / sqrt( Sgg(l) Shh(l) )

        where Sgh, Shh and Sgg are the cross-power and power spectra of the
        functions g (self) and h (input).
        """
        if not isinstance(hlm, _SHCoeffs):
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
        sgg = _spectrum(self.coeffs, normalization=self.normalization,
                        lmax=lmax)

        with _np.errstate(invalid='ignore', divide='ignore'):
            admit = sgh / shh
            corr = sgh / _np.sqrt(sgg * shh)
            if errors:
                sigma = (sgg / shh) * (1. - corr**2) / _np.arange(lmax+1) / 2.
                admit = _np.column_stack((admit, _np.sqrt(sigma)))

        if function == 'geoid':
            admit *= 1000. * self.r0
        else:
            degrees = _np.arange(lmax+1)
            if errors:
                admit *= 1.e8 * self.gm * (degrees[:, None] + 1) / self.r0**2
            else:
                admit *= 1.e8 * self.gm * (degrees + 1) / self.r0**2

        return admit, corr

    # ---- Operations that return a new SHGravCoeffs class instance ----
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
            raise ValueError('convention must be a string. Input type is {:s}.'
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

        rot = self._rotate(angles, dj_matrix, gm=self.gm, r0=self.r0,
                           omega=self.omega, backend=backend,
                           nthreads=nthreads)
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
            return SHGravCoeffs.from_array(
                coeffs, gm=self.gm, r0=self.r0, omega=self.omega,
                errors=errors, error_kind=self.error_kind,
                normalization=normalization.lower(),
                csphase=csphase, epoch=self.epoch, copy=False)
        else:
            coeffs = self.to_array(normalization=normalization.lower(),
                                   csphase=csphase, lmax=lmax)
            return SHGravCoeffs.from_array(
                coeffs, gm=self.gm, r0=self.r0, omega=self.omega,
                normalization=normalization.lower(), csphase=csphase,
                epoch=self.epoch, copy=False)

    def pad(self, lmax, copy=True):
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

        Notes
        -----
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
    def expand(self, a=None, f=None, colat=None, lat=None, lon=None,
               degrees=True, lmax=None, lmax_calc=None, normal_gravity=True,
               sampling=2, extend=True):
        """
        Create 2D cylindrical maps on a flattened and rotating ellipsoid of the
        three components of the gravity vector, the gravity disturbance, and
        the gravity potential. Alternatively, compute the gravity vector at
        specified coordinates.

        Usage
        -----
        grids = x.expand([a, f, lmax, lmax_calc, normal_gravity, sampling,
                          extend])
        g = x.expand(lat, lon, [a, f, lmax, lmax_calc, degrees])
        g = x.expand(colat, lon, [a, f, lmax, lmax_calc, degrees])

        Returns
        -------
        grids : SHGravGrid class instance.
        g     : (r, theta, phi) components of the gravity vector at the
                specified points.

        Parameters
        ----------
        a : optional, float, default = self.r0
            The semi-major axis of the flattened ellipsoid on which the field
            is computed.
        f : optional, float, default = 0
            The flattening of the reference ellipsoid: f=(a-b)/a.
        lat : int, float, ndarray, or list, optional, default = None
            Latitude coordinates where the gravity is to be evaluated.
        colat : int, float, ndarray, or list, optional, default = None
            Colatitude coordinates where the gravity is to be evaluated.
        lon : int, float, ndarray, or list, optional, default = None
            Longitude coordinates where the gravity is to be evaluated.
        degrees : bool, optional, default = True
            True if lat, colat and lon are in degrees, False if in radians.
        lmax : optional, integer, default = self.lmax
            The maximum spherical harmonic degree, which determines the number
            of samples of the output grids, n=2lmax+2, and the latitudinal
            sampling interval, 90/(lmax+1).
        lmax_calc : optional, integer, default = lmax
            The maximum spherical harmonic degree used in evaluating the
            functions. This must be less than or equal to lmax.
        normal_gravity : optional, bool, default = True
            If True (and if a, f and x.omega are set explicitly), the normal
            gravity (the gravitational acceleration on the rotating ellipsoid)
            will be subtracted from the total gravitational acceleration,
            yielding the "gravity disturbance." This is done using Somigliana's
            formula (after converting geocentric to geodetic coordinates).
        sampling : optional, integer, default = 2
            If 1 the output grids are equally sampled (n by n). If 2 (default),
            the grids are equally spaced in degrees.
        extend : bool, optional, default = True
            If True, compute the longitudinal band for 360 E and the
            latitudinal band for 90 S.

        Notes
        -----
        This method will create 2-dimensional cylindrical maps of the three
        components of the gravity vector (gravitational force + centrifugal
        force), the magnitude of the gravity vector, and the gravity
        potential, and return these as an SHGravGrid class instance. Each map
        is stored as an SHGrid class instance using Driscoll and Healy grids
        that are either equally sampled (n by n) or equally spaced in degrees
        latitude and longitude. All grids use geocentric coordinates, the
        output is in SI units, and the sign of the radial components is
        positive when directed upwards. If latitude and longitude coordinates
        are specified, this method will instead return the gravity vector.

        If the angular rotation rate omega is specified in the SHGravCoeffs
        instance, both the potential and gravity vectors will be calculated in
        a body-fixed rotating reference frame and will include the contribution
        from the centrifugal force. If normal_gravity is set to True, and a, f,
        and omega are all set explicitly, the normal gravity will be removed
        from the magnitude of the gravity vector, yielding the gravity
        disturbance.

        The gravitational potential is given by

            V = GM/r Sum_{l=0}^lmax (r0/r)^l Sum_{m=-l}^l C_{lm} Y_{lm},

        and the gravitational acceleration is

            B = Grad V.

        The coefficients are referenced to the radius r0, and the function is
        computed on a flattened ellipsoid with semi-major axis a (i.e., the
        mean equatorial radius) and flattening f. To convert m/s^2 to mGals,
        multiply the gravity grids by 10^5.
        """
        if lat is not None and colat is not None:
            raise ValueError('lat and colat can not both be specified.')

        if a is not None and f is not None and self.omega is not None \
                and normal_gravity is True:
            ng = 1
        else:
            ng = 0
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
                                        degrees=degrees, lmax_calc=lmax_calc,
                                        omega=self.omega)
            return values

        else:
            if lmax is None:
                lmax = self.lmax
            if lmax_calc is None:
                lmax_calc = lmax

            coeffs = self.to_array(normalization='4pi', csphase=1,
                                   errors=False)
            rad, theta, phi, total, pot = _MakeGravGridDH(
                coeffs, self.gm, self.r0, a=a, f=f, lmax=lmax,
                lmax_calc=lmax_calc, sampling=sampling, omega=self.omega,
                normal_gravity=ng, extend=extend)

            return _SHGravGrid(rad, theta, phi, total, pot, self.gm, a, f,
                               self.omega, normal_gravity, lmax, lmax_calc,
                               units='m/s2', pot_units='m2/s2',
                               epoch=self.epoch)

    def tensor(self, a=None, f=None, lmax=None, lmax_calc=None, degree0=False,
               sampling=2, extend=True):
        """
        Create 2D cylindrical maps on a flattened ellipsoid of the 9
        components of the gravity "gradient" tensor in a local north-oriented
        reference frame, and return an SHGravTensor class instance.

        Usage
        -----
        tensor = x.tensor([a, f, lmax, lmax_calc, sampling, extend])

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
            the grids are equally spaced in degrees.
        extend : bool, optional, default = True
            If True, compute the longitudinal band for 360 E and the
            latitudinal band for 90 S.

        Notes
        -----
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

        coeffs = self.to_array(normalization='4pi', csphase=1, errors=False)

        if degree0 is False:
            coeffs[0, 0, 0] = 0.

        vxx, vyy, vzz, vxy, vxz, vyz = _MakeGravGradGridDH(
            coeffs, self.gm, self.r0, a=a, f=f, lmax=lmax,
            lmax_calc=lmax_calc, sampling=sampling, extend=extend)

        return _SHGravTensor(1.e9*vxx, 1.e9*vyy, 1.e9*vzz, 1.e9*vxy, 1.e9*vxz,
                             1.e9*vyz, self.gm, a, f, lmax, lmax_calc,
                             units='Eötvös', epoch=self.epoch)

    def geoid(self, potref, a=None, f=None, r=None, omega=None, order=2,
              lmax=None, lmax_calc=None, grid='DH2', extend=True):
        """
        Create a global map of the height of the geoid and return an SHGeoid
        class instance.

        Usage
        -----
        geoid = x.geoid(potref, [a, f, r, omega, order, lmax, lmax_calc, grid,
                                 extend])

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
            'DH' or 'DH1' for an equally sampled grid with nlat=nlon, or
            'DH2' for an equally spaced grid in degrees latitude and longitude.
        extend : bool, optional, default = True
            If True, compute the longitudinal band for 360 E and the
            latitudinal band for 90 S.

        Notes
        -----
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
                    "Input value is {:s}.".format(repr(grid)))

        coeffs = self.to_array(normalization='4pi', csphase=1, errors=False)

        if omega is None:
            omega = self.omega

        geoid = _MakeGeoidGridDH(coeffs, self.r0, self.gm, potref, lmax=lmax,
                                 omega=omega, r=r, order=order,
                                 lmax_calc=lmax_calc, a=a, f=f,
                                 sampling=sampling, extend=extend)

        return _SHGeoid(geoid, self.gm, potref, a, f, omega, r, order,
                        lmax, lmax_calc, units='m', epoch=self.epoch)

    # ---- Plotting routines ----
    def plot_spectrum(self, function='geoid', unit='per_l', base=10.,
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
            axes.set_xscale('log', base=base)
        if yscale == 'log':
            axes.set_yscale('log', base=base)

        if self.errors is not None:
            axes.plot(ls[2:lmax + 1], spectrum[2:lmax + 1], label=legend,
                      **kwargs)
            axes.plot(ls[2:lmax + 1], error_spectrum[2:lmax + 1],
                      label=legend_error, **kwargs)
        else:
            axes.plot(ls[2:lmax + 1], spectrum[2: lmax + 1], label=legend,
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

    def plot_spectrum2d(self, function='geoid', ticks='WSen',
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
        x.plot_spectrum2d([function, ticks, tick_interval, minor_tick_interval,
                           degree_label, order_label, title, colorbar, origin,
                           cmap, cmap_limits, cmap_rlimits, cmap_reverse,
                           cmap_scale, cb_triangles, cb_label, cb_offset,
                           cb_width, lmax, errors, xscale, yscale, grid,
                           titlesize, axes_labelsize, tick_labelsize, ax,
                           show, fname])

        Parameters
        ----------
        function : str, optional, default = 'geoid'
            The type of power spectrum to calculate: 'potential' for the
            gravitational potential, 'geoid' for the geoid, 'radial' for
            the radial gravity, or 'total' for the total gravitational field.
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
        the parameter function: 'potential' for the gravitational potential,
        'geoid' for the geoid, 'radial' for the radial gravity, or 'total' for
        the total gravitational field. In all cases, the total power of the
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
                if function == 'geoid':
                    cb_label = 'Power, m$^2$'
                elif function == 'potential':
                    cb_label = 'Power, m$^4$ s$^{-4}$'
                elif function == 'radial':
                    cb_label = 'Power, mGal$^2$'
                elif function == 'total':
                    cb_label = 'Power, mGal$^2$'

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

    def plot_admitcorr(self, hlm, errors=True, function='radial',
                       style='separate', lmax=None, grid=True, legend=None,
                       legend_loc='best', axes_labelsize=None,
                       tick_labelsize=None, elinewidth=0.75, ax=None, ax2=None,
                       show=True, fname=None, **kwargs):
        """
        Plot the admittance and/or correlation with another function.

        Usage
        -----
        x.plot_admitcorr(hlm, [errors, function, style, lmax, grid, legend,
                               legend_loc, axes_labelsize, tick_labelsize,
                               elinewidth, ax, ax2, show, fname, **kwargs])

        Parameters
        ----------
        hlm : SHCoeffs class instance.
            The second function used in computing the admittance and
            correlation.
        errors : bool, optional, default = True
            Plot the uncertainty of the admittance.
        function : str, optional, default = 'radial'
            The type of admittance to return: 'geoid' for using the geoid, in
            units of m/km, or 'radial' for using the radial gravity in units
            of mGal/km.
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
        ax : matplotlib axes object, optional, default = None
            A single matplotlib axes object where the plot will appear.
        ax2 : matplotlib axes object, optional, default = None
            A single matplotlib axes object where the second plot will appear
            when style is 'separate'.
        show : bool, optional, default = True
            If True, plot to the screen.
        fname : str, optional, default = None
            If present, and if ax is not specified, save the image to the
            specified file.
        **kwargs : keyword arguments, optional
            Keyword arguments for pyplot.plot() and pyplot.errorbar().

        Notes
        -----
        If gravity g and topography h are related by the equation

            glm = Z(l) hlm + nlm

        where nlm is a zero-mean random variable, the admittance and spectral
        correlation gamma(l) can be estimated using

            Z(l) = Sgh(l) / Shh(l)
            gamma(l) = Sgh(l) / sqrt( Sgg(l) Shh(l) )

        where Sgh, Shh and Sgg are the cross-power and power spectra of g
        (self) and h (input).
        """
        if lmax is None:
            lmax = min(self.lmax, hlm.lmax)

        if style in ('combined', 'separate'):
            admit, corr = self.admitcorr(hlm, errors=errors, function=function,
                                         lmax=lmax)
        elif style == 'corr':
            corr = self.correlation(hlm, lmax=lmax)
        elif style == 'admit':
            admit = self.admittance(hlm, errors=errors, function=function,
                                    lmax=lmax)
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
            if function == 'radial':
                admitax.set_ylabel('Admittance, mGal/km',
                                   fontsize=axes_labelsize)
            else:
                admitax.set_ylabel('Admittance, m/km',
                                   fontsize=axes_labelsize)
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

    def plot_admittance(self, hlm, errors=True, function='radial',
                        lmax=None, grid=True, legend=None,
                        legend_loc='best', axes_labelsize=None,
                        tick_labelsize=None, elinewidth=0.75, ax=None,
                        show=True, fname=None, **kwargs):
        """
        Plot the admittance with another function.

        Usage
        -----
        x.plot_admittance(hlm, [errors, function, lmax, grid, legend,
                                legend_loc, axes_labelsize, tick_labelsize,
                                elinewidth, ax, show, fname, **kwargs])

        Parameters
        ----------
        hlm : SHCoeffs class instance.
            The second function used in computing the admittance.
        errors : bool, optional, default = True
            Plot the uncertainty of the admittance.
        function : str, optional, default = 'radial'
            The type of admittance to return: 'geoid' for using the geoid, in
            units of m/km, or 'radial' for using the radial gravity in units
            of mGal/km.
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
        If gravity g and topography h are related by the equation

            glm = Z(l) hlm + nlm

        where nlm is a zero-mean random variable, the admittance can be
        estimated using

            Z(l) = Sgh(l) / Shh(l)

        where Sgh and Shh are the cross-power and power spectra of the
        g (self) and h (input).
        """
        return self.plot_admitcorr(hlm, errors=errors, function=function,
                                   style='admit', lmax=lmax, grid=grid,
                                   legend=legend, legend_loc=legend_loc,
                                   axes_labelsize=axes_labelsize,
                                   tick_labelsize=tick_labelsize,
                                   elinewidth=elinewidth, show=True,
                                   fname=fname, ax=ax, **kwargs)

    def plot_correlation(self, hlm, lmax=None, grid=True, legend=None,
                         legend_loc='best', axes_labelsize=None,
                         tick_labelsize=None, elinewidth=0.75, ax=None,
                         show=True, fname=None, **kwargs):
        """
        Plot the correlation with another function.

        Usage
        -----
        x.plot_correlation(hlm, [lmax, grid, legend, legend_loc,
                                 axes_labelsize, tick_labelsize, elinewidth,
                                 ax, show, fname, **kwargs])

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
        return self.plot_admitcorr(hlm, style='corr', lmax=lmax, grid=grid,
                                   legend=legend, legend_loc=legend_loc,
                                   axes_labelsize=axes_labelsize,
                                   tick_labelsize=tick_labelsize,
                                   show=True, fname=fname, ax=ax, **kwargs)


class SHGravRealCoeffs(SHGravCoeffs):
    """
    Real spherical harmonic coefficient class for the gravitational potential.
    """

    def __init__(self, coeffs, gm=None, r0=None, omega=None, errors=None,
                 error_kind=None, normalization='4pi', csphase=1, copy=True,
                 header=None, header2=None, name=None, epoch=None):
        """Initialize real gravitational potential coefficients class."""
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
        self.gm = gm
        self.r0 = r0
        self.omega = omega
        self.name = name
        self.epoch = epoch
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
                'GM (m3 / s2) = {:s}\n'
                'r0 (m) = {:s}\n'
                'Omega (rad / s) = {:s}\n'
                'error_kind = {:s}\n'
                'header = {:s}\n'
                'header2 = {:s}\n'
                'name = {:s}\n'
                'epoch = {:s}'
                .format(repr(self.kind), repr(self.normalization),
                        self.csphase, self.lmax, repr(self.gm), repr(self.r0),
                        repr(self.omega), repr(self.error_kind),
                        repr(self.header), repr(self.header2),
                        repr(self.name), repr(self.epoch)))

    def _rotate(self, angles, dj_matrix, gm=None, r0=None, omega=None,
                backend=None, nthreads=None):
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
            return SHGravCoeffs.from_array(
                temp, errors=self.errors, normalization=self.normalization,
                csphase=self.csphase, copy=False, gm=gm, r0=r0, omega=omega,
                epoch=self.epoch)
        else:
            return SHGravCoeffs.from_array(coeffs, errors=self.errors,
                                           gm=gm, r0=r0, omega=omega,
                                           epoch=self.epoch, copy=False)

    def _expand_coord(self, a, f, lat, lon, degrees, lmax_calc, omega):
        """Evaluate the gravity at the coordinates lat and lon."""
        coeffs = self.to_array(normalization='4pi', csphase=1, errors=False)

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

            return _MakeGravGridPoint(coeffs, gm=self.gm, r0=self.r0,
                                      r=r, lat=latin, lon=lonin,
                                      lmax=lmax_calc, omega=self.omega)
        elif type(lat) is _np.ndarray:
            values = _np.empty((len(lat), 3), dtype=_np.float64)
            for i, (latitude, longitude) in enumerate(zip(latin, lonin)):
                if f == 0.:
                    r = a
                else:
                    r = _np.cos(_np.deg2rad(latitude))**2 + \
                        _np.sin(_np.deg2rad(latitude))**2 / (1.0 - f)**2
                    r = a * _np.sqrt(1. / r)

                values[i, :] = _MakeGravGridPoint(coeffs, gm=self.gm,
                                                  r0=self.r0, r=r,
                                                  lat=latitude, lon=longitude,
                                                  lmax=lmax_calc,
                                                  omega=self.omega)
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
                    _MakeGravGridPoint(coeffs, gm=self.gm, r0=self.r0,
                                       r=r, lat=latitude, lon=longitude,
                                       lmax=lmax_calc, omega=self.omega))
            return values
        else:
            raise ValueError('lat and lon must be either an int, float, '
                             'ndarray, or list. Input types are {:s} and {:s}.'
                             .format(repr(type(lat)), repr(type(lon))))
