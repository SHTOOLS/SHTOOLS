"""
    Class for spherical harmonic coefficients of gravitational fields.

    The class SHGravCoeffs is used solely to initialize an instance of the
    class SHGravRealCoeffs, which is a subclass of SHRealCoeffs.
"""
from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

import numpy as _np
import matplotlib as _mpl
import matplotlib.pyplot as _plt

from .shcoeffsgrid import SHCoeffs as _SHCoeffs
from .shcoeffsgrid import SHRealCoeffs as _SHRealCoeffs
from .shcoeffsgrid import DHRealGrid as _DHRealGrid
from ..shtools import CilmPlusRhoHDH as _CilmPlusRhoHDH
from ..shtools import CilmPlusDH as _CilmPlusDH
from ..spectralanalysis import spectrum as _spectrum
from ..shio import shread as _shread
from ..constant import grav_constant as _G


# =============================================================================
# =========    SHGravCoeffs class    =========================================
# =============================================================================

class SHGravCoeffs(object):
    """
    Spherical harmonic coefficients class for the gravitational potential.

    The coefficients of this class can be initialized using one of the
    four constructor methods:

        x = SHGravCoeffs.from_array(numpy.zeros((2, lmax+1, lmax+1)))
        x = SHGravCoeffs.from_random(powerspectrum[0:lmax+1])
        x = SHGravCoeffs.from_zeros(lmax)
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
    gm            : The mass times the gravitational constant that is
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
    kind          : The coefficient data type: either 'complex' or 'real'.
    header        : A list of values (of type str) from the header line of the
                    input file used to initialize the class (for 'shtools'
                    formatted files).

    Each class instance provides the following methods:

    set_gm()              : Set the gm of the class instance.
    set_r0()              : Set the reference radius of the class instance.
    set_omega()           : Set the angular rotation rate of the body.
    to_array()            : Return an array of spherical harmonic coefficients
                            with a different normalization convention.
    to_file()             : Save raw spherical harmonic coefficients as a file.
    degrees()             : Return an array listing the spherical harmonic
                            degrees from 0 to lmax.
    spectrum()            : Return the spectrum of the function as a function
                            of spherical harmonic degree.
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
    copy()                : Return a copy of the class instance.
    plot_spectrum()       : Plot the  spectrum as a function of spherical
                            harmonic degree.
    plot_spectrum2d()     : Plot the 2D spectrum of all spherical harmonic
                            coefficients.
    info()                : Print a summary of the data stored in the SHCoeffs
                            instance.
    """

    def __init__(self):
        """Unused constructor of the super class."""
        print('Initialize the class using one of the class methods:\n'
              '>>> SHGravCoeffs.from_array\n'
              '>>> SHGravCoeffs.from_random\n'
              '>>> SHGravCoeffs.from_zeros\n'
              '>>> SHGravCoeffs.from_file\n'
              '>>> SHGravCoeffs.from_shape\n')

    # ---- Factory methods ----
    @classmethod
    def from_array(self, coeffs, gm, r0, omega=None,
                   normalization='4pi', csphase=1, lmax=None, copy=True):
        """
        Initialize the class with spherical harmonic coefficients from an input
        array.

        Usage
        -----
        x = SHGravCoeffs.from_array(array, gm, r0, [omega, normalization,
                                                    csphase, lmax, copy])

        Returns
        -------
        x : SHGravCoeffs class instance.

        Parameters
        ----------
        array : ndarray, shape (2, lmaxin+1, lmaxin+1).
            The input spherical harmonic coefficients.
        gm : float
            The mass times the gravitational constant that is associated with
            the gravitational potential coefficients.
        r0 : float
            The reference radius of the spherical harmonic coefficients.
        omega : float, optional, default = None
            The angular rotation rate of the body.
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
        """
        if _np.iscomplexobj(coeffs):
            raise TypeError('The input array must be real.')

        if type(normalization) != str:
            raise ValueError('normalization must be a string. ' +
                             'Input type was {:s}'
                             .format(str(type(normalization))))

        if normalization.lower() not in ('4pi', 'ortho', 'schmidt', 'unnorm'):
            raise ValueError(
                "The normalization must be '4pi', 'ortho', 'schmidt', " +
                "or 'unnorm'. Input value was {:s}."
                .format(repr(normalization))
                )

        if csphase != 1 and csphase != -1:
            raise ValueError(
                "csphase must be either 1 or -1. Input value was {:s}."
                .format(repr(csphase))
                )

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
                           "85. Input value was {:d}.".format(lmax),
                           category=RuntimeWarning)
            lmax = 85

        clm = SHGravRealCoeffs(coeffs[:, 0:lmax+1, 0:lmax+1],
                               normalization=normalization.lower(),
                               csphase=csphase, copy=copy, gm=gm, r0=r0,
                               omega=omega)
        return clm

    @classmethod
    def from_zeros(self, lmax, gm, r0, omega=None,
                   normalization='4pi', csphase=1):
        """
        Initialize the class with spherical harmonic coefficients set to zero
        from degree 0 to lmax.

        Usage
        -----
        x = SHGravCoeffs.from_zeros(lmax, gm, r0, [omega, normalization,
                                                   csphase])

        Returns
        -------
        x : SHGravCoeffs class instance.

        Parameters
        ----------
        lmax : int
            The highest spherical harmonic degree l of the coefficients.
        gm : float
            The mass times the gravitational constant that is associated with
            the gravitational potential coefficients.
        r0 : float
            The reference radius of the spherical harmonic coefficients.
        omega : float, optional, default = None
            The angular rotation rate of the body.
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
                "The normalization must be '4pi', 'ortho', 'schmidt', " +
                "or 'unnorm'. Input value was {:s}."
                .format(repr(normalization))
                )

        if csphase != 1 and csphase != -1:
            raise ValueError(
                "csphase must be either 1 or -1. Input value was {:s}."
                .format(repr(csphase))
                )

        if normalization.lower() == 'unnorm' and lmax > 85:
            _warnings.warn("Calculations using unnormalized coefficients " +
                           "are stable only for degrees less than or equal " +
                           "to 85. lmax for the coefficients will be set to " +
                           "85. Input value was {:d}.".format(lmax),
                           category=RuntimeWarning)
            lmax = 85

        coeffs = _np.zeros((2, lmax + 1, lmax + 1))
        clm = SHGravRealCoeffs(coeffs, normalization=normalization.lower(),
                               csphase=csphase, gm=gm, r0=r0, omega=omega)
        return clm

    @classmethod
    def from_file(self, fname, lmax=None, format='shtools',
                  normalization='4pi', skip=0, header=True, csphase=1,
                  r0_index=0, gm_index=1, omega_index=None, gm=None, r0=None,
                  omega=None, **kwargs):
        """
        Initialize the class with spherical harmonic coefficients from a file.

        Usage
        -----
        x = SHGravCoeffs.from_file(filename, [format='shtools', lmax,
                                              normalization, csphase, skip,
                                              header, gm_index, r0_index,
                                              omega_index, gm, r0, omega])
        x = SHGravCoeffs.from_file(filename, [format='npy', normalization,
                                              csphase, gm, r0, omega,
                                              **kwargs])

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
            The mass times the gravitational constant that is associated with
            the gravitational potential coefficients.
        r0 : float, optional, default = None
            The reference radius of the spherical harmonic coefficients.
        omega : float, optional, default = None
            The angular rotation rate of the body.
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
        omega.

        For shtools formatted files, all lines that do not start with 2
        integers and that are less than 3 words long will be treated as
        comments and ignored. For this format, each line of the file must
        contain

        l, m, coeffs[0, l, m], coeffs[1, l, m]

        where l and m are the spherical harmonic degree and order,
        respectively. The terms coeffs[1, l, 0] can be neglected as they are
        zero. For more information, see `shio.shread()`.

        If format='npy', a binary numpy 'npy' file will be read using
        numpy.load().
        """
        if type(normalization) != str:
            raise ValueError('normalization must be a string. ' +
                             'Input type was {:s}'
                             .format(str(type(normalization))))

        if normalization.lower() not in ('4pi', 'ortho', 'schmidt', 'unnorm'):
            raise ValueError(
                "The input normalization must be '4pi', 'ortho', 'schmidt', " +
                "or 'unnorm'. Provided value was {:s}"
                .format(repr(normalization))
                )
        if csphase != 1 and csphase != -1:
            raise ValueError(
                "csphase must be 1 or -1. Input value was {:s}"
                .format(repr(csphase))
                )

        if r0_index is not None and r0 is not None:
            raise ValueError('Can not specify both r0_index and r0')
        if gm_index is not None and gm is not None:
            raise ValueError('Can not specify both gm_index and gm')
        if omega_index is not None and omega is not None:
            raise ValueError('Can not specify both omega_index and omega')

        header_list = None
        if format.lower() == 'shtools':
            if header is True:
                coeffs, lmaxout, header_list = _shread(fname, lmax=lmax,
                                                       skip=skip, header=True)
            else:
                coeffs, lmaxout = _shread(fname, lmax=lmax, skip=skip)
        elif format.lower() == 'npy':
            coeffs = _np.load(fname, **kwargs)
        else:
            raise NotImplementedError(
                'format={:s} not implemented'.format(repr(format)))

        if _np.iscomplexobj(coeffs):
            raise TypeError('The input coefficients must be real.')

        if normalization.lower() == 'unnorm' and lmaxout > 85:
            _warnings.warn("Calculations using unnormalized coefficients " +
                           "are stable only for degrees less than or equal " +
                           "to 85. lmax for the coefficients will be set to " +
                           "85. Input value was {:d}.".format(lmaxout),
                           category=RuntimeWarning)
            lmaxout = 85

        clm = SHGravRealCoeffs(coeffs, normalization=normalization.lower(),
                               csphase=csphase, header=header_list)

        if gm is not None:
            clm.gm = gm
        if r0 is not None:
            clm.r0 = r0
        if omega is not None:
            clm.omega = omega

        if format.lower() == 'shtools' and header is True:
            if r0_index is not None:
                clm.r0 = float(header_list[r0_index])
            if gm_index is not None:
                clm.gm = float(header_list[gm_index])
            if omega_index is not None:
                clm.omega = float(header_list[omega_index])

        return clm

    @classmethod
    def from_random(self, power, gm, r0, omega=None, spectrum='potential',
                    lmax=None, normalization='4pi', csphase=1,
                    exact_power=False):
        """
        Initialize the class of gravitational potential spherical harmonic
        coefficients as random variables.

        Usage
        -----
        x = SHGravCoeffs.from_random(power, gm, r0, [omega, spectrum, lmax,
                                                     normalization,
                                                     csphase, exact_power])

        Returns
        -------
        x : SHGravCoeffs class instance.

        Parameters
        ----------
        power : ndarray, shape (L+1)
            numpy array of shape (L+1) that specifies the expected power per
            degree l of the gravitational potential coefficients, where L is
            the maximum spherical harmonic bandwidth.
        gm : float
            The mass times the gravitational constant that is associated with
            the gravitational potential coefficients.
        r0 : float
            The reference radius of the spherical harmonic coefficients.
        omega : float, optional, default = None
            The angular rotation rate of the body.
        spectrum : str, optional, default = 'potential'
            The type of input power spectrum: 'potential' for the gravitational
            potential, 'geoid' for the geoid, or 'gravity' for the radial
            gravity.
        lmax : int, optional, default = len(power) - 1
            The highest spherical harmonic degree l of the output coefficients.
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
            power. This means that the distribution of power at degree l
            amongst the angular orders is random, but the total power is fixed.

        Description
        -----------
        This routine returns a random realization of spherical harmonic
        gravitational potential coefficients obtained from a normal
        distribution. The variance of each coefficient at degree l is equal to
        the total power at degree l divided by the number of coefficients at
        that degree (2l+1). These coefficients are then divided by a prefactor
        that depends on the type of the input spectrum: gm/r0 for 'potential',
        r0 for 'geoid', and gm/(r**2)*(l+1) for (radial) 'gravity'. The power
        spectrum of the random realization can be fixed exactly to the input
        spectrum using the keyword exact_power.
        """
        if type(normalization) != str:
            raise ValueError('normalization must be a string. ' +
                             'Input type was {:s}'
                             .format(str(type(normalization))))

        if spectrum.lower() not in ('potential', 'geoid', 'gravity'):
            raise ValueError(
                "Spectrum must be of type 'potential', "
                "'geoid', or 'gravity'. Provided value was {:s}"
                .format(repr(spectrum))
                )

        if normalization.lower() not in ('4pi', 'ortho', 'schmidt', 'unnorm'):
            raise ValueError(
                "The input normalization must be '4pi', 'ortho', 'schmidt', " +
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
            _warnings.warn("Calculations using unnormalized coefficients " +
                           "are stable only for degrees less than or equal " +
                           "to 85. lmax for the coefficients will be set to " +
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
                power[0:nl] / (2. * degrees + 1.))[_np.newaxis, :, _np.newaxis]

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

        if spectrum.lower() == 'potential':
            coeffs /= gm / r0
        elif spectrum.lower() == 'geoid':
            coeffs /= r0
        elif spectrum.lower() == 'gravity':
            for l in degrees:
                coeffs[:, l, :l+1] /= gm * (l + 1) / r0**2

        if lmax > nl - 1:
            coeffs = _np.pad(coeffs, ((0, 0), (0, lmax - nl + 1),
                             (0, lmax - nl + 1)), 'constant')

        clm = SHGravRealCoeffs(coeffs, normalization=normalization.lower(),
                               csphase=csphase, gm=gm, r0=r0, omega=omega)
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
            The mass times the gravitational constant that is associated with
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
        Initialize an SHGravCoeffs class instance by calculatin the spherical
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
        performed on a grid that can resolves degrees up to lmax_grid, where
        only the first lmax_calc coefficients are used. It should be emphasized
        that the input shape must correspond to absolute radii as the degree 0
        term determines the reference radius of the coefficients.

        As an intermediate step, this routine calculates the spherical harmonic
        coefficients of the interface referenced raised to the nth power, i.e.,
        (shape-r0)**n, where r0 is the mean radius of shape. If the input shape
        is bandlimited to degree L, the resulting function will thus be
        bandlimited to degree L*nmax. This subroutine assumes implicitly that
        the input shape has an effective spherical harmonic bandwidth
        greater or equal to this value. If this is not the case, aliasing will
        occur. In practice, for accurate results, the effective bandwidth needs
        only to be about three times the size of L, though this should be
        verified for each application.
        """
        mass = gm / _G

        if type(shape) is not _SHRealCoeffs and type(shape) is not _DHRealGrid:
            raise ValueError('shape must be of type SHRealCoeffs '
                             'or DHRealGrid. Input type is {:s}'
                             .format(repr(type(shape))))

        if (type(rho) is not float and type(rho) is not int
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

        clm = SHGravRealCoeffs(cilm, normalization='4pi', csphase=1, gm=gm,
                               r0=d, omega=omega)
        return clm

    # ---- Define gravity specific routines ----
    def set_gm(self, gm):
        """
        Set the value of GM for the class instance.

        Usage
        -----
        x.set_gm(gm)

        Parameters
        ----------
        gm : float
            The value of GM associatiated with the spherical harmonic
            coefficients of the gravitational potential.
        """
        self.gm = gm

    def set_r0(self, r0):
        """
        Set the reference radius of the class instance.

        Usage
        -----
        x.set_r0(r0)

        Parameters
        ----------
        r0 : float
            The reference radius of the gravitational potential spherical
            harmonic coefficients.
        """
        self.r0 = r0

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

    # ---- IO routines ----
    def to_file(self, filename, format='shtools', header=None, **kwargs):
        """
        Save raw spherical harmonic coefficients to a file.

        Usage
        -----
        x.to_file(filename, [format='shtools', header])
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
        **kwargs : keyword argument list, optional for format = 'npy'
            Keyword arguments of numpy.save().

        Description
        -----------
        If format='shtools', the coefficients and meta-data will be written to
        and ascii formatted file. The first line is an optional user provided
        header line, and the following line provides the reference radius, GM,
        angular rotation rate of the body, and the maximum spherical harmonic
        degree. The spherical harmonic coefficients are then listed, with
        increasing degree and order, with the format

        l, m, coeffs[0, l, m], coeffs[1, l, m]

        where l and m are the spherical harmonic degree and order,
        respectively.

        If format='npy', the spherical harmonic coefficients (but not the
        meta-data) will be saved to a binary numpy 'npy' file using
        numpy.save().
        """
        if format is 'shtools':
            with open(filename, mode='w') as file:
                if header is not None:
                    file.write(header + '\n')
                file.write('{:e}, {:e}, {:e}, {:d}'.format(self.r0, self.gm,
                                                           self.omega,
                                                           self.lmax))
                for l in range(self.lmax+1):
                    for m in range(l+1):
                        file.write('{:d}, {:d}, {:e}, {:e}\n'
                                   .format(l, m, self.coeffs[0, l, m],
                                           self.coeffs[1, l, m]))
        elif format is 'npy':
            _np.save(filename, self.coeffs, **kwargs)
        else:
            raise NotImplementedError(
                'format={:s} not implemented'.format(repr(format)))

    def copy(self):
        """Return a deep copy of the class instance."""
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
    # ---- Mathematical operators
    # ----
    # ---- Operations that involve a change of units are not permitted, such as
    # ---- SHGravCoeffs*SHGravCoeffs, SHGravCoeffs/SHGravCoeffs, and
    # ---- SHGravCoeffs+SHCoeffs.
    # -------------------------------------------------------------------------
    def __add__(self, other):
        """
        Add two similar sets of gravitational potential coefficients:
        self + other.
        """
        if isinstance(other, SHGravCoeffs):
            if (self.gm == other.gm and self.r0 == other.r0 and
                    self.normalization == other.normalization and
                    self.csphase == other.csphase and self.kind == other.kind):
                coeffs = _np.empty([2, self.lmax+1, self.lmax+1],
                                   dtype=self.coeffs.dtype)
                coeffs[self.mask] = (self.coeffs[self.mask] +
                                     other.coeffs[self.mask])
                return SHGravCoeffs.from_array(
                    coeffs, gm=self.gm, r0=self.r0, omega=self.omega,
                    csphase=self.csphase, normalization=self.normalization)
            else:
                raise ValueError('Addition is permitted only for two sets of '
                                 'coefficients must be of the same kind, have '
                                 'the same normalization and csphase, and '
                                 'have the same gm and r0.')
        else:
            raise TypeError('Addition is permitted only for two SHGravCoeffs '
                            'instances. Type of other is {:s}'
                            .format(repr(type(other))))

    def __radd__(self, other):
        """
        Add two similar sets of coefficients or coefficients and a scalar:
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
                    self.csphase == other.csphase and self.kind == other.kind):
                coeffs = _np.empty([2, self.lmax+1, self.lmax+1],
                                   dtype=self.coeffs.dtype)
                coeffs[self.mask] = (self.coeffs[self.mask] -
                                     other.coeffs[self.mask])
                return SHGravCoeffs.from_array(
                    coeffs, gm=self.gm, r0=self.r0, omega=self.omega,
                    csphase=self.csphase, normalization=self.normalization)
            else:
                raise ValueError('Subtraction is permitted only for two sets '
                                 'of coefficients must be of the same kind, '
                                 'have the same normalization and csphase, '
                                 'and have the same gm and r0.')
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
                    self.csphase == other.csphase and self.kind == other.kind):
                coeffs = _np.empty([2, self.lmax+1, self.lmax+1],
                                   dtype=self.coeffs.dtype)
                coeffs[self.mask] = (other.coeffs[self.mask] -
                                     self.coeffs[self.mask])
                return SHGravCoeffs.from_array(
                    coeffs, gm=self.gm, r0=self.r0, omega=self.omega,
                    csphase=self.csphase, normalization=self.normalization)
            else:
                raise ValueError('Subtraction is permitted only for two sets '
                                 'of coefficients must be of the same kind, '
                                 'have the same normalization and csphase, '
                                 'and have the same gm and r0.')
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
                    self.csphase == other.csphase and self.kind == other.kind):
                coeffs = _np.empty([2, self.lmax+1, self.lmax+1],
                                   dtype=self.coeffs.dtype)
                coeffs[self.mask] = (self.coeffs[self.mask] *
                                     other.coeffs[self.mask])
                return SHGravCoeffs.from_array(
                    coeffs, gm=self.gm, r0=self.r0, omega=self.omega,
                    csphase=self.csphase, normalization=self.normalization)
            else:
                raise ValueError('The two sets of coefficients must be of ' +
                                 'the same kind and have the same ' +
                                 'normalization and csphase.')
        elif _np.isscalar(other) is True:
            coeffs = _np.empty([2, self.lmax+1, self.lmax+1],
                               dtype=self.coeffs.dtype)
            if self.kind == 'real' and _np.iscomplexobj(other):
                raise ValueError('Can not multiply real coefficients by ' +
                                 'a complex constant.')
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
                    self.csphase == other.csphase and self.kind == other.kind):
                coeffs = _np.empty([2, self.lmax+1, self.lmax+1],
                                   dtype=self.coeffs.dtype)
                coeffs[self.mask] = (self.coeffs[self.mask] /
                                     other.coeffs[self.mask])
                return SHGravCoeffs.from_array(
                    coeffs, gm=self.gm, r0=self.r0, omega=self.omega,
                    csphase=self.csphase, normalization=self.normalization)
            else:
                raise ValueError('The two sets of coefficients must be of ' +
                                 'the same kind and have the same ' +
                                 'normalization and csphase.')
        elif _np.isscalar(other) is True:
            coeffs = _np.empty([2, self.lmax+1, self.lmax+1],
                               dtype=self.coeffs.dtype)
            if self.kind == 'real' and _np.iscomplexobj(other):
                raise ValueError('Can not divide real coefficients by ' +
                                 'a complex constant.')
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
                    self.csphase == other.csphase and self.kind == other.kind):
                coeffs = _np.empty([2, self.lmax+1, self.lmax+1],
                                   dtype=self.coeffs.dtype)
                coeffs[self.mask] = (self.coeffs[self.mask] /
                                     other.coeffs[self.mask])
                return SHGravCoeffs.from_array(
                    coeffs, gm=self.gm, r0=self.r0, omega=self.omega,
                    csphase=self.csphase, normalization=self.normalization)
            else:
                raise ValueError('The two sets of coefficients must be of ' +
                                 'the same kind and have the same ' +
                                 'normalization and csphase.')
        elif _np.isscalar(other) is True:
            coeffs = _np.empty([2, self.lmax+1, self.lmax+1],
                               dtype=self.coeffs.dtype)
            if self.kind == 'real' and _np.iscomplexobj(other):
                raise ValueError('Can not multiply real coefficients by ' +
                                 'a complex constant.')
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

    def spectrum(self, spectrum='geoid', lmax=None, convention='power',
                 unit='per_l', base=10.):
        """
        Return the spectrum as a function of spherical harmonic degree.

        Usage
        -----
        power = x.spectrum([spectrum, lmax, convention, unit, base])

        Returns
        -------
        power : ndarray, shape (lmax+1)
            1-D numpy ndarray of the spectrum, where lmax is the maximum
            spherical harmonic degree.

        Parameters
        ----------
        spectrum : str, optional, default = 'geoid'
            The type of power spectrum to return: 'potential' for the
            gravitational potential, 'geoid' for the geoid, or 'gravity' for
            the radial gravity.
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

        Description
        -----------
        This function returns either the power spectrum, energy spectrum, or
        l2-norm spectrum. The output spectrum can be one of three types as
        defined by the optional parameter 'spectrum': 'potential' for the
        gravitational potential, 'geoid' for the geoid, or 'gravity' for the
        radial gravity.

        Total power is defined as the integral of the function squared over all
        space, divided by the area the function spans. If the mean of the
        function is zero, this is equivalent to the variance of the function.
        The total energy is the integral of the function squared over all space
        and is 4pi times the total power. For normalized coefficients ('4pi',
        'ortho', or 'schmidt'), the l2-norm is the sum of the magnitude of the
        coefficients squared.

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
        if spectrum.lower() not in ('potential', 'geoid', 'gravity'):
            raise ValueError(
                "Spectrum must be of type 'potential', "
                "'geoid', or 'gravity'. Provided value was {:s}"
                .format(repr(spectrum))
                )

        s = _spectrum(self.coeffs, normalization=self.normalization,
                      convention=convention, unit=unit, base=base, lmax=lmax)

        if spectrum.lower() == 'potential':
            s *= (self.gm / self.r0)**2
        elif spectrum.lower() == 'geoid':
            s *= self.r0**2
        elif spectrum.lower() == 'gravity':
            degrees = _np.arange(len(s))
            s *= (self.gm * (degrees + 1) / self.r0**2)**2

        return s


class SHGravRealCoeffs(SHGravCoeffs):
    """
    Real spherical harmonic coefficients class for the gravitational potential.
    """

    def __init__(self, coeffs, gm=None, r0=None, omega=None,
                 normalization='4pi', csphase=1, copy=True, header=None, ):
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
        self.gm = gm
        self.r0 = r0
        self.omega = omega

        if copy:
            self.coeffs = _np.copy(coeffs)
            self.coeffs[~mask] = 0.
        else:
            self.coeffs = coeffs

    def __repr__(self):
        return ('kind = {:s}\n'
                'normalization = {:s}\n'
                'csphase = {:d}\n'
                'lmax = {:d}\n'
                'GM (m3 / s2) = {:s}\n'
                'r0 (m) = {:s}\n'
                'Omega (rad / s) = {:s}\n'
                'header = {:s}\n'
                .format(repr(self.kind), repr(self.normalization),
                        self.csphase, self.lmax, repr(self.gm), repr(self.r0),
                        repr(self.omega), repr(self.header),))

# modify to_file to be to_file = to_file(with header)
# del unsed functions

# need routine to upward/downward continue coefficients, change_r0,
# change_reference_radius()
