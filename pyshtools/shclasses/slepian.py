"""
    Class for Slepian functions on the sphere.

        Slepian: SlepianCap, SlepianMask
"""
import numpy as _np
import matplotlib as _mpl
import matplotlib.pyplot as _plt
from mpl_toolkits.axes_grid1 import make_axes_locatable as _make_axes_locatable
import copy as _copy

from ..backends import shtools as _shtools
from ..spectralanalysis import spectrum as _spectrum

from .shcoeffs import SHCoeffs
from .shgrid import SHGrid
from .slepiancoeffs import SlepianCoeffs


__all__ = ['Slepian', 'SlepianCap', 'SlepianMask']


class Slepian(object):
    """
    Class for Slepian functions on the sphere.

    The Slepian class can be initialized from:

    >>>  x = Slepian.from_cap()
    >>>  x = Slepian.from_mask()

    Each class instance defines the following class attributes:

    kind            : Either 'cap' or 'mask'.
    tapers          : Matrix containing the spherical harmonic coefficients
                      (in packed form) of either the unrotated spherical cap
                      Slepian functions or the Slepian functions corresponding
                      to the input mask.
    coeffs          : Array of spherical harmonic coefficients of the rotated
                      spherical cap Slepian functions. These are '4pi'
                      normalized and do not use the Condon-Shortley phase
                      factor.
    shannon         : The Shannon number, which approximates the number of
                      well localized Slepian functions.
    area            : Area of the concentration domain, in radians.
    eigenvalues     : Concentration factors of the Slepian functions.
    orders          : The angular orders for each of the spherical cap Slepian
                      functions.
    lmax            : Spherical harmonic bandwidth of the Slepian functions.
    theta           : Angular radius of the spherical cap localization domain
                      (default in degrees).
    theta_degrees   : True (default) if theta is in degrees.
    nmax            : The number of Slepian functions. Default is (lmax+1)**2.
    nrot            : The number of best-concentrated spherical cap Slepian
                      functions that were rotated and whose coefficients are
                      stored in coeffs.
    clat, clon      : Latitude and longitude of the center of the rotated
                      spherical cap Slepian functions (default in degrees).
    coord_degrees   : True (default) if clat and clon are in degrees.
    slepian_degrees : Boolean or int array defining which spherical harmonic
                      degrees were used to construct the Slepian functions.

    Each class instance provides the following methods:

    expand()               : Expand the input function in Slepian functions.
    coupling_matrix()      : Compute the spherical harmonic coupling matrix.
    number_concentrated()  : Return the number of functions that have
                             concentration factors greater or equal to a
                             specified value.
    to_array()             : Return an array of the spherical harmonic
                             coefficients for function alpha, where alpha=0 is
                             the best concentrated, optionally using a
                             different normalization convention.
    to_shcoeffs()          : Return the spherical harmonic coefficients of
                             function alpha, where alpha=0 is the best
                             concentrated, as a new SHCoeffs class instance,
                             optionally using a different normalization
                             convention.
    to_shgrid()            : Return as a new SHGrid instance a grid of function
                             alpha, where alpha=0 is the best concentrated.
    degrees()              : Return an array containing the spherical harmonic
                             degrees of the Slepian functions, from 0 to lmax.
    spectra()              : Return the spectra of one or more Slepian
                             function.
    rotate()               : Rotate the spherical cap Slepian functions,
                             originally located at the North pole, to clat and
                             clon and save the spherical harmonic coefficients
                             in the attribute coeffs.
    variance()             : Calculate the theoretical variance of the power of
                             a function expanded in spherical-cap Slepian
                             functions.
    copy()                 : Return a copy of the class instance.
    plot()                 : Plot the best concentrated Slepian functions using
                             a simple cylindrical projection.
    plot_spectra()         : Plot the spectra of the best-concentrated Slepian
                             functions.
    plot_coupling_matrix() : Plot the spherical harmonic coupling matrix.
    info()                 : Print a summary of the data stored in the Slepian
                             instance.
"""

    def __init__(self):
        """Initialize with a factory method."""
        print('Initialize the class using one of the class methods:\n'
              '>>> pyshtools.Slepian.from_cap\n'
              '>>> pyshtools.Slepian.from_mask')

    # ---- factory methods:
    @classmethod
    def from_cap(cls, theta, lmax, clat=None, clon=None, nmax=None,
                 theta_degrees=True, coord_degrees=True, dj_matrix=None,
                 slepian_degrees=None):
        """
        Construct spherical cap Slepian functions.

        Usage
        -----
        x = Slepian.from_cap(theta, lmax, [clat, clon, nmax, theta_degrees,
                                           coord_degrees, dj_matrix,
                                           slepian_degrees])

        Returns
        -------
        x : Slepian class instance

        Parameters
        ----------
        theta : float
            Angular radius of the spherical-cap localization domain (default
            in degrees).
        lmax : int
            Spherical harmonic bandwidth of the Slepian functions.
        clat, clon : float, optional, default = None
            Latitude and longitude of the center of the rotated spherical-cap
            Slepian functions (default in degrees).
        nmax : int, optional, default (lmax+1)**2
            Number of Slepian functions to compute.
        theta_degrees : bool, optional, default = True
            True if theta is in degrees.
        coord_degrees : bool, optional, default = True
            True if clat and clon are in degrees.
        dj_matrix : ndarray, optional, default = None
            The djpi2 rotation matrix computed by a call to djpi2.
        slepian_degrees : bool or int, optional, dimension (lmax+1),
                          default = None
            Boolean or int array defining which spherical harmonic degrees were
            used (True or 1) to construct the Slepian functions.
        """
        if theta_degrees:
            tapers, eigenvalues, taper_order = _shtools.SHReturnTapers(
                _np.radians(theta), lmax, degrees=slepian_degrees)
        else:
            tapers, eigenvalues, taper_order = _shtools.SHReturnTapers(
                theta, lmax, degrees=slepian_degrees)

        return SlepianCap(theta, tapers, eigenvalues, taper_order, clat, clon,
                          nmax, theta_degrees, coord_degrees, dj_matrix,
                          slepian_degrees, copy=False)

    @classmethod
    def from_mask(cls, dh_mask, lmax, nmax=None, slepian_degrees=None):
        """
        Construct Slepian functions that are optimally concentrated within
        the region specified by a mask.

        Usage
        -----
        x = Slepian.from_mask(dh_mask, lmax, [nmax, slepian_degrees])

        Returns
        -------
        x : Slepian class instance

        Parameters
        ----------
        dh_mask :ndarray or SHGrid class instance, shape (nlat, nlon)
            A Driscoll and Healy sampled grid describing the concentration
            region R. All elements should either be 1 or 0 for inside or
            outside of the concentration region, respectively. The grid must
            have dimensions nlon=nlat, nlon=2*nlat, or nlon=2*nlat-1.
        lmax : int
            The spherical harmonic bandwidth of the Slepian functions.
        nmax : int, optional, default = (lmax+1)**2
            The number of best-concentrated eigenvalues and eigenfunctions to
            return.
        slepian_degrees : bool or int, optional, dimension (lmax+1),
                          default = None
            Boolean or int array defining which spherical harmonic degrees were
            used (True or 1) to construct the Slepian functions.
        """
        if nmax is None:
            nmax = (lmax + 1)**2
        else:
            if nmax > (lmax + 1)**2:
                raise ValueError('nmax must be less than or equal to ' +
                                 '(lmax + 1)**2. lmax = {:d} and nmax = {:d}.'
                                 .format(lmax, nmax))

        if isinstance(dh_mask, _np.ndarray):
            mask = SHGrid.from_array(dh_mask, grid='DH', copy=False)
            data = mask.data[:mask.nlat-mask.extend, :mask.nlon-mask.extend]
            area = 4 * _np.pi * mask.expand(lmax_calc=0).coeffs[0, 0, 0]

        elif isinstance(dh_mask, SHGrid):
            if dh_mask.grid != 'DH':
                raise ValueError("The grid type of dh_mask must be 'DH'. "
                                 'Input grid is {:s}.'
                                 .format(repr(dh_mask.grid)))
            data = dh_mask.data[:dh_mask.nlat-dh_mask.extend,
                                :dh_mask.nlon-dh_mask.extend]
            area = 4 * _np.pi * dh_mask.expand(lmax_calc=0).coeffs[0, 0, 0]

        else:
            raise ValueError('dh_mask must be an numpy.ndarrary or '
                             'pyshtools.SHGrid class instance. '
                             'Input type is {:s}.'
                             .format(str(type(dh_mask))))

        tapers, eigenvalues = _shtools.SHReturnTapersMap(
            data, lmax, ntapers=nmax, degrees=slepian_degrees)

        return SlepianMask(tapers, eigenvalues, area, slepian_degrees,
                           copy=False)

    def copy(self):
        """Return a deep copy of the class instance."""
        return _copy.deepcopy(self)

    def degrees(self):
        """
        Return a numpy array listing the spherical harmonic degrees of the
        Slepian functions from 0 to lmax.

        Usage
        -----
        degrees = x.degrees()

        Returns
        -------
        degrees : ndarray, shape (lmax+1)
            numpy ndarray containing a list of the spherical harmonic degrees.
        """
        return _np.arange(self.lmax + 1)

    def number_concentrated(self, concentration):
        """
        Return the number of Slepian functions that have concentration factors
        greater or equal to lambda.

        Usage
        -----
        k = x.number_concentrated(lambda)

        Returns
        -------
        k : int
            The number of Slepian functions with concentration factors greater
            or equal to lambda.

        Parameters
        ----------
        lambda : float
            The concentration factor, which is the power of the function within
            the concentration region divided by the total power.
        """
        return len(self.eigenvalues[self.eigenvalues >= concentration])

    def expand(self, flm, nmax=None):
        """
        Return the Slepian expansion coefficients of the input function.

        Usage
        -----
        s = x.expand(flm, [nmax])

        Returns
        -------
        s : SlepianCoeff class instance
            The Slepian expansion coefficients of the input function.

        Parameters
        ----------
        flm : SHCoeffs class instance
            The input function to expand in Slepian functions.
        nmax : int, optional, default = (x.lmax+1)**2
            The number of Slepian expansion coefficients to compute.

        Notes
        -----
        The global function f is input using its spherical harmonic
        expansion coefficients flm. The expansion coefficients of the function
        f using Slepian functions g is given by

        f_alpha = sum_{lm}^{lmax} f_lm g(alpha)_lm
        """
        if nmax is None:
            nmax = (self.lmax+1)**2
        elif nmax is not None and nmax > (self.lmax+1)**2:
            raise ValueError(
                "nmax must be less than or equal to (lmax+1)**2, "
                "where lmax is {:s}. Input value is {:s}."
                .format(repr(self.lmax), repr(nmax))
                )

        coeffsin = flm.to_array(normalization='4pi', csphase=1, lmax=self.lmax)

        return self._expand(coeffsin, nmax)

    def coupling_matrix(self, nmax=None):
        """
        Return the spherical harmonic coupling matrix. This matrix relates the
        power spectrum expectation of the function expressed in a subset of the
        best-localized Slepian functions to the expectation of the global
        power spectrum.

        Usage
        -----
        kij = x.coupling_matrix([nmax])

        Returns
        -------
        kij : ndarray, shape (lmax+1, lmax+1)
            The coupling matrix that relates the power spectrum expectation of
            the function expressed in a subset of the best-localized Slepian
            functions to the expectation of the global power spectrum.

        Parameters
        ----------
        nmax : int, optional, default = x.nmax
            The number of Slepian functions used in reconstructing the
            function.
        """
        if nmax is None:
            nmax = self.nmax

        return self._coupling_matrix(nmax=nmax)

    def to_array(self, alpha, normalization='4pi', csphase=1):
        """
        Return the spherical harmonic coefficients of Slepian function i as a
        numpy array.

        Usage
        -----
        coeffs = x.to_array(alpha, [normalization, csphase])

        Returns
        -------
        coeffs : ndarray, shape (2, lmax+1, lmax+11)
            3-D numpy ndarray of the spherical harmonic coefficients
            of the Slepian function.

        Parameters
        ----------
        alpha : int
            Function number, where alpha=0 is the best concentrated Slepian
            function.
        normalization : str, optional, default = '4pi'
            Normalization of the output coefficients: '4pi', 'ortho' or
            'schmidt' for geodesy 4pi normalized, orthonormalized, or Schmidt
            semi-normalized coefficients, respectively.
        csphase : int, optional, default = 1
            Condon-Shortley phase convention: 1 to exclude the phase factor,
            or -1 to include it.
        """
        if type(normalization) != str:
            raise ValueError('normalization must be a string. ' +
                             'Input type is {:s}.'
                             .format(str(type(normalization))))

        if normalization.lower() not in ('4pi', 'ortho', 'schmidt'):
            raise ValueError(
                "normalization must be '4pi', 'ortho' " +
                "or 'schmidt'. Provided value is {:s}."
                .format(repr(normalization))
                )
        if csphase != 1 and csphase != -1:
            raise ValueError(
                'csphase must be 1 or -1. Input value is {:s}.'
                .format(repr(csphase))
                )

        return self._to_array(
            alpha, normalization=normalization.lower(), csphase=csphase)

    def to_shcoeffs(self, alpha, normalization='4pi', csphase=1):
        """
        Return the spherical harmonic coefficients of Slepian function i as a
        SHCoeffs class instance.

        Usage
        -----
        clm = x.to_shcoeffs(alpha, [normalization, csphase])

        Returns
        -------
        clm : SHCoeffs class instance

        Parameters
        ----------
        alpha : int
            Function number, where alpha=0 is the best concentrated Slepian
            function.
        normalization : str, optional, default = '4pi'
            Normalization of the output class: '4pi', 'ortho' or 'schmidt' for
            geodesy 4pi-normalized, orthonormalized, or Schmidt semi-normalized
            coefficients, respectively.
        csphase : int, optional, default = 1
            Condon-Shortley phase convention: 1 to exclude the phase factor,
            or -1 to include it.
        """
        if type(normalization) != str:
            raise ValueError('normalization must be a string. ' +
                             'Input type is {:s}.'
                             .format(str(type(normalization))))

        if normalization.lower() not in set(['4pi', 'ortho', 'schmidt']):
            raise ValueError(
                "normalization must be '4pi', 'ortho' " +
                "or 'schmidt'. Provided value is {:s}."
                .format(repr(normalization))
                )
        if csphase != 1 and csphase != -1:
            raise ValueError(
                "csphase must be 1 or -1. Input value is {:s}."
                .format(repr(csphase))
                )

        coeffs = self.to_array(alpha, normalization=normalization.lower(),
                               csphase=csphase)
        return SHCoeffs.from_array(coeffs, normalization=normalization.lower(),
                                   csphase=csphase, copy=False)

    def to_shgrid(self, alpha, grid='DH2', zeros=None, extend=True):
        """
        Evaluate the coefficients of Slepian function i on a grid and return
        an SHGrid class instance.

        Usage
        -----
        f = x.to_shgrid(alpha, [grid, zeros])

        Returns
        -------
        f : SHGrid class instance

        Parameters
        ----------
        alpha : int
            Function number, where alpha=0 is the best concentrated Slepian
            function.
        grid : str, optional, default = 'DH2'
            'DH' or 'DH1' for an equisampled lat/lon grid with nlat=nlon, 'DH2'
            for an equidistant lat/lon grid with nlon=2*nlat, or 'GLQ' for a
            Gauss-Legendre quadrature grid.
        zeros : ndarray, optional, default = None
            The cos(colatitude) nodes used in the Gauss-Legendre Quadrature
            grids.
        extend : bool, optional, default = True
            If True, compute the longitudinal band for 360 E (DH and GLQ grids)
            and the latitudinal band for 90 S (DH grids only).

        Notes
        -----
        For more information concerning the spherical harmonic expansions and
        the properties of the output grids, see the documentation for
        SHExpandDH and SHExpandGLQ.
        """
        if type(grid) != str:
            raise ValueError('grid must be a string. Input type is {:s}.'
                             .format(str(type(grid))))

        if grid.upper() in ('DH', 'DH1'):
            gridout = _shtools.MakeGridDH(self.to_array(alpha), sampling=1,
                                          norm=1, csphase=1, extend=extend)
            return SHGrid.from_array(gridout, grid='DH', copy=False)
        elif grid.upper() == 'DH2':
            gridout = _shtools.MakeGridDH(self.to_array(alpha), sampling=2,
                                          norm=1, csphase=1, extend=extend)
            return SHGrid.from_array(gridout, grid='DH', copy=False)
        elif grid.upper() == 'GLQ':
            if zeros is None:
                zeros, weights = _shtools.SHGLQ(self.lmax)
            gridout = _shtools.MakeGridGLQ(self.to_array(alpha), zeros,
                                           norm=1, csphase=1, extend=extend)
            return SHGrid.from_array(gridout, grid='GLQ', copy=False)
        else:
            raise ValueError(
                "grid must be 'DH', 'DH1', 'DH2', or 'GLQ'. " +
                "Input value is {:s}.".format(repr(grid)))

    def spectra(self, alpha=None, nmax=None, convention='power', unit='per_l',
                base=10.):
        """
        Return the spectra of one or more Slepian functions.

        Usage
        -----
        spectra = x.spectra([alpha, nmax, convention, unit, base])

        Returns
        -------
        spectra : ndarray, shape (lmax+1, nmax)
             A matrix with each column containing the spectrum of a Slepian
             function, and where the functions are arranged with increasing
             concentration factors. If alpha is set, only a single vector is
             returned, whereas if nmax is set, the first nmax spectra are
             returned.

        Parameters
        ----------
        alpha : int, optional, default = None
            The function number of the output spectrum, where alpha=0
            corresponds to the best concentrated Slepian function.
        nmax : int, optional, default = 1
            The number of best concentrated Slepian function power spectra
            to return.
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
        This function returns either the power spectrum, energy spectrum, or
        l2-norm spectrum of one or more of the Slepian funtions. Total power
        is defined as the integral of the function squared over all space,
        divided by the area the function spans. If the mean of the function is
        zero, this is equivalent to the variance of the function. The total
        energy is the integral of the function squared over all space and is
        4pi times the total power. The l2-norm is the sum of the magnitude of
        the coefficients squared.

        The output spectrum can be expresed using one of three units. 'per_l'
        returns the contribution to the total spectrum from all angular orders
        at degree l. 'per_lm' returns the average contribution to the total
        spectrum from a single coefficient at degree l. The 'per_lm' spectrum
        is equal to the 'per_l' spectrum divided by (2l+1). 'per_dlogl' returns
        the contribution to the total spectrum from all angular orders over an
        infinitessimal logarithmic degree band. The contrubution in the band
        dlog_a(l) is spectrum(l, 'per_dlogl')*dlog_a(l), where a is the base,
        and where spectrum(l, 'per_dlogl) is equal to
        spectrum(l, 'per_l')*l*log(a).

         """
        if alpha is None:
            if nmax is None:
                nmax = self.nmax
            spectra = _np.zeros((self.lmax+1, nmax))

            for iwin in range(nmax):
                coeffs = self.to_array(iwin)
                spectra[:, iwin] = _spectrum(coeffs, normalization='4pi',
                                             convention=convention, unit=unit,
                                             base=base)
        else:
            coeffs = self.to_array(alpha)
            spectra = _spectrum(coeffs, normalization='4pi',
                                convention=convention, unit=unit, base=base)

        return spectra

    def variance(self, power, k, lmax=None):
        """
        Calculate the theoretical variance of the power of a function expanded
        in spherical-cap Slepian functions.

        Usage
        -----
        variance = x.variance(power, k, [lmax])

        Returns
        -------
        variance : ndarray, shape (lmax+1)
            The theoretical variance of the spectrum estimate.

        Parameters
        ----------
        power : ndarray, dimension (lmax_in+1)
            The input global power spectrum.
        k : int
            The number of Slepian functions used to represent the function.
        lmax : int, optional, default = min(lmax_in, self.lmax)
            The maximum spherical harmonic degree of the variance to compute.
        """
        if lmax is None:
            lmax = min(len(power) - 1, self.lmax)
        else:
            if lmax > self.lmax:
                raise ValueError('lmax must be less than or equal to '
                                 'self.lmax. Input value is {:s}, and '
                                 'self.lmax is {:s}.'.format(repr(lmax),
                                                             repr(self.lmax)))

        return self._variance(power, k, lmax=lmax)

    def plot(self, nmax, projection=None, lmax=None, maxcolumns=3,
             ticks='WSen', tick_interval=[60, 45],
             minor_tick_interval=[None, None], xlabel='Longitude',
             ylabel='Latitude', title=True, colorbar=None, cmap='viridis',
             cmap_limits=None, cmap_reverse=False, cb_triangles='neither',
             cb_label=None, cb_ylabel=None, cb_tick_interval=None,
             cb_minor_tick_interval=None, grid=False, loss=False,
             axes_labelsize=None, tick_labelsize=None, titlesize=8,
             show=True, ax=None, fname=None, cb_offset=None, cb_width=None):
        """
        Plot the best-concentrated Slepian functions.

        Usage
        -----
        x.plot(nmax, [projections, lmax, maxcolumns, tick_interval,
                      minor_tick_interval, ticks, xlabel, ylabel, title,
                      titlesize, colorbar, cmap, cmap_limits, cmap_reverse,
                      cb_triangles, cb_label, cb_ylabel, cb_tick_interval,
                      cb_minor_tick_interval, cb_offset, cb_width, grid, loss,
                      axes_labelsize, tick_labelsize, ax, show, fname])

        Parameters
        ----------
        nmax : int
            The number of Slepian functions to plot.
        projection : Cartopy projection class, optional, default = None
            The Cartopy projection class used to project the gridded data,
            for Driscoll and Healy sampled grids only.
        lmax : int, optional, default = self.lmax
            The maximum degree to use when plotting the Slepian function, which
            controls the number of samples in latitude and longitude.
        maxcolumns : int, optional, default = 3
            The maximum number of columns to use when plotting multiple Slepian
            functions.
        tick_interval : list or tuple, optional, default = [60, 45]
            Intervals to use when plotting the x and y ticks. If set to None,
            ticks will not be plotted.
        minor_tick_interval : list or tuple, optional, default = [None, None]
            Intervals to use when plotting the minor x and y ticks. If set to
            None, minor ticks will not be plotted.
        ticks : str, optional, default = 'WSen'
            Specify which axes should have ticks drawn and annotated. Capital
            letters plot the ticks and annotations, whereas small letters plot
            only the ticks. 'W', 'S', 'E', and 'N' denote the west, south, east
            and north boundaries of the plot.
        xlabel : str, optional, default = 'longitude'
            Label for the longitude axis.
        ylabel : str, optional, default = 'latitude'
            Label for the latitude axis.
        title : bool, optional, default = True
            If True, plot a title on top of each subplot providing the taper
            number and 1 minus the concentration factor.
        colorbar : str, optional, default = None
            Plot a colorbar along the 'top', 'right', 'bottom', or 'left' axis.
        cmap : str, optional, default = 'viridis'
            The color map to use when plotting the data.
        cmap_limits : list, optional, default = [self.min(), self.max()]
            Set the lower and upper limits of the data used by the colormap,
            and optionally an interval for each color band. If the
            interval is specified, the number of discrete colors will be
            (cmap_limits[1]-cmap_limits[0])/cmap_limits[2].
        cmap_reverse : bool, optional, default = False
            Set to True to reverse the sense of the color progression in the
            color table.
        cb_triangles : str, optional, default = 'neither'
            Add triangles to the edges of the colorbar for minimum and maximum
            values. Can be 'neither', 'both', 'min', or 'max'.
        cb_label : str, optional, default = None
            Text label for the colorbar.
        cb_ylabel : str, optional, default = None
            Text label for the y axis of the colorbar
        cb_tick_interval : float, optional, default = None
            Colorbar major tick and annotation interval.
        cb_minor_tick_interval : float, optional, default = None
            Colorbar minor tick interval.
        cb_offset : float or int, optional, default = None
            Offset of the colorbar from the map edge in points. If None,
            the offset will be calculated automatically.
        cb_width : float, optional, default = None
            Width of the colorbar in percent with respect to the width of the
            respective image axis. Defaults are 2.5 and 5 for vertical and
            horizontal colorbars, respectively.
        grid : bool, optional, default = False
            If True, plot major grid lines.
        loss : bool, optional, default = False
            When plotting titles, provide the loss factor instead of the
            concentration factor (loss=1-concentration).
        titlesize : int, optional, default = 8
            The font size for the subplot titles.
        axes_labelsize : int, optional, default = None
            The font size for the x and y axes labels.
        tick_labelsize : int, optional, default = None
            The font size for the x and y tick labels.
        ax : matplotlib axes object, optional, default = None
            An array of matplotlib axes objects where the plots will appear.
        show : bool, optional, default = True
            If True, plot the image to the screen.
        fname : str, optional, default = None
            If present, save the image to the specified file.
        """
        if self.kind == 'cap':
            if self.nrot is not None and self.nrot <= nmax:
                nmax = self.nrot

        ncolumns = min(maxcolumns, nmax)
        nrows = _np.ceil(nmax / ncolumns).astype(int)
        figsize = (_mpl.rcParams['figure.figsize'][0],
                   _mpl.rcParams['figure.figsize'][0]
                   * 0.55 * nrows / ncolumns + 0.41)

        if ax is None:
            fig, axes = _plt.subplots(nrows, ncolumns, figsize=figsize,
                                      sharex='all', sharey='all')
        else:
            if hasattr(ax, 'flatten') and ax.size < nmax:
                raise ValueError('ax.size must be greater or equal to nmax. ' +
                                 'nmax = {:s}'.format(repr(nmax)) +
                                 ' and ax.size = {:s}.'.format(repr(ax.size)))
            axes = ax

        for alpha in range(min(self.nmax, nmax)):
            evalue = self.eigenvalues[alpha]
            if min(self.nmax, nmax) == 1 and ax is None:
                axtemp = axes
            elif hasattr(axes, 'flatten'):
                axtemp = axes.flatten()[alpha]
            else:
                axtemp = axes[alpha]
            coeffs = self.to_shcoeffs(alpha)
            if lmax is not None:
                coeffs = coeffs.pad(lmax=lmax, copy=False)
            grid_temp = coeffs.expand()

            if title:
                if loss:
                    title_str = '#{:d} [loss={:2.2g}]'.format(alpha, 1-evalue)
                else:
                    title_str = '#{:d} [concentration={:2.2g}]'.format(
                        alpha, evalue)
            else:
                title_str = None

            grid_temp.plot(projection=projection, tick_interval=tick_interval,
                           minor_tick_interval=minor_tick_interval,
                           title=title_str, ticks=ticks,
                           xlabel=xlabel, ylabel=ylabel, grid=grid,
                           cmap=cmap, cmap_reverse=cmap_reverse,
                           axes_labelsize=axes_labelsize,
                           tick_labelsize=tick_labelsize,
                           colorbar=colorbar, cmap_limits=cmap_limits,
                           cb_triangles=cb_triangles, cb_label=cb_label,
                           cb_ylabel=cb_ylabel, cb_offset=cb_offset,
                           cb_tick_interval=cb_tick_interval,
                           cb_width=cb_width,
                           cb_minor_tick_interval=cb_minor_tick_interval,
                           titlesize=titlesize, ax=axtemp)

        if ax is None and projection is None:
            if nrows > 1:
                for axtemp in axes[:-1, :].flatten():
                    for xlabel_i in axtemp.get_xticklabels():
                        xlabel_i.set_visible(False)
                    axtemp.set_xlabel('', visible=False)
                for axtemp in axes[:, 1:].flatten():
                    for ylabel_i in axtemp.get_yticklabels():
                        ylabel_i.set_visible(False)
                    axtemp.set_ylabel('', visible=False)
            elif nmax > 1:
                for axtemp in axes[1:].flatten():
                    for ylabel_i in axtemp.get_yticklabels():
                        ylabel_i.set_visible(False)
                    axtemp.set_ylabel('', visible=False)

        if ax is None:
            fig.tight_layout(pad=0.5)
            if show:
                fig.show()
            if fname is not None:
                fig.savefig(fname)
            return fig, axes

    def plot_spectra(self, nmax, convention='power', unit='per_l', base=10.,
                     maxcolumns=3, xscale='lin', yscale='log', grid=True,
                     xlim=(None, None), ylim=(None, None), show=True,
                     title=True, axes_labelsize=None, tick_labelsize=None,
                     titlesize=None, ax=None, fname=None):
        """
        Plot the spectra of the best-concentrated Slepian functions.

        Usage
        -----
        x.plot_spectra(nmax, [convention, unit, base, maxcolumns, xscale,
                              yscale, grid, xlim, ylim, show, title,
                              axes_labelsize, tick_labelsize, titlesize,
                              ax, fname])

        Parameters
        ----------
        nmax : int
            The number of Slepian functions to plot.
        convention : str, optional, default = 'power'
            The type of spectra to plot: 'power' for power spectrum, and
            'energy' for energy spectrum.
        unit : str, optional, default = 'per_l'
            If 'per_l', return the total contribution to the spectrum for each
            spherical harmonic degree l. If 'per_lm', return the average
            contribution to the spectrum for each coefficient at spherical
            harmonic degree l. If 'per_dlogl', return the spectrum per log
            interval dlog_a(l).
        base : float, optional, default = 10.
            The logarithm base when calculating the 'per_dlogl' spectrum.
        maxcolumns : int, optional, default = 3
            The maximum number of columns to use when plotting the spectra
            of multiple localization windows.
        xscale : str, optional, default = 'lin'
            Scale of the x axis: 'lin' for linear or 'log' for logarithmic.
        yscale : str, optional, default = 'log'
            Scale of the y axis: 'lin' for linear or 'log' for logarithmic.
        grid : bool, optional, default = True
            If True, plot grid lines.
        xlim : tuple, optional, default = (None, None)
            The upper and lower limits used for the x axis.
        ylim : tuple, optional, default = (None, None)
            The lower and upper limits used for the y axis.
        show : bool, optional, default = True
            If True, plot the image to the screen.
        title : bool, optional, default = True
            If True, plot a legend on top of each subplot providing the taper
            number and 1 minus the concentration factor.
        axes_labelsize : int, optional, default = None
            The font size for the x and y axes labels.
        tick_labelsize : int, optional, default = None
            The font size for the x and y tick labels.
        titlesize : int, optional, default = None
            The font size for the subplot titles.
        ax : matplotlib axes object, optional, default = None
            An array of matplotlib axes objects where the plots will appear.
        fname : str, optional, default = None
            If present, save the image to the file.
        """
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

        degrees = self.degrees()
        spectrum = self.spectra(nmax=nmax, convention=convention, unit=unit,
                                base=base)

        ncolumns = min(maxcolumns, nmax)
        nrows = _np.ceil(nmax / ncolumns).astype(int)
        figsize = (_mpl.rcParams['figure.figsize'][0],
                   _mpl.rcParams['figure.figsize'][0]
                   * 0.7 * nrows / ncolumns + 0.41)

        if ax is None:
            fig, axes = _plt.subplots(nrows, ncolumns, figsize=figsize,
                                      sharex='all', sharey='all')
        else:
            if hasattr(ax, 'flatten') and ax.size < nmax:
                raise ValueError('ax.size must be greater or equal to nmax. ' +
                                 'nmax = {:s}'.format(repr(nmax)) +
                                 ' and ax.size = {:s}.'.format(repr(ax.size)))
            axes = ax

        if ax is None:
            if nrows > 1:
                for axtemp in axes[:-1, :].flatten():
                    for xlabel_i in axtemp.get_xticklabels():
                        xlabel_i.set_visible(False)
                    axtemp.set_xlabel('', visible=False)
                for axtemp in axes[:, 1:].flatten():
                    for ylabel_i in axtemp.get_yticklabels():
                        ylabel_i.set_visible(False)
                    axtemp.set_ylabel('', visible=False)
            elif nmax > 1:
                for axtemp in axes[1:].flatten():
                    for ylabel_i in axtemp.get_yticklabels():
                        ylabel_i.set_visible(False)
                    axtemp.set_ylabel('', visible=False)

        if ylim == (None, None):
            upper = spectrum[:, :min(self.nmax, nmax)].max()
            lower = upper * 1.e-6
            ylim = (lower, 5 * upper)
        if xlim == (None, None):
            if xscale == 'lin':
                xlim = (degrees[0], degrees[-1])

        for alpha in range(min(self.nmax, nmax)):
            evalue = self.eigenvalues[alpha]
            if min(self.nmax, nmax) == 1 and ax is None:
                axtemp = axes
            elif hasattr(axes, 'flatten'):
                axtemp = axes.flatten()[alpha]
            else:
                axtemp = axes[alpha]
            if (convention == 'power'):
                axtemp.set_ylabel('Power', fontsize=axes_labelsize)
            else:
                axtemp.set_ylabel('Energy', fontsize=axes_labelsize)

            if yscale == 'log':
                axtemp.set_yscale('log', base=base)

            if xscale == 'log':
                axtemp.set_xscale('log', base=base)
                axtemp.plot(degrees[1:], spectrum[1:, alpha],
                            label='#{:d} [loss={:2.2g}]'
                            .format(alpha, 1-evalue))
            else:
                axtemp.plot(degrees[0:], spectrum[0:, alpha],
                            label='#{:d} [loss={:2.2g}]'
                            .format(alpha, 1-evalue))
            axtemp.set_xlabel('Spherical harmonic degree',
                              fontsize=axes_labelsize)
            axtemp.set(xlim=xlim, ylim=ylim)
            axtemp.minorticks_on()
            axtemp.grid(grid, which='major')
            axtemp.tick_params(labelsize=tick_labelsize)
            if title is True:
                axtemp.set_title('#{:d} [loss={:2.2g}]'
                                 .format(alpha, 1-evalue),
                                 fontsize=titlesize)

        if ax is None:
            fig.tight_layout(pad=0.5)
            if show:
                fig.show()
            if fname is not None:
                fig.savefig(fname)
            return fig, axes

    def plot_coupling_matrix(self, nmax=None, vmin=None, vmax=None,
                             xlabel='Input degree', ylabel='Output degree',
                             title=None, axes_labelsize=None,
                             tick_labelsize=None, titlesize=None,
                             colorbar=None, cb_label=None, normalize=False,
                             show=True, ax=None, fname=None, **kwargs):
        """
        Plot the spherical harmonic coupling matrix. This matrix relates the
        power spectrum expectation of the function expressed in a subset of the
        best-localized Slepian functions to the expectation of the global
        power spectrum.

        Usage
        -----
        x.plot_coupling_matrix([nmax, vmin, vmax, xlabel, ylabel, title
                                axes_labelsize, tick_labelsize,
                                titlesize, colorbar, cb_label, normalize,
                                show, ax, fname, **kwargs])

        Parameters
        ----------
        nmax : int, optional, default = x.nmax
            The number of Slepian functions used in reconstructing the
            function.
        vmin : float, optional, default=None
            The minmum range of the colormap. If None, the minimum value of the
            spectrum will be used.
        vmax : float, optional, default=None
            The maximum range of the colormap. If None, the maximum value of
            the spectrum will be used.
        xlabel : str, optional, default = 'Input degree'
            Label for the x axis.
        ylabel : str, optional, default = 'Output degree'
            Label for the y axis.
        title : str, optional, default = None
            Add a title to the plot.
        axes_labelsize : int, optional, default = None
            The font size for the x and y axes labels.
        tick_labelsize : int, optional, default = None
            The font size for the x and y tick labels.
        titlesize : int, optional, default = None
            The font size for the title.
        colorbar : str, optional, default = None
            Plot a colorbar that is either 'horizontal' or 'vertical'.
        cb_label : str, optional, default = None
            Text label for the colorbar.
        normalize : bool, optional, default = False
            Normalize the coupling maxtrix such that the maximum value is 1.
        show : bool, optional, default = True
            If True, plot the image to the screen.
        ax : matplotlib axes object, optional, default = None
            An array of matplotlib axes objects where the plots will appear.
        fname : str, optional, default = None
            If present, save the image to the specified file.
        kwargs : optional
            Keyword arguements that will be sent to plt.imshow(), such as cmap.
        """
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

        if ax is None:
            if colorbar is not None:
                if colorbar.lower()[0] == 'h':
                    scale = 1.1
                elif colorbar.lower()[0] == 'v':
                    scale = 0.85
                else:
                    raise ValueError("colorbar must be either 'horizontal' or "
                                     "'vertical'. Input value is {:s}."
                                     .format(repr(colorbar)))
            else:
                scale = 1

            figsize = (_mpl.rcParams['figure.figsize'][0],
                       _mpl.rcParams['figure.figsize'][0] * scale)
            fig, axes = _plt.subplots(1, 1, figsize=figsize)
        else:
            axes = ax

        kll = self.coupling_matrix(nmax=nmax)
        if normalize:
            kll = kll / kll.max()

        cim = axes.imshow(kll, aspect='equal', vmin=vmin, vmax=vmax, **kwargs)
        if xlabel is not None:
            axes.set_xlabel(xlabel, fontsize=axes_labelsize)
        if ylabel is not None:
            axes.set_ylabel(ylabel, fontsize=axes_labelsize)
        if title is not None:
            axes.set_title(title, fontsize=titlesize)
        axes.tick_params(labelsize=tick_labelsize)
        axes.minorticks_on()

        if colorbar is not None:
            if colorbar.lower()[0] == 'v':
                divider = _make_axes_locatable(axes)
                cax = divider.append_axes("right", size="2.5%", pad=0.15)
                cbar = _plt.colorbar(cim, cax=cax, orientation='vertical')
            elif colorbar.lower()[0] == 'h':
                divider = _make_axes_locatable(axes)
                cax = divider.append_axes("bottom", size="2.5%", pad=0.5)
                cbar = _plt.colorbar(cim, cax=cax, orientation='horizontal')
            if cb_label is not None:
                cbar.set_label(cb_label, fontsize=axes_labelsize)
            cbar.ax.tick_params(labelsize=tick_labelsize)

        if ax is None:
            fig.tight_layout(pad=0.5)
            if show:
                fig.show()
            if fname is not None:
                fig.savefig(fname)
            return fig, axes

    def info(self):
        """
        Print a summary of the data stored in the Slepian instance.

        Usage
        -----
        x.info()
        """
        self._info()


class SlepianCap(Slepian):
    """Class for Slepian functions concentrated within a spherical cap."""

    @staticmethod
    def istype(kind):
        return kind == 'cap'

    def __init__(self, theta, tapers, eigenvalues, taper_order, clat, clon,
                 nmax, theta_degrees, coord_degrees, dj_matrix,
                 slepian_degrees, copy=True):
        self.kind = 'cap'
        self.theta = theta
        self.clat = clat
        self.clon = clon
        self.lmax = tapers.shape[0] - 1
        self.theta_degrees = theta_degrees
        self.coord_degrees = coord_degrees
        self.dj_matrix = dj_matrix
        self.nrot = None
        self.slepian_degrees = slepian_degrees

        if (self.theta_degrees):
            self.area = 2 * _np.pi * (1 - _np.cos(_np.radians(self.theta)))
        else:
            self.area = 2 * _np.pi * (1 - _np.cos(self.theta))

        if nmax is not None:
            self.nmax = nmax
        else:
            self.nmax = tapers.shape[1]

        if self.nmax > (self.lmax + 1)**2:
            raise ValueError('nmax must be less than or equal to ' +
                             '(lmax+1)**2. nmax = {:s} and lmax = {:s}.'
                             .format(repr(self.nmax), repr(self.lmax)))

        if copy:
            self.tapers = _np.copy(tapers[:, :self.nmax])
            self.eigenvalues = _np.copy(eigenvalues[:self.nmax])
            self.orders = _np.copy(taper_order[:self.nmax])
        else:
            self.tapers = tapers[:, :self.nmax]
            self.eigenvalues = eigenvalues[:self.nmax]
            self.orders = taper_order[:self.nmax]

        if self.slepian_degrees is None:
            self.shannon = (self.lmax + 1)**2 / (4 * _np.pi) * self.area
        else:
            self.shannon = sum(self.eigenvalues)

        # If the windows aren't rotated, don't store them.
        if self.clat is None and self.clon is None:
            self.coeffs = None
        else:
            self.rotate(clat=self.clat, clon=self.clon,
                        coord_degrees=self.coord_degrees,
                        dj_matrix=self.dj_matrix)

    def _expand(self, coeffsin, nmax):
        """
        Determine the Slepian expansion coefficients of a function.
        """
        if self.coeffs is None:
            self.rotate(clat=90., clon=0., nrot=nmax)

        if (nmax > self.nrot):
            raise ValueError('nmax must be less than or equal to ' +
                             'nrot. nmax = {:d}, nrot = {:d}.'
                             .format(nmax, self.nrot))
        falpha = _shtools.SlepianCoeffs(self.coeffs, coeffsin, nmax)

        return SlepianCoeffs(falpha, self)

    def _coupling_matrix(self, nmax):
        """Return the coupling matrix."""
        return _shtools.SHSCouplingMatrixCap(self.tapers, self.orders, nmax)

    def _taper2coeffs(self, alpha):
        """
        Return the spherical harmonic coefficients of the unrotated Slepian
        function i as an array, where i = 0 is the best concentrated function.
        """
        taperm = self.orders[alpha]
        coeffs = _np.zeros((2, self.lmax + 1, self.lmax + 1))
        if taperm < 0:
            coeffs[1, :, abs(taperm)] = self.tapers[:, alpha]
        else:
            coeffs[0, :, abs(taperm)] = self.tapers[:, alpha]

        return coeffs

    def _to_array(self, alpha, normalization='4pi', csphase=1):
        """
        Return the spherical harmonic coefficients of Slepian function i as an
        array, where i = 0 is the best concentrated function.
        """
        if self.coeffs is None:
            coeffs = _np.copy(self._taper2coeffs(alpha))
        else:
            if alpha > self.nrot - 1:
                raise ValueError('alpha must be less than or equal to ' +
                                 'nrot - 1. alpha = {:d}, nrot = {:d}.'
                                 .format(alpha, self.nrot))
            coeffs = _shtools.SHVectorToCilm(self.coeffs[:, alpha])

        if normalization == 'schmidt':
            for l in range(self.lmax + 1):
                coeffs[:, l, :l+1] *= _np.sqrt(2.0 * l + 1.0)
        elif normalization == 'ortho':
            coeffs *= _np.sqrt(4.0 * _np.pi)

        if csphase == -1:
            for m in range(self.lmax + 1):
                if m % 2 == 1:
                    coeffs[:, :, m] = - coeffs[:, :, m]

        return coeffs

    def rotate(self, clat, clon, coord_degrees=True, dj_matrix=None,
               nrot=None):
        """"
        Rotate the spherical-cap Slepian functions centered on the North pole
        to clat and clon, and save the spherical harmonic coefficients in the
        attribute coeffs.

        Usage
        -----
        x.rotate(clat, clon [coord_degrees, dj_matrix, nrot])

        Parameters
        ----------
        clat, clon : float
            Latitude and longitude of the center of the rotated spherical-cap
            Slepian functions (default in degrees).
        coord_degrees : bool, optional, default = True
            True if clat and clon are in degrees.
        dj_matrix : ndarray, optional, default = None
            The djpi2 rotation matrix computed by a call to djpi2.
        nrot : int, optional, default = (lmax+1)**2
            The number of best-concentrated Slepian functions to rotate, where
            lmax is the spherical harmonic bandwidth of the functions.

        Notes
        -----
        This function will take the spherical-cap Slepian functions centered at
        the North pole (and saved in the attributes tapers and orders), rotate
        each function to the coordinate (clat, clon), and save the spherical
        harmonic coefficients in the attribute coeffs. Each column of coeffs
        contains a single Slepian function, and the coefficients are ordered
        according to the convention in SHCilmToVector.
        """
        self.coeffs = _np.zeros(((self.lmax + 1)**2, self.nmax))
        self.clat = clat
        self.clon = clon
        self.coord_degrees = coord_degrees

        if nrot is not None:
            self.nrot = nrot
        else:
            self.nrot = self.nmax

        if self.coord_degrees:
            angles = _np.radians(_np.array([0., -(90. - clat), -clon]))
        else:
            angles = _np.array([0., -(_np.pi/2. - clat), -clon])

        if dj_matrix is None:
            if self.dj_matrix is None:
                self.dj_matrix = _shtools.djpi2(self.lmax + 1)
                dj_matrix = self.dj_matrix
            else:
                dj_matrix = self.dj_matrix

        if ((coord_degrees is True and clat == 90. and clon == 0.) or
                (coord_degrees is False and clat == _np.pi/2. and clon == 0.)):
            for i in range(self.nrot):
                coeffs = self._taper2coeffs(i)
                self.coeffs[:, i] = _shtools.SHCilmToVector(coeffs)

        else:
            coeffs = _shtools.SHRotateTapers(self.tapers, self.orders,
                                             self.nrot, angles, dj_matrix)
            self.coeffs = coeffs

    def _variance(self, power, k, lmax=None):
        """
        Compute the theoretical variance of the power of a function expanded in
        Slepian functions.
        """
        var = _np.zeros(lmax+1)
        for l in range(lmax+1):
            var[l] = _shtools.SHSlepianVar(l, self.tapers, self.orders, power,
                                           kmax=k)

        return var

    def _info(self):
        """Print a summary of the data in the Slepian instance."""
        print(repr(self))

    def __repr__(self):
        str = 'kind = {:s}\n'.format(repr(self.kind))

        if self.theta_degrees:
            str += 'theta = {:f} degrees\n'.format(self.theta)
        else:
            str += 'theta = {:f} radians'.format(self.theta)

        str += ('lmax = {:d}\n'
                'nmax = {:d}\n'
                'nrot = {:s}\n'
                'shannon = {:f}\n'
                'area (radians) = {:e}\n'
                .format(self.lmax, self.nmax, repr(self.nrot), self.shannon,
                        self.area))

        if self.clat is not None:
            if self.coord_degrees:
                str += 'clat = {:f} degrees\n'.format(self.clat)
            else:
                str += 'clat = {:f} radians\n'.format(self.clat)
        else:
            str += 'clat is not specified\n'

        if self.clon is not None:
            if self.coord_degrees:
                str += 'clon = {:f} degrees\n'.format(self.clon)
            else:
                str += 'clon = {:f} radians\n'.format(self.clon)
        else:
            str += 'clon is not specified\n'

        if self.dj_matrix is not None:
            str += 'dj_matrix is stored\n'
        else:
            str += 'dj_matrix is not stored\n'

        return str


class SlepianMask(Slepian):
    """
    Class for Slepian functions concentrated within a specified mask and
    for a given spherical harmonic bandwidth.
    """

    @staticmethod
    def istype(kind):
        return kind == 'mask'

    def __init__(self, tapers, eigenvalues, area, slepian_degrees, copy=True):
        self.kind = 'mask'
        self.lmax = _np.sqrt(tapers.shape[0]).astype(int) - 1
        self.nmax = tapers.shape[1]
        self.slepian_degrees = slepian_degrees

        if copy:
            self.tapers = _np.copy(tapers)
            self.eigenvalues = _np.copy(eigenvalues)
        else:
            self.tapers = tapers
            self.eigenvalues = eigenvalues

        self.area = area

        if self.slepian_degrees is None:
            self.shannon = (self.lmax + 1)**2 / (4 * _np.pi) * self.area
        else:
            self.shannon = sum(self.eigenvalues)

    def _expand(self, coeffsin, nmax):
        """
        Determine the Slepian expansion coefficients of a function.
        """
        falpha = _shtools.SlepianCoeffs(self.tapers, coeffsin, nmax)

        return SlepianCoeffs(falpha, self)

    def _coupling_matrix(self, nmax):
        """Return the coupling matrix."""
        return _shtools.SHSCouplingMatrix(self.tapers, nmax)

    def _to_array(self, alpha, normalization='4pi', csphase=1):
        """
        Return the spherical harmonic coefficients of Slepian function i as an
        array, where i=0 is the best concentrated function.
        """
        coeffs = _shtools.SHVectorToCilm(self.tapers[:, alpha])

        if normalization == 'schmidt':
            for l in range(self.lmax + 1):
                coeffs[:, l, :l+1] *= _np.sqrt(2.0 * l + 1.0)
        elif normalization == 'ortho':
            coeffs *= _np.sqrt(4.0 * _np.pi)

        if csphase == -1:
            for m in range(self.lmax + 1):
                if m % 2 == 1:
                    coeffs[:, :, m] = - coeffs[:, :, m]

        return coeffs

    def _variance(self, power, nwin, lmax=None):
        """
        Compute the theoretical variance of the power of a function expanded in
        Slepian functions.
        """
        raise RuntimeError('Computation of the theoretical variance is '
                           'not yet implemented for arbitrary windows.')

    def _info(self):
        """Print a summary of the data in the Slepian instance."""
        print(repr(self))

    def __repr__(self):
        str = ('kind = {:s}\n'
               'lmax = {:d}\n'
               'nmax = {:d}\n'
               'shannon = {:f}\n'
               'area (radians) = {:e}\n'.format(repr(self.kind), self.lmax,
                                                self.nmax, self.shannon,
                                                self.area))

        return str
