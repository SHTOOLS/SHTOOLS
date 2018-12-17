"""
    Class for Slepian functions on the sphere.

        Slepian: SlepianCap, SlepianMask
"""

from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

import numpy as _np
import matplotlib as _mpl
import matplotlib.pyplot as _plt
import copy as _copy

from .. import shtools as _shtools
from ..spectralanalysis import spectrum as _spectrum

from .shcoeffsgrid import SHCoeffs
from .shcoeffsgrid import SHGrid
from .slepiancoeffs import SlepianCoeffs


__all__ = ['Slepian', 'SlepianCap', 'SlepianMask']


class Slepian(object):
    """
    Class for Slepian functions on the sphere.

    The Slepian class can be initialized from:

    >>>  x = Slepian.from_cap(theta, lmax, [clat, clon, nmax])
    >>>  x = Slepian.from_mask(SHGrid)

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

    Each class instance provides the following methods:

    expand()              : Expand the input function in Slepian functions.
    to_array()            : Return an array of the spherical harmonic
                            coefficients for function alpha, where alpha=0 is
                            the best concentrated, optionally using a different
                            normalization convention.
    to_shcoeffs()         : Return the spherical harmonic coefficients of
                            function alpha, where alpha=0 is the best
                            concentrated, as a new SHCoeffs class instance,
                            optionally using a different normalization
                            convention.
    to_shgrid()           : Return as a new SHGrid instance a grid of function
                            alpha, where alpha=0 is the best concentrated.
    number_concentrated() : Return the number of functions that have
                            concentration factors greater or equal to a
                            specified value.
    degrees()             : Return an array containing the spherical harmonic
                            degrees of the Slepian functions, from 0 to lmax.
    spectra()             : Return the spectra of one or more Slepian function.
    rotate()              : Rotate the spherical cap Slepian functions,
                            originally located at the North pole, to clat and
                            clon and save the spherical harmonic coefficients
                            in the attribute coeffs.
    copy()                : Return a copy of the class instance.
    plot()                : Plot the best concentrated Slepian functions using
                            a simple cylindrical projection.
    plot_spectra()        : Plot the spectra of the best-concentrated Slepian
                            functions.
    info()                : Print a summary of the data stored in the Slepian
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
                 theta_degrees=True, coord_degrees=True, dj_matrix=None):
        """
        Construct spherical cap Slepian functions.

        Usage
        -----
        x = Slepian.from_cap(theta, lmax, [clat, clon, nmax, theta_degrees,
                                           coord_degrees, dj_matrix])

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
        """
        if theta_degrees:
            tapers, eigenvalues, taper_order = _shtools.SHReturnTapers(
                _np.radians(theta), lmax)
        else:
            tapers, eigenvalues, taper_order = _shtools.SHReturnTapers(
                theta, lmax)

        return SlepianCap(theta, tapers, eigenvalues, taper_order, clat, clon,
                          nmax, theta_degrees, coord_degrees, dj_matrix,
                          copy=False)

    @classmethod
    def from_mask(cls, dh_mask, lmax, nmax=None):
        """
        Construct Slepian functions that are optimally concentrated within
        the region specified by a mask.

        Usage
        -----
        x = Slepian.from_mask(dh_mask, lmax, [nmax])

        Returns
        -------
        x : Slepian class instance

        Parameters
        ----------
        dh_mask :ndarray, shape (nlat, nlon)
            A Driscoll and Healy (1994) sampled grid describing the
            concentration region R. All elements should either be 1 (for inside
            the concentration region) or 0 (for outside the concentration
            region). The grid must have dimensions nlon=nlat or nlon=2*nlat,
            where nlat is even.
        lmax : int
            The spherical harmonic bandwidth of the Slepian functions.
        nmax : int, optional, default = (lmax+1)**2
            The number of best-concentrated eigenvalues and eigenfunctions to
            return.
        """
        if nmax is None:
            nmax = (lmax + 1)**2
        else:
            if nmax > (lmax + 1)**2:
                raise ValueError('nmax must be less than or equal to ' +
                                 '(lmax + 1)**2. lmax = {:d} and nmax = {:d}'
                                 .format(lmax, nmax))

        if dh_mask.shape[0] % 2 != 0:
            raise ValueError('The number of latitude bands in dh_mask ' +
                             'must be even. nlat = {:d}'
                             .format(dh_mask.shape[0]))

        if dh_mask.shape[1] == dh_mask.shape[0]:
            _sampling = 1
        elif dh_mask.shape[1] == 2 * dh_mask.shape[0]:
            _sampling = 2
        else:
            raise ValueError('dh_mask must be dimensioned as (n, n) or ' +
                             '(n, 2 * n). Input shape is ({:d}, {:d})'
                             .format(dh_mask.shape[0], dh_mask.shape[1]))

        mask_lm = _shtools.SHExpandDH(dh_mask, sampling=_sampling, lmax_calc=0)
        area = mask_lm[0, 0, 0] * 4 * _np.pi

        tapers, eigenvalues = _shtools.SHReturnTapersMap(dh_mask, lmax,
                                                         ntapers=nmax)

        return SlepianMask(tapers, eigenvalues, area, copy=False)

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

        Description
        -----------
        The global function f is input using its spherical harmonic
        expansion coefficients flm. The expansion coefficients of the function
        f using Slepian functions g is given by
        
        f_alpha = sum_{lm}^{lmax} f_lm g(alpha)_lm
        """
        if nmax is None:
            nmax = (self.lmax+1)**2
        elif nmax is not None and nmax > (self.lmax+1)**2:
            raise ValueError(
                "nmax must be less than or equal to (lmax+1)**2 " +
                "where lmax is {:s}. Input value is {:s}"
                .format(repr(self.lmax), repr(nmax))
                )

        coeffsin = flm.to_array(normalization='4pi', csphase=1, lmax=self.lmax)

        return self._expand(coeffsin, nmax)

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
                             'Input type was {:s}'
                             .format(str(type(normalization))))

        if normalization.lower() not in ('4pi', 'ortho', 'schmidt'):
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

        coeffs = self.to_array(alpha, normalization=normalization.lower(),
                               csphase=csphase)
        return SHCoeffs.from_array(coeffs, normalization=normalization.lower(),
                                   csphase=csphase, copy=False)

    def to_shgrid(self, alpha, grid='DH2', zeros=None):
        """
        Evaluate the coefficients of Slepian function i on a spherical grid and
        return a SHGrid class instance.

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

        Description
        -----------
        For more information concerning the spherical harmonic expansions and
        the properties of the output grids, see the documentation for
        SHExpandDH and SHExpandGLQ.
        """
        if type(grid) != str:
            raise ValueError('grid must be a string. ' +
                             'Input type was {:s}'
                             .format(str(type(grid))))

        if grid.upper() in ('DH', 'DH1'):
            gridout = _shtools.MakeGridDH(self.to_array(alpha), sampling=1,
                                          norm=1, csphase=1)
            return SHGrid.from_array(gridout, grid='DH', copy=False)
        elif grid.upper() == 'DH2':
            gridout = _shtools.MakeGridDH(self.to_array(alpha), sampling=2,
                                          norm=1, csphase=1)
            return SHGrid.from_array(gridout, grid='DH', copy=False)
        elif grid.upper() == 'GLQ':
            if zeros is None:
                zeros, weights = _shtools.SHGLQ(self.lmax)
            gridout = _shtools.MakeGridGLQ(self.to_array(alpha), zeros,
                                           norm=1, csphase=1)
            return SHGrid.from_array(gridout, grid='GLQ', copy=False)
        else:
            raise ValueError(
                "grid must be 'DH', 'DH1', 'DH2', or 'GLQ'. " +
                "Input value was {:s}".format(repr(grid)))

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

        Description
        -----------
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

    def plot(self, nmax, lmax=None, maxcolumns=3,
             tick_interval=[60, 45], minor_tick_interval=None,
             xlabel='Longitude', ylabel='Latitude',
             axes_labelsize=None, tick_labelsize=None,
             title_labelsize=None, grid=False, show=True, title=True,
             ax=None, fname=None):
        """
        Plot the best-concentrated Slepian functions.

        Usage
        -----
        x.plot(nmax, [lmax, maxcolumns, tick_interval, minor_tick_interval,
                      xlabel, ylabel, grid, show, title, axes_labelsize,
                      tick_labelsize, title_labelsize, ax, fname])

        Parameters
        ----------
        nmax : int
            The number of Slepian functions to plot.
        lmax : int, optional, default = self.lmax
            The maximum degree to use when plotting the Slepian function, which
            controls the number of samples in latitude and longitude.
        maxcolumns : int, optional, default = 3
            The maximum number of columns to use when plotting multiple Slepian
            functions.
        tick_interval : list or tuple, optional, default = [60, 45]
            Intervals to use when plotting the x and y ticks. If set to None,
            ticks will not be plotted.
        minor_tick_interval : list or tuple, optional, default = None
            Intervals to use when plotting the minor x and y ticks. If set to
            None, minor ticks will not be plotted.
        xlabel : str, optional, default = 'longitude'
            Label for the longitude axis.
        ylabel : str, optional, default = 'latitude'
            Label for the latitude axis.
        grid : bool, optional, default = False
            If True, plot grid lines.
        show : bool, optional, default = True
            If True, plot the image to the screen.
        title : bool, optional, default = True
            If True, plot a title on top of each subplot providing the taper
            number and 1 minus the concentration factor.
        axes_labelsize : int, optional, default = None
            The font size for the x and y axes labels.
        tick_labelsize : int, optional, default = None
            The font size for the x and y tick labels.
        title_labelsize : int, optional, default = None
            The font size for the subplot titles.
        ax : matplotlib axes object, optional, default = None
            An array of matplotlib axes objects where the plots will appear.
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

        if tick_interval is None:
            xticks = []
            yticks = []
        else:
            xticks = _np.linspace(0, 360, num=360//tick_interval[0]+1,
                                  endpoint=True)
            yticks = _np.linspace(-90, 90, num=180//tick_interval[1]+1,
                                  endpoint=True)

        if axes_labelsize is None:
            axes_labelsize = _mpl.rcParams['axes.labelsize']
        if tick_labelsize is None:
            tick_labelsize = _mpl.rcParams['xtick.labelsize']
        if title_labelsize is None:
            title_labelsize = _mpl.rcParams['axes.titlesize']

        if minor_tick_interval is None:
            minor_xticks = []
            minor_yticks = []
        else:
            minor_xticks = _np.linspace(
                0, 360, num=360//minor_tick_interval[0]+1, endpoint=True)
            minor_yticks = _np.linspace(
                -90, 90, num=180//minor_tick_interval[1]+1, endpoint=True)

        deg = '$^{\circ}$'
        xticklabels = [str(int(y)) + deg for y in xticks]
        yticklabels = [str(int(y)) + deg for y in yticks]

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

        for alpha in range(min(self.nmax, nmax)):
            evalue = self.eigenvalues[alpha]
            if min(self.nmax, nmax) == 1 and ax is None:
                axtemp = axes
            elif hasattr(axes, 'flatten'):
                axtemp = axes.flatten()[alpha]
            else:
                axtemp = axes[alpha]
            gridout = _shtools.MakeGridDH(self.to_array(alpha), sampling=2,
                                          lmax=lmax, norm=1, csphase=1)
            axtemp.imshow(gridout, origin='upper',
                          extent=(0., 360., -90., 90.))
            axtemp.set(xticks=xticks, yticks=yticks)
            axtemp.set_xlabel(xlabel, fontsize=axes_labelsize)
            axtemp.set_ylabel(ylabel, fontsize=axes_labelsize)
            axtemp.set_xticklabels(xticklabels, fontsize=tick_labelsize)
            axtemp.set_yticklabels(yticklabels, fontsize=tick_labelsize)
            axtemp.set_xticks(minor_xticks, minor=True)
            axtemp.set_yticks(minor_yticks, minor=True)
            axtemp.grid(grid, which='major')
            if title is True:
                axtemp.set_title('#{:d} [loss={:2.2g}]'
                                 .format(alpha, 1-evalue),
                                 fontsize=title_labelsize)

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
                     title_labelsize=None, ax=None, fname=None):
        """
        Plot the spectra of the best-concentrated Slepian functions.

        Usage
        -----
        x.plot_spectra(nmax, [convention, unit, base, maxcolumns, xscale,
                              yscale, grid, xlim, ylim, show, title,
                              axes_labelsize, tick_labelsize, title_labelsize,
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
        title_labelsize : int, optional, default = None
            The font size for the subplot titles.
        ax : matplotlib axes object, optional, default = None
            An array of matplotlib axes objects where the plots will appear.
        fname : str, optional, default = None
            If present, save the image to the file.
        """
        if axes_labelsize is None:
            axes_labelsize = _mpl.rcParams['axes.labelsize']
        if tick_labelsize is None:
            tick_labelsize = _mpl.rcParams['xtick.labelsize']
        if title_labelsize is None:
            title_labelsize = _mpl.rcParams['axes.titlesize']

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
                axtemp.set_yscale('log', basey=base)

            if xscale == 'log':
                axtemp.set_xscale('log', basex=base)
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
                                 fontsize=title_labelsize)

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
                 nmax, theta_degrees, coord_degrees, dj_matrix, copy=True):
        self.kind = 'cap'
        self.theta = theta
        self.clat = clat
        self.clon = clon
        self.lmax = tapers.shape[0] - 1
        self.theta_degrees = theta_degrees
        self.coord_degrees = coord_degrees
        self.dj_matrix = dj_matrix
        self.nrot = None

        if (self.theta_degrees):
            self.area = 2 * _np.pi * (1 - _np.cos(_np.radians(self.theta)))
        else:
            self.area = 2 * _np.pi * (1 - _np.cos(self.theta))

        self.shannon = (self.lmax + 1)**2 / (4 * _np.pi) * self.area

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
            falpha = _shtools.SlepianCoeffs(self.coeffs, coeffsin, self.nrot)
        else:
            falpha = _shtools.SlepianCoeffs(self.coeffs, coeffsin, self.nrot)

        return SlepianCoeffs(falpha, self)

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
                                 'nrot - 1. alpha = {:d}, nrot = {:d}'
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

        Description
        -----------
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
                'shannon = {:e}\n'
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

    def __init__(self, tapers, eigenvalues, area, copy=True):
        self.kind = 'mask'
        self.lmax = _np.sqrt(tapers.shape[0]).astype(int) - 1
        self.nmax = tapers.shape[1]
        if copy:
            self.tapers = _np.copy(tapers)
            self.eigenvalues = _np.copy(eigenvalues)
        else:
            self.tapers = tapers
            self.eigenvalues = eigenvalues

        self.area = area
        self.shannon = (self.lmax + 1)**2 / (4 * _np.pi) * self.area

    def _expand(self, coeffsin, nmax):
        """
        Determine the Slepian expansion coefficients of a function.
        """
        falpha = _shtools.SlepianCoeffs(self.tapers, coeffsin, nmax)

        return SlepianCoeffs(falpha, self)

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

    def _info(self):
        """Print a summary of the data in the Slepian instance."""
        print(repr(self))

    def __repr__(self):
        str = ('kind = {:s}\n'
               'lmax = {:d}\n'
               'nmax = {:d}\n'
               'shannon = {:e}\n'
               'area (radians) = {:e}\n'.format(repr(self.kind), self.lmax,
                                                self.nmax, self.shannon,
                                                self.area))

        return str
