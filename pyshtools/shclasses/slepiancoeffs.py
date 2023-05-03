"""
    Class for Slepian expansion coefficients.
"""
import numpy as _np
import copy as _copy

from ..backends import backend_module
from ..backends import preferred_backend
from ..backends import shtools as _shtools

from .shcoeffs import SHCoeffs
from .shgrid import SHGrid


__all__ = ['SlepianCoeffs']


class SlepianCoeffs(object):
    """
    Class for Slepian expansion coefficients.

    The SlepianCoeffs class is initialized by:

    >>>  x = Slepian.expand(flm)

    Each class instance defines the following class attributes:

    falpha          : Array of the Slepian expansion coefficients.
    galpha          : A Slepian class instance that contains the associated
                      Slepian functions.
    nmax            : The number of Slepian expansion coefficients.
    name            : The name of the dataset.

    Each class instance provides the following methods:

    expand()        : Expand the function on a grid an return an SHGrid class
                      instance.
    to_shcoeffs()   : Return the spherical harmonic coefficients of the
                      function as an SHCoeffs class instance.
    plot_spectrum() : Plot the spectrum as a function of spherical harmonic
                      degree.
    copy()          : Return a copy of the class instance.
    info()          : Print a summary of the data stored in the SlepianCoeffs
                      instance.
"""
    def __init__(self, falpha, galpha, name=None, copy=True):
        """
        Initialize the SlepianCoeffs class.
        """
        if copy:
            self.falpha = _np.copy(falpha)
            self.galpha = _copy.deepcopy(galpha)
        else:
            self.falpha = falpha
            self.galpha = galpha

        self.nmax = len(self.falpha)
        self.name = name

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
        Print a summary of the data stored in the SlepianCoeffs class instance.

        Usage
        -----
        x.info()
        """
        print(repr(self))

    def __repr__(self):
        str = ('nmax = {:d}\n'
               'lmax = {:d}\n'
               'name = {:s}\n'
               .format(self.nmax, self.galpha.lmax, repr(self.name)))
        str += '\nSlepian functions:\n' + self.galpha.__repr__()
        return str

    def expand(self, nmax=None, grid='DH2', zeros=None, extend=True,
               backend=None, nthreads=None):
        """
        Expand the function on a grid using the first n Slepian coefficients.

        Usage
        -----
        f = x.expand([nmax, grid, zeros, extend])

        Returns
        -------
        f : SHGrid class instance

        Parameters
        ----------
        nmax : int, optional, default = x.nmax
            The number of expansion coefficients to use when calculating the
            spherical harmonic coefficients.
        grid : str, optional, default = 'DH2'
            'DH' or 'DH1' for an equally sampled grid with nlat=nlon, 'DH2'
            for an equally spaced grid in degrees latitude and longitude, or
            'GLQ' for a Gauss-Legendre quadrature grid.
        zeros : ndarray, optional, default = None
            The cos(colatitude) nodes used in the Gauss-Legendre Quadrature
            grids.
        extend : bool, optional, default = True
            If True, compute the longitudinal band for 360 E (DH and GLQ grids)
            and the latitudinal band for 90 S (DH grids only).
        backend : str, optional, default = preferred_backend()
            Name of the preferred backend, either 'shtools' or 'ducc'.
        nthreads : int, optional, default = 1
            Number of threads to use for the 'ducc' backend. Setting this
            parameter to 0 will use as many threads as there are hardware
            threads on the system.
        """
        if type(grid) != str:
            raise ValueError('grid must be a string. ' +
                             'Input type was {:s}'
                             .format(str(type(grid))))

        if nmax is None:
            nmax = self.nmax
        if backend is None:
            backend = preferred_backend()

        if self.galpha.kind == 'cap':
            shcoeffs = _shtools.SlepianCoeffsToSH(self.falpha,
                                                  self.galpha.coeffs, nmax)
        else:
            shcoeffs = _shtools.SlepianCoeffsToSH(self.falpha,
                                                  self.galpha.tapers, nmax)

        if grid.upper() in ('DH', 'DH1'):
            gridout = backend_module(
                backend=backend, nthreads=nthreads).MakeGridDH(
                    shcoeffs, sampling=1, norm=1, csphase=1, extend=extend)
            return SHGrid.from_array(gridout, grid='DH', copy=False)
        elif grid.upper() == 'DH2':
            gridout = backend_module(
                backend=backend, nthreads=nthreads).MakeGridDH(
                    shcoeffs, sampling=2, norm=1, csphase=1, extend=extend)
            return SHGrid.from_array(gridout, grid='DH', copy=False)
        elif grid.upper() == 'GLQ':
            if backend == "shtools" and zeros is None:
                zeros, weights = _shtools.SHGLQ(self.galpha.lmax)
            gridout = backend_module(
                backend=backend, nthreads=nthreads).MakeGridGLQ(
                    shcoeffs, zeros=zeros, norm=1, csphase=1, extend=extend)
            return SHGrid.from_array(gridout, grid='GLQ', copy=False)
        else:
            raise ValueError(
                "grid must be 'DH', 'DH1', 'DH2', or 'GLQ'. " +
                "Input value was {:s}".format(repr(grid)))

    def to_shcoeffs(self, nmax=None, normalization='4pi', csphase=1):
        """
        Return the spherical harmonic coefficients using the first n Slepian
        coefficients.

        Usage
        -----

        s = x.to_shcoeffs([nmax])

        Returns
        -------
        s : SHCoeffs class instance
            The spherical harmonic coefficients obtained from using the first
            n Slepian expansion coefficients.

        Parameters
        ----------
        nmax : int, optional, default = x.nmax
            The maximum number of expansion coefficients to use when
            calculating the spherical harmonic coefficients.
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

        if nmax is None:
            nmax = self.nmax

        if self.galpha.kind == 'cap':
            shcoeffs = _shtools.SlepianCoeffsToSH(self.falpha,
                                                  self.galpha.coeffs, nmax)
        else:
            shcoeffs = _shtools.SlepianCoeffsToSH(self.falpha,
                                                  self.galpha.tapers, nmax)

        temp = SHCoeffs.from_array(shcoeffs, normalization='4pi', csphase=1)

        if normalization != '4pi' or csphase != 1:
            return temp.convert(normalization=normalization, csphase=csphase)
        else:
            return temp

    def plot_spectrum(self, nmax=None, convention='power', unit='per_l',
                      base=10., lmax=None, xscale='lin', yscale='log',
                      grid=True, legend=None, legend_loc='best',
                      axes_labelsize=None, tick_labelsize=None, show=True,
                      ax=None, fname=None, **kwargs):
        """
        Plot the spectrum as a function of spherical harmonic degree.

        Usage
        -----
        x.plot_spectrum([nmax, convention, unit, base, lmax, xscale, yscale,
                         grid, legend, legend_loc, axes_labelsize,
                         tick_labelsize, legend, show, ax, fname, **kwargs])

        Parameters
        ----------
        nmax : int, optional, default = x.nmax
            The maximum number of expansion coefficients to use when
            calculating the spherical harmonic coefficients.
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
            If present, and if ax is not specified, save the image to the
            specified file.
        **kwargs : keyword arguments, optional
            Keyword arguments for pyplot.plot().

        Notes
        -----
        This method plots either the power spectrum, energy spectrum, or
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
        temp = self.to_shcoeffs(nmax=nmax)

        if ax is None:
            fig, axes = temp.plot_spectrum(convention=convention, unit=unit,
                                           base=base, lmax=lmax, xscale=xscale,
                                           yscale=yscale, grid=grid,
                                           axes_labelsize=axes_labelsize,
                                           tick_labelsize=tick_labelsize,
                                           legend=legend,
                                           legend_loc=legend_loc, show=show,
                                           ax=ax, fname=fname, **kwargs)
            return fig, axes
        else:
            temp.plot_spectrum(convention=convention, unit=unit,
                               base=base, lmax=lmax, xscale=xscale,
                               yscale=yscale, grid=grid,
                               axes_labelsize=axes_labelsize,
                               tick_labelsize=tick_labelsize,
                               legend=legend, legend_loc=legend_loc, show=show,
                               ax=ax, fname=fname, **kwargs)
