"""
    Spherical Harmonic Window classes

        SHWindow: SHWindowCap, SHWindowMask
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


__all__ = ['SHWindow', 'SHWindowCap', 'SHWindowMask']


class SHWindow(object):
    """
    Class for spatio-spectral localization windows on the sphere.

    The windows can be initialized from:

    >>>  x = SHWindow.from_cap(theta, lwin, [clat, clon, nwin])
    >>>  x = SHWindow.from_mask(SHGrid)

    Each class instance defines the following class attributes:

    kind            : Either 'cap' or 'mask'.
    tapers          : Matrix containing the spherical harmonic coefficients
                      (in packed form) of either the unrotated spherical cap
                      localization windows or the localization windows
                      corresponding to the input mask.
    coeffs          : Array of spherical harmonic coefficients of the
                      rotated spherical cap localization windows. These are
                      '4pi' normalized and do not use the Condon-Shortley phase
                      factor.
    eigenvalues     : Concentration factors of the localization windows.
    orders          : The angular orders for each of the spherical cap
                      localization windows.
    weights         : Taper weights used with the multitaper spectral analyses.
                      Defaut is None.
    lwin            : Spherical harmonic bandwidth of the localization windows.
    theta           : Angular radius of the spherical cap localization domain
                      (default in degrees).
    theta_degrees   : True (default) if theta is in degrees.
    nwin            : Number of localization windows. Default is (lwin+1)**2.
    nwinrot         : The number of best concentrated windows that were rotated
                      and whose coefficients are stored in coeffs.
    clat, clon      : Latitude and longitude of the center of the rotated
                      spherical cap localization windows (default in degrees).
    coord_degrees   : True (default) if clat and clon are in degrees.

    Each class instance provides the following methods:

    to_array()            : Return an array of the spherical harmonic
                            coefficients for taper i, where i=0 is the best
                            concentrated, optionally using a different
                            normalization convention.
    to_shcoeffs()         : Return the spherical harmonic coefficients of taper
                            i, where i=0 is the best concentrated, as a new
                            SHCoeffs class instance, optionally using a
                            different normalization convention.
    to_shgrid()           : Return as a new SHGrid instance a grid of taper i,
                            where i=0 is the best concentrated window.
    number_concentrated() : Return the number of windows that have
                            concentration factors greater or equal to a
                            specified value.
    degrees()             : Return an array containing the spherical harmonic
                            degrees of the localization windows, from 0 to
                            lwin.
    spectra()             : Return the spectra of one or more localization
                            windows.
    rotate()              : Rotate the spherical cap tapers, originally located
                            at the north pole, to clat and clon and save the
                            spherical harmonic coefficients in the attribute
                            coeffs.
    coupling_matrix()     : Return the coupling matrix of the first nwin
                            localization windows.
    biased_spectrum()     : Calculate the multitaper (cross-) spectrum
                            expectation of a localized function.
    multitaper_spectrum()       : Return the multitaper power spectrum
                                  estimate and uncertainty for the input
                                  SHCoeffs class instance.
    multitaper_cross_spectrum() : Return the multitaper cross-power
                                  spectrum estimate and uncertainty for
                                  two input SHCoeffs class instances.
    copy()                 : Return a copy of the class instance.
    plot_windows()         : Plot the best concentrated localization windows
                             using a simple cylindrical projection.
    plot_spectra()         : Plot the spectra of the best-concentrated
                             localization windows.
    plot_coupling_matrix() : Plot the multitaper coupling matrix.
    info()                 : Print a summary of the data stored in the SHWindow
                             instance.
"""

    def __init__(self):
        """Initialize with a factory method."""
        pass

    # ---- factory methods:
    @classmethod
    def from_cap(cls, theta, lwin, clat=None, clon=None, nwin=None,
                 theta_degrees=True, coord_degrees=True, dj_matrix=None,
                 weights=None):
        """
        Construct spherical cap localization windows.

        Usage
        -----
        x = SHWindow.from_cap(theta, lwin, [clat, clon, nwin, theta_degrees,
                                            coord_degrees, dj_matrix, weights])

        Returns
        -------
        x : SHWindow class instance

        Parameters
        ----------
        theta : float
            Angular radius of the spherical cap localization domain (default
            in degrees).
        lwin : int
            Spherical harmonic bandwidth of the localization windows.
        clat, clon : float, optional, default = None
            Latitude and longitude of the center of the rotated spherical cap
            localization windows (default in degrees).
        nwin : int, optional, default (lwin+1)**2
            Number of localization windows.
        theta_degrees : bool, optional, default = True
            True if theta is in degrees.
        coord_degrees : bool, optional, default = True
            True if clat and clon are in degrees.
        dj_matrix : ndarray, optional, default = None
            The djpi2 rotation matrix computed by a call to djpi2.
        weights : ndarray, optional, default = None
            Taper weights used with the multitaper spectral analyses.
        """
        if theta_degrees:
            tapers, eigenvalues, taper_order = _shtools.SHReturnTapers(
                _np.radians(theta), lwin)
        else:
            tapers, eigenvalues, taper_order = _shtools.SHReturnTapers(
                theta, lwin)

        return SHWindowCap(theta, tapers, eigenvalues, taper_order,
                           clat, clon, nwin, theta_degrees, coord_degrees,
                           dj_matrix, weights, copy=False)

    @classmethod
    def from_mask(cls, dh_mask, lwin, nwin=None, weights=None):
        """
        Construct localization windows that are optimally concentrated within
        the region specified by a mask.

        Usage
        -----
        x = SHWindow.from_mask(dh_mask, lwin, [nwin, weights])

        Returns
        -------
        x : SHWindow class instance

        Parameters
        ----------
        dh_mask :ndarray, shape (nlat, nlon)
            A Driscoll and Healy (1994) sampled grid describing the
            concentration region R. All elements should either be 1 (for inside
            the concentration region) or 0 (for outside the concentration
            region). The grid must have dimensions nlon=nlat or nlon=2*nlat,
            where nlat is even.
        lwin : int
            The spherical harmonic bandwidth of the localization windows.
        nwin : int, optional, default = (lwin+1)**2
            The number of best concentrated eigenvalues and eigenfunctions to
            return.
        weights ndarray, optional, default = None
            Taper weights used with the multitaper spectral analyses.
        """
        if nwin is None:
            nwin = (lwin + 1)**2
        else:
            if nwin > (lwin + 1)**2:
                raise ValueError('nwin must be less than or equal to ' +
                                 '(lwin + 1)**2. lwin = {:d} and nwin = {:d}'
                                 .format(lwin, nwin))

        if dh_mask.shape[0] % 2 != 0:
            raise ValueError('The number of latitude bands in dh_mask ' +
                             'must be even. nlat = {:d}'
                             .format(dh_mask.shape[0]))
        if (dh_mask.shape[1] != dh_mask.shape[0] and
                dh_mask.shape[1] != 2 * dh_mask.shape[0]):
            raise ValueError('dh_mask must be dimensioned as (n, n) or ' +
                             '(n, 2 * n). Input shape is ({:d}, {:d})'
                             .format(dh_mask.shape[0], dh_mask.shape[1]))

        tapers, eigenvalues = _shtools.SHReturnTapersMap(dh_mask, lwin,
                                                         ntapers=nwin)
        return SHWindowMask(tapers, eigenvalues, weights, copy=False)

    def copy(self):
        """Return a deep copy of the class instance."""
        return _copy.deepcopy(self)

    def degrees(self):
        """
        Return a numpy array listing the spherical harmonic degrees of the
        localization windows from 0 to lwin.

        Usage
        -----
        degrees = x.degrees()

        Returns
        -------
        degrees : ndarray, shape (lwin+1)
            numpy ndarray containing a list of the spherical harmonic degrees.
        """
        return _np.arange(self.lwin + 1)

    def number_concentrated(self, alpha):
        """
        Return the number of localization windows that have concentration
        factors greater or equal to alpha.

        Usage
        -----
        k = x.number_concentrated(alpha)

        Returns
        -------
        k : int
            The number of windows with concentration factors greater or equal
            to alpha.

        Parameters
        ----------
        alpha : float
            The concentration factor, which is the power of the window within
            the concentration region divided by the total power.
        """
        return len(self.eigenvalues[self.eigenvalues >= alpha])

    def to_array(self, itaper, normalization='4pi', csphase=1):
        """
        Return the spherical harmonic coefficients of taper i as a numpy
        array.

        Usage
        -----
        coeffs = x.to_array(itaper, [normalization, csphase])

        Returns
        -------
        coeffs : ndarray, shape (2, lwin+1, lwin+11)
            3-D numpy ndarray of size of the spherical harmonic coefficients
            of the window.

        Parameters
        ----------
        itaper : int
            Taper number, where itaper=0 is the best concentrated.
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
            itaper, normalization=normalization.lower(), csphase=csphase)

    def to_shcoeffs(self, itaper, normalization='4pi', csphase=1):
        """
        Return the spherical harmonic coefficients of taper i as a SHCoeffs
        class instance.

        Usage
        -----
        clm = x.to_shcoeffs(itaper, [normalization, csphase])

        Returns
        -------
        clm : SHCoeffs class instance

        Parameters
        ----------
        itaper : int
            Taper number, where itaper=0 is the best concentrated.
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

        coeffs = self.to_array(itaper, normalization=normalization.lower(),
                               csphase=csphase)
        return SHCoeffs.from_array(coeffs, normalization=normalization.lower(),
                                   csphase=csphase, copy=False)

    def to_shgrid(self, itaper, grid='DH2', zeros=None):
        """
        Evaluate the coefficients of taper i on a spherical grid and return
        a SHGrid class instance.

        Usage
        -----
        f = x.to_shgrid(itaper, [grid, zeros])

        Returns
        -------
        f : SHGrid class instance

        Parameters
        ----------
        itaper : int
            Taper number, where itaper=0 is the best concentrated.
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
            gridout = _shtools.MakeGridDH(self.to_array(itaper), sampling=1,
                                          norm=1, csphase=1)
            return SHGrid.from_array(gridout, grid='DH', copy=False)
        elif grid.upper() == 'DH2':
            gridout = _shtools.MakeGridDH(self.to_array(itaper), sampling=2,
                                          norm=1, csphase=1)
            return SHGrid.from_array(gridout, grid='DH', copy=False)
        elif grid.upper() == 'GLQ':
            if zeros is None:
                zeros, weights = _shtools.SHGLQ(self.lwin)
            gridout = _shtools.MakeGridGLQ(self.to_array(itaper), zeros,
                                           norm=1, csphase=1)
            return SHGrid.from_array(gridout, grid='GLQ', copy=False)
        else:
            raise ValueError(
                "grid must be 'DH', 'DH1', 'DH2', or 'GLQ'. " +
                "Input value was {:s}".format(repr(grid)))

    def multitaper_spectrum(self, clm, k, convention='power', unit='per_l',
                            **kwargs):
        """
        Return the multitaper spectrum estimate and standard error.

        Usage
        -----
        mtse, sd = x.multitaper_spectrum(clm, k, [convention, unit, lmax,
                                                  taper_wt, clat, clon,
                                                  coord_degrees])

        Returns
        -------
        mtse : ndarray, shape (lmax-lwin+1)
            The localized multitaper spectrum estimate, where lmax is the
            spherical-harmonic bandwidth of clm, and lwin is the
            spherical-harmonic bandwidth of the localization windows.
        sd : ndarray, shape (lmax-lwin+1)
            The standard error of the localized multitaper spectrum
            estimate.

        Parameters
        ----------
        clm : SHCoeffs class instance
            SHCoeffs class instance containing the spherical harmonic
            coefficients of the global field to analyze.
        k : int
            The number of tapers to be utilized in performing the multitaper
            spectral analysis.
        convention : str, optional, default = 'power'
            The type of output spectra: 'power' for power spectra, and
            'energy' for energy spectra.
        unit : str, optional, default = 'per_l'
            The units of the output spectra. If 'per_l', the spectra contain
            the total contribution for each spherical harmonic degree l. If
            'per_lm', the spectra contain the average contribution for each
            coefficient at spherical harmonic degree l.
        lmax : int, optional, default = clm.lmax
            The maximum spherical-harmonic degree of clm to use.
        taper_wt : ndarray, optional, default = None
            1-D numpy array of the weights used in calculating the multitaper
            spectral estimates and standard error.
        clat, clon : float, optional, default = 90., 0.
            Latitude and longitude of the center of the spherical-cap
            localization windows.
        coord_degrees : bool, optional, default = True
            True if clat and clon are in degrees.
        """
        return self._multitaper_spectrum(clm, k, convention=convention,
                                         unit=unit, **kwargs)

    def multitaper_cross_spectrum(self, clm, slm, k, convention='power',
                                  unit='per_l', **kwargs):
        """
        Return the multitaper cross-spectrum estimate and standard error.

        Usage
        -----
        mtse, sd = x.multitaper_cross_spectrum(clm, slm, k, [convention, unit,
                                                             lmax, taper_wt,
                                                             clat, clon,
                                                             coord_degrees])

        Returns
        -------
        mtse : ndarray, shape (lmax-lwin+1)
            The localized multitaper cross-spectrum estimate, where lmax is the
            smaller of the two spherical-harmonic bandwidths of clm and slm,
            and lwin is the spherical-harmonic bandwidth of the localization
            windows.
        sd : ndarray, shape (lmax-lwin+1)
            The standard error of the localized multitaper cross-spectrum
            estimate.

        Parameters
        ----------
        clm : SHCoeffs class instance
            SHCoeffs class instance containing the spherical harmonic
            coefficients of the first global field to analyze.
        slm : SHCoeffs class instance
            SHCoeffs class instance containing the spherical harmonic
            coefficients of the second global field to analyze.
        k : int
            The number of tapers to be utilized in performing the multitaper
            spectral analysis.
        convention : str, optional, default = 'power'
            The type of output spectra: 'power' for power spectra, and
            'energy' for energy spectra.
        unit : str, optional, default = 'per_l'
            The units of the output spectra. If 'per_l', the spectra contain
            the total contribution for each spherical harmonic degree l. If
            'per_lm', the spectra contain the average contribution for each
            coefficient at spherical harmonic degree l.
        lmax : int, optional, default = min(clm.lmax, slm.lmax)
            The maximum spherical-harmonic degree of the input coefficients
            to use.
        taper_wt : ndarray, optional, default = None
            The weights used in calculating the multitaper cross-spectral
            estimates and standard error.
        clat, clon : float, optional, default = 90., 0.
            Latitude and longitude of the center of the spherical-cap
            localization windows.
        coord_degrees : bool, optional, default = True
            True if clat and clon are in degrees.
        """
        return self._multitaper_cross_spectrum(clm, slm, k,
                                               convention=convention,
                                               unit=unit, **kwargs)

    def biased_spectrum(self, power, k, convention='power', unit='per_l',
                        **kwargs):
        """
        Calculate the multitaper (cross-)spectrum expectation of a
        localized function.

        Usage
        -----
        outspectrum = x.biased_spectrum(spectrum, k, [unit, power, taper_wt,
                                                      save_cg, ldata])

        Returns
        -------
        outspectrum : ndarray, shape (ldata+lwin+1)
            The expectation of the windowed spectrum, where ldata is the
            spherical-harmonic bandwidth of the input spectrum, and lwin is the
            spherical-harmonic bandwidth of the localization windows.

        Parameters
        ----------
        spectrum : ndarray, shape (ldata+1)
            The global spectrum.
        k : int
            The number of best-concentrated localization windows to use in
            constructing the windowed spectrum.
        convention : str, optional, default = 'power'
            The type of input global and output biased spectra: 'power' for
            power spectra, and 'energy' for energy spectra.
        unit : str, optional, default = 'per_l'
            The units of the input global and output biased spectra. If
            'per_l', the spectra contain the total contribution for each
            spherical harmonic degree l. If 'per_lm', the spectra contain the
            average contribution for each coefficient at spherical harmonic
            degree l.
        taper_wt : ndarray, optional, default = None
            The weights used in calculating the multitaper spectral estimates
            and standard error.
        save_cg : int, optional, default = 0
            If 1, the Clebsch-Gordon coefficients will be precomputed and saved
            for future use. If 0, the Clebsch-Gordon coefficients will be
            recomputed for each call.
        ldata : int, optional, default = len(power)-1
            The maximum degree of the global unwindowed spectrum.
        """
        return self._biased_spectrum(power, k, convention=convention,
                                     unit=unit, **kwargs)

    def spectra(self, itaper=None, nwin=None, convention='power', unit='per_l',
                base=10.):
        """
        Return the spectra of one or more localization windows.

        Usage
        -----
        spectra = x.spectra([itaper, nwin, convention, unit, base])

        Returns
        -------
        spectra : ndarray, shape (lwin+1, nwin)
             A matrix with each column containing the spectrum of a
             localization window, and where the windows are arranged with
             increasing concentration factors. If itaper is set, only a single
             vector is returned, whereas if nwin is set, the first nwin spectra
             are returned.

        Parameters
        ----------
        itaper : int, optional, default = None
            The taper number of the output spectrum, where itaper=0
            corresponds to the best concentrated taper.
        nwin : int, optional, default = 1
            The number of best concentrated localization window power spectra
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
        l2-norm spectrum of one or more of the localization windows.
        Total power is defined as the integral of the function squared over all
        space, divided by the area the function spans. If the mean of the
        function is zero, this is equivalent to the variance of the function.
        The total energy is the integral of the function squared over all space
        and is 4pi times the total power. The l2-norm is the sum of the
        magnitude of the coefficients squared.

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
        if itaper is None:
            if nwin is None:
                nwin = self.nwin
            spectra = _np.zeros((self.lwin+1, nwin))

            for iwin in range(nwin):
                coeffs = self.to_array(iwin)
                spectra[:, iwin] = _spectrum(coeffs, normalization='4pi',
                                             convention=convention, unit=unit,
                                             base=base)
        else:
            coeffs = self.to_array(itaper)
            spectra = _spectrum(coeffs, normalization='4pi',
                                convention=convention, unit=unit, base=base)

        return spectra

    def coupling_matrix(self, lmax, nwin=None, weights=None, mode='full'):
        """
        Return the coupling matrix of the first nwin tapers. This matrix
        relates the global power spectrum to the expectation of the localized
        multitaper spectrum.

        Usage
        -----
        Mmt = x.coupling_matrix(lmax, [nwin, weights, mode])

        Returns
        -------
        Mmt : ndarray, shape (lmax+lwin+1, lmax+1) or (lmax+1, lmax+1) or
              (lmax-lwin+1, lmax+1)

        Parameters
        ----------
        lmax : int
            Spherical harmonic bandwidth of the global power spectrum.
        nwin : int, optional, default = x.nwin
            Number of tapers used in the mutlitaper spectral analysis.
        weights : ndarray, optional, default = x.weights
            Taper weights used with the multitaper spectral analyses.
        mode : str, opitonal, default = 'full'
            'full' returns a biased output spectrum of size lmax+lwin+1. The
            input spectrum is assumed to be zero for degrees l>lmax.
            'same' returns a biased output spectrum with the same size
            (lmax+1) as the input spectrum. The input spectrum is assumed to be
            zero for degrees l>lmax.
            'valid' returns a biased spectrum with size lmax-lwin+1. This
            returns only that part of the biased spectrum that is not
            influenced by the input spectrum beyond degree lmax.
        """
        if weights is not None:
            if nwin is not None:
                if len(weights) != nwin:
                    raise ValueError(
                        'Length of weights must be equal to nwin. ' +
                        'len(weights) = {:d}, nwin = {:d}'.format(len(weights),
                                                                  nwin))
            else:
                if len(weights) != self.nwin:
                    raise ValueError(
                        'Length of weights must be equal to nwin. ' +
                        'len(weights) = {:d}, nwin = {:d}'.format(len(weights),
                                                                  self.nwin))

        if mode == 'full':
            return self._coupling_matrix(lmax, nwin=nwin, weights=weights)
        elif mode == 'same':
            cmatrix = self._coupling_matrix(lmax, nwin=nwin,
                                            weights=weights)
            return cmatrix[:lmax+1, :]
        elif mode == 'valid':
            cmatrix = self._coupling_matrix(lmax, nwin=nwin,
                                            weights=weights)
            return cmatrix[:lmax - self.lwin+1, :]
        else:
            raise ValueError("mode has to be 'full', 'same' or 'valid', not "
                             "{}".format(mode))

    def plot_windows(self, nwin, lmax=None, maxcolumns=5,
                     tick_interval=[60, 45], xlabel='longitude',
                     ylabel='latitude', show=True, ax=None, legend=True,
                     fname=None):
        """
        Plot the best-concentrated localization windows.

        Usage
        -----
        x.plot_windows(nwin, [lmax, maxcolumns, tick_interval, xlabel, ylabel,
                              show, ax, legend, fname])

        Parameters
        ----------
        nwin : int
            The number of localization windows to plot.
        lmax : int, optional, default = self.lwin
            The maximum degree to use when plotting the windows, which controls
            the number of samples in latitude and longitude.
        maxcolumns : int, optional, default = 5
            The maximum number of columns to use when plotting multiple
            localization windows.
        tick_interval : list or tuple, optional, default = [60, 45]
            Intervals to use when plotting the x and y ticks. If set to None,
            ticks will not be plotted.
        xlabel : str, optional, default = 'longitude'
            Label for the longitude axis.
        ylabel : str, optional, default = 'latitude'
            Label for the latitude axis.
        show : bool, optional, default = True
            If True, plot the image to the screen.
        ax : matplotlib axes object, optional, default = None
            An array of matplotlib axes objects where the plots will appear.
        legend : bool, optional, default = True
            If True, plot a legend in each subplot providing the taper number
            and 1 minus the concentration factor.
        fname : str, optional, default = None
            If present, save the image to the specified file.
        """
        if self.kind == 'cap':
            if self.nwinrot is not None and self.nwinrot <= nwin:
                nwin = self.nwinrot

        ncolumns = min(maxcolumns, nwin)
        nrows = _np.ceil(nwin / ncolumns).astype(int)
        figsize = ncolumns * 2.4, nrows * 1.2 + 0.5

        if ax is None:
            fig, axes = _plt.subplots(nrows, ncolumns, figsize=figsize)
        else:
            if hasattr(ax, 'flatten') and ax.size < nwin:
                raise ValueError('ax.size must be greater or equal to nwin. ' +
                                 'nwin = {:s}'.format(repr(nwin)) +
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

        if ax is None:
            if nrows > 1:
                for axtemp in axes[:-1, :].flatten():
                    for xlabel_i in axtemp.get_xticklabels():
                        xlabel_i.set_visible(False)
                for axtemp in axes[:, 1:].flatten():
                    for ylabel_i in axtemp.get_yticklabels():
                        ylabel_i.set_visible(False)
            elif nwin > 1:
                for axtemp in axes[1:].flatten():
                    for ylabel_i in axtemp.get_yticklabels():
                        ylabel_i.set_visible(False)

        for itaper in range(min(self.nwin, nwin)):
            evalue = self.eigenvalues[itaper]
            if min(self.nwin, nwin) == 1 and ax is None:
                axtemp = axes
            elif hasattr(axes, 'flatten'):
                axtemp = axes.flatten()[itaper]
            else:
                axtemp = axes[itaper]
            gridout = _shtools.MakeGridDH(self.to_array(itaper), sampling=2,
                                          lmax=lmax, norm=1, csphase=1)
            axtemp.imshow(gridout, origin='upper',
                          extent=(0., 360., -90., 90.))
            axtemp.set(xlabel=xlabel, ylabel=ylabel, xticks=xticks,
                       yticks=yticks)
            if legend is True:
                axtemp.text(0.02, 0.95,
                            '#{:d} [loss={:2.2g}]'.format(itaper, 1-evalue),
                            transform=axtemp.transAxes,
                            va='top', ha='left', color='black',
                            bbox={'boxstyle': 'round',  'edgecolor': 'none',
                                  'facecolor': 'white'})
        if ax is None:
            if show:
                _plt.show()

            fig.tight_layout(pad=0.5)
            if fname is not None:
                fig.savefig(fname)
            return fig, axes

    def plot_spectra(self, nwin, convention='power', unit='per_l', base=10.,
                     maxcolumns=5, xscale='lin', yscale='log',
                     xlim=(None, None), ylim=(None, None), show=True, ax=None,
                     fname=None):
        """
        Plot the spectra of the best-concentrated localization windows.

        Usage
        -----
        x.plot_spectra(nwin, [convention, unit, base, maxcolumns, xscale,
                              yscale, xlim, ylim, show, ax, fname])

        Parameters
        ----------
        nwin : int
            The number of localization windows to plot.
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
        maxcolumns : int, optional, default = 5
            The maximum number of columns to use when plotting the spectra
            of multiple localization windows.
        xscale : str, optional, default = 'lin'
            Scale of the x axis: 'lin' for linear or 'log' for logarithmic.
        yscale : str, optional, default = 'log'
            Scale of the y axis: 'lin' for linear or 'log' for logarithmic.
        xlim : tuple, optional, default = (None, None)
            The upper and lower limits used for the x axis.
        ylim : tuple, optional, default = (None, None)
            The lower and upper limits used for the y axis.
        show : bool, optional, default = True
            If True, plot the image to the screen.
        ax : matplotlib axes object, optional, default = None
            An array of matplotlib axes objects where the plots will appear.
        fname : str, optional, default = None
            If present, save the image to the file.
        """
        degrees = self.degrees()
        spectrum = self.spectra(nwin=nwin, convention=convention, unit=unit,
                                base=base)

        ncolumns = min(maxcolumns, nwin)
        nrows = _np.ceil(nwin / ncolumns).astype(int)
        figsize = ncolumns * 2.4, nrows * 1.2 + 0.5

        if ylim == (None, None):
            upper = spectrum[:, :min(self.nwin, nwin)].max()
            lower = upper * 1.e-6
            ylim = (lower, 5 * upper)

        if ax is None:
            fig, axes = _plt.subplots(nrows, ncolumns, figsize=figsize)
        else:
            if hasattr(ax, 'flatten') and ax.size < nwin:
                raise ValueError('ax.size must be greater or equal to nwin. ' +
                                 'nwin = {:s}'.format(repr(nwin)) +
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
            elif nwin > 1:
                for axtemp in axes[1:].flatten():
                    for ylabel_i in axtemp.get_yticklabels():
                        ylabel_i.set_visible(False)
                    axtemp.set_ylabel('', visible=False)

        for itaper in range(min(self.nwin, nwin)):
            evalue = self.eigenvalues[itaper]
            if min(self.nwin, nwin) == 1 and ax is None:
                axtemp = axes
            elif hasattr(axes, 'flatten'):
                axtemp = axes.flatten()[itaper]
            else:
                axtemp = axes[itaper]
            if (convention == 'power'):
                axtemp.set_ylabel('power')
            else:
                axtemp.set_ylabel('energy')

            if yscale == 'log':
                axtemp.set_yscale('log', basey=base)

            if xscale == 'log':
                axtemp.set_xscale('log', basex=base)
                axtemp.plot(degrees[1:], spectrum[1:, itaper])
            else:
                axtemp.plot(degrees[0:], spectrum[0:, itaper],
                            label='#{:d}'.format(itaper))
            axtemp.set(xlabel='degree l', xlim=xlim, ylim=ylim)
            axtemp.grid(True, which='both')
            axtemp.set_title('#{:d} [loss={:2.2g}]'.format(itaper, 1-evalue))

        if ax is None:
            if show:
                _plt.show()

            fig.tight_layout(pad=0.5)
            if fname is not None:
                fig.savefig(fname)
            return fig, axes

    def plot_coupling_matrix(self, lmax, nwin=None, weights=None, mode='full',
                             show=True, ax=None, fname=None):
        """
        Plot the multitaper coupling matrix.

        This matrix relates the global power spectrum to the expectation of
        the localized multitaper spectrum.

        Usage
        -----
        x.plot_coupling_matrix(lmax, [nwin, weights, mode, show, ax, fname])

        Parameters
        ----------
        lmax : int
            Spherical harmonic bandwidth of the global power spectrum.
        nwin : int, optional, default = x.nwin
            Number of tapers used in the mutlitaper spectral analysis.
        weights : ndarray, optional, default = x.weights
            Taper weights used with the multitaper spectral analyses.
        mode : str, opitonal, default = 'full'
            'full' returns a biased output spectrum of size lmax+lwin+1. The
            input spectrum is assumed to be zero for degrees l>lmax.
            'same' returns a biased output spectrum with the same size
            (lmax+1) as the input spectrum. The input spectrum is assumed to be
            zero for degrees l>lmax.
            'valid' returns a biased spectrum with size lmax-lwin+1. This
            returns only that part of the biased spectrum that is not
            influenced by the input spectrum beyond degree lmax.
        show : bool, optional, default = True
            If True, plot the image to the screen.
        ax : matplotlib axes object, optional, default = None
            An array of matplotlib axes objects where the plots will appear.
        fname : str, optional, default = None
            If present, save the image to the specified file.
        """
        figsize = _mpl.rcParams['figure.figsize']
        figsize[0] = figsize[1]

        if ax is None:
            fig = _plt.figure(figsize=figsize)
            axes = fig.add_subplot(111)
        else:
            axes = ax

        axes.imshow(self.coupling_matrix(lmax, nwin=nwin, weights=weights,
                                         mode=mode), aspect='auto')
        axes.set_xlabel('input power')  # matrix index 1 (columns)
        axes.set_ylabel('output power')  # matrix index 0 (rows)

        if show:
            _plt.show()

        if ax is None:
            fig.tight_layout(pad=0.1)
            if fname is not None:
                fig.savefig(fname)
            return fig, axes

    def info(self):
        """
        Print a summary of the data stored in the SHWindow instance.

        Usage
        -----
        x.info()
        """
        self._info()


class SHWindowCap(SHWindow):
    """Class for localization windows concentrated within a spherical cap."""

    @staticmethod
    def istype(kind):
        return kind == 'cap'

    def __init__(self, theta, tapers, eigenvalues, taper_order,
                 clat, clon, nwin, theta_degrees, coord_degrees, dj_matrix,
                 weights, copy=True):
        self.kind = 'cap'
        self.theta = theta
        self.clat = clat
        self.clon = clon
        self.lwin = tapers.shape[0] - 1
        self.theta_degrees = theta_degrees
        self.coord_degrees = coord_degrees
        self.dj_matrix = dj_matrix
        self.weights = weights
        self.nwinrot = None

        if nwin is not None:
            self.nwin = nwin
        else:
            self.nwin = tapers.shape[1]

        if self.nwin > (self.lwin + 1)**2:
            raise ValueError('nwin must be less than or equal to ' +
                             '(lwin+1)**2. nwin = {:s} and lwin = {:s}.'
                             .format(repr(self.nwin), repr(self.lwin)))

        if copy:
            self.tapers = _np.copy(tapers[:, :self.nwin])
            self.eigenvalues = _np.copy(eigenvalues[:self.nwin])
            self.orders = _np.copy(taper_order[:self.nwin])
        else:
            self.tapers = tapers[:, :self.nwin]
            self.eigenvalues = eigenvalues[:self.nwin]
            self.orders = taper_order[:self.nwin]

        # If the windows aren't rotated, don't store them.
        if self.clat is None and self.clon is None:
            self.coeffs = None
        else:
            self.rotate(clat=self.clat, clon=self.clon,
                        coord_degrees=self.coord_degrees,
                        dj_matrix=self.dj_matrix)

    def _taper2coeffs(self, itaper):
        """
        Return the spherical harmonic coefficients of the unrotated taper i
        as an array, where i = 0 is the best concentrated.
        """
        taperm = self.orders[itaper]
        coeffs = _np.zeros((2, self.lwin + 1, self.lwin + 1))
        if taperm < 0:
            coeffs[1, :, abs(taperm)] = self.tapers[:, itaper]
        else:
            coeffs[0, :, abs(taperm)] = self.tapers[:, itaper]

        return coeffs

    def _to_array(self, itaper, normalization='4pi', csphase=1):
        """
        Return the spherical harmonic coefficients of taper i as an
        array, where i = 0 is the best concentrated.
        """
        if self.coeffs is None:
            coeffs = _np.copy(self._taper2coeffs(itaper))
        else:
            if itaper > self.nwinrot - 1:
                raise ValueError('itaper must be less than or equal to ' +
                                 'nwinrot - 1. itaper = {:d}, nwinrot = {:d}'
                                 .format(itaper, self.nwinrot))
            coeffs = _shtools.SHVectorToCilm(self.coeffs[:, itaper])

        if normalization == 'schmidt':
            for l in range(self.lwin + 1):
                coeffs[:, l, :l+1] *= _np.sqrt(2.0 * l + 1.0)
        elif normalization == 'ortho':
            coeffs *= _np.sqrt(4.0 * _np.pi)

        if csphase == -1:
            for m in range(self.lwin + 1):
                if m % 2 == 1:
                    coeffs[:, :, m] = - coeffs[:, :, m]

        return coeffs

    def rotate(self, clat, clon, coord_degrees=True, dj_matrix=None,
               nwinrot=None):
        """"
        Rotate the spherical-cap windows centered on the north pole to clat
        and clon, and save the spherical harmonic coefficients in the
        attribute coeffs.

        Usage
        -----
        x.rotate(clat, clon [coord_degrees, dj_matrix, nwinrot])

        Parameters
        ----------
        clat, clon : float
            Latitude and longitude of the center of the rotated spherical-cap
            localization windows (default in degrees).
        coord_degrees : bool, optional, default = True
            True if clat and clon are in degrees.
        dj_matrix : ndarray, optional, default = None
            The djpi2 rotation matrix computed by a call to djpi2.
        nwinrot : int, optional, default = (lwin+1)**2
            The number of best concentrated windows to rotate, where lwin is
            the spherical harmonic bandwidth of the localization windows.

        Description
        -----------
        This function will take the spherical-cap localization windows
        centered at the north pole (and saved in the attributes tapers and
        orders), rotate each function to the coordinate (clat, clon), and save
        the spherical harmonic coefficients in the attribute coeffs. Each
        column of coeffs contains a single window, and is ordered according to
        the convention in SHCilmToVector.
        """
        self.coeffs = _np.zeros(((self.lwin + 1)**2, self.nwin))
        self.clat = clat
        self.clon = clon
        self.coord_degrees = coord_degrees

        if nwinrot is not None:
            self.nwinrot = nwinrot
        else:
            self.nwinrot = self.nwin

        if self.coord_degrees:
            angles = _np.radians(_np.array([0., -(90. - clat), -clon]))
        else:
            angles = _np.array([0., -(_np.pi/2. - clat), -clon])

        if dj_matrix is None:
            if self.dj_matrix is None:
                self.dj_matrix = _shtools.djpi2(self.lwin + 1)
                dj_matrix = self.dj_matrix
            else:
                dj_matrix = self.dj_matrix

        for i in range(self.nwinrot):
            if ((coord_degrees is True and clat == 90. and clon == 0.) or
                    (coord_degrees is False and clat == _np.pi/2. and
                     clon == 0.)):
                coeffs = self._taper2coeffs(i)
                self.coeffs[:, i] = _shtools.SHCilmToVector(coeffs)

            else:
                coeffs = _shtools.SHRotateRealCoef(self._taper2coeffs(i),
                                                   angles, dj_matrix)
                self.coeffs[:, i] = _shtools.SHCilmToVector(coeffs)

    def _coupling_matrix(self, lmax, nwin=None, weights=None):
        """Return the coupling matrix of the first nwin tapers."""
        if nwin is None:
            nwin = self.nwin

        if weights is None:
            weights = self.weights

        if weights is None:
            return _shtools.SHMTCouplingMatrix(lmax, self.tapers**2, k=nwin)
        else:
            return _shtools.SHMTCouplingMatrix(lmax, self.tapers**2, k=nwin,
                                               taper_wt=self.weights)

    def _multitaper_spectrum(self, clm, k, convention='power', unit='per_l',
                             clat=None, clon=None, coord_degrees=True,
                             lmax=None, taper_wt=None):
        """
        Return the multitaper spectrum estimate and standard error for an
        input SHCoeffs class instance.
        """
        if lmax is None:
            lmax = clm.lmax

        if (clat is not None and clon is not None and clat == self.clat and
                clon == self.clon and coord_degrees is self.coord_degrees and
                k <= self.nwinrot):
            # use the already stored coeffs
            pass
        elif (clat is None and clon is None) and \
                (self.clat is not None and self.clon is not None and
                 k <= self.nwinrot):
            # use the already stored coeffs
            pass
        else:
            if clat is None:
                clat = self.clat
            if clon is None:
                clon = self.clon
            if (clat is None and clon is not None) or \
                    (clat is not None and clon is None):
                raise ValueError('clat and clon must both be input. ' +
                                 'clat = {:s}, clon = {:s}'
                                 .format(repr(clat), repr(clon)))
            if clat is None and clon is None:
                self.rotate(clat=90., clon=0., coord_degrees=True, nwinrot=k)
            else:
                self.rotate(clat=clat, clon=clon, coord_degrees=coord_degrees,
                            nwinrot=k)

        sh = clm.to_array(normalization='4pi', csphase=1, lmax=lmax)

        if taper_wt is None:
            mtse, sd = _shtools.SHMultiTaperMaskSE(sh, self.coeffs,
                                                   lmax=lmax, k=k)
        else:
            mtse, sd = _shtools.SHMultiTaperMaskSE(sh, self.coeffs, lmax=lmax,
                                                   k=k, taper_wt=taper_wt)

        if (unit == 'per_l'):
            pass
        elif (unit == 'per_lm'):
            degree_l = _np.arange(len(mtse))
            mtse /= (2.0 * degree_l + 1.0)
            sd /= (2.0 * degree_l + 1.0)
        else:
            raise ValueError(
                "unit must be 'per_l' or 'per_lm'." +
                "Input value was {:s}".format(repr(unit)))

        if (convention == 'power'):
            return mtse, sd
        elif (convention == 'energy'):
            return mtse * 4.0 * _np.pi, sd * 4.0 * _np.pi
        else:
            raise ValueError(
                "convention must be 'power' or 'energy'." +
                "Input value was {:s}".format(repr(convention)))

    def _multitaper_cross_spectrum(self, clm, slm, k, convention='power',
                                   unit='per_l', clat=None, clon=None,
                                   coord_degrees=True, lmax=None,
                                   taper_wt=None):
        """
        Return the multitaper cross-spectrum estimate and standard error for
        two input SHCoeffs class instances.
        """
        if lmax is None:
            lmax = min(clm.lmax, slm.lmax)

        if (clat is not None and clon is not None and clat == self.clat and
                clon == self.clon and coord_degrees is self.coord_degrees and
                k <= self.nwinrot):
            # use the already stored coeffs
            pass
        elif (clat is None and clon is None) and \
                (self.clat is not None and self.clon is not None and
                 k <= self.nwinrot):
            # use the already stored coeffs
            pass
        else:
            if clat is None:
                clat = self.clat
            if clon is None:
                clon = self.clon
            if (clat is None and clon is not None) or \
                    (clat is not None and clon is None):
                raise ValueError('clat and clon must both be input. ' +
                                 'clat = {:s}, clon = {:s}'
                                 .format(repr(clat), repr(clon)))
            if clat is None and clon is None:
                self.rotate(clat=90., clon=0., coord_degrees=True, nwinrot=k)
            else:
                self.rotate(clat=clat, clon=clon, coord_degrees=coord_degrees,
                            nwinrot=k)

        sh1 = clm.to_array(normalization='4pi', csphase=1, lmax=lmax)
        sh2 = slm.to_array(normalization='4pi', csphase=1, lmax=lmax)

        if taper_wt is None:
            mtse, sd = _shtools.SHMultiTaperMaskCSE(sh1, sh2, self.coeffs,
                                                    lmax1=lmax, lmax2=lmax,
                                                    k=k)
        else:
            mtse, sd = _shtools.SHMultiTaperMaskCSE(sh1, sh2, self.coeffs,
                                                    lmax1=lmax, lmax2=lmax,
                                                    k=k, taper_wt=taper_wt)

        if (unit == 'per_l'):
            pass
        elif (unit == 'per_lm'):
            degree_l = _np.arange(len(mtse))
            mtse /= (2.0 * degree_l + 1.0)
            sd /= (2.0 * degree_l + 1.0)
        else:
            raise ValueError(
                "unit must be 'per_l' or 'per_lm'." +
                "Input value was {:s}".format(repr(unit)))

        if (convention == 'power'):
            return mtse, sd
        elif (convention == 'energy'):
            return mtse * 4.0 * _np.pi, sd * 4.0 * _np.pi
        else:
            raise ValueError(
                "convention must be 'power' or 'energy'." +
                "Input value was {:s}".format(repr(convention)))

    def _biased_spectrum(self, spectrum, k, convention='power', unit='per_l',
                         **kwargs):
        """
        Calculate the multitaper (cross-) spectrum expectation of a function
        localized by spherical cap windows.
        """
        # The equation is not modified if the in- and out- spectra are power
        # or energy. However, the convention can not be l2norm, which depends
        # upon the normalization of the coefficients.
        if (convention != 'power' and convention != 'energy'):
            raise ValueError(
                "convention must be 'power' or 'energy'." +
                "Input value was {:s}".format(repr(convention)))

        if (unit == 'per_l'):
            outspectrum = _shtools.SHBiasK(self.tapers, spectrum, k=k,
                                           **kwargs)
        elif (unit == 'per_lm'):
            degree_l = _np.arange(len(spectrum))
            temp = spectrum * (2.0 * degree_l + 1.0)
            outspectrum = _shtools.SHBiasK(self.tapers, temp, k=k,
                                           **kwargs)
            outspectrum /= (2.0 * degree_l + 1.0)
        else:
            raise ValueError(
                "unit must be 'per_l' or 'per_lm'." +
                "Input value was {:s}".format(repr(unit)))

        return outspectrum

    def _info(self):
        """Print a summary of the data in the SHWindow instance."""
        print('kind = {:s}\n'.format(repr(self.kind)), end='')

        if self.theta_degrees:
            print('theta = {:f} degrees\n'.format(self.theta), end='')
        else:
            print('theta = {:f} radians'.format(self.theta), end='')

        print('lwin = {:d}\n'.format(self.lwin), end='')
        print('nwin = {:d}\n'.format(self.nwin), end='')

        if self.clat is not None:
            if self.coord_degrees:
                print('clat = {:f} degrees\n'.format(self.clat), end='')
            else:
                print('clat = {:f} radians\n'.format(self.clat), end='')
        else:
            print('clat is not specified')

        if self.clon is not None:
            if self.coord_degrees:
                print('clon = {:f} degrees\n'.format(self.clon), end='')
            else:
                print('clon = {:f} radians\n'.format(self.clon), end='')
        else:
            print('clon is not specified')

        print('nwinrot = {:s}'.format(repr(self.nwinrot)))

        if self.dj_matrix is not None:
            print('dj_matrix is stored')
        else:
            print('dj_matrix is not stored')

        if self.weights is None:
            print('Taper weights are not set.')
        else:
            print('Taper weights are set.')


class SHWindowMask(SHWindow):
    """
    Class for localization windows concentrated within a specified mask and
    for a given spherical harmonic bandwidth.
    """

    @staticmethod
    def istype(kind):
        return kind == 'mask'

    def __init__(self, tapers, eigenvalues, weights, copy=True):
        self.kind = 'mask'
        self.lwin = _np.sqrt(tapers.shape[0]).astype(int) - 1
        self.nwin = tapers.shape[1]
        if copy:
            self.weights = weights
            self.tapers = _np.copy(tapers)
            self.eigenvalues = _np.copy(eigenvalues)
        else:
            self.weights = weights
            self.tapers = tapers
            self.eigenvalues = eigenvalues

    def _to_array(self, itaper, normalization='4pi', csphase=1):
        """
        Return the spherical harmonic coefficients of taper i as an
        array, where i=0 is the best concentrated.
        """
        coeffs = _shtools.SHVectorToCilm(self.tapers[:, itaper])

        if normalization == 'schmidt':
            for l in range(self.lwin + 1):
                coeffs[:, l, :l+1] *= _np.sqrt(2.0 * l + 1.0)
        elif normalization == 'ortho':
            coeffs *= _np.sqrt(4.0 * _np.pi)

        if csphase == -1:
            for m in range(self.lwin + 1):
                if m % 2 == 1:
                    coeffs[:, :, m] = - coeffs[:, :, m]

        return coeffs

    def _coupling_matrix(self, lmax, nwin=None, weights=None):
        """Return the coupling matrix of the first nwin tapers."""
        if nwin is None:
            nwin = self.nwin

        if weights is None:
            weights = self.weights

        tapers_power = _np.zeros((self.lwin+1, nwin))
        for i in range(nwin):
            tapers_power[:, i] = _spectrum(self.to_array(i),
                                           normalization='4pi',
                                           convention='power', unit='per_l')

        if weights is None:
            return _shtools.SHMTCouplingMatrix(lmax, tapers_power, k=nwin)
        else:
            return _shtools.SHMTCouplingMatrix(lmax, tapers_power, k=nwin,
                                               taper_wt=self.weights)

    def _multitaper_spectrum(self, clm, k, convention='power', unit='per_l',
                             lmax=None, taper_wt=None):
        """
        Return the multitaper spectrum estimate and standard error for an
        input SHCoeffs class instance.
        """
        if lmax is None:
            lmax = clm.lmax

        sh = clm.to_array(normalization='4pi', csphase=1, lmax=lmax)

        if taper_wt is None:
            mtse, sd = _shtools.SHMultiTaperMaskSE(sh, self.tapers, lmax=lmax,
                                                   k=k)
        else:
            mtse, sd = _shtools.SHMultiTaperMaskSE(sh, self.tapers, lmax=lmax,
                                                   k=k, taper_wt=taper_wt)

        if (unit == 'per_l'):
            pass
        elif (unit == 'per_lm'):
            degree_l = _np.arange(len(mtse))
            mtse /= (2.0 * degree_l + 1.0)
            sd /= (2.0 * degree_l + 1.0)
        else:
            raise ValueError(
                "unit must be 'per_l' or 'per_lm'." +
                "Input value was {:s}".format(repr(unit)))

        if (convention == 'power'):
            return mtse, sd
        elif (convention == 'energy'):
            return mtse * 4.0 * _np.pi, sd * 4.0 * _np.pi
        else:
            raise ValueError(
                "convention must be 'power' or 'energy'." +
                "Input value was {:s}".format(repr(convention)))

    def _multitaper_cross_spectrum(self, clm, slm, k, convention='power',
                                   unit='per_l', lmax=None, taper_wt=None):
        """
        Return the multitaper cross-spectrum estimate and standard error for
        two input SHCoeffs class instances.
        """
        if lmax is None:
            lmax = min(clm.lmax, slm.lmax)

        sh1 = clm.to_array(normalization='4pi', csphase=1, lmax=lmax)
        sh2 = slm.to_array(normalization='4pi', csphase=1, lmax=lmax)

        if taper_wt is None:
            mtse, sd = _shtools.SHMultiTaperMaskCSE(sh1, sh2, self.tapers,
                                                    lmax=lmax, k=k)
        else:
            mtse, sd = _shtools.SHMultiTaperMaskCSE(sh1, sh2, self.tapers,
                                                    lmax=lmax, k=k,
                                                    taper_wt=taper_wt)

        if (unit == 'per_l'):
            pass
        elif (unit == 'per_lm'):
            degree_l = _np.arange(len(mtse))
            mtse /= (2.0 * degree_l + 1.0)
            sd /= (2.0 * degree_l + 1.0)
        else:
            raise ValueError(
                "unit must be 'per_l' or 'per_lm'." +
                "Input value was {:s}".format(repr(unit)))

        if (convention == 'power'):
            return mtse, sd
        elif (convention == 'energy'):
            return mtse * 4.0 * _np.pi, sd * 4.0 * _np.pi
        else:
            raise ValueError(
                "convention must be 'power' or 'energy'." +
                "Input value was {:s}".format(repr(convention)))

    def _biased_spectrum(self, spectrum, k, convention='power', unit='per_l',
                         **kwargs):
        """
        Calculate the multitaper (cross-) spectrum expectation of a function
        localized by arbitary windows.
        """
        # The equation is not modified if the in- and out- spectra are power
        # or energy. However, the convention can not be l2norm, which depends
        # upon the normalization of the coefficients.
        if (convention != 'power' and convention != 'energy'):
            raise ValueError(
                "convention must be 'power' or 'energy'." +
                "Input value was {:s}".format(repr(convention)))

        if (unit == 'per_l'):
            outspectrum = _shtools.SHBiasKMask(self.tapers, spectrum, k=k,
                                               **kwargs)
        elif (unit == 'per_lm'):
            degree_l = _np.arange(len(spectrum))
            temp = spectrum * (2.0 * degree_l + 1.0)
            outspectrum = _shtools.SHBiasKMask(self.tapers, temp, k=k,
                                               **kwargs)
            outspectrum /= (2.0 * degree_l + 1.0)
        else:
            raise ValueError(
                "unit must be 'per_l' or 'per_lm'." +
                "Input value was {:s}".format(repr(unit)))

        return outspectrum

    def _info(self):
        """Print a summary of the data in the SHWindow instance."""
        print('kind = {:s}\n'.format(repr(self.kind)), end='')

        print('lwin = {:d}\n'.format(self.lwin), end='')
        print('nwin = {:d}\n'.format(self.nwin), end='')

        if self.weights is None:
            print('Taper weights are not set.')
        else:
            print('Taper weights are set.')
