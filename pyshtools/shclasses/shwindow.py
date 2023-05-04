"""
    Class for localized spectral analyses on the sphere.

        SHWindow: SHWindowCap, SHWindowMask
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


__all__ = ['SHWindow', 'SHWindowCap', 'SHWindowMask']


class SHWindow(object):
    """
    Class for localized spectral analyses on the sphere.

    The windows can be initialized from:

    >>>  x = SHWindow.from_cap()
    >>>  x = SHWindow.from_mask()

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
    shannon         : The Shannon number, which approximates the number of
                      well localized windows.
    area            : Area of the concentration domain, in radians.
    eigenvalues     : Concentration factors of the localization windows.
    orders          : The angular orders for each of the spherical cap
                      localization windows.
    weights         : Taper weights used with the multitaper spectral analyses.
                      Defaut is None.
    lwin            : Spherical harmonic bandwidth of the localization windows.
    theta           : Angular radius of the spherical cap localization domain
                      (default in degrees).
    theta_degrees   : True (default) if theta is in degrees.
    nwin            : The number of localization windows. Default is
                      (lwin+1)**2.
    nwinrot         : The number of best-concentrated spherical cap
                      localization windows that were rotated
                      and whose coefficients are stored in coeffs.
    clat, clon      : Latitude and longitude of the center of the rotated
                      spherical cap localization windows (default in degrees).
    coord_degrees   : True (default) if clat and clon are in degrees.
    taper_degrees   : Boolean or int array defining which spherical harmonic
                      degrees were used to construct the windows.

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
    variance()             : Compute the theoretical variance of a windowed
                             function for a given input power spectrum.
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
        print('Initialize the class using one of the class methods:\n'
              '>>> pyshtools.SHWindow.from_cap\n'
              '>>> pyshtools.SHWindow.from_mask')

    # ---- factory methods:
    @classmethod
    def from_cap(cls, theta, lwin, clat=None, clon=None, nwin=None,
                 theta_degrees=True, coord_degrees=True, dj_matrix=None,
                 weights=None, taper_degrees=None):
        """
        Construct spherical cap localization windows.

        Usage
        -----
        x = SHWindow.from_cap(theta, lwin, [clat, clon, nwin, theta_degrees,
                                            coord_degrees, dj_matrix, weights,
                                            taper_degrees])

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
        nwin : int, optional, default = (lwin+1)**2
            Number of localization windows.
        theta_degrees : bool, optional, default = True
            True if theta is in degrees.
        coord_degrees : bool, optional, default = True
            True if clat and clon are in degrees.
        dj_matrix : ndarray, optional, default = None
            The djpi2 rotation matrix computed by a call to djpi2.
        weights : ndarray, optional, default = None
            Taper weights used with the multitaper spectral analyses.
        taper_degrees : bool or int, optional, dimension (lmax+1),
                        default = None
            Boolean or int array defining which spherical harmonic degrees were
            used (True or 1) to construct the windows.
        """
        if theta_degrees:
            tapers, eigenvalues, taper_order = _shtools.SHReturnTapers(
                _np.radians(theta), lwin, degrees=taper_degrees)
        else:
            tapers, eigenvalues, taper_order = _shtools.SHReturnTapers(
                theta, lwin, degrees=taper_degrees)

        return SHWindowCap(theta, tapers, eigenvalues, taper_order,
                           clat, clon, nwin, theta_degrees, coord_degrees,
                           dj_matrix, weights, taper_degrees, copy=False)

    @classmethod
    def from_mask(cls, dh_mask, lwin, nwin=None, weights=None,
                  taper_degrees=None):
        """
        Construct localization windows that are optimally concentrated within
        the region specified by a mask.

        Usage
        -----
        x = SHWindow.from_mask(dh_mask, lwin, [nwin, weights, taper_degrees])

        Returns
        -------
        x : SHWindow class instance

        Parameters
        ----------
        dh_mask :ndarray or SHGrid class instance, shape (nlat, nlon)
            A Driscoll and Healy sampled grid describing the concentration
            region R. All elements should either be 1 or 0 for inside or
            outside of the concentration region, respectively. The grid must
            have dimensions nlon=nlat, nlon=2*nlat, or nlon=2*nlat-1.
        lwin : int
            The spherical harmonic bandwidth of the localization windows.
        nwin : int, optional, default = (lwin+1)**2
            The number of best-concentrated eigenvalues and eigenfunctions to
            return.
        weights : ndarray, optional, default = None
            Taper weights used with the multitaper spectral analyses.
        taper_degrees : bool or int, optional, dimension (lmax+1),
                        default = None
            Boolean or int array defining which spherical harmonic degrees were
            used (True or 1) to construct the windows.
        """
        if nwin is None:
            nwin = (lwin + 1)**2
        else:
            if nwin > (lwin + 1)**2:
                raise ValueError('nwin must be less than or equal to ' +
                                 '(lwin + 1)**2. lwin = {:d} and nwin = {:d}.'
                                 .format(lwin, nwin))

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
            data, lwin, ntapers=nwin, degrees=taper_degrees)

        return SHWindowMask(tapers, eigenvalues, weights, area, taper_degrees,
                            copy=False)

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

    def number_concentrated(self, concentration):
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
        return len(self.eigenvalues[self.eigenvalues >= concentration])

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
            3-D numpy ndarray of the spherical harmonic coefficients of the
            window.

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
                "csphase must be 1 or -1. Input value is {:s}."
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

        coeffs = self.to_array(itaper, normalization=normalization.lower(),
                               csphase=csphase)
        return SHCoeffs.from_array(coeffs, normalization=normalization.lower(),
                                   csphase=csphase, copy=False)

    def to_shgrid(self, itaper, grid='DH2', lmax=None, zeros=None,
                  extend=True):
        """
        Evaluate the coefficients of taper i on a spherical grid and return
        a SHGrid class instance.

        Usage
        -----
        f = x.to_shgrid(itaper, [grid, lmax, zeros, extend])

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
        lmax : int, optional, default = x.lwin
            The maximum spherical harmonic degree, which determines the grid
            spacing of the output grid.
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
        if lmax is None:
            lmax = self.lwin
        if type(grid) != str:
            raise ValueError('grid must be a string. Input type is {:s}.'
                             .format(str(type(grid))))

        if grid.upper() in ('DH', 'DH1'):
            gridout = _shtools.MakeGridDH(self.to_array(itaper), sampling=1,
                                          lmax=lmax, norm=1, csphase=1,
                                          extend=extend)
            return SHGrid.from_array(gridout, grid='DH', copy=False)
        elif grid.upper() == 'DH2':
            gridout = _shtools.MakeGridDH(self.to_array(itaper), sampling=2,
                                          lmax=lmax, norm=1, csphase=1,
                                          extend=extend)
            return SHGrid.from_array(gridout, grid='DH', copy=False)
        elif grid.upper() == 'GLQ':
            if zeros is None:
                zeros, weights = _shtools.SHGLQ(self.lmax)
            gridout = _shtools.MakeGridGLQ(self.to_array(itaper), zeros,
                                           lmax=lmax, norm=1, csphase=1,
                                           extend=extend)
            return SHGrid.from_array(gridout, grid='GLQ', copy=False)
        else:
            raise ValueError(
                "grid must be 'DH', 'DH1', 'DH2', or 'GLQ'. " +
                "Input value is {:s}.".format(repr(grid)))

    def multitaper_spectrum(self, clm, k, convention='power', unit='per_l',
                            lmax=None, weights=None, clat=None, clon=None,
                            coord_degrees=True):
        """
        Return the multitaper spectrum estimate and standard error.

        Usage
        -----
        mtse, sd = x.multitaper_spectrum(clm, k, [convention, unit, lmax,
                                                  clat, clon, weights
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
        weights : ndarray, optional, dimension (k), default = x.weights
            1-D numpy array of the weights used in calculating the multitaper
            spectral estimates and standard error.
        clat, clon : float, optional, default = 90., 0.
            When using spherical-cap localization windows, rotate the center of
            the localization windows to the latitude and longitude coordinates
            clat and clon, respectively.
        coord_degrees : bool, optional, default = True
            True if clat and clon are in degrees.
        """
        if weights is not None:
            if len(weights) != k:
                raise ValueError('Length of weights must be equal to k. '
                                 'len(weights) = {:d}, k = {:d}.'
                                 .format(len(weights), k))
        else:
            weights = self.weights

        return self._multitaper_spectrum(clm, k, convention=convention,
                                         unit=unit, lmax=lmax, weights=weights,
                                         clat=clat, clon=clon,
                                         coord_degrees=coord_degrees)

    def multitaper_cross_spectrum(self, clm, slm, k, convention='power',
                                  unit='per_l', lmax=None, weights=None,
                                  clat=None, clon=None, coord_degrees=True):
        """
        Return the multitaper cross-spectrum estimate and standard error.

        Usage
        -----
        mtse, sd = x.multitaper_cross_spectrum(clm, slm, k, [convention, unit,
                                                             lmax, weights,
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
        weights : ndarray, optional, dimension (k), default = x.weights
            The weights used in calculating the multitaper cross-spectral
            estimates and standard error.
        clat, clon : float, optional, default = 90., 0.
            When using spherical-cap localization windows, rotate the center of
            the localization windows to the latitude and longitude coordinates
            clat and clon, respectively.
        coord_degrees : bool, optional, default = True
            True if clat and clon are in degrees.
        """
        if weights is not None:
            if len(weights) != k:
                raise ValueError('Length of weights must be equal to k. '
                                 'len(weights) = {:d}, k = {:d}.'
                                 .format(len(weights), k))
        else:
            weights = self.weights

        return self._multitaper_cross_spectrum(clm, slm, k,
                                               convention=convention,
                                               unit=unit, lmax=lmax,
                                               weights=weights,
                                               clat=clat, clon=clon,
                                               coord_degrees=coord_degrees)

    def biased_spectrum(self, power, k, convention='power', unit='per_l',
                        weights=None, save_cg=None, ldata=None):
        """
        Calculate the multitaper (cross-)spectrum expectation of a
        localized function.

        Usage
        -----
        outspectrum = x.biased_spectrum(spectrum, k, [convention, unit,
                                                      weights, save_cg, ldata])

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
        weights : ndarray, optional, dimension (k), default = x.weights
            The weights used in calculating the multitaper spectral estimates
            and standard error.
        save_cg : int, optional, default = 0
            If 1, the Clebsch-Gordon coefficients will be precomputed and saved
            for future use. If 0, the Clebsch-Gordon coefficients will be
            recomputed for each call.
        ldata : int, optional, default = len(power)-1
            The maximum degree of the global unwindowed spectrum.
        """
        if weights is not None:
            if len(weights) != k:
                raise ValueError('Length of weights must be equal to k. '
                                 'len(weights) = {:d}, k = {:d}.'
                                 .format(len(weights), k))
        else:
            weights = self.weights

        return self._biased_spectrum(power, k, convention=convention,
                                     unit=unit, weights=weights,
                                     save_cg=save_cg, ldata=ldata)

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

        Notes
        -----
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

    def coupling_matrix(self, lmax, k=None, weights=None, mode='full'):
        """
        Return the coupling matrix of the first nwin tapers. This matrix
        relates the global power spectrum to the expectation of the localized
        multitaper spectrum.

        Usage
        -----
        Mmt = x.coupling_matrix(lmax, [k, weights, mode])

        Returns
        -------
        Mmt : ndarray, shape (lmax+lwin+1, lmax+1) or (lmax+1, lmax+1) or
              (lmax-lwin+1, lmax+1)

        Parameters
        ----------
        lmax : int
            Spherical harmonic bandwidth of the global power spectrum.
        k : int, optional, default = x.nwin
            Number of tapers used in the mutlitaper spectral analysis.
        weights : ndarray, optional, dimension (k), default = x.weights
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
            if k is not None:
                if len(weights) != k:
                    raise ValueError(
                        'Length of weights must be equal to k. '
                        'len(weights) = {:d}, k = {:d}.'
                        .format(len(weights), k))
            else:
                if len(weights) != self.nwin:
                    raise ValueError(
                        'Length of weights must be equal to nwin when k is '
                        'not specified. len(weights) = {:d}, nwin = {:d}.'
                        .format(len(weights), self.nwin))
        else:
            weights = self.weights

        if mode == 'full':
            return self._coupling_matrix(lmax, k=k, weights=weights)
        elif mode == 'same':
            cmatrix = self._coupling_matrix(lmax, k=k, weights=weights)
            return cmatrix[:lmax+1, :]
        elif mode == 'valid':
            cmatrix = self._coupling_matrix(lmax, k=k, weights=weights)
            return cmatrix[:lmax - self.lwin+1, :]
        else:
            raise ValueError("mode has to be 'full', 'same' or 'valid', not "
                             "{}.".format(mode))

    def variance(self, power, k, lmax=None, weights=None):
        """
        Compute the theoretical variance of a windowed function for a given
        input power spectrum (using spherical-cap localization windows only).

        Usage
        -----
        variance = x.variance(power, k, [lmax, weights])

        Returns
        -------
        variance : ndarray, shape (min(lmax+1, lmax_in+1-lwin))
            The theoretical variance of the windowed function.

        Parameters
        ----------
        power : ndarray, dimension (lmax_in+1)
            The input global power spectrum.
        k : int
            The number of tapers to be utilized in performing the multitaper
            spectral analysis.
        lmax : int, optional, default = lmax_in
            The maximum spherical harmonic degree of the variance to compute.
        weights : ndarray, optional, dimension (k), default = x.weights
            Taper weights used with the multitaper spectral analyses.
        """
        if weights is not None:
            if len(weights) != k:
                raise ValueError(
                    'Length of weights must be equal to k. '
                    'len(weights) = {:d}, k = {:d}.'
                    .format(len(weights), k))
        else:
            weights = self.weights

        if lmax is None:
            lmax = len(power) - 1 - self.lwin
        else:
            if lmax > len(power) - 1 - self.lwin:
                raise ValueError('lmax must be less than or equal to '
                                 'len(power) - 1 - lwin.')

        return self._variance(power, k, lmax=lmax, weights=weights)

    def plot_windows(self, nwin, projection=None, lmax=None, maxcolumns=3,
                     tick_interval=[60, 45], minor_tick_interval=[None, None],
                     ticks='WSen', xlabel='Longitude', ylabel='Latitude',
                     title=True, colorbar=None, cmap='viridis',
                     cmap_limits=None, cmap_reverse=False, cb_offset=None,
                     cb_triangles='neither', cb_label=None, cb_ylabel=None,
                     cb_tick_interval=None, cb_minor_tick_interval=None,
                     grid=False, loss=False, axes_labelsize=None,
                     tick_labelsize=None, titlesize=9, show=True, ax=None,
                     cb_width=None, fname=None):
        """
        Plot the best-concentrated localization windows.

        Usage
        -----
        x.plot_windows(nwin, [projection, lmax, maxcolumns, tick_interval,
                              minor_tick_interval, ticks, xlabel, ylabel,
                              title, colorbar, cmap, cmap_limits, cmap_reverse,
                              cb_triangles, cb_label, cb_ylabel,
                              cb_tick_interval, cb_minor_tick_interval,
                              cb_offset, cb_width, grid, loss, titlesize,
                              axes_labelsize, tick_labelsize, ax, show, fname])

        Parameters
        ----------
        nwin : int
            The number of localization windows to plot.
        projection : Cartopy projection class, optional, default = None
            The Cartopy projection class used to project the gridded data,
            for Driscoll and Healy sampled grids only.
        lmax : int, optional, default = self.lwin
            The maximum degree to use when plotting the windows, which controls
            the number of samples in latitude and longitude.
        maxcolumns : int, optional, default = 3
            The maximum number of columns to use when plotting multiple
            localization windows.
        tick_interval : list or tuple, optional, default = [60, 45]
            Intervals to use when plotting the x and y ticks. If set to None,
            ticks will not be plotted.
        minor_tick_interval : list or tuple, optional, default = [None, None]
            Intervals to use when plotting the minor x and y ticks. If set to
            None, minor ticks will not be plotted.
        ticks : str, optional, default = 'WSen'
            Specify which axes should have ticks drawn and annotated. Capital
            letters will plot the ticks and annotations, whereas small letters
            will plot only the ticks. 'W', 'S', 'E', and 'N' denote the west,
            south, east and north boundaries of the plot.
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
            If True, plot grid lines.
        loss : bool, optional, default = False
            When plotting titles, provide the loss factor instead of the
            concentration factor (loss=1-concentration).
        titlesize : int, optional, default = 9
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
            if self.nwinrot is not None and self.nwinrot <= nwin:
                nwin = self.nwinrot

        ncolumns = min(maxcolumns, nwin)
        nrows = _np.ceil(nwin / ncolumns).astype(int)
        figsize = (_mpl.rcParams['figure.figsize'][0],
                   _mpl.rcParams['figure.figsize'][0]
                   * 0.6 * nrows / ncolumns + 0.41)

        if ax is None:
            fig, axes = _plt.subplots(nrows, ncolumns, figsize=figsize,
                                      sharex='all', sharey='all')
        else:
            if hasattr(ax, 'flatten') and ax.size < nwin:
                raise ValueError('ax.size must be greater or equal to nwin. ' +
                                 'nwin = {:s}'.format(repr(nwin)) +
                                 ' and ax.size = {:s}.'.format(repr(ax.size)))
            axes = ax

        for itaper in range(min(self.nwin, nwin)):
            evalue = self.eigenvalues[itaper]
            if min(self.nwin, nwin) == 1 and ax is None:
                axtemp = axes
            elif hasattr(axes, 'flatten'):
                axtemp = axes.flatten()[itaper]
            else:
                axtemp = axes[itaper]
            coeffs = self.to_shcoeffs(itaper)
            if lmax is not None:
                coeffs = coeffs.pad(lmax=lmax, copy=False)
            grid_temp = coeffs.expand()

            if title:
                if loss:
                    title_str = '#{:d} [loss={:2.2g}]'.format(itaper, 1-evalue)
                else:
                    title_str = '#{:d} [concentration={:2.2g}]'.format(
                        itaper, evalue)
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
            elif nwin > 1:
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

    def plot_spectra(self, nwin, convention='power', unit='per_l', base=10.,
                     maxcolumns=3, xscale='lin', yscale='log', grid=True,
                     xlim=(None, None), ylim=(None, None), show=True,
                     title=True, axes_labelsize=None, tick_labelsize=None,
                     loss=False, titlesize=9, ax=None, fname=None):
        """
        Plot the spectra of the best-concentrated localization windows.

        Usage
        -----
        x.plot_spectra(nwin, [convention, unit, base, maxcolumns, xscale,
                              yscale, grid, xlim, ylim, show, title, loss,
                              axes_labelsize, tick_labelsize, titlesize,
                              ax, fname])

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
        loss : bool, optional, default = False
            When plotting titles, provide the loss factor instead of the
            concentration factor (loss=1-concentration).
        axes_labelsize : int, optional, default = None
            The font size for the x and y axes labels.
        tick_labelsize : int, optional, default = None
            The font size for the x and y tick labels.
        titlesize : int, optional, default = 9
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
        spectrum = self.spectra(nwin=nwin, convention=convention, unit=unit,
                                base=base)

        ncolumns = min(maxcolumns, nwin)
        nrows = _np.ceil(nwin / ncolumns).astype(int)
        figsize = (_mpl.rcParams['figure.figsize'][0],
                   _mpl.rcParams['figure.figsize'][0]
                   * 0.8 * nrows / ncolumns + 0.41)

        if ax is None:
            fig, axes = _plt.subplots(nrows, ncolumns, figsize=figsize,
                                      sharex='all', sharey='all')
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

        if ylim == (None, None):
            upper = spectrum[:, :min(self.nwin, nwin)].max()
            lower = upper * 1.e-6
            ylim = (lower, 5 * upper)
        if xlim == (None, None):
            if xscale == 'lin':
                xlim = (degrees[0], degrees[-1])

        for itaper in range(min(self.nwin, nwin)):
            evalue = self.eigenvalues[itaper]
            if min(self.nwin, nwin) == 1 and ax is None:
                axtemp = axes
            elif hasattr(axes, 'flatten'):
                axtemp = axes.flatten()[itaper]
            else:
                axtemp = axes[itaper]
            if (convention == 'power'):
                axtemp.set_ylabel('Power', fontsize=axes_labelsize)
            else:
                axtemp.set_ylabel('Energy', fontsize=axes_labelsize)

            if yscale == 'log':
                axtemp.set_yscale('log', base=base)

            if xscale == 'log':
                axtemp.set_xscale('log', base=base)
                axtemp.plot(degrees[1:], spectrum[1:, itaper],
                            label='#{:d} [loss={:2.2g}]'
                            .format(itaper, 1-evalue))
            else:
                axtemp.plot(degrees[0:], spectrum[0:, itaper],
                            label='#{:d} [loss={:2.2g}]'
                            .format(itaper, 1-evalue))
            axtemp.set_xlabel('Spherical harmonic degree',
                              fontsize=axes_labelsize)
            axtemp.set(xlim=xlim, ylim=ylim)
            axtemp.minorticks_on()
            axtemp.grid(grid, which='major')
            axtemp.tick_params(labelsize=tick_labelsize)

            if title:
                if loss:
                    title_str = '#{:d} [loss={:2.2g}]'.format(itaper, 1-evalue)
                else:
                    title_str = '#{:d} [concentration={:2.2g}]'.format(
                        itaper, evalue)
                axtemp.set_title(title_str, fontsize=titlesize)

        if ax is None:
            fig.tight_layout(pad=0.5)
            if show:
                fig.show()
            if fname is not None:
                fig.savefig(fname)
            return fig, axes

    def plot_coupling_matrix(self, lmax, k=None, weights=None, mode='full',
                             vmin=None, vmax=None, xlabel='Input degree',
                             ylabel='Output degree', title=None,
                             axes_labelsize=None, tick_labelsize=None,
                             titlesize=None, colorbar=None,
                             cb_label=None, normalize=False, show=True,
                             ax=None, fname=None, **kwargs):
        """
        Plot the multitaper coupling matrix.

        This matrix relates the global power spectrum to the expectation of
        the localized multitaper spectrum.

        Usage
        -----
        x.plot_coupling_matrix(lmax, [k, weights, mode, vmin, vmax, xlabel,
                                      ylabel, title, axes_labelsize,
                                      tick_labelsize, titlesize,
                                      colorbar, cb_label, normalize, show, ax,
                                      fname, weights, **kwargs])

        Parameters
        ----------
        lmax : int
            Spherical harmonic bandwidth of the global power spectrum.
        k : int, optional, default = x.nwin
            Number of tapers used in the mutlitaper spectral analysis.
        weights : ndarray, optional, dimension (k), default = x.weights
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
        if weights is not None:
            if k is not None:
                if len(weights) != k:
                    raise ValueError(
                        'Length of weights must be equal to k. '
                        'len(weights) = {:d}, k = {:d}.'
                        .format(len(weights), k))
            else:
                if len(weights) != self.nwin:
                    raise ValueError(
                        'Length of weights must be equal to nwin when k is '
                        'not specified. len(weights) = {:d}, nwin = {:d}.'
                        .format(len(weights), self.nwin))
        else:
            weights = self.weights

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

        mmt = self.coupling_matrix(lmax, k=k, weights=weights, mode=mode)
        if normalize:
            mmt = mmt / mmt.max()

        cim = axes.imshow(mmt, aspect='equal', vmin=vmin, vmax=vmax, **kwargs)
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
            else:
                raise ValueError("colorbar must be either 'horizontal' or "
                                 "'vertical'. Input value is {:s}."
                                 .format(repr(colorbar)))
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
                 weights, taper_degrees, copy=True):
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
        self.taper_degrees = taper_degrees

        if (self.theta_degrees):
            self.area = 2 * _np.pi * (1 - _np.cos(_np.radians(self.theta)))
        else:
            self.area = 2 * _np.pi * (1 - _np.cos(self.theta))

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

        if self.taper_degrees is None:
            self.shannon = (self.lwin + 1)**2 / (4 * _np.pi) * self.area
        else:
            self.shannon = sum(self.eigenvalues)

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
                                 'nwinrot - 1. itaper = {:d}, nwinrot = {:d}.'
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
        Rotate the spherical-cap windows centered on the North pole to clat
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

        Notes
        -----
        This function will take the spherical-cap localization windows
        centered at the North pole (and saved in the attributes tapers and
        orders), rotate each function to the coordinate (clat, clon), and save
        the spherical harmonic coefficients in the attribute coeffs. Each
        column of coeffs contains a single window, and the coefficients are
        ordered according to the convention in SHCilmToVector.
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

        if ((coord_degrees is True and clat == 90. and clon == 0.) or
                (coord_degrees is False and clat == _np.pi/2. and clon == 0.)):
            for i in range(self.nwinrot):
                coeffs = self._taper2coeffs(i)
                self.coeffs[:, i] = _shtools.SHCilmToVector(coeffs)

        else:
            coeffs = _shtools.SHRotateTapers(self.tapers, self.orders,
                                             self.nwinrot, angles, dj_matrix)
            self.coeffs = coeffs

    def _coupling_matrix(self, lmax, k=None, weights=None):
        """Return the coupling matrix of the first nwin tapers."""
        if k is None:
            k = self.nwin

        if weights is None:
            weights = self.weights

        return _shtools.SHMTCouplingMatrix(lmax, self.tapers**2, k=k,
                                           taper_wt=weights)

    def _multitaper_spectrum(self, clm, k, convention='power', unit='per_l',
                             clat=None, clon=None, coord_degrees=True,
                             lmax=None, weights=None):
        """
        Return the multitaper spectrum estimate and standard error for an
        input SHCoeffs class instance.
        """
        if lmax is None:
            lmax = clm.lmax

        if weights is None:
            weights = self.weights

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
                                 'clat = {:s}, clon = {:s}.'
                                 .format(repr(clat), repr(clon)))
            if clat is None and clon is None:
                self.rotate(clat=90., clon=0., coord_degrees=True, nwinrot=k)
            else:
                self.rotate(clat=clat, clon=clon, coord_degrees=coord_degrees,
                            nwinrot=k)

        sh = clm.to_array(normalization='4pi', csphase=1, lmax=lmax)

        mtse, sd = _shtools.SHMultiTaperMaskSE(sh, self.coeffs, lmax=lmax,
                                               k=k, taper_wt=weights)

        if (unit == 'per_l'):
            pass
        elif (unit == 'per_lm'):
            degree_l = _np.arange(len(mtse))
            mtse /= (2.0 * degree_l + 1.0)
            sd /= (2.0 * degree_l + 1.0)
        else:
            raise ValueError(
                "unit must be 'per_l' or 'per_lm'." +
                "Input value is {:s}.".format(repr(unit)))

        if (convention == 'power'):
            return mtse, sd
        elif (convention == 'energy'):
            return mtse * 4.0 * _np.pi, sd * 4.0 * _np.pi
        else:
            raise ValueError(
                "convention must be 'power' or 'energy'." +
                "Input value is {:s}.".format(repr(convention)))

    def _multitaper_cross_spectrum(self, clm, slm, k, convention='power',
                                   unit='per_l', clat=None, clon=None,
                                   coord_degrees=True, lmax=None,
                                   weights=None):
        """
        Return the multitaper cross-spectrum estimate and standard error for
        two input SHCoeffs class instances.
        """
        if lmax is None:
            lmax = min(clm.lmax, slm.lmax)

        if weights is None:
            weights = self.weights

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
                                 'clat = {:s}, clon = {:s}.'
                                 .format(repr(clat), repr(clon)))
            if clat is None and clon is None:
                self.rotate(clat=90., clon=0., coord_degrees=True, nwinrot=k)
            else:
                self.rotate(clat=clat, clon=clon, coord_degrees=coord_degrees,
                            nwinrot=k)

        sh1 = clm.to_array(normalization='4pi', csphase=1, lmax=lmax)
        sh2 = slm.to_array(normalization='4pi', csphase=1, lmax=lmax)

        mtse, sd = _shtools.SHMultiTaperMaskCSE(sh1, sh2, self.coeffs,
                                                lmax1=lmax, lmax2=lmax, k=k,
                                                taper_wt=weights)

        if (unit == 'per_l'):
            pass
        elif (unit == 'per_lm'):
            degree_l = _np.arange(len(mtse))
            mtse /= (2.0 * degree_l + 1.0)
            sd /= (2.0 * degree_l + 1.0)
        else:
            raise ValueError(
                "unit must be 'per_l' or 'per_lm'." +
                "Input value is {:s}.".format(repr(unit)))

        if (convention == 'power'):
            return mtse, sd
        elif (convention == 'energy'):
            return mtse * 4.0 * _np.pi, sd * 4.0 * _np.pi
        else:
            raise ValueError(
                "convention must be 'power' or 'energy'." +
                "Input value is {:s}.".format(repr(convention)))

    def _biased_spectrum(self, spectrum, k, convention='power', unit='per_l',
                         weights=None, save_cg=None, ldata=None):
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
                "Input value is {:s}.".format(repr(convention)))

        if (unit == 'per_l'):
            outspectrum = _shtools.SHBiasK(self.tapers, spectrum, k=k,
                                           taper_wt=weights, save_cg=save_cg,
                                           ldata=ldata)
        elif (unit == 'per_lm'):
            degree_l = _np.arange(len(spectrum))
            temp = spectrum * (2.0 * degree_l + 1.0)
            outspectrum = _shtools.SHBiasK(self.tapers, temp, k=k,
                                           taper_wt=weights, save_cg=save_cg,
                                           ldata=ldata)
            outspectrum /= (2.0 * degree_l + 1.0)
        else:
            raise ValueError(
                "unit must be 'per_l' or 'per_lm'." +
                "Input value is {:s}.".format(repr(unit)))

        return outspectrum

    def _variance(self, power, k, lmax=None, weights=None):
        """
        Compute the theoretical variance of a windowed function.
        """
        var = _np.zeros(lmax+1)
        for l in range(lmax+1):
            var[l] = _shtools.SHMTVar(l, self.tapers, self.orders, power,
                                      kmax=k, taper_wt=weights)

        return var

    def _info(self):
        """Print a summary of the data in the SHWindow instance."""
        print(repr(self))

    def __repr__(self):
        str = 'kind = {:s}\n'.format(repr(self.kind))

        if self.theta_degrees:
            str += 'theta = {:f} degrees\n'.format(self.theta)
        else:
            str += 'theta = {:f} radians'.format(self.theta)

        str += ('lwin = {:d}\n'
                'nwin = {:d}\n'
                'nwinrot = {:s}\n'
                'shannon = {:f}\n'
                'area (radians) = {:e}\n'
                .format(self.lwin, self.nwin, repr(self.nwinrot), self.shannon,
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

        if self.weights is None:
            str += 'Taper weights are not set'
        else:
            str += 'Taper weights are set'

        return str


class SHWindowMask(SHWindow):
    """
    Class for localization windows concentrated within a specified mask and
    for a given spherical harmonic bandwidth.
    """

    @staticmethod
    def istype(kind):
        return kind == 'mask'

    def __init__(self, tapers, eigenvalues, weights, area, taper_degrees,
                 copy=True):
        self.kind = 'mask'
        self.lwin = _np.sqrt(tapers.shape[0]).astype(int) - 1
        self.nwin = tapers.shape[1]
        self.taper_degrees = taper_degrees

        if copy:
            self.weights = weights
            self.tapers = _np.copy(tapers)
            self.eigenvalues = _np.copy(eigenvalues)
        else:
            self.weights = weights
            self.tapers = tapers
            self.eigenvalues = eigenvalues

        self.area = area

        if self.taper_degrees is None:
            self.shannon = (self.lwin + 1)**2 / (4 * _np.pi) * self.area
        else:
            self.shannon = sum(self.eigenvalues)

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

    def _coupling_matrix(self, lmax, k=None, weights=None):
        """Return the coupling matrix of the first nwin tapers."""
        if k is None:
            k = self.nwin

        if weights is None:
            weights = self.weights

        tapers_power = _np.zeros((self.lwin+1, k))
        for i in range(k):
            tapers_power[:, i] = _spectrum(self.to_array(i),
                                           normalization='4pi',
                                           convention='power', unit='per_l')

        return _shtools.SHMTCouplingMatrix(lmax, tapers_power, k=k,
                                           taper_wt=weights)

    def _multitaper_spectrum(self, clm, k, convention='power', unit='per_l',
                             lmax=None, weights=None, clat=None, clon=None,
                             coord_degrees=True):
        """
        Return the multitaper spectrum estimate and standard error for an
        input SHCoeffs class instance.
        """
        if lmax is None:
            lmax = clm.lmax

        if weights is None:
            weights = self.weights

        sh = clm.to_array(normalization='4pi', csphase=1, lmax=lmax)

        mtse, sd = _shtools.SHMultiTaperMaskSE(sh, self.tapers, k=k,
                                               taper_wt=weights)

        if (unit == 'per_l'):
            pass
        elif (unit == 'per_lm'):
            degree_l = _np.arange(len(mtse))
            mtse /= (2.0 * degree_l + 1.0)
            sd /= (2.0 * degree_l + 1.0)
        else:
            raise ValueError(
                "unit must be 'per_l' or 'per_lm'." +
                "Input value is {:s}.".format(repr(unit)))

        if (convention == 'power'):
            return mtse, sd
        elif (convention == 'energy'):
            return mtse * 4.0 * _np.pi, sd * 4.0 * _np.pi
        else:
            raise ValueError(
                "convention must be 'power' or 'energy'." +
                "Input value is {:s}.".format(repr(convention)))

    def _multitaper_cross_spectrum(self, clm, slm, k, convention='power',
                                   unit='per_l', lmax=None, weights=None,
                                   clat=None, clon=None, coord_degrees=True):
        """
        Return the multitaper cross-spectrum estimate and standard error for
        two input SHCoeffs class instances.
        """
        if lmax is None:
            lmax = min(clm.lmax, slm.lmax)

        if weights is None:
            weights = self.weights

        sh1 = clm.to_array(normalization='4pi', csphase=1, lmax=lmax)
        sh2 = slm.to_array(normalization='4pi', csphase=1, lmax=lmax)

        mtse, sd = _shtools.SHMultiTaperMaskCSE(sh1, sh2, self.tapers,
                                                k=k, taper_wt=weights)

        if (unit == 'per_l'):
            pass
        elif (unit == 'per_lm'):
            degree_l = _np.arange(len(mtse))
            mtse /= (2.0 * degree_l + 1.0)
            sd /= (2.0 * degree_l + 1.0)
        else:
            raise ValueError(
                "unit must be 'per_l' or 'per_lm'." +
                "Input value is {:s}.".format(repr(unit)))

        if (convention == 'power'):
            return mtse, sd
        elif (convention == 'energy'):
            return mtse * 4.0 * _np.pi, sd * 4.0 * _np.pi
        else:
            raise ValueError(
                "convention must be 'power' or 'energy'." +
                "Input value is {:s}.".format(repr(convention)))

    def _biased_spectrum(self, spectrum, k, convention='power', unit='per_l',
                         weights=None, save_cg=None, ldata=None):
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
                "Input value is {:s}.".format(repr(convention)))

        if (unit == 'per_l'):
            outspectrum = _shtools.SHBiasKMask(self.tapers, spectrum, k=k,
                                               taper_wt=weights,
                                               save_cg=save_cg, ldata=ldata)
        elif (unit == 'per_lm'):
            degree_l = _np.arange(len(spectrum))
            temp = spectrum * (2.0 * degree_l + 1.0)
            outspectrum = _shtools.SHBiasKMask(self.tapers, temp, k=k,
                                               taper_wt=weights,
                                               save_cg=save_cg, ldata=ldata)
            outspectrum /= (2.0 * degree_l + 1.0)
        else:
            raise ValueError(
                "unit must be 'per_l' or 'per_lm'." +
                "Input value is {:s}.".format(repr(unit)))

        return outspectrum

    def _variance(self, power, nwin, lmax=None, weights=None):
        """
        Compute the theoretical variance of a windowed function.
        """
        raise RuntimeError('Computation of the theoretical variance is '
                           'not yet implemented for arbitrary windows.')

    def _info(self):
        """Print a summary of the data in the SHWindow instance."""
        print(repr(self))

    def __repr__(self):
        str = ('kind = {:s}\n'
               'lwin = {:d}\n'
               'nwin = {:d}\n'
               'shannon = {:f}\n'
               'area (radians) = {:e}\n'.format(repr(self.kind), self.lwin,
                                                self.nwin, self.shannon,
                                                self.area))

        if self.weights is None:
            str += 'Taper weights are not set'
        else:
            str += 'Taper weights are set'

        return str
