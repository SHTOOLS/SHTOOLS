# ==== SPHERICAL HARMONICS WINDOW FUNCTION CLASS ====

from __future__ import absolute_import, division, print_function

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from .. import _SHTOOLS as shtools


class SHWindow(object):
    """
    EXPERIMENTAL:
    This class contains collections of spherical harmonics windows that
    provide spectral estimates about a specific region
    """
    def __init__(self):
        print("use one of the following constructors: [...]")

    @classmethod
    def from_cap(self, lmax, nwins, theta, clat=0., clon=0., degrees=True):
        """
        constructs a spherical cap window
        """
        if degrees:
            theta = np.radians(theta)

        tapers, eigenvalues, taper_order = shtools.SHReturnTapers(theta, lmax)
        return SHSymmetricWindow(tapers, eigenvalues, taper_order,
                                 clat=clat, clon=clon)

    @classmethod
    def from_mask(self, lmax, nwins, dh_mask, sampling=1):
        """
        constructs optimal window functions in a masked region (needs dh grid)
        """
        tapers, eigenvalues = shtools.SHReturnTapersMap(
            dh_mask, lmax, sampling=sampling, Ntapers=nwins)
        return SHAsymmetricWindow(tapers, eigenvalues)

    def plot(self, nwins, show=True, fname=None):
        """
        plots the best concentrated spherical harmonics taper functions
        """
        # ---- setup figure and axes
        maxcolumns = 5
        ncolumns = min(maxcolumns, nwins)
        nrows = np.ceil(nwins / ncolumns).astype(int)
        figsize = ncolumns * 1.2, nrows * 1.2 + 0.5
        fig, axes = plt.subplots(nrows, ncolumns, figsize=figsize)
        for ax in axes[:-1, :].flatten():
            for xlabel_i in ax.get_xticklabels():
                xlabel_i.set_visible(False)
        for ax in axes[:, 1:].flatten():
            for ylabel_i in ax.get_yticklabels():
                ylabel_i.set_visible(False)

        # loop through tapers and plot them
        for itaper in range(min(self.nwins, nwins)):
            evalue = self.eigenvalues[itaper]
            coeffs = self._coeffs(itaper)
            ax = axes.flatten()[itaper]
            grid = shtools.MakeGridDH(coeffs)
            ax.imshow(grid)
            ax.set_title('concentration: {:2.2f}'.format(evalue))
        fig.tight_layout(pad=0.5)

        if show:
            plt.show()
        if fname is not None:
            fig.savefig(fname)

    def get_spectrum(self, shcoeffs, nwins):
        """Returns the regional spherical harmonics spectrum"""
        for itaper in range(nwins):
            tapercoeffs = self._coeffs(itaper)
            modelcoeffs = shcoeffs.get_coeffs(normalization='4pi', kind='real')
            coeffs = shtools.SHMultiply(tapercoeffs, modelcoeffs)

    def get_couplingmatrix(self, lmax, nwins):
        """returns the coupling matrix of the first nwins tapers"""
        # store sqrt of taper power in 'tapers' array:
        if nwins > self.nwins:
            nwins = self.nwins
        tapers = np.zeros((self.nl, nwins))
        for itaper in range(nwins):
            tapers[:, itaper] = np.sqrt(shtools.SHPowerSpectrum(
                self._coeffs(itaper)))

        # compute coupling matrix of the first nwins tapers:
        coupling_matrix = shtools.SHMTCouplingMatrix(lmax, tapers[:, :nwins])
        return coupling_matrix

    def plot_couplingmatrix(self, lmax, nwins, show=True, fname=None):
        """plots the window's coupling strength"""
        figsize = mpl.rcParams['figure.figsize']
        figsize[0] = figsize[1]
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        coupling_matrix = self.get_couplingmatrix(lmax, nwins)
        ax.imshow(coupling_matrix)
        ax.set_xlabel('output power')
        ax.set_ylabel('input power')
        fig.tight_layout(pad=0.1)

        if show:
            plt.show()
        if fname is not None:
            fig.savefig(fname)

    def info(self):
        """print meta information about the tapers"""
        self._info()


class SHSymmetricWindow(SHWindow):
    """
    This class saves a symmetric spherical window function. It needs to
    save only the m=0 coefficients
    """
    @staticmethod
    def istype(kind):
        return kind == 'Symmetric'

    def __init__(self, tapers, eigenvalues, orders, clat=0., clon=0.):
        # center of cap window
        self.clat, self.clon = clat, clon
        # nl: number of degrees, nwins: number of windows
        self.nl, self.nwins = tapers.shape
        # lmax: maximum degree
        self.lmax = self.nl - 1
        # tapers[nl,nwins]: ith window coefs with m=orders[iwin]
        self.tapers = tapers
        # concentration factor of the ith taper
        self.eigenvalues = eigenvalues
        # order m of the ith taper
        self.orders = orders

    def _coeffs(self, itaper):
        taperm = self.orders[itaper]
        coeffs = np.zeros((2, self.nl, self.nl))
        if taperm < 0:
            coeffs[1, :, abs(taperm)] = self.tapers[:, itaper]
        else:
            coeffs[0, :, abs(taperm)] = self.tapers[:, itaper]
        return coeffs

    def _info(self):
        print('Cap window with {:d} tapers'.format(self.nwins))


class SHAsymmetricWindow(SHWindow):
    """
    This class saves a asymmetric spherical window function and is much like a
    set of real sherical harmonics. It could maybe be merged at some point...
    """
    @staticmethod
    def istype(kind):
        return kind == 'Asymmetric'

    def __init__(self, tapers, eigenvalues):
        ncoeffs, self.nwins = tapers.shape
        self.nl = np.sqrt(ncoeffs).astype(int)
        self.lmax = self.nl-1
        self.tapers = tapers
        self.eigenvalues = eigenvalues

    def _coeffs(self, itaper):
        return shtools.SHVectorToCilm(self.tapers[:, itaper], self.lmax)

    def _info(self):
        print('Asymmetric window with {:d} tapers'.format(self.nwins))
