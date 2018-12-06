"""
    Class for Slepian expansion coefficients.
"""

from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

import numpy as _np
import matplotlib as _mpl
import matplotlib.pyplot as _plt
import copy as _copy

from .. import shtools as _shtools

from .shcoeffsgrid import SHCoeffs
from .shcoeffsgrid import SHGrid


__all__ = ['SlepianCoeffs']


class SlepianCoeffs(object):
    """
    Class for Slepian expansion coefficients.

    The SlepianCoeffs class is initialized by:

    >>>  x = Slepian.expand(flm)

    Each class instance defines the following class attributes:

    salpha         : Array of the Slepian expansion coefficients.
    galpha         : A Slepian class instance that contains the associated
                     Slepian functions.
    nalpha         : The number of Slepian expansion coefficients

    Each class instance provides the following methods:

    expand()              : Expand the function of a grid an return an SHGrid
                            class instance.
    to_shcoeffs()         : Return the spherical harmonic coefficients of the
                            function.
    copy()                : Return a copy of the class instance.
    info()                : Print a summary of the data stored in the
                            SlepianCoeffs instance.
"""
    def __init__(self, salpha, galpha, copy=True):
        """
        Initialize the SlepianCoeffs class.
        """
        if copy:
            self.salpha = _np.copy(salpha)
            self.galpha = _copy.deepcopy(galpha)
        else:
            self.salpha = salpha
            self.galpha = galpha

        if self.galpha.kind is 'cap':
            self.nalpha = self.galpha.nwinrot
        else:
            self.nalpha = self.galpha.nwin

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
        str = ('Number of coefficients = {:d}\n'
               'L = {:d}\n'
               .format(self.nalpha, self.galpha.lwin))
        return str

    def expand(self, n=None, grid='DH2', zeros=None):
        """
        Expand the function on a grid using the first n Slepian coefficients.

        Usage
        -----
        f = x.expand([n, grid, zeros])

        Returns
        -------
        f : SHGrid class instance

        Parameters
        ----------
        n : int, optional, default = x.nalpha
            The number of expansion coefficients to use when calculating the
            spherical harmonic coefficients.
        grid : str, optional, default = 'DH2'
            'DH' or 'DH1' for an equisampled lat/lon grid with nlat=nlon, 'DH2'
            for an equidistant lat/lon grid with nlon=2*nlat, or 'GLQ' for a
            Gauss-Legendre quadrature grid.
        zeros : ndarray, optional, default = None
            The cos(colatitude) nodes used in the Gauss-Legendre Quadrature
            grids.
        """
        if type(grid) != str:
            raise ValueError('grid must be a string. ' +
                             'Input type was {:s}'
                             .format(str(type(grid))))

        if n is None:
            n = self.nalpha

        if self.galpha.kind == 'cap':
            shcoeffs = _shtools.SlepianCoeffsToSH(self.salpha,
                                                  self.galpha.coeffs, n)
        else:
            shcoeffs = _shtools.SlepianCoeffsToSH(self.salpha,
                                                  self.galpha.tapers, n)

        if grid.upper() in ('DH', 'DH1'):
            gridout = _shtools.MakeGridDH(shcoeffs, sampling=1,
                                          norm=1, csphase=1)
            return SHGrid.from_array(gridout, grid='DH', copy=False)
        elif grid.upper() == 'DH2':
            gridout = _shtools.MakeGridDH(shcoeffs, sampling=2,
                                          norm=1, csphase=1)
            return SHGrid.from_array(gridout, grid='DH', copy=False)
        elif grid.upper() == 'GLQ':
            if zeros is None:
                zeros, weights = _shtools.SHGLQ(self.galpha.lwin)
            gridout = _shtools.MakeGridGLQ(shcoeffs, zeros,
                                           norm=1, csphase=1)
            return SHGrid.from_array(gridout, grid='GLQ', copy=False)
        else:
            raise ValueError(
                "grid must be 'DH', 'DH1', 'DH2', or 'GLQ'. " +
                "Input value was {:s}".format(repr(grid)))

    def to_shcoeffs(self, n=None, normalization='4pi', csphase=1):
        """
        Return the spherical harmonic coefficients using the first n Slepian
        coefficients.

        Usage
        -----

        s = x.to_shcoeffs([n])

        Returns
        -------
        s : SHCoeffs class instance
            The spherical harmonic coefficients obtained from using the first
            n Slepian expansion coefficients.

        Parameters
        ----------
        n : int, optional, default = x.nalpha
            The number of expansion coefficients to use when calculating the
            spherical harmonic coefficients.
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

        if n is None:
            n = self.nalpha

        if self.galpha.kind == 'cap':
            shcoeffs = _shtools.SlepianCoeffsToSH(self.salpha,
                                                  self.galpha.coeffs, n)
        else:
            shcoeffs = _shtools.SlepianCoeffsToSH(self.salpha,
                                                  self.galpha.tapers, n)

        temp = SHCoeffs.from_array(shcoeffs, normalization='4pi', csphase=1)

        if normalization != '4pi' or csphase != 1:
            return temp.convert(normalization=normalization, csphase=csphase)
        else:
            return temp
