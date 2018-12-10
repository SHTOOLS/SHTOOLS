"""
    Convenience functions for computing the spherical harmonic functions.

    spharm        Compute all the spherical harmonic functions up to a maximum
                  degree and order.
    spharm_lm     Compute the spherical harmonic function for a specific degree
                  l and order m.
"""
import numpy as _np

from ..legendre import legendre as _legendre
from ..legendre import legendre_lm as _legendre_lm


def spharm(lmax, theta, phi, normalization='4pi', kind='real', csphase=1,
           packed=False, degrees=True):
    """

    """

    if normalization.lower() not in ('4pi', 'ortho', 'schmidt', 'unnorm'):
        raise ValueError(
            "The normalization must be '4pi', 'ortho', 'schmidt', " +
            "or 'unnorm'. Input value was {:s}."
            .format(repr(normalization))
            )

    if kind.lower() not in ('real', 'complex'):
        raise ValueError(
            "kind must be 'real' or 'complex'. " +
            "Input value was {:s}.".format(repr(kind))
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

    if degrees is True:
        theta = _np.deg2rad(theta)
        phi = _np.deg2rad(phi)

    if kind.lower() == 'real':
        p = _legendre(lmax, _np.cos(theta), normalization=normalization,
                      csphase=csphase, cnorm=0, packed=packed)
    else:
        p = _legendre(lmax, _np.cos(theta), normalization=normalization,
                      csphase=csphase, cnorm=1, packed=packed)

    if packed is False:
        if kind.lower() == 'real':
            ylm = _np.zeros((2, lmax+1, lmax+1), dtype=float_)
            ylm[0, :, :] = p
            ylm[1, :, :] = p
            for m in range(lmax+1):
                ylm[0, m:lmax+1, m] *= _np.cos(m*phi)
                ylm[1, m:lmax+1, m] *= _np.sin(m*phi)
        else:
            ylm = _np.zeros((2, lmax+1, lmax+1), dtype=np.complex_)
            ylm[0, :, :] = p
            ylm[1, :, :] = p
            for m in range(lmax+1):
                ylm[0, m:lmax+1, m] *= (_np.cos(m*phi) + 1j*j_np.sin(m*phi))
                ylm[1, m:lmax+1, m] = ylm[0, m:lmax+1, m].conj() * (-1)**m

    else:
        if kind.lower() == 'real':

        else:

    return ylm


def spharm_lm(l, m, theta, phi, normalization='4pi', csphase=1, cnorm=0,
              degrees=True):
    """

    """

    if m > l:
        raise ValueError(
            "The degree l must be greater or equal than the order m. " +
            "Input degree and order are {:s} and {:s}."
            .format(repr(l), repr(m))
            )

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

    if cnorm != 0 and cnorm != 1:
        raise ValueError(
            "cnorm must be either 0 or 1. Input value was {:s}."
            .format(repr(cnorm))
            )

    if normalization.lower() == 'unnorm' and lmax > 85:
        _warnings.warn("Calculations using unnormalized coefficients " +
                       "are stable only for degrees less than or equal " +
                       "to 85. lmax for the coefficients will be set to " +
                       "85. Input value was {:d}.".format(lmax),
                       category=RuntimeWarning)
        lmax = 85

    if normalization == '4pi':
        p = _PlmBar(l, z, csphase=csphase, cnorm=cnorm)
    elif normalization == 'ortho':
        p = _PlmON(l, z, csphase=csphase, cnorm=cnorm)
    elif normalization == 'schmidt':
        p = _PlmSchmidt(l, z, csphase=csphase, cnorm=cnorm)
    elif normalization == 'unnorm':
        p = _PLegendreA(l, z, csphase=csphase, cnorm=cnorm)

    return p[(l*(l+1))//2+m]
