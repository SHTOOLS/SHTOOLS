"""
    Convenience functions for computing the spherical harmonic functions.

    spharm        Compute all the spherical harmonic functions up to a maximum
                  degree and order.
    spharm_lm     Compute the spherical harmonic function for specific degrees
                  l and orders m.
"""
import numpy as _np
import warnings as _warnings

from ..legendre import legendre as _legendre


def spharm(lmax, theta, phi, normalization='4pi', kind='real', csphase=1,
           packed=False, degrees=True):
    """
    Compute all the spherical harmonic functions up to a maximum degree.

    Usage
    -----
    ylm = spharm (lmax, theta, phi, [normalization, kind, csphase, packed,
                                     degrees])

    Returns
    -------
    ylm : float or complex, dimension (2, lmax+1, lmax+1) or
                                      (2, (lmax+1)*(lmax+2)/2)
        An array of spherical harmonic functions, ylm[i, l, m], where l and m
        are the spherical harmonic degree and (positive) order, respectively.
        The index i provides the positive (0) and negative (1) order. If packed
        is True, the array is 2-dimensional with the index of the second column
        corresponding to l*(l+1)/2+m.

    Parameters
    ----------
    lmax : integer
        The maximum degree of the spherical harmonic functions to be computed.
    theta : float
        The colatitude in degrees. Use radians if 'degrees' is set to False.
    phi : float
        The longitude in degrees. Use radians if 'degrees' is set to False.
    normalization : str, optional, default = '4pi'
        '4pi', 'ortho', 'schmidt', or 'unnorm' for geodesy 4pi normalized,
        orthonormalized, Schmidt semi-normalized, or unnormalized spherical
        harmonic functions, respectively.
    kind : str, optional, default = 'real'
        'real' or 'complex' spherical harmonic coefficients.
    csphase : integer, optional, default = 1
        If 1 (default), the Condon-Shortley phase will be excluded. If -1, the
        Condon-Shortley phase of (-1)^m will be appended to the spherical
        harmonic functions.
    packed : bool, optional, default = False
        If True, return a 2-dimensional packed array where the index of the
        second column corresponds to l*(l+1)/2+m, where l and m are
        respectively the degree and order.
    degrees : bool, optional, default = True
        If True, `theta` and `phi` are expressed in degrees.

    Notes
    -----
    spharm will calculate all of the spherical harmonic functions up to degree
    lmax for a given colatitude theta and longitude phi. Three parameters
    determine how the spherical harmonic functions are defined. normalization
    can be either '4pi' (default), 'ortho', 'schmidt', or 'unnorm' for 4pi
    normalized, orthonormalized, Schmidt semi-normalized, or unnormalized
    spherical harmonic functions, respectively. kind can be either 'real' or
    'complex', and csphase determines whether to include or exclude (default)
    the Condon-Shortley phase factor.

    By default, the routine will return a 3-dimensional array, ylm[i, l, m],
    where l and m are the spherical harmonic degree and (positive) order,
    respectively. The index i=0 corresponds to the positive orders, whereas i=1
    corresponds to the negative orders. If the optional parameter packed is set
    to True, the output will instead be a 2-dimensional array where the indices
    of the second column correspond to l*(l+1)/2+m.

    The spherical harmonic functions are calculated using the standard three-
    term recursion formula, and in order to prevent overflows, the scaling
    approach of Holmes and Featherstone (2002) is utilized. The resulting
    functions are accurate to about degree 2800. See Wieczorek and Meschede
    (2018) for exact definitions on how the spherical harmonic functions are
    defined.

    References
    ----------
    Holmes, S. A., and W. E. Featherstone, A unified approach to the Clenshaw
    summation and the recursive computation of very high degree and order
    normalised associated Legendre functions, J. Geodesy, 76, 279-299,
    doi:10.1007/s00190-002-0216-2, 2002.

    Wieczorek, M. A., and M. Meschede. SHTools — Tools for working with
    spherical harmonics, Geochem., Geophys., Geosyst., 19, 2574-2592,
    doi:10.1029/2018GC007529, 2018.
    """
    if lmax < 0:
        raise ValueError(
            "lmax must be greater or equal to 0. Input value was {:s}."
            .format(repr(lmax))
            )

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
            ylm = _np.zeros((2, lmax+1, lmax+1), dtype=_np.float64)
            ylm[0, :, :] = p[:, :]
            ylm[1, :, :] = p[:, :]
            for m in range(lmax+1):
                ylm[0, m:lmax+1, m] *= _np.cos(m*phi)
                ylm[1, m:lmax+1, m] *= _np.sin(m*phi)
        else:
            ylm = _np.zeros((2, lmax+1, lmax+1), dtype=_np.complex128)
            ylm[0, :, :] = p[:, :]
            for m in range(lmax+1):
                ylm[0, m:lmax+1, m] *= (_np.cos(m*phi) + 1j * _np.sin(m*phi))
                ylm[1, m:lmax+1, m] = ylm[0, m:lmax+1, m].conj()
                if _np.mod(m, 2) == 1:
                    ylm[1, m:lmax+1, m] = - ylm[1, m:lmax+1, m]

    else:
        if kind.lower() == 'real':
            ylm = _np.zeros((2, (lmax+1)*(lmax+2)//2), dtype=_np.float64)
            ylm[0, :] = p[:]
            ylm[1, :] = p[:]
            for m in range(lmax+1):
                cos = _np.cos(m*phi)
                sin = _np.sin(m*phi)
                for l in range(m, lmax+1):
                    ind = l*(l+1)//2+m
                    ylm[0, ind] *= cos
                    ylm[1, ind] *= sin
        else:
            ylm = _np.zeros((2, (lmax+1)*(lmax+2)//2), dtype=_np.complex128)
            ylm[0, :] = p[:]
            ylm[1, :] = p[:]
            for m in range(lmax+1):
                eimphi = (_np.cos(m*phi) + 1j * _np.sin(m*phi))
                for l in range(m, lmax+1):
                    ind = l*(l+1)//2+m
                    ylm[0, ind] *= eimphi
                    ylm[1, ind] = ylm[0, ind].conj()
                    if _np.mod(m, 2) == 1:
                        ylm[1, ind] = - ylm[1, ind]

    return ylm


@_np.vectorize
def spharm_lm(l, m, theta, phi, normalization='4pi', kind='real', csphase=1,
              degrees=True):
    """
    Compute the spherical harmonic function for specific degrees and orders.

    Usage
    -----
    ylm = spharm_lm (l, m, theta, phi, [normalization, kind, csphase, degrees])

    Returns
    -------
    ylm : float or complex, ndarray
        The spherical harmonic function ylm, where l and m are the spherical
        harmonic degree and order, respectively.

    Parameters
    ----------
    l : integer, array_like
        The spherical harmonic degree.
    m : integer, array_like
        The spherical harmonic order.
    theta : float, array_like
        The colatitude in degrees. Use radians if 'degrees' is set to False.
    phi : float, array_like
        The longitude in degrees. Use radians if 'degrees' is set to False.
    normalization : str, array_like, optional, default = '4pi'
        '4pi', 'ortho', 'schmidt', or 'unnorm' for geodesy 4pi normalized,
        orthonormalized, Schmidt semi-normalized, or unnormalized spherical
        harmonic functions, respectively.
    kind : str, array_like, optional, default = 'real'
        'real' or 'complex' spherical harmonic coefficients.
    csphase : integer, array_like, optional, default = 1
        If 1 (default), the Condon-Shortley phase will be excluded. If -1, the
        Condon-Shortley phase of (-1)^m will be appended to the spherical
        harmonic functions.
    degrees : bool, array_like, optional, default = True
        If True, `theta` and `phi` are expressed in degrees.

    Notes
    -----
    spharm_lm will calculate the spherical harmonic function for specific
    degrees l, orders m, colatitudes theta and longitudes phi. Three parameters
    determine how the spherical harmonic functions are defined. normalization
    can be either '4pi' (default), 'ortho', 'schmidt', or 'unnorm' for 4pi
    normalized, orthonormalized, Schmidt semi-normalized, or unnormalized
    spherical harmonic functions, respectively. kind can be either 'real' or
    'complex', and csphase determines whether to include or exclude (default)
    the Condon-Shortley phase factor.

    The spherical harmonic functions are calculated using the standard
    three-term recursion formula, and in order to prevent overflows, the
    scaling approach of Holmes and Featherstone (2002) is utilized.
    The resulting functions are accurate to about degree 2800. See Wieczorek
    and Meschede (2018) for exact definitions on how the spherical harmonic
    functions are defined.

    References
    ----------
    Holmes, S. A., and W. E. Featherstone, A unified approach to the Clenshaw
    summation and the recursive computation of very high degree and order
    normalised associated Legendre functions, J. Geodesy, 76, 279-299,
    doi:10.1007/s00190-002-0216-2, 2002.

    Wieczorek, M. A., and M. Meschede. SHTools — Tools for working with
    spherical harmonics, Geochem., Geophys., Geosyst., 19, 2574-2592,
    doi:10.1029/2018GC007529, 2018.
    """
    if l < 0:
        raise ValueError(
            "The degree l must be greater or equal than 0. " +
            "Input value was {:s}.".format(repr(l))
            )

    if m > l:
        raise ValueError(
            "The order m must be less than or equal to the degree l. " +
            "Input values were l={:s} and m={:s}.".format(repr(l), repr(m))
            )

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

    if normalization.lower() == 'unnorm' and l > 85:
        raise ValueError("Calculations using unnormalized coefficients " +
                         "are stable only for degrees less than or equal " +
                         "to 85. Input value was {:d}.".format(l))

    ind = (l*(l+1))//2 + abs(m)

    if degrees is True:
        theta = _np.deg2rad(theta)
        phi = _np.deg2rad(phi)

    if kind.lower() == 'real':
        p = _legendre(l, _np.cos(theta), normalization=normalization,
                      csphase=csphase, cnorm=0, packed=True)
        if m >= 0:
            ylm = p[ind] * _np.cos(m*phi)
        else:
            ylm = p[ind] * _np.sin(abs(m)*phi)

    else:
        p = _legendre(l, _np.cos(theta), normalization=normalization,
                      csphase=csphase, cnorm=1, packed=True)
        ylm = p[ind] * (_np.cos(m*phi) + 1j * _np.sin(abs(m)*phi))  # Yl|m|

        if m < 0:
            ylm = ylm.conj()
            if _np.mod(m, 2) == 1:
                ylm = - ylm

    return ylm
