"""
    Convenience functions for computing the associated Legendre functions.

    legendre      Compute all the associated Legendre functions up to a maximum
                  degree and order.
    legendre_lm   Compute the associated Legendre function for specific
                  degrees l and orders m.
"""
import numpy as _np
import warnings as _warnings

from . import PlmBar as _PlmBar
from . import PlmON as _PlmON
from . import PlmSchmidt as _PlmSchmidt
from . import PLegendreA as _PLegendreA


def legendre(lmax, z, normalization='4pi', csphase=1, cnorm=0, packed=False):
    """
    Compute all the associated Legendre functions up to a maximum degree and
    order.

    Usage
    -----
    plm = legendre (lmax, z, [normalization, csphase, cnorm, packed])

    Returns
    -------
    plm : float, dimension (lmax+1, lmax+1) or ((lmax+1)*(lmax+2)/2)
        An array of associated Legendre functions, plm[l, m], where l and m
        are the degree and order, respectively. If packed is True, the array
        is 1-dimensional with the index corresponding to l*(l+1)/2+m.

    Parameters
    ----------
    lmax : integer
        The maximum degree of the associated Legendre functions to be computed.
    z : float
        The argument of the associated Legendre functions.
    normalization : str, optional, default = '4pi'
        '4pi', 'ortho', 'schmidt', or 'unnorm' for use with geodesy 4pi
        normalized, orthonormalized, Schmidt semi-normalized, or unnormalized
        spherical harmonic functions, respectively.
    csphase : integer, optional, default = 1
        If 1 (default), the Condon-Shortley phase will be excluded. If -1, the
        Condon-Shortley phase of (-1)^m will be appended to the associated
        Legendre functions.
    cnorm : integer, optional, default = 0
        If 1, the complex normalization of the associated Legendre functions
        will be used. The default is to use the real normalization.
    packed : bool, optional, default = False
        If True, return a 1-dimensional packed array with the index
        corresponding to l*(l+1)/2+m, where l and m are respectively the
        degree and order.

    Notes
    -----
    legendre will calculate all of the associated Legendre functions up to
    degree lmax for a given argument. The Legendre functions are used typically
    as a part of the spherical harmonic functions, and three parameters
    determine how they are defined. normalization can be either '4pi'
    (default), 'ortho', 'schmidt', or 'unnorm' for use with 4pi normalized,
    orthonormalized, Schmidt semi-normalized, or unnormalized spherical
    harmonic functions, respectively. csphase determines whether to include
    or exclude (default) the Condon-Shortley phase factor. cnorm determines
    whether to normalize the Legendre functions for use with real (default)
    or complex spherical harmonic functions.

    By default, the routine will return a 2-dimensional array, p[l, m]. If the
    optional parameter packed is set to True, the output will instead be a
    1-dimensional array where the indices correspond to l*(l+1)/2+m. The
    Legendre functions are calculated using the standard three-term recursion
    formula, and in order to prevent overflows, the scaling approach of Holmes
    and Featherstone (2002) is utilized. The resulting functions are accurate
    to about degree 2800. See Wieczorek and Meschede (2018) for exact
    definitions on how the Legendre functions are defined.

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
        p = _PlmBar(lmax, z, csphase=csphase, cnorm=cnorm)
    elif normalization == 'ortho':
        p = _PlmON(lmax, z, csphase=csphase, cnorm=cnorm)
    elif normalization == 'schmidt':
        p = _PlmSchmidt(lmax, z, csphase=csphase, cnorm=cnorm)
    elif normalization == 'unnorm':
        p = _PLegendreA(lmax, z, csphase=csphase)

    if packed is True:
        return p
    else:
        plm = _np.zeros((lmax+1, lmax+1))
        for l in range(lmax+1):
            for m in range(l+1):
                plm[l, m] = p[(l*(l+1))//2+m]

        return plm


@_np.vectorize
def legendre_lm(l, m, z, normalization='4pi', csphase=1, cnorm=0):
    """
    Compute the associated Legendre function for specific degrees and orders.

    Usage
    -----
    plm = legendre_lm (l, m, z, [normalization, csphase, cnorm])

    Returns
    -------
    plm : float, ndarray
        The associated Legendre functions for degree l and order m.

    Parameters
    ----------
    l : integer, array_like
        The spherical harmonic degree.
    m : integer, array_like
        The spherical harmonic order.
    z : float, array_like
        The argument of the associated Legendre functions.
    normalization : str, array_like, optional, default = '4pi'
        '4pi', 'ortho', 'schmidt', or 'unnorm' for use with geodesy 4pi
        normalized, orthonormalized, Schmidt semi-normalized, or unnormalized
        spherical harmonic functions, respectively.
    csphase : integer, array_like, optional, default = 1
        If 1 (default), the Condon-Shortley phase will be excluded. If -1, the
        Condon-Shortley phase of (-1)^m will be appended to the associated
        Legendre functions.
    cnorm : integer, array_like, optional, default = 0
        If 1, the complex normalization of the associated Legendre functions
        will be used. The default is to use the real normalization.

    Notes
    -----
    legendre_lm will calculate the associated Legendre function for specific
    degrees l and orders m. The Legendre functions are used typically as a part
    of the spherical harmonic functions, and three parameters determine how
    they are defined. normalization can be either '4pi' (default), 'ortho',
    'schmidt', or 'unnorm' for use with 4pi normalized, orthonormalized,
    Schmidt semi-normalized, or unnormalized spherical harmonic functions,
    respectively. csphase determines whether to include or exclude (default)
    the Condon-Shortley phase factor. cnorm determines whether to normalize
    the Legendre functions for use with real (default) or complex spherical
    harmonic functions.

    The Legendre functions are calculated using the standard three-term
    recursion formula, and in order to prevent overflows, the scaling approach
    of Holmes and Featherstone (2002) is utilized. The resulting functions are
    accurate to about degree 2800. See Wieczorek and Meschede (2018) for exact
    definitions on how the Legendre functions are defined.

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
            "The degree l must be greater or equal to 0. Input value was {:s}."
            .format(repr(l))
            )

    if m < 0:
        raise ValueError(
            "The order m must be greater or equal to 0. Input value was {:s}."
            .format(repr(m))
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

    if normalization.lower() == 'unnorm' and l > 85:
        raise ValueError("Calculations using unnormalized coefficients " +
                         "are stable only for degrees less than or equal " +
                         "to 85. Input value was {:d}.".format(l))

    if normalization == '4pi':
        p = _PlmBar(l, z, csphase=csphase, cnorm=cnorm)
    elif normalization == 'ortho':
        p = _PlmON(l, z, csphase=csphase, cnorm=cnorm)
    elif normalization == 'schmidt':
        p = _PlmSchmidt(l, z, csphase=csphase, cnorm=cnorm)
    elif normalization == 'unnorm':
        p = _PLegendreA(l, z, csphase=csphase)

    return p[(l*(l+1))//2+m]
