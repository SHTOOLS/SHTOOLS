"""
    Routines for computing spherical harmonic coefficients by least squares
    inversion.

    shlsq         Generic function for computing the spherical harmonic
                  coefficients using a (weighted) least squares inversion.
"""
import numpy as _np
import warnings as _warnings

from . import SHExpandLSQ as _SHExpandLSQ
from . import SHExpandLSQ_G as _SHExpandLSQ_G
from . import SHExpandWLSQ as _SHExpandWLSQ
from . import SHExpandWLSQ_G as _SHExpandWLSQ_G


def shlsq(data, latitude, longitude, lmax, weights=None, g=None,
          normalization='4pi', kind='real', csphase=1, degrees=True):
    """
    Determine the spherical harmonic coefficients of an irregularly sampled
    function using a (weighted) least squares inversion, optionally with a
    precomputed data kernel matrix.

    Usage
    -----
    coeffs, chi2 = shlsq (data, latitude, longitude, lmax, [weights, g,
                          normalization, csphase, degrees])

    Returns
    -------
    coeffs : float, dimension (2, lmax+1, lmax+1)
        The real spherical harmonic coefficients of the function, where
        coeffs[0, :, :] and coeffs[1, :, :] refer to the cosine and sine terms,
        respectively.
    chi2 : float
        The (weighted) residual sum of squares misfit.

    Parameters
    ----------
    data : float, dimension (nmax)
        The value of the function at the coordinates (lat, lon).
    latitude : float, dimension (nmax)
        The latitude in degrees corresponding to the value in d.
    longitude : float, dimension (nmax)
        The longitude in degrees corresponding to the value in d.
    lmax : integer
        The maximum spherical harmonic degree of the output coefficients.
    weights : float, dimension (nmax)
        The weights used for a weighted least squares inversion.
    g : float, dimension(nmax, (lmax+1)**2)
        The precomputed data kernel matrix G obtained from LSQ_G.
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
    degrees : bool, optional, default = True
        If True, latitude and longitude are in degrees, otherwise they are in
        radians.

    Notes
    -----
    If this routine is used several times with the same latitude and longitude
    coordinates, the data kernel matrix G can be precomputed using LSG_G. When
    the number of data points is greater or equal to the number of spherical
    harmonic coefficients (i.e., nmax>=(lmax+1)**2), the solution of the
    overdetermined system will be determined. If there are more coefficients
    than data points, then the solution of the underdetermined system that
    minimizes the solution norm will be determined. The inversions are
    performed using the LAPACK routine DGELS.

    When weigths are present, they should be set equal to the inverse of the
    data variance. Here, it is assumed explicitly that each measurement is
    statistically independent (i.e., the weighting matrix is diagonal). The
    weighted least squares inversion must be overdetermined (i.e.,
    nmax>(lmax+1)**2), and the inversion is performed using the same LAPACK
    routine after scaling the data vector and inversion matrix.
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

    if kind == 'complex':
        raise NotImplementedError(
            'Complex coefficients are not implemented in shlsq.')

    if degrees is False:
        latitude = _np.rad2deg(latitude)
        longitude = _np.dad2deg(longitude)

    if normalization.lower() == '4pi':
        norm = 1
    elif normalization.lower() == 'schmidt':
        norm = 2
    elif normalization.lower() == 'unnorm':
        norm = 3
    else:
        norm = 4

    if weights is not None:
        if g is not None:
            coeffs, chi2 = _SHExpandWLSQ_G(data, weights, latitude, longitude,
                                           lmax, g, norm=norm, csphase=csphase)
        else:
            coeffs, chi2 = _SHExpandWLSQ(data, weights, latitude, longitude,
                                         lmax, norm=norm, csphase=csphase)
    else:
        if g is not None:
            coeffs, chi2 = _SHExpandLSQ_G(data, latitude, longitude, lmax, g,
                                          norm=norm, csphase=csphase)
        else:
            coeffs, chi2 = _SHExpandLSQ(data, latitude, longitude, lmax,
                                        norm=norm, csphase=csphase)

    return coeffs, chi2
