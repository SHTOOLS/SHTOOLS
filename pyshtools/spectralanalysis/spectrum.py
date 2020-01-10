import numpy as _np
from scipy.special import factorial as _factorial


def spectrum(clm, normalization='4pi', degrees=None, lmax=None,
             convention='power', unit='per_l', base=10.):
    """
    Return the spectrum of the spherical harmonic coefficients as a function
    of spherical harmonic degree.

    Usage
    -----
    array = spectrum(clm, [normalization, degrees, lmax, convention,
                           unit, base])

    Returns
    -------
    array : ndarray, shape (len(degrees))
        1-D ndarray of the spectrum.

    Parameters
    ----------
    clm : ndarray, shape (2, lmax + 1, lmax + 1)
        ndarray containing the spherical harmonic coefficients.
    normalization : str, optional, default = '4pi'
        '4pi', 'ortho', 'schmidt', or 'unnorm' for geodesy 4pi normalized,
        orthonormalized, Schmidt semi-normalized, or unnormalized coefficients,
        respectively.
    lmax : int, optional, default = len(clm[0,:,0]) - 1.
        Maximum spherical harmonic degree to output.
    degrees : ndarray, optional, default = numpy.arange(lmax+1)
        Array containing the spherical harmonic degrees where the spectrum
        is computed.
    convention : str, optional, default = 'power'
        The type of spectrum to return: 'power' for power spectrum, 'energy'
        for energy spectrum, and 'l2norm' for the l2-norm spectrum.
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
    l2-norm spectrum. Total power is defined as the integral of the
    function squared over all space, divided by the area the function
    spans. If the mean of the function is zero, this is equivalent to the
    variance of the function. The total energy is the integral of the
    function squared over all space and is 4pi times the total power. The
    l2-norm is the sum of the magnitude of the coefficients squared.

    The output spectrum can be expresed using one of three units. 'per_l'
    returns the contribution to the total spectrum from all angular orders
    at degree l. 'per_lm' returns the average contribution to the total
    spectrum from a single coefficient at degree l, and is equal to the
    'per_l' spectrum divided by (2l+1). 'per_dlogl' returns the contribution to
    the total spectrum from all angular orders over an infinitessimal
    logarithmic degree band. The contrubution in the band dlog_a(l) is
    spectrum(l, 'per_dlogl')*dlog_a(l), where a is the base, and where
    spectrum(l, 'per_dlogl) is equal to spectrum(l, 'per_l')*l*log(a).
    """
    if normalization.lower() not in ('4pi', 'ortho', 'schmidt', 'unnorm'):
        raise ValueError("The normalization must be '4pi', 'ortho', " +
                         "'schmidt', or 'unnorm'. Input value was {:s}."
                         .format(repr(normalization)))

    if convention.lower() not in ('power', 'energy', 'l2norm'):
        raise ValueError("convention must be 'power', 'energy', or " +
                         "'l2norm'. Input value was {:s}"
                         .format(repr(convention)))

    if unit.lower() not in ('per_l', 'per_lm', 'per_dlogl'):
        raise ValueError("unit must be 'per_l', 'per_lm', or 'per_dlogl'." +
                         "Input value was {:s}".format(repr(unit)))

    if lmax is None:
        lmax = len(clm[0, :, 0]) - 1

    if degrees is None:
        degrees = _np.arange(lmax+1)

    array = _np.empty(len(degrees))

    if normalization.lower() == 'unnorm':
        if convention.lower() == 'l2norm':
            raise ValueError("convention can not be set to 'l2norm' when " +
                             "using unnormalized harmonics.")

        for i, l in enumerate(degrees):
            ms = _np.arange(l+1)
            conv = _factorial(l+ms) / (2. * l + 1.) / _factorial(l-ms)

            if _np.iscomplexobj(clm):
                array[i] = (conv[0:l + 1] * clm[0, l, 0:l + 1] *
                            clm[0, l, 0:l + 1].conjugate()).real.sum() + \
                           (conv[1:l + 1] * clm[1, l, 1:l + 1] *
                            clm[1, l, 1:l + 1].conjugate()).real.sum()
            else:
                conv[1:l + 1] = conv[1:l + 1] / 2.
                array[i] = (conv[0:l + 1] * clm[0, l, 0:l+1]**2).sum() + \
                           (conv[1:l + 1] * clm[1, l, 1:l+1]**2).sum()

    else:
        for i, l in enumerate(degrees):
            if _np.iscomplexobj(clm):
                array[i] = (clm[0, l, 0:l + 1] *
                            clm[0, l, 0:l + 1].conjugate()).real.sum() + \
                           (clm[1, l, 1:l + 1] *
                            clm[1, l, 1:l + 1].conjugate()).real.sum()
            else:
                array[i] = (clm[0, l, 0:l+1]**2).sum() + \
                           (clm[1, l, 1:l+1]**2).sum()

        if convention.lower() == 'l2norm':
            return array
        else:
            if normalization.lower() == '4pi':
                pass
            elif normalization.lower() == 'schmidt':
                array /= (2. * degrees + 1.)
            elif normalization.lower() == 'ortho':
                array /= (4. * _np.pi)

    if convention.lower() == 'energy':
        array *= 4. * _np.pi

    if unit.lower() == 'per_l':
        pass
    elif unit.lower() == 'per_lm':
        array /= (2. * degrees + 1.)
    elif unit.lower() == 'per_dlogl':
        array *= degrees * _np.log(base)

    return array
