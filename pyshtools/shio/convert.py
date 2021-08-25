"""
    Functions for converting spherical harmonic coefficients to a different
    normalization conventions.
"""
import numpy as _np
import warnings as _warnings
from scipy.special import factorial as _factorial


def convert(coeffs_in, normalization_in=None, normalization_out=None,
            csphase_in=None, csphase_out=None, lmax=None):
    """
    Convert an array of spherical harmonic coefficients to a different
    normalization convention.

    Usage
    -----
    coeffs_out = convert(coeffs_in, [normalization_in, normalization_out,
                                         csphase_in, csphase_out, lmax])

    Returns
    -------
    coeffs_out : ndarray, size (2, lmax+1, lmax+1)
        An array of spherical harmonic coefficients with the new
        normalization convention.

    Parameters
    ----------
    coeffs_in : ndarray
        The array of imput spherical harmonic coefficients.
    normalization_in : str, optional, default = None
        Normalization of the output coefficients: '4pi', 'ortho'
        'schmidt', or 'unnorm', for geodesy 4pi normalized,
        orthonormalized, Schmidt semi-normalized, or unnormalized
        coefficients, respectively.
    normalization_out : str, optional, default = None
        Normalization of the output coefficients: '4pi', 'ortho'
        'schmidt', or 'unnorm', for geodesy 4pi normalized,
        orthonormalized, Schmidt semi-normalized, or unnormalized
        coefficients, respectively.
    csphase_in : int, optional, default = None
        Condon-Shortley phase convention of the input coefficients: 1 to
        exclude the phase factor, or -1 to include it.
    csphase_out : int, optional, default = None
        Condon-Shortley phase convention of the output coefficients: 1 to
        exclude the phase factor, or -1 to include it.
    lmax : int, optional, default = coeffs.shape[1] - 1
        Maximum spherical harmonic degree to output. If lmax is larger than
        that of the input coefficients, the output array will be zero
        padded.

    Notes
    -----
    This routine will convert an array of spherical harmonic coefficients
    to a different normalization convention and different Condon-Shortley
    phase convention. Optionally, a different maximum spherical harmonic
    degree can be specified. If this degree is smaller than that of the
    input coefficients, the input coefficients will be truncated. If this
    degree is larger than the input coefficients, then the output
    coefficients will be zero padded.
    """

    # check argument consistency
    if normalization_in is not None:
        if type(normalization_in) != str:
            raise ValueError('normalization_in must be a string. ' +
                             'Input type was {:s}'
                             .format(str(type(normalization_in))))
        if normalization_in.lower() not in ('4pi', 'ortho', 'schmidt',
                                            'unnorm'):
            raise ValueError(
                "normalization_in must be '4pi', 'ortho', 'schmidt', or " +
                "'unnorm'. Provided value was {:s}"
                .format(repr(normalization_in))
                )
        if normalization_out is None:
            raise ValueError("normalization_in and normalization_out " +
                             "must both be specified.")
    if normalization_out is not None:
        if type(normalization_out) != str:
            raise ValueError('normalization_out must be a string. ' +
                             'Input type was {:s}'
                             .format(str(type(normalization_out))))
        if normalization_out.lower() not in ('4pi', 'ortho', 'schmidt',
                                             'unnorm'):
            raise ValueError(
                "normalization_out must be '4pi', 'ortho', 'schmidt', or" +
                " 'unnorm'. Provided value was {:s}"
                .format(repr(normalization_out))
                )
        if normalization_in is None:
            raise ValueError("normalization_in and normalization_out " +
                             "must both be specified.")
    if csphase_in is not None:
        if csphase_in != 1 and csphase_in != -1:
            raise ValueError(
                "csphase_in must be 1 or -1. Input value was {:s}"
                .format(repr(csphase_in)))
        if csphase_out is None:
            raise ValueError("csphase_in and csphase_out must both be " +
                             "specified.")
    if csphase_out is not None:
        if csphase_out != 1 and csphase_out != -1:
            raise ValueError(
                "csphase_out must be 1 or -1. Input value was {:s}"
                .format(repr(csphase_out)))
        if csphase_in is None:
            raise ValueError("csphase_in and csphase_out must both be " +
                             "specified.")

    lmaxin = coeffs_in.shape[1] - 1

    if lmax is None:
        lmaxout = lmaxin
    else:
        lmaxout = lmax

    lconv = min(lmaxin, lmaxout)

    if ((normalization_in == 'unnorm' or normalization_out ==
            'unnorm') and lconv > 85):
        _warnings.warn("Conversions with unnormalized coefficients are " +
                       "stable only for degrees less than or equal to " +
                       "85. lmax of the output coefficients will be " +
                       "truncated after degree 85. The spherical " +
                       "harmonic degree of the input coefficients was " +
                       "{:d}.".format(lmaxin), category=RuntimeWarning)
        lconv = 85

    degrees = _np.arange(lconv + 1)

    if _np.iscomplexobj(coeffs_in):
        coeffs = _np.zeros((2, lmaxout+1, lmaxout+1), dtype=_np.complex128)
    else:
        coeffs = _np.zeros((2, lmaxout+1, lmaxout+1))

    coeffs[:, :lconv+1, :lconv+1] = coeffs_in[:, :lconv+1, :lconv+1]

    if normalization_in == normalization_out:
        pass
    elif normalization_in == '4pi' and normalization_out == 'schmidt':
        for l in degrees:
            coeffs[:, l, :l+1] *= _np.sqrt(2. * l + 1.)
    elif normalization_in == '4pi' and normalization_out == 'ortho':
        coeffs *= _np.sqrt(4. * _np.pi)
    elif normalization_in == '4pi' and normalization_out == 'unnorm':
        for l in degrees:
            ms = _np.arange(l+1)
            conv = (2. * l + 1.) * _factorial(l-ms) / _factorial(l+ms)
            if not _np.iscomplexobj(coeffs):
                conv[1:] *= 2.
            coeffs[:, l, :l+1] *= _np.sqrt(conv)
    elif normalization_in == 'schmidt' and normalization_out == '4pi':
        for l in degrees:
            coeffs[:, l, :l+1] /= _np.sqrt(2. * l + 1.)
    elif normalization_in == 'schmidt' and normalization_out == 'ortho':
        for l in degrees:
            coeffs[:, l, :l+1] *= _np.sqrt(4. * _np.pi / (2. * l + 1.))
    elif normalization_in == 'schmidt' and normalization_out == 'unnorm':
        for l in degrees:
            ms = _np.arange(l+1)
            conv = _factorial(l-ms) / _factorial(l+ms)
            if not _np.iscomplexobj(coeffs):
                conv[1:] *= 2.
            coeffs[:, l, :l+1] *= _np.sqrt(conv)
    elif normalization_in == 'ortho' and normalization_out == '4pi':
        coeffs /= _np.sqrt(4. * _np.pi)
    elif normalization_in == 'ortho' and normalization_out == 'schmidt':
        for l in degrees:
            coeffs[:, l, :l+1] *= _np.sqrt((2. * l + 1.) / (4. * _np.pi))
    elif normalization_in == 'ortho' and normalization_out == 'unnorm':
        for l in degrees:
            ms = _np.arange(l+1)
            conv = (2. * l + 1.) * _factorial(l-ms) \
                / 4. / _np.pi / _factorial(l+ms)
            if not _np.iscomplexobj(coeffs):
                conv[1:] *= 2.
            coeffs[:, l, :l+1] *= _np.sqrt(conv)
    elif normalization_in == 'unnorm' and normalization_out == '4pi':
        for l in degrees:
            ms = _np.arange(l+1)
            conv = _factorial(l+ms) / (2. * l + 1.) / _factorial(l-ms)
            if not _np.iscomplexobj(coeffs):
                conv[1:] /= 2.
            coeffs[:, l, :l+1] *= _np.sqrt(conv)
    elif normalization_in == 'unnorm' and normalization_out == 'schmidt':
        for l in degrees:
            ms = _np.arange(l+1)
            conv = _factorial(l+ms) / _factorial(l-ms)
            if not _np.iscomplexobj(coeffs):
                conv[1:] /= 2.
            coeffs[:, l, :l+1] *= _np.sqrt(conv)
    elif normalization_in == 'unnorm' and normalization_out == 'ortho':
        for l in degrees:
            ms = _np.arange(l+1)
            conv = 4. * _np.pi * _factorial(l+ms) / (2. * l + 1.) / \
                _factorial(l-ms)
            if not _np.iscomplexobj(coeffs):
                conv[1:] /= 2.
            coeffs[:, l, :l+1] *= _np.sqrt(conv)

    if csphase_in != csphase_out:
        for m in degrees:
            if m % 2 == 1:
                coeffs[:, m:lconv+1, m] = - coeffs[:, m:lconv+1, m]

    return coeffs
