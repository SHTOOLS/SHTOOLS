def YilmIndexVector(i, l, m):
    """
    Compute the index of a 1D array of spherical harmonic coefficients
    corresponding to i, l, and m.

    Usage
    -----
    index = YilmIndexVector (i, l, m)

    Returns
    -------
    index : integer
        Index of a 1D array of spherical harmonic coefficients corresponding
        to i, l, and m.

    Parameters
    ----------
    i : integer
        1 corresponds to the cosine coefficient Ylm = cilm[0,:,:], and 2
        corresponds to the sine coefficient Yl,-m = cilm[1,:,:].
    l : integer
        The spherical harmonic degree.
    m : integer
        The angular order, which must be greater or equal to zero.

    Notes
    -----
    YilmIndexVector will calculate the index of a 1D vector of spherical
    harmonic coefficients corresponding to degree l, (positive) angular order
    m and i (1 = cosine, 2 = sine). The index is given by l**2+(i-1)*l+m.
    """
    if l < 0:
        raise ValueError('The spherical harmonic degree must be positive. '
                         'Input value is {:s}'.format(repr(l)))
    if m < 0:
        raise ValueError('The angular order must be positive. '
                         'Input value is {:s}'.format(repr(m)))
    if m >= l:
        raise ValueError('The angular order must be less than or equal to '
                         'the spherical harmonic degree. Input degree is {:s}.'
                         ' Input order is {:s}.'.format(repr(l), repr(m)))

    return l**2 + (i - 1) * l + m
