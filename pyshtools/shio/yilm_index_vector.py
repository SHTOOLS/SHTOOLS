def YilmIndexVector(i, l, m):
    """
    Compute the index of an 1D array of spherical harmonic coefficients
    corresponding to i, l, and m.

    Usage
    -----
    index = YilmIndexVector (i, l, m)

    Returns
    -------
    index : integer
        Index of an 1D array of spherical harmonic coefficients corresponding
        to i, l, and m.

    Parameters
    ----------
    i : integer
        1 corresponds to the cosine coefficient cilm[0,:,:], and 2 corresponds
        to the sine coefficient cilm[1,:,:].
    l : integer
        The spherical harmonic degree.
    m : integer
        The angular order.

    Description
    -----------
    YilmIndexVector will calculate the index of a 1D vector of spherical
    harmonic coefficients corresponding to degree l, angular order m and i
    (1 = cosine, 2 = sine). The index is given by l**2+(i-1)*l+m.
    """

    return l**2 + (i - 1) * l + m
