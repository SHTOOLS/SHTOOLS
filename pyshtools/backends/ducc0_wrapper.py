"""
    Distinctly Useful Code Collection (DUCC)

    DUCC wrapper functions for use in pyshtools.
"""
import numpy as _np

try:
    import ducc0

    major, minor, patch = ducc0.__version__.split(".")
    if int(major) < 1 and int(minor) < 15:
        raise RuntimeError
except:
    ducc0 = None

# setup a few required variables
if ducc0 is not None:
    import os as _os

    try:
        nthreads = int(_os.environ["OMP_NUM_THREADS"])
    except:
        nthreads = 0


def _fixdtype(arr):
    return arr.astype(_np.float64, copy=False)


def set_nthreads(ntnew):
    global nthreads
    nthreads = ntnew


def available():
    return ducc0 is not None


def _nalm(lmax, mmax):
    return ((mmax + 1) * (mmax + 2)) // 2 + (mmax + 1) * (lmax - mmax)


# ducc0's only accepted conventions are
#   normalization = 'ortho'
#   csphase = 1
# so we need to make the required adjustments


def _get_norm(lmax, norm):
    if norm == 1:
        return _np.full(lmax + 1, _np.sqrt(4 * _np.pi))
    if norm == 2:
        return _np.sqrt(4 * _np.pi / (2 * _np.arange(lmax + 1) + 1.0))
    if norm == 3:
        return _np.sqrt(2 * _np.pi / (2 * _np.arange(lmax + 1) + 1.0))
    if norm == 4:
        return _np.ones(lmax + 1)
    raise RuntimeError("unsupported normalization")


def _rcilm2alm(cilm, lmax):
    alm = _np.empty((_nalm(lmax, lmax),), dtype=_np.complex128)
    alm[0: lmax + 1] = cilm[0, :, 0]
    ofs = lmax + 1
    for m in range(1, lmax + 1):
        alm[ofs: ofs + lmax + 1 - m].real = cilm[0, m:, m]
        alm[ofs: ofs + lmax + 1 - m].imag = cilm[1, m:, m]
        ofs += lmax + 1 - m
    return alm


def _ralm2cilm(alm, lmax):
    cilm = _np.zeros((2, lmax + 1, lmax + 1), dtype=_np.float64)
    cilm[0, :, 0] = alm[0: lmax + 1].real
    ofs = lmax + 1
    for m in range(1, lmax + 1):
        cilm[0, m:, m] = alm[ofs: ofs + lmax + 1 - m].real
        cilm[1, m:, m] = alm[ofs: ofs + lmax + 1 - m].imag
        ofs += lmax + 1 - m
    return cilm


def _apply_norm(alm, lmax, norm, csphase, reverse):
    lnorm = _get_norm(lmax, norm)
    if reverse:
        lnorm = 1.0 / lnorm
    alm[0: lmax + 1] *= lnorm[0: lmax + 1]
    lnorm *= _np.sqrt(2.0) if reverse else (1.0 / _np.sqrt(2.0))
    mlnorm = -lnorm
    ofs = lmax + 1
    for m in range(1, lmax + 1):
        if csphase == 1:
            if m & 1:
                alm[ofs: ofs + lmax + 1 - m].real *= mlnorm[m:]
                alm[ofs: ofs + lmax + 1 - m].imag *= lnorm[m:]
            else:
                alm[ofs: ofs + lmax + 1 - m].real *= lnorm[m:]
                alm[ofs: ofs + lmax + 1 - m].imag *= mlnorm[m:]
        else:
            alm[ofs: ofs + lmax + 1 - m].real *= lnorm[m:]
            alm[ofs: ofs + lmax + 1 - m].imag *= mlnorm[m:]
        ofs += lmax + 1 - m
    if norm == 3:  # special treatment for unnormalized a_lm
        r = _np.arange(lmax + 1)
        fct = _np.ones(lmax + 1)
        ofs = lmax + 1
        if reverse:
            alm[0: lmax + 1] /= _np.sqrt(2)
            for m in range(1, lmax + 1):
                fct[m:] *= _np.sqrt((r[m:] + m) * (r[m:] - m + 1))
                alm[ofs: ofs + lmax + 1 - m] /= fct[m:]
                ofs += lmax + 1 - m
        else:
            alm[0: lmax + 1] *= _np.sqrt(2)
            for m in range(1, lmax + 1):
                fct[m:] *= _np.sqrt((r[m:] + m) * (r[m:] - m + 1))
                alm[ofs: ofs + lmax + 1 - m] *= fct[m:]
                ofs += lmax + 1 - m
    return alm


def _make_alm(cilm, lmax, norm, csphase):
    alm = _rcilm2alm(cilm, lmax)
    return _apply_norm(alm, lmax, norm, csphase, False)


def _extract_alm(alm, lmax, norm, csphase):
    _apply_norm(alm, lmax, norm, csphase, True)
    return _ralm2cilm(alm, lmax)


def _synthesize_DH(alm, lmax, extend, out):
    ducc0.sht.experimental.synthesis_2d(
        alm=alm.reshape((1, -1)),
        map=out[:, : out.shape[1] - extend].reshape(
            (1, out.shape[0], out.shape[1] - extend)
        ),
        spin=0,
        lmax=lmax,
        geometry="CC" if extend else "DH",
        nthreads=nthreads,
    )
    if extend:
        out[:, -1] = out[:, 0]
    return out


def _synthesize_DH_deriv1(alm, lmax, extend, out):
    ducc0.sht.experimental.synthesis_2d_deriv1(
        alm=alm.reshape((1, -1)),
        map=out[:, :, : out.shape[2] - extend],
        lmax=lmax,
        geometry="CC" if extend else "DH",
        nthreads=nthreads,
    )
    out[:, 0, :] = 0.0
    if extend:
        out[:, -1, :] = 0.0
        out[:, :, -1] = out[:, :, 0]
    return out


def _synthesize_GLQ(alm, lmax, extend, out):
    ducc0.sht.experimental.synthesis_2d(
        alm=alm.reshape((1, -1)),
        map=out[:, : out.shape[1] - extend].reshape(
            (1, out.shape[0], out.shape[1] - extend)
        ),
        spin=0,
        lmax=lmax,
        geometry="GL",
        nthreads=nthreads,
    )
    if extend:
        out[:, -1] = out[:, 0]
    return out


def _analyze_DH(map, lmax):
    alm = ducc0.sht.experimental.analysis_2d(
        map=map.reshape((1, map.shape[0], map.shape[1])),
        spin=0,
        lmax=lmax,
        geometry="DH",
        nthreads=nthreads,
    )
    return alm[0]


def _analyze_GLQ(map, lmax):
    alm = ducc0.sht.experimental.analysis_2d(
        map=map.reshape((1, map.shape[0], map.shape[1])),
        spin=0,
        lmax=lmax,
        geometry="GL",
        nthreads=nthreads,
    )
    return alm[0]


def _ccilm2almr(cilm):
    lmax = cilm.shape[1] - 1
    alm = _np.empty((_nalm(lmax, lmax),), dtype=_np.complex128)
    fct = (-1) ** _np.arange(lmax + 1)
    alm[0: lmax + 1] = cilm[0, :, 0].real
    ofs = lmax + 1
    for m in range(1, lmax + 1):
        tmp = _np.conj(cilm[1, m:, m])
        tmp *= fct[m]
        tmp += cilm[0, m:, m]
        tmp *= 1.0 / _np.sqrt(2.0)
        alm[ofs: ofs + lmax + 1 - m] = _np.conj(tmp)
        ofs += lmax + 1 - m
    return alm


def _ccilm2almi(cilm):
    lmax = cilm.shape[1] - 1
    alm = _np.empty((_nalm(lmax, lmax),), dtype=_np.complex128)
    fct = (-1) ** _np.arange(lmax + 1)
    alm[0: lmax + 1] = cilm[0, :, 0].imag
    ofs = lmax + 1
    for m in range(1, lmax + 1):
        tmp = _np.conj(cilm[1, m:, m])
        tmp *= -fct[m]
        tmp += cilm[0, m:, m]
        tmp *= 1.0 / _np.sqrt(2.0)
        alm[ofs: ofs + lmax + 1 - m] = tmp.imag + 1j * tmp.real
        ofs += lmax + 1 - m
    return alm


def _addRealpart(cilm, alm):
    lmax = cilm.shape[1] - 1
    cilm[0, :, 0].real += alm[0: lmax + 1].real
    ofs = lmax + 1
    for m in range(1, lmax + 1):
        tmp = alm[ofs: ofs + lmax + 1 - m] / _np.sqrt(2.0)
        cilm[0, m:, m].real += tmp.real
        cilm[0, m:, m].imag -= tmp.imag
        if m & 1:
            cilm[1, m:, m] -= tmp
        else:
            cilm[1, m:, m] += tmp
        ofs += lmax + 1 - m
    return cilm


def _addImagpart(cilm, alm):
    lmax = cilm.shape[1] - 1
    cilm[0, :, 0].imag += alm[0: lmax + 1].real
    ofs = lmax + 1
    for m in range(1, lmax + 1):
        tmp = alm[ofs: ofs + lmax + 1 - m] / _np.sqrt(2.0)
        cilm[0, m:, m].real += tmp.imag
        cilm[0, m:, m].imag += tmp.real
        if m & 1:
            cilm[1, m:, m].real += tmp.imag
            cilm[1, m:, m].imag -= tmp.real
        else:
            cilm[1, m:, m].real -= tmp.imag
            cilm[1, m:, m].imag += tmp.real
        ofs += lmax + 1 - m
    return cilm


def _prep_lmax(lmax, lmax_calc, cilm):
    if lmax is None:
        lmax = cilm.shape[1] - 1
    if lmax_calc is None:
        lmax_calc = cilm.shape[1] - 1
    if lmax_calc > lmax:
        raise RuntimeError(
            "lmax_calc ({}) must be less than or equal to lmax ({})".format(
                lmax_calc, lmax
            )
        )
    # lmax_calc need must not be higher than cilm.shape[1] - 1.
    lmax_calc = min(cilm.shape[1] - 1, lmax_calc)
    return lmax, lmax_calc, cilm[:, : lmax_calc + 1, : lmax_calc + 1]


def SHRotateRealCoef(cilm, x, dj=None):
    """Determine the spherical harmonic coefficients of a real function rotated by
    three Euler angles.

    Usage
    -----
    cilmrot = SHRotateRealCoef (cilm, x, dj, [lmax])

    Returns
    -------
    cilmrot : float, dimension (2, lmax+1, lmax+1)
        The spherical harmonic coefficients of the rotated function, normalized for
        use with the geodesy 4-pi spherical harmonics.

    Parameters
    ----------
    cilm : float, dimension (2, lmaxin+1, lmaxin+1)
        The input real spherical harmonic coefficients. The coefficients must
        correspond to geodesy 4-pi normalized spherical harmonics that do not
        possess the Condon-Shortley phase convention.
    x : float, dimension(3)
        The three Euler angles, alpha, beta, and gamma, in radians.
    dj : optional, ignored
        This parameter only exists to maintain interface compatibility with the
        "shtools" backend.
    lmax : optional, integer, default = lmaxin
        The maximum spherical harmonic degree of the input and output coefficients.

    Description
    -----------
    SHRotateRealCoef will take the real spherical harmonic coefficients of a
    function, rotate it according to the three Euler anlges in x, and output the
    spherical harmonic coefficients of the rotated function. The input and output
    coefficients must correspond to geodesy 4-pi normalized spherical harmonics that
    do not possess the Condon-Shortley phase convention.

    The rotation of a coordinate system or body can be viewed in two complementary
    ways involving three successive rotations. Both methods have the same initial
    and final configurations, and the angles listed in both schemes are the same.
    This routine uses the 'y convention', where the second rotation axis corresponds
    to the y axis.

    Scheme A:

    (I) Rotation about the z axis by alpha.
    (II) Rotation about the new y axis by beta.
    (III) Rotation about the new z axis by gamma.

    Scheme B:

    (I) Rotation about the z axis by gamma.
    (II) Rotation about the initial y axis by beta.
    (III) Rotation about the initial z axis by alpha.

    The rotations can further be viewed either as a rotation of the coordinate
    system or the physical body. For a rotation of the coordinate system without
    rotation of the physical body, use

    x(alpha, beta, gamma).

    For a rotation of the physical body without rotation of the coordinate system,
    use

    x(-gamma, -beta, -alpha).

    The inverse transform of x(alpha, beta, gamma) is x(-gamma, -beta, -alpha).

    Note that this routine uses the "y convention", where the second rotation is
    with respect to the new y axis. If alpha, beta, and gamma were originally
    defined in terms of the "x convention", where the second rotation was with
    respect to the new x axis, the Euler angles according to the y convention would
    be alpha_y=alpha_x-pi/2, beta_x=beta_y, and gamma_y=gamma_x+pi/2.
    """ # noqa
    lmax = cilm.shape[1] - 1
    alm = _make_alm(cilm, lmax, 1, 1)
    alm = ducc0.sht.rotate_alm(alm, lmax, -x[0], -x[1], -x[2],
                               nthreads=nthreads)
    return _extract_alm(alm, lmax, 1, 1)


def SHRotateComplexCoef(cilm, x, dj=None):
    """Determine the spherical harmonic coefficients of a complex-valued function
    rotated by three Euler angles.

    Usage
    -----
    cilmrot = SHRotateComplexCoef (cilm, x, dj, [lmax])

    Returns
    -------
    cilmrot : complex, dimension (2, lmax+1, lmax+1)
        The spherical harmonic coefficients of the rotated function, normalized for
        use with the geodesy 4-pi spherical harmonics.

    Parameters
    ----------
    cilm : complex, dimension (2, lmaxin+1, lmaxin+1)
        The input complex spherical harmonic coefficients. The coefficients must
        correspond to geodesy 4-pi normalized spherical harmonics that do not
        possess the Condon-Shortley phase convention.
    x : float, dimension(3)
        The three Euler angles, alpha, beta, and gamma, in radians.
    dj : optional, ignored
        This parameter only exists to maintain interface compatibility with the
        "shtools" backend.
    lmax : optional, integer, default = lmaxin
        The maximum spherical harmonic degree of the input and output coefficients.

    Description
    -----------
    SHRotateCoplexCoef will take the complex spherical harmonic coefficients of a
    function, rotate it according to the three Euler anlges in x, and output the
    spherical harmonic coefficients of the rotated function. The input and output
    coefficients must correspond to geodesy 4-pi normalized spherical harmonics that
    do not possess the Condon-Shortley phase convention.

    The rotation of a coordinate system or body can be viewed in two complementary
    ways involving three successive rotations. Both methods have the same initial
    and final configurations, and the angles listed in both schemes are the same.
    This routine uses the 'y convention', where the second rotation axis corresponds
    to the y axis.

    Scheme A:

    (I) Rotation about the z axis by alpha.
    (II) Rotation about the new y axis by beta.
    (III) Rotation about the new z axis by gamma.

    Scheme B:

    (I) Rotation about the z axis by gamma.
    (II) Rotation about the initial y axis by beta.
    (III) Rotation about the initial z axis by alpha.

    The rotations can further be viewed either as a rotation of the coordinate
    system or the physical body. For a rotation of the coordinate system without
    rotation of the physical body, use

    x(alpha, beta, gamma).

    For a rotation of the physical body without rotation of the coordinate system,
    use

    x(-gamma, -beta, -alpha).

    The inverse transform of x(alpha, beta, gamma) is x(-gamma, -beta, -alpha).

    Note that this routine uses the "y convention", where the second rotation is
    with respect to the new y axis. If alpha, beta, and gamma were originally
    defined in terms of the "x convention", where the second rotation was with
    respect to the new x axis, the Euler angles according to the y convention would
    be alpha_y=alpha_x-pi/2, beta_x=beta_y, and gamma_y=gamma_x+pi/2.
    """ # noqa
    lmax = cilm.shape[1] - 1
    alm = _ccilm2almr(cilm)
    alm = _apply_norm(alm, lmax, 1, 1, False)
    alm = ducc0.sht.rotate_alm(alm, lmax, -x[0], -x[1], -x[2],
                               nthreads=nthreads)
    alm = _apply_norm(alm, lmax, 1, 1, True)
    res = _np.zeros((2, lmax + 1, lmax + 1), dtype=_np.complex128)
    _addRealpart(res, alm)
    alm = _ccilm2almi(cilm)
    alm = _apply_norm(alm, lmax, 1, 1, False)
    alm = ducc0.sht.rotate_alm(alm, lmax, -x[0], -x[1], -x[2],
                               nthreads=nthreads)
    alm = _apply_norm(alm, lmax, 1, 1, True)
    _addImagpart(res, alm)
    return res


def MakeGridDH(
    cilm,
    lmax=None,
    norm=1,
    sampling=1,
    csphase=1,
    lmax_calc=None,
    extend=False,
):
    """Create a 2D map from a set of spherical harmonic coefficients using the Driscoll
    and Healy (1994) sampling theorem.

    Usage
    -----
    griddh = MakeGridDH (cilm, [lmax, norm, sampling, csphase, lmax_calc, extend])

    Returns
    -------
    griddh : float, dimension (nlat, nlong)
        A 2D map of the input spherical harmonic coefficients cilm that conforms to
        the sampling theorem of Driscoll and Healy (1994). If sampling is 1, the
        grid is equally sampled and is dimensioned as (n by n), where n is 2lmax+2.
        If sampling is 2, the grid is equally spaced and is dimensioned as (n by
        2n). The first latitudinal band of the grid corresponds to 90 N, the
        latitudinal sampling interval is 180/n degrees, and the default behavior is
        to exclude the latitudinal band for 90 S. The first longitudinal band of the
        grid is 0 E, by default the longitudinal band for 360 E is not included, and
        the longitudinal sampling interval is 360/n for an equally sampled and 180/n
        for an equally spaced grid, respectively. If extend is 1, the longitudinal
        band for 360 E and the latitudinal band for 90 S will be included, which
        increases each of the dimensions of the grid by 1.

    Parameters
    ----------
    cilm : float, dimension (2, lmaxin+1, lmaxin+1)
        The real spherical harmonic coefficients of the function. The coefficients
        cilm[0,l,m] and cilm[1,l,m] refer to the "cosine" (Clm) and "sine" (Slm)
        coefficients, respectively.
    lmax : optional, integer, default = lmaxin
        The maximum spherical harmonic degree of the function, which determines the
        sampling n of the output grid.
    norm : optional, integer, default = 1
        1 = 4-pi (geodesy) normalized harmonics; 2 = Schmidt semi-normalized
        harmonics; 3 = unnormalized harmonics;  4 = orthonormal harmonics.
    sampling : optional, integer, default = 1
        If 1 (default) the input grid is equally sampled (n by n). If 2, the grid is
        equally spaced (n by 2n).
    csphase : optional, integer, default = 1
        1 (default) = do not apply the Condon-Shortley phase factor to the
        associated Legendre functions; -1 = append the Condon-Shortley phase factor
        of (-1)^m to the associated Legendre functions.
    lmax_calc : optional, integer, default = lmax
        The maximum spherical harmonic degree used in evaluating the  function. This
        must be less than or equal to lmax, and does not affect the number of
        samples of the output grid.
    extend : input, optional, bool, default = False
        If True, compute the longitudinal band for 360 E and the latitudinal band
        for 90 S. This increases each of the dimensions of griddh by 1.

    Description
    -----------
    MakeGridDH will create a 2-dimensional map equally sampled or equally spaced in
    latitude and longitude from a set of input spherical harmonic coefficients. This
    grid conforms with the sampling theorem of Driscoll and Healy (1994) and this
    routine is the inverse of SHExpandDH. The function is evaluated at each
    longitudinal band by inverse Fourier transforming the sin and cos terms for each
    degree l, and then summing over all degrees. When evaluating the function, the
    maximum spherical harmonic degree that is considered is the minimum of lmaxin,
    lmax, and lmax_calc (if specified).

    The default is to use an input grid that is equally sampled (n by n), but this
    can be changed to use an equally spaced grid (n by 2n) by the optional argument
    sampling. The redundant longitudinal band for 360 E and the latitudinal band for
    90 S are excluded by default, but these can be computed by specifying the
    optional argument extend. The employed spherical harmonic normalization and
    Condon-Shortley phase convention can be set by the optional arguments norm and
    csphase; if not set, the default is to use geodesy 4-pi normalized harmonics
    that exclude the Condon-Shortley phase of (-1)^m.

    The normalized legendre functions are calculated in
    this routine using the recurrence given by Ishioka (2018), which are accurate
    to at least degree 100000. The unnormalized functions are accurate
    only to about degree 15.

    References
    ----------
    Driscoll, J.R. and D.M. Healy, Computing Fourier transforms and convolutions on
    the 2-sphere, Adv. Appl. Math., 15, 202-250, 1994.

    Ishioka, K.: Journal of the Meteorological Society of Japan, 96, 241−249, 2018
    """ # noqa
    lmax, lmax_calc, cilm = _prep_lmax(lmax, lmax_calc, cilm)
    alm = _make_alm(cilm, lmax_calc, norm, csphase)
    out = _np.empty([2 * lmax + 2 + extend,
                     sampling * (2 * lmax + 2) + extend])
    return _synthesize_DH(alm, lmax_calc, extend, out)


def MakeGridDHC(
    cilm,
    lmax=None,
    norm=1,
    sampling=1,
    csphase=1,
    lmax_calc=None,
    extend=False,
):
    """Create a 2D complex map from a set of complex spherical harmonic coefficients
    that conforms with Driscoll and Healy's (1994) sampling theorem.

    Usage
    -----
    griddh = MakeGridDHC (cilm, [lmax, norm, sampling, csphase, lmax_calc, extend])

    Returns
    -------
    griddh : complex, dimension (nlat, nlong)
        A 2D complex map of the input spherical harmonic coefficients cilm that
        conforms to the sampling theorem of Driscoll and Healy (1994). If sampling
        is 1, the grid is equally sampled and is dimensioned as (n by n), where n is
        2lmax+2. If sampling is 2, the grid is equally spaced and is dimensioned as
        (n by 2n). The first latitudinal band of the grid corresponds to 90 N, the
        latitudinal sampling interval is 180/n degrees, and the default behavior is
        to exclude the latitudinal band for 90 S. The first longitudinal band of the
        grid is 0 E, by default the longitudinal band for 360 E is not included, and
        the longitudinal sampling interval is 360/n for an equally sampled and 180/n
        for an equally spaced grid, respectively. If extend is 1, the longitudinal
        band for 360 E and the latitudinal band for 90 S will be included, which
        increases each of the dimensions of the grid by 1.

    Parameters
    ----------
    cilm : complex, dimension (2, lmaxin+1, lmaxin+1)
        The complex spherical harmonic coefficients of the function.  The first
        index specifies the coefficient corresponding to the positive and negative
        order of m, respectively, with Clm=cilm[0,l,m] and Cl,-m=cilm[1,l,m)].
    lmax : optional, integer, default = lmaxin
        The maximum spherical harmonic degree of the function, which determines the
        sampling n of the output grid.
    norm : optional, integer, default = 1
        1 = 4-pi (geodesy) normalized harmonics; 2 = Schmidt semi-normalized
        harmonics; 3 = unnormalized harmonics;  4 = orthonormal harmonics.
    sampling : optional, integer, default = 1
        If 1 (default) the input grid is equally sampled (n by n). If 2, the grid is
        equally spaced (n by 2n).
    csphase : optional, integer, default = 1
        1 (default) = do not apply the Condon-Shortley phase factor to the
        associated Legendre functions; -1 = append the Condon-Shortley phase factor
        of (-1)^m to the associated Legendre functions.
    lmax_calc : optional, integer, default = lmax
        The maximum spherical harmonic degree used in evaluating the  function. This
        must be less than or equal to lmax, and does not affect the number of
        samples of the output grid.
    extend : input, optional, bool, default = False
        If True, compute the longitudinal band for 360 E and the latitudinal band
        for 90 S. This increases each of the dimensions of griddh by 1.

    Description
    -----------
    MakeGridDHC will create a 2-dimensional complex map equally sampled (n by n) or
    equally spaced (n by 2n) in latitude and longitude from a set of input complex
    spherical harmonic coefficients, where N is 2lmax+2. This grid conforms with the
    sampling theorem of Driscoll and Healy (1994) and this routine is the inverse of
    SHExpandDHC. The function is evaluated at each longitudinal band by inverse
    Fourier transforming the exponential terms for each degree l, and then summing
    over all degrees. When evaluating the function, the maximum spherical harmonic
    degree that is considered is the minimum of lmax, the size of cilm-1, or
    lmax_calc (if specified).

    The default is to use an input grid that is equally sampled (n by n), but this
    can be changed to use an equally spaced grid (n by 2n) by the optional argument
    sampling. The redundant longitudinal band for 360 E and the latitudinal band for
    90 S are excluded by default, but these can be computed by specifying the
    optional argument extend. The employed spherical harmonic normalization and
    Condon-Shortley phase convention can be set by the optional arguments norm and
    csphase; if not set, the default is to use geodesy 4-pi normalized harmonics
    that exclude the Condon-Shortley phase of (-1)^m.

    The normalized legendre functions are calculated in
    this routine using the recurrence given by Ishioka (2018), which are accurate
    to at least degree 100000. The unnormalized functions are accurate
    only to about degree 15.


    References
    ----------
    Driscoll, J.R. and D.M. Healy, Computing Fourier transforms and convolutions on
    the 2-sphere, Adv. Appl. Math., 15, 202-250, 1994.

    Ishioka, K.: Journal of the Meteorological Society of Japan, 96, 241−249, 2018
    """ # noqa
    lmax, lmax_calc, cilm = _prep_lmax(lmax, lmax_calc, cilm)
    alm = _ccilm2almi(cilm)
    alm = _apply_norm(alm, lmax, norm, csphase, False)
    res = _np.empty(
        [2 * lmax + 2 + extend, sampling * (2 * lmax + 2) + extend],
        dtype=_np.complex128,
    )
    _synthesize_DH(alm, lmax_calc, extend, res.imag)
    alm = _ccilm2almr(cilm)
    alm = _apply_norm(alm, lmax, norm, csphase, False)
    _synthesize_DH(alm, lmax_calc, extend, res.real)
    return res


def SHExpandDH(griddh, norm=1, sampling=1, csphase=1, lmax_calc=None):
    """Expand an equally sampled or equally spaced grid into spherical harmonics using
    Driscoll and Healy's (1994) sampling theorem.

    Usage
    -----
    cilm = SHExpandDH (griddh, [norm, sampling, csphase, lmax_calc])

    Returns
    -------
    cilm : float, dimension (2, n/2, n/2) or (2, lmax_calc+1, lmax_calc+1)
        The real spherical harmonic coefficients of the function. These will be
        exact if the function is bandlimited to degree lmax=n/2-1. The coefficients
        c1lm and c2lm refer to the cosine (clm) and sine (slm) coefficients,
        respectively, with clm=cilm[0,l,m] and slm=cilm[1,l,m].

    Parameters
    ----------
    griddh : float, dimension (n, n) or (n, 2*n)
        A 2D equally sampled (default) or equally spaced grid that conforms to the
        sampling theorem of Driscoll and Healy (1994). The first latitudinal band
        corresponds to 90 N, the latitudinal band for 90 S is not included, and the
        latitudinal sampling interval is 180/n degrees. The first longitudinal band
        is 0 E, the longitude band for 360 E is not included, and the longitudinal
        sampling interval is 360/n for an equally and 180/n for an equally spaced
        grid, respectively.
    norm : optional, integer, default = 1
        1 (default) = 4-pi (geodesy) normalized harmonics; 2 = Schmidt semi-
        normalized harmonics; 3 = unnormalized harmonics; 4 = orthonormal harmonics.
    sampling : optional, integer, default = 1
        If 1 (default) the input grid is equally sampled (n by n). If 2, the grid is
        equally spaced (n by 2n).
    csphase : optional, integer, default = 1
        1 (default) = do not apply the Condon-Shortley phase factor to the
        associated Legendre functions; -1 = append the Condon-Shortley phase factor
        of (-1)^m to the associated Legendre functions.
    lmax_calc : optional, integer, default = n/2-1
        The maximum spherical harmonic degree calculated in the spherical harmonic
        expansion.

    Description
    -----------
    SHExpandDH will expand an equally sampled (n by n) or equally spaced grid (n by
    2n) into spherical harmonics using the sampling theorem of Driscoll and Healy
    (1994). The number of latitudinal samples, n, must be even, and the transform is
    exact if the function is bandlimited to spherical harmonic degree n/2-1. The
    inverse transform is given by the routine MakeGridDH. If the optional parameter
    lmax_calc is specified, the spherical harmonic coefficients will only be
    calculated to this degree instead of n/2-1. The algorithm is based on performing
    FFTs in longitude and then integrating over latitude using an exact quadrature
    rule.

    The default is to use an input grid that is equally sampled (n by n), but this
    can be changed to use an equally spaced grid (n by 2n) by the optional argument
    sampling.  When using an equally spaced grid, the Fourier components
    corresponding to degrees greater than n/2-1 are simply discarded; this is done
    to prevent aliasing that would occur if an equally sampled grid was constructed
    from an equally spaced grid by discarding every other column of the input grid.

    The employed spherical harmonic normalization and Condon-Shortley phase
    convention can be set by the optional arguments norm and csphase; if not set,
    the default is to use geodesy 4-pi normalized harmonics that exclude the Condon-
    Shortley phase of (-1)^m. The normalized legendre functions are calculated in
    this routine using the recurrence given by Ishioka (2018), which are accurate
    to at least degree 100000. The unnormalized functions are accurate
    only to about degree 15.

    References
    ----------
    Driscoll, J.R. and D.M. Healy, Computing Fourier transforms and convolutions on
    the 2-sphere, Adv. Appl. Math., 15, 202-250, 1994.

    Ishioka, K.: Journal of the Meteorological Society of Japan, 96, 241−249, 2018
    """ # noqa

    griddh = _fixdtype(griddh)
    if griddh.shape[1] != sampling * griddh.shape[0]:
        raise RuntimeError("grid resolution mismatch")
    if lmax_calc is None:
        lmax_calc = griddh.shape[0] // 2 - 1
    if lmax_calc > (griddh.shape[0] // 2 - 1):
        raise RuntimeError("lmax_calc too high")
    alm = _analyze_DH(griddh, lmax_calc)
    return _extract_alm(alm, lmax_calc, norm, csphase)


def SHExpandDHC(griddh, norm=1, sampling=1, csphase=1, lmax_calc=None):
    """Expand an equally sampled or equally spaced complex grid into complex spherical
    harmonics using Driscoll and Healy's (1994) sampling theorem.

    Usage
    -----
    cilm = SHExpandDHC (griddh, [norm, sampling, csphase, lmax_calc])

    Returns
    -------
    cilm : complex, dimension (2, n/2, n/2) or (2, lmax_calc+1, lmax_calc+1)
        The complex spherical harmonic coefficients of the function. These will be
        exact if the function is bandlimited to degree lmax=n/2-1. The first index
        specifies the coefficient corresponding to the positive and negative order
        of m, respectively, with Clm=cilm[0,l,m] and Cl,-m=cilm[1,l,m].

    Parameters
    ----------
    griddh : complex, dimension (n, n) or (n, 2*n)
        A 2D equally sampled (default) or equally spaced complex grid that conforms
        to the sampling theorem of Driscoll and Healy (1994). The first latitudinal
        band corresponds to 90 N, the latitudinal band for 90 S is not included, and
        the latitudinal sampling interval is 180/n degrees. The first longitudinal
        band is 0 E, the longitude band for 360 E is not included, and the
        longitudinal sampling interval is 360/n for an equally and 180/n for an
        equally spaced grid, respectively.
    norm : optional, integer, default = 1
        1 (default) = 4-pi (geodesy) normalized harmonics; 2 = Schmidt semi-
        normalized harmonics; 3 = unnormalized harmonics; 4 = orthonormal harmonics.
    sampling : optional, integer, default = 1
        If 1 (default) the input grid is equally sampled (n by n). If 2, the grid is
        equally spaced (n by 2n).
    csphase : optional, integer, default = 1
        1 (default) = do not apply the Condon-Shortley phase factor to the
        associated Legendre functions; -1 = append the Condon-Shortley phase factor
        of (-1)^m to the associated Legendre functions.
    lmax_calc : optional, integer, default = n/2-1
        The maximum spherical harmonic degree calculated in the spherical harmonic
        expansion.

    Description
    -----------
    SHExpandDHC will expand an equally sampled (n by n) or equally spaced complex
    grid (n by 2n) into complex spherical harmonics using the sampling theorem of
    Driscoll and Healy (1994). The number of latitudinal samples n must be even, and
    the transform is exact if the function is bandlimited to spherical harmonic
    degree n/2 - 1. The inverse transform is given by the routine MakeGridDHC. If
    the optional parameter lmax_calc is specified, the spherical harmonic
    coefficients will only be calculated to this degree instead of n/2 - 1. The
    algorithm is based on performing FFTs in longitude and then integrating over
    latitude using an exact quadrature rule.

    The default is to use an input grid that is equally sampled (n by n), but this
    can be changed to use an equally spaced grid (n by 2n) by the optional argument
    sampling.  When using an equally spaced grid, the Fourier components
    corresponding to degrees greater than n/2 - 1 are simply discarded; this is done
    to prevent aliasing that would occur if an equally sampled grid was constructed
    from an equally spaced grid by discarding every other column of the input grid.

    The employed spherical harmonic normalization and Condon-Shortley phase
    convention can be set by the optional arguments norm and csphase; if not set,
    the default is to use geodesy 4-pi normalized harmonics that exclude the Condon-
    Shortley phase of (-1)^m. The normalized legendre functions are calculated in
    this routine using the recurrence given by Ishioka (2018), which are accurate
    to at least degree 100000. The unnormalized functions are accurate
    only to about degree 15.

    References
    ----------
    Driscoll, J.R. and D.M. Healy, Computing Fourier transforms and convolutions on
    the 2-sphere, Adv. Appl. Math., 15, 202-250, 1994.

    Ishioka, K.: Journal of the Meteorological Society of Japan, 96, 241−249, 2018
    """ # noqa
    if griddh.shape[1] != sampling * griddh.shape[0]:
        raise RuntimeError("grid resolution mismatch")
    if lmax_calc is None:
        lmax_calc = griddh.shape[0] // 2 - 1
    if lmax_calc > (griddh.shape[0] // 2 - 1):
        raise RuntimeError("lmax_calc too high")
    lmax = griddh.shape[0] // 2 - 1 if lmax_calc is None else lmax_calc
    res = _np.zeros((2, lmax + 1, lmax + 1), dtype=_np.complex128)
    alm = _analyze_DH(_fixdtype(griddh.real), lmax_calc)
    alm = _apply_norm(alm, lmax, norm, csphase, True)
    _addRealpart(res, alm)
    alm = _analyze_DH(_fixdtype(griddh.imag), lmax_calc)
    alm = _apply_norm(alm, lmax, norm, csphase, True)
    _addImagpart(res, alm)
    return res


# zero is ignored (they are computed internally)
def MakeGridGLQ(
    cilm, zero=None, lmax=None, norm=1, csphase=1, lmax_calc=None, extend=False
):
    """Create a 2D map from a set of spherical harmonic coefficients sampled on the
    Gauss-Legendre quadrature nodes.

    Usage
    -----
    gridglq = MakeGridGLQ (cilm, zero, [lmax,  norm, csphase, lmax_calc, extend])

    Returns
    -------
    gridglq : float, dimension (nlat, nlong)
        A 2D map of the function sampled on the Gauss-Legendre quadrature nodes,
        dimensioned as (lmax+1, 2*lmax+1) if extend is 0 or (lmax+1, 2*lmax+2) if
        extend is 1.

    Parameters
    ----------
    cilm : float, dimension (2, lmaxin+1, lmaxin+1)
        The real spherical harmonic coefficients of the function. When evaluating
        the function, the maximum spherical harmonic degree considered is the
        minimum of lmax, lmaxin, or lmax_calc (if specified). The first index
        specifies the coefficient corresponding to the positive and negative order
        of m, respectively, with Clm=cilm[0,l,m+] and Cl,-m=cilm[1,l,m].
    zero : optional, ignored
        This parameter only exists to maintain interface compatibility with the
        "shtools" backend.
    lmax : optional, integer, default = lamxin
        The maximum spherical harmonic bandwidth of the function. This determines
        the sampling nodes and dimensions of the output grid.
    norm : optional, integer, default = 1
        1 (default) = Geodesy 4-pi normalized harmonics; 2 = Schmidt semi-normalized
        harmonics; 3 = unnormalized harmonics; 4 = orthonormal harmonics.
    csphase : optional, integer, default = 1
        1 (default) = do not apply the Condon-Shortley phase factor to the
        associated Legendre functions; -1 = append the Condon-Shortley phase factor
        of (-1)^m to the associated Legendre functions.
    lmax_calc : optional, integer, default = lmax
        The maximum spherical harmonic degree used in evaluating the function. This
        must be less than or equal to lmax.
    extend : input, optional, bool, default = False
        If True, compute the longitudinal band for 360 E.

    Description
    -----------
    MakeGridGLQ will create a 2-dimensional map from a set of input spherical
    harmonic coefficients sampled on the Gauss-Legendre quadrature nodes. This is
    the inverse of the routine SHExpandGLQ. The latitudinal nodes correspond to the
    zeros of the Legendre polynomial of degree lmax+1, and the longitudinal nodes
    are equally spaced with an interval of 360/(2*lmax+1) degrees. When evaluating
    the function, the maximum spherical harmonic degree that is considered is the
    minimum of lmax, the size of cilm-1, or lmax_calc (if specified).

    The redundant longitudinal band for 360 E is excluded from the grid by default,
    but this can be computed by specifying the optional argument extend. The
    employed spherical harmonic normalization and Condon-Shortley phase convention
    can be set by the optional arguments norm and csphase; if not set, the default
    is to use geodesy 4-pi normalized harmonics that exclude the Condon-Shortley
    phase of (-1)^m. The normalized legendre functions are calculated in
    this routine using the recurrence given by Ishioka (2018), which are accurate to
    at least degree 100000. The unnormalized functions are accurate
    only to about degree 15.

    The zeros of the Legendre polynomials and the quadrature weights are computed
    using the method described by Bogaert (2014).

    References
    ----------
    Bogaert, I.: SIAM Journal on Scientific Computing, 36, A1008-A1026, 2014

    Ishioka, K.: Journal of the Meteorological Society of Japan, 96, 241−249, 2018
    """ # noqa
    lmax, lmax_calc, cilm = _prep_lmax(lmax, lmax_calc, cilm)
    alm = _make_alm(cilm, lmax_calc, norm, csphase)
    out = _np.empty([lmax + 1, (2 * lmax + 1) + extend])
    return _synthesize_GLQ(alm, lmax_calc, extend, out)


# zero is ignored (they are computed internally)
def MakeGridGLQC(
    cilm, zero=None, lmax=None, norm=1, csphase=1, lmax_calc=None, extend=False
):
    """Create a 2D complex map from a set of complex spherical harmonic coefficients
    sampled on the Gauss-Legendre quadrature nodes.

    Usage
    -----
    gridglq = MakeGridGLQC (cilm, zero, [lmax, norm, csphase, lmax_calc, extend])

    Returns
    -------
    gridglq : complex, dimension (nlat, nlong)
        A 2D complex map of the function sampled on the Gauss-Legendre quadrature
        nodes, dimensioned as (lmax+1, 2*lmax+1) if extend is 0 or (lmax+1,
        2*lmax+2) if extend is 1.


    Parameters
    ----------
    cilm : complex, dimension (2, lmaxin+1, lmaxin+1)
        The complex spherical harmonic coefficients of the function. When evaluating
        the function, the maximum spherical harmonic degree considered is the
        minimum of lmax, lmaxin, or lmax_calc (if specified). The first index
        specifies the coefficient corresponding to the positive and negative order
        of m, respectively, with Clm=cilm[0,l,m+] and Cl,-m=cilm[1,l,m].
    zero : optional, ignored
        This parameter only exists to maintain interface compatibility with the
        "shtools" backend.
    lmax : optional, integer, default = lmaxin
        The maximum spherical harmonic bandwidth of the function. This determines
        the sampling nodes and dimensions of the output grid.
    norm : optional, integer, default = 1
        1 (default) = Geodesy 4-pi normalized harmonics; 2 = Schmidt semi-normalized
        harmonics; 3 = unnormalized harmonics; 4 = orthonormal harmonics.
    csphase : optional, integer, default = 1
        1 (default) = do not apply the Condon-Shortley phase factor to the
        associated Legendre functions; -1 = append the Condon-Shortley phase factor
        of (-1)^m to the associated Legendre functions.
    lmax_calc : optional, integer, default = lmax
        The maximum spherical harmonic degree used in evaluating the function. This
        must be less than or equal to lmax.
    extend : input, optional, bool, default = False
        If True, compute the longitudinal band for 360 E.

    Description
    -----------
    MakeGridGLQC will create a 2-dimensional complex map from a set of input complex
    spherical harmonic coefficients sampled on the Gauss-Legendre quadrature nodes.
    This is the inverse of the routine SHExpandGLQC. The latitudinal nodes
    correspond to the zeros of the Legendre polynomial of degree lmax+1, and the
    longitudinal nodes are equally spaced with an interval of 360/(2*lmax+1)
    degrees. When evaluating the function, the maximum spherical harmonic degree
    that is considered is the minimum of lmax, lmaxin, or lmax_calc (if specified).

    The redundant longitudinal band for 360 E is excluded from the grid by default,
    but this can be computed by specifying the optional argument extend. The
    employed spherical harmonic normalization and Condon-Shortley phase convention
    can be set by the optional arguments norm and csphase; if not set, the default
    is to use geodesy 4-pi normalized harmonics that exclude the Condon-Shortley
    phase of (-1)^m. The normalized legendre functions are calculated in
    this routine using the recurrence given by Ishioka (2018), which are accurate to
    at least degree 100000. The unnormalized functions are accurate
    only to about degree 15.

    The zeros of the Legendre polynomials and the quadrature weights are computed
    using the method described by Bogaert (2014).

    References
    ----------
    Bogaert, I.: SIAM Journal on Scientific Computing, 36, A1008-A1026, 2014

    Ishioka, K.: Journal of the Meteorological Society of Japan, 96, 241−249, 2018
    """ # noqa
    lmax, lmax_calc, cilm = _prep_lmax(lmax, lmax_calc, cilm)
    alm = _ccilm2almi(cilm)
    alm = _apply_norm(alm, lmax, norm, csphase, False)
    res = _np.empty([lmax + 1, 2 * lmax + 1 + extend], dtype=_np.complex128)
    _synthesize_GLQ(alm, lmax_calc, extend, res.imag)
    alm = _ccilm2almr(cilm)
    alm = _apply_norm(alm, lmax, norm, csphase, False)
    _synthesize_GLQ(alm, lmax_calc, extend, res.real)
    return res


# weights and zeros are ignored (they are computed internally)
def SHExpandGLQ(gridglq, w=None, zero=None, norm=1, csphase=1, lmax_calc=None):
    """
    Expand a 2D grid sampled on the Gauss-Legendre quadrature nodes into spherical
    harmonics.

    Usage
    -----
    cilm = SHExpandGLQ (gridglq, [w, zero, norm, csphase, lmax_calc])

    Returns
    -------
    cilm : float, dimension (2, lmax+1, lmax+1) or (2, lmax_calc+1, lmax_calc+1)
        The real spherical harmonic coefficients of the function. The coefficients
        C0lm and Cilm refer to the "cosine" (Clm) and "sine" (Slm) coefficients,
        respectively, with Clm=cilm[0,l,m] and Slm=cilm[1,l,m].

    Parameters
    ----------
    gridglq : float, dimension (lmax+1, 2*lmax+1)
        A 2D grid of data sampled on the Gauss-Legendre quadrature nodes. The
        latitudinal nodes correspond to the zeros of the Legendre polynomial of
        degree lmax+1, and the longitudinal nodes are equally spaced with an
        interval of 360/(2*lmax+1) degrees. See also GLQGridCoord.
    w : optional, ignored
        This parameter only exists to maintain interface compatibility with the
        "shtools" backend.
    zero : optional, ignored
        This parameter only exists to maintain interface compatibility with the
        "shtools" backend.
    norm : optional, integer, default = 1
        1 (default) = Geodesy 4-pi normalized harmonics; 2 = Schmidt semi-normalized
        harmonics; 3 = unnormalized harmonics; 4 = orthonormal harmonics.
    csphase : optional, integer, default = 1
        1 (default) = do not apply the Condon-Shortley phase factor to the
        associated Legendre functions; -1 = append the Condon-Shortley phase factor
        of (-1)^m to the associated Legendre functions.
    lmax_calc : optional, integer, default = lmax
        The maximum spherical harmonic degree calculated in the spherical harmonic
        expansion.

    Description
    -----------
    SHExpandGLQ will expand a 2-dimensional grid of data sampled on the Gauss-
    Legendre quadrature nodes into spherical harmonics. This is the inverse of the
    routine MakeGridGLQ. The latitudinal nodes of the input grid correspond to the
    zeros of the Legendre polynomial of degree lmax+1, and the longitudinal nodes
    are equally spaced with an interval of 360/(2*lmax+1) degrees. It is implicitly
    assumed that the function is bandlimited to degree lmax. If the optional
    parameter lmax_calc is specified, the spherical harmonic coefficients will be
    calculated up to this degree, instead of lmax.

    The employed spherical harmonic normalization and Condon-Shortley phase
    convention can be set by the optional arguments norm and csphase; if not set,
    the default is to use geodesy 4-pi normalized harmonics that exclude the Condon-
    Shortley phase of (-1)^m. The normalized legendre functions are calculated in
    this routine using the recurrence given by Ishioka (2018), which are accurate to
    at least degree 100000. The unnormalized functions are accurate
    only to about degree 15.

    The zeros of the Legendre polynomials and the quadrature weights are computed
    using the method described by Bogaert (2014).

    References
    ----------
    Bogaert, I.: SIAM Journal on Scientific Computing, 36, A1008-A1026, 2014

    Ishioka, K.: Journal of the Meteorological Society of Japan, 96, 241−249, 2018
    """ # noqa
    gridglq = _fixdtype(gridglq)
    if lmax_calc is None:
        lmax_calc = gridglq.shape[0] - 1
    if lmax_calc > (gridglq.shape[0] - 1):
        raise RuntimeError("lmax_calc too high")
    alm = _analyze_GLQ(gridglq, lmax_calc)
    return _extract_alm(alm, lmax_calc, norm, csphase)


# weights and zeros are ignored (they are computed internally)
def SHExpandGLQC(gridglq, w=None, zero=None, norm=1, csphase=1,
                 lmax_calc=None):
    """Expand a 2D grid sampled on the Gauss-Legendre quadrature nodes into spherical
    harmonics.

    Usage
    -----
    cilm = SHExpandGLQC (gridglq, w, zero, [norm, csphase, lmax_calc])

    Returns
    -------
    cilm : complex, dimension (2, lmax+1, lmax+1) or (2, lmax_calc+1, lmax_calc+1)
        The complex spherical harmonic coefficients of the complex function. The
        first index specifies the coefficient corresponding to the positive and
        negative order of m, respectively, with Clm=cilm[0,l,m] and Cl,-m
        =cilm[1,l,m].

    Parameters
    ----------
    gridglq : complex, dimension (lmax+1, 2*lmax+1)
        A 2D grid of complex data sampled on the Gauss-Legendre quadrature nodes.
        The latitudinal nodes correspond to the zeros of the Legendre polynomial of
        degree lmax+1, and the longitudinal nodes are equally spaced with an
        interval of 360/(2*lmax+1) degrees. See also GLQGridCoord.
    w : optional, ignored
        This parameter only exists to maintain interface compatibility with the
        "shtools" backend.
    zero : optional, ignored
        This parameter only exists to maintain interface compatibility with the
        "shtools" backend.
    norm : optional, integer, default = 1
        1 (default) = 4-pi (geodesy) normalized harmonics; 2 = Schmidt semi-
        normalized harmonics; 3 = unnormalized harmonics; 4 = orthonormal harmonics.
    csphase : optional, integer, default = 1
        1 (default) = do not apply the Condon-Shortley phase factor to the
        associated Legendre functions; -1 = append the Condon-Shortley phase factor
        of (-1)^m to the associated Legendre functions.
    lmax_calc : optional, integer, default = lmax
        The maximum spherical harmonic degree calculated in the spherical harmonic
        expansion.

    Description
    -----------
    SHExpandGLQC will expand a 2-dimensional grid of complex data sampled on the
    Gauss-Legendre quadrature nodes into complex spherical harmonics. This is the
    inverse of the routine MakeGridGLQC. The latitudinal nodes of the input grid
    correspond to the zeros of the Legendre polynomial of degree lmax+1, and the
    longitudinal nodes are equally spaced with an interval of 360/(2*lmax+1)
    degrees. It is implicitly assumed that the function is bandlimited to degree
    lmax. If the optional parameter lmax_calc is specified, the spherical harmonic
    coefficients will be calculated up to this degree, instead of lmax.

    The employed spherical harmonic normalization and Condon-Shortley phase
    convention can be set by the optional arguments norm and csphase; if not set,
    the default is to use geodesy 4-pi normalized harmonics that exclude the Condon-
    Shortley phase of (-1)^m. The normalized legendre functions are calculated in
    this routine using the recurrence given by Ishioka (2018), which are accurate
    to at least degree 100000. The unnormalized functions are accurate
    only to about degree 15.

    The zeros of the Legendre polynomials and the quadrature weights are computed
    using the method described by Bogaert (2014).

    References
    ----------
    Bogaert, I.: SIAM Journal on Scientific Computing, 36, A1008-A1026, 2014

    Ishioka, K.: Journal of the Meteorological Society of Japan, 96, 241−249, 2018
    """ # noqa
    if lmax_calc is None:
        lmax_calc = gridglq.shape[0] - 1
    if lmax_calc > (gridglq.shape[0] - 1):
        raise RuntimeError("lmax_calc too high")
    res = _np.zeros((2, lmax_calc + 1, lmax_calc + 1), dtype=_np.complex128)
    alm = _analyze_GLQ(_fixdtype(gridglq.real), lmax_calc)
    alm = _apply_norm(alm, lmax_calc, norm, csphase, True)
    _addRealpart(res, alm)
    alm = _analyze_GLQ(_fixdtype(gridglq.imag), lmax_calc)
    alm = _apply_norm(alm, lmax_calc, norm, csphase, True)
    _addImagpart(res, alm)
    return res


def MakeGradientDH(
    cilm, lmax=None, sampling=1, lmax_calc=None, extend=False, radius=None
):
    """Compute the gradient of a scalar function and return grids of the two horizontal
    components that conform with Driscoll and Healy's (1994) sampling theorem.

    Usage
    -----
    theta, phi = MakeGradientDH (cilm, [lmax, sampling, lmax_calc, extend, radius])

    Returns
    -------
    theta : float, dimension (nlat, nlong)
        A 2D map of the theta component of the horizontal gradient that conforms to
        the sampling theorem of Driscoll and Healy (1994). If sampling is 1, the
        grid is equally sampled and is dimensioned as (n by n), where n is 2lmax+2.
        If sampling is 2, the grid is equally spaced and is dimensioned as (n by
        2n). The first latitudinal band of the grid corresponds to 90 N, the
        latitudinal sampling interval is 180/n degrees, and the default behavior is
        to exclude the latitudinal band for 90 S. The first longitudinal band of the
        grid is 0 E, by default the longitudinal band for 360 E is not included, and
        the longitudinal sampling interval is 360/n for an equally sampled and 180/n
        for an equally spaced grid, respectively. If extend is 1, the longitudinal
        band for 360 E and the latitudinal band for 90 S will be included, which
        increases each of the dimensions of the grid by 1.
    phi : float, dimension (nlat, nlong)
        A 2D equally sampled or equally spaced grid of the phi component of the
        horizontal gradient.

    Parameters
    ----------
    cilm : float, dimension (2, lmaxin+1, lmaxin+1)
        The real 4-pi normalized spherical harmonic coefficients of a scalar
        function. The coefficients c1lm and c2lm refer to the cosine and sine
        coefficients, respectively, with c1lm=cilm[0,l,m] and c2lm=cilm[1,l,m].
    lmax : optional, integer, default = lmaxin
        The maximum spherical harmonic degree of the coefficients cilm. This
        determines the number of samples of the output grids, n=2lmax+2, and the
        latitudinal sampling interval, 90/(lmax+1).
    sampling : optional, integer, default = 2
        If 1 (default) the output grids are equally sampled (n by n). If 2, the
        grids are equally spaced (n by 2n).
    lmax_calc : optional, integer, default = lmax
        The maximum spherical harmonic degree used in evaluating the functions. This
        must be less than or equal to lmax.
    extend : optional, bool, default = False
        If True, compute the longitudinal band for 360 E and the latitudinal band
        for 90 S. This increases each of the dimensions of griddh by 1.
    radius : optional, float, default = 1.0
        The radius of the sphere used when computing the gradient of the function.


    Description
    -----------
    MakeGradientDH will compute the horizontal gradient of a scalar function on a
    sphere defined by the spherical harmonic coefficients cilm. The output grids of
    the theta and phi components of the gradient are either equally sampled (n by n)
    or equally spaced (n by 2n) in latitude and longitude. The gradient is given by
    the formula

    Grad F = 1/r dF/theta theta-hat + 1/(r sin theta) dF/dphi phi-hat.

    where theta is colatitude and phi is longitude. The radius r is by default set
    to 1, but this can be modified by use of the optional parameter radius.

    The default is to use an input grid that is equally sampled (n by n), but this
    can be changed to use an equally spaced grid (n by 2n) by the optional argument
    sampling. The redundant longitudinal band for 360 E and the latitudinal band for
    90 S are excluded by default, but these can be computed by specifying the
    optional argument extend.

    Reference
    ---------
    Driscoll, J.R. and D.M. Healy, Computing Fourier transforms and convolutions on
    the 2-sphere, Adv. Appl. Math., 15, 202-250, 1994.
    """ # noqa
    lmax, lmax_calc, cilm = _prep_lmax(lmax, lmax_calc, cilm)
    alm = _make_alm(cilm, lmax_calc, norm=1, csphase=1)
    res = _np.empty((2, 2 * lmax + 2 + extend,
                     sampling * (2 * lmax + 2) + extend))
    res = _synthesize_DH_deriv1(alm, lmax_calc, extend, res)
    if radius is not None:
        res *= 1.0 / radius
    return res
