import numpy as np

try:
    import ducc0
except:
    ducc0 = None

# setup a few required variables
if ducc0 is not None:
    import os

    try:
        nthreads = int(os.environ["OMP_NUM_THREADS"])
    except:
        nthreads = 0


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
        return np.full(lmax + 1, np.sqrt(4 * np.pi))
    if norm == 2:
        return np.sqrt(4 * np.pi / (2 * np.arange(lmax + 1) + 1.0))
    if norm == 3:
        return np.sqrt(2 * np.pi / (2 * np.arange(lmax + 1) + 1.0))
    if norm == 4:
        return np.ones(lmax + 1)
    raise RuntimeError("unsupported normalization")


def _make_alm(cilm, lmax, norm, csphase):
    alm = np.empty((_nalm(lmax, lmax),), dtype=np.complex128)
    lnorm = _get_norm(lmax, norm)
    alm[0 : lmax + 1] = cilm[0, 0 : lmax + 1, 0]
    alm[0 : lmax + 1] *= lnorm[0 : lmax + 1]
    lnorm /= np.sqrt(2.0)
    mlnorm = -lnorm
    ofs = lmax + 1
    for m in range(1, lmax + 1):
        alm[ofs : ofs + lmax + 1 - m].real = cilm[0, m : lmax + 1, m]
        alm[ofs : ofs + lmax + 1 - m].imag = cilm[1, m : lmax + 1, m]
        if csphase == 1:
            if m & 1:
                alm[ofs : ofs + lmax + 1 - m].real *= mlnorm[m : lmax + 1]
                alm[ofs : ofs + lmax + 1 - m].imag *= lnorm[m : lmax + 1]
            else:
                alm[ofs : ofs + lmax + 1 - m].real *= lnorm[m : lmax + 1]
                alm[ofs : ofs + lmax + 1 - m].imag *= mlnorm[m : lmax + 1]
        else:
            alm[ofs : ofs + lmax + 1 - m].real *= lnorm[m : lmax + 1]
            alm[ofs : ofs + lmax + 1 - m].imag *= mlnorm[m : lmax + 1]
        ofs += lmax + 1 - m
    if norm == 3:  # special treatment for unnormalized a_lm
        r = np.arange(lmax + 1)
        fct = np.ones(lmax + 1)
        ofs = lmax + 1
        for m in range(1, lmax + 1):
            fct[m:] *= np.sqrt((r[m:] + m) * (r[m:] - m + 1))
            alm[ofs : ofs + lmax + 1 - m] *= fct[m:]
            ofs += lmax + 1 - m
        alm[0 : lmax + 1] *= np.sqrt(2)
    return alm


def _extract_alm(alm, lmax, norm, csphase):
    almx = alm.reshape((-1,))
    cilm = np.zeros((2, lmax + 1, lmax + 1), dtype=np.float64)
    lnorm = 1.0 / _get_norm(lmax, norm)
    cilm[0, 0 : lmax + 1, 0] = almx[0 : lmax + 1].real * lnorm[0 : lmax + 1]
    lnorm *= np.sqrt(2.0)
    mlnorm = -lnorm
    tmp = np.empty(lmax + 1, dtype=np.complex128)
    ofs = lmax + 1
    for m in range(1, lmax + 1):
        tmp[m:] = almx[ofs : ofs + lmax + 1 - m]
        if csphase == 1:
            if m & 1:
                tmp[m:].real *= mlnorm[m : lmax + 1]
                tmp[m:].imag *= lnorm[m : lmax + 1]
            else:
                tmp[m:].real *= lnorm[m : lmax + 1]
                tmp[m:].imag *= mlnorm[m : lmax + 1]
        else:
            tmp[m:].real *= lnorm[m : lmax + 1]
            tmp[m:].imag *= mlnorm[m : lmax + 1]

        cilm[0, m : lmax + 1, m] = tmp[m:].real
        cilm[1, m : lmax + 1, m] = tmp[m:].imag
        ofs += lmax + 1 - m
    if norm == 3:  # special treatment for unnormalized a_lm
        r = np.arange(lmax + 1)
        fct = np.ones(lmax + 1)
        for m in range(1, lmax + 1):
            fct[m:] /= np.sqrt((r[m:] + m) * (r[m:] - m + 1))
            cilm[:, m:, m] *= fct[m:]
        cilm[:, 0] /= np.sqrt(2)
    return cilm


def _synthesize_DH(alm, lmax, shp, extend):
    ntheta, nphi = shp
    map_out = np.empty((ntheta + extend, nphi + extend), dtype=np.float64)
    ducc0.sht.experimental.synthesis_2d(
        alm=alm.reshape((1,-1)),
        map=map_out[:,:nphi].reshape((1,ntheta+extend,nphi)),
        spin=0,
        lmax=lmax,
        geometry = "CC" if extend else "DH",
        nthreads=nthreads)
    if extend:
        map_out[:, -1] = map_out[:, 0]
    return map_out


def _synthesize_DH_deriv1(alm, lmax, shp, extend):
    ntheta, nphi = shp
    map_out = np.empty((2, ntheta + extend, nphi + extend), dtype=np.float64)
    ducc0.sht.experimental.synthesis_2d_deriv1(
        alm=alm.reshape((1,-1)),
        map=map_out[:,:,:nphi],
        lmax=lmax,
        geometry = "CC" if extend else "DH",
        nthreads=nthreads)
    map_out[:, 0, :] = 0.0
    if extend:
        map_out[:, -1, :] = 0.0
        map_out[:, :, -1] = map_out[:, :, 0]
    return map_out


def _synthesize_GLQ(alm, lmax, shp, extend):
    ntheta, nphi = shp
    map_out = np.empty((ntheta, nphi + extend), dtype=np.float64)
    ducc0.sht.experimental.synthesis_2d(
        alm=alm.reshape((1,-1)),
        map=map_out[:,:nphi].reshape((1,ntheta,nphi)),
        spin=0,
        lmax=lmax,
        geometry = "GL",
        nthreads=nthreads)
    if extend:
        map_out[:, -1] = map_out[:, 0]
    return map_out


def _analyze_DH(map, lmax):
    ntheta, nphi = map.shape
    alm = np.empty(_nalm(lmax,lmax), dtype=np.complex128)
    ducc0.sht.experimental.analysis_2d(
        alm=alm.reshape((1,-1)),
        map=map.reshape((1,ntheta,nphi)),
        spin=0,
        lmax=lmax,
        geometry = "DH",
        nthreads=nthreads)
    return alm


def _analyze_GLQ(map, lmax):
    ntheta, nphi = map.shape
    alm = np.empty(_nalm(lmax,lmax), dtype=np.complex128)
    ducc0.sht.experimental.analysis_2d(
        alm=alm.reshape((1,-1)),
        map=map.reshape((1,ntheta,nphi)),
        spin=0,
        lmax=lmax,
        geometry = "GL",
        nthreads=nthreads)
    return alm


def _extractRealPart(ccoeffs):
    fct = (-1) ** np.arange(ccoeffs.shape[2]).reshape((1, -1))
    tmp = np.conj(ccoeffs[1])
    tmp *= fct
    tmp += ccoeffs[0]
    tmp *= 1.0 / np.sqrt(2.0)
    rpart = np.empty(ccoeffs.shape, dtype=np.float64)
    rpart[0] = tmp.real
    rpart[1] = -tmp.imag
    rpart[0, :, 0] = ccoeffs[0, :, 0].real
    rpart[1, :, 0] = 0.0
    return rpart


def _extractImagPart(ccoeffs):
    fct = (-1) ** np.arange(ccoeffs.shape[2]).reshape((1, -1))
    tmp = np.conj(ccoeffs[1])
    tmp *= -fct
    tmp += ccoeffs[0]
    tmp *= 1.0 / np.sqrt(2.0)
    ipart = np.empty(ccoeffs.shape, dtype=np.float64)
    ipart[0] = tmp.imag
    ipart[1] = tmp.real
    ipart[0, :, 0] = ccoeffs[0, :, 0].imag
    ipart[1, :, 0] = 0.0
    return ipart


def _combineComplexCoef(rpart, ipart):
    fct = (-1) ** np.arange(rpart.shape[2]).reshape((1, 1, -1))
    res = np.empty(rpart.shape, dtype=np.complex128)
    res[0].real = rpart[0] + ipart[1]
    res[0].imag = -rpart[1] + ipart[0]
    res[1].real = fct * (rpart[0] - ipart[1])
    res[1].imag = fct * (rpart[1] + ipart[0])
    res *= 1.0 / np.sqrt(2.0)
    res[0, :, 0] = rpart[0, :, 0] + 1j * ipart[0, :, 0]
    res[1, :, 0] = 0.0
    return res


def SHRotateRealCoef(rcoeffs, angles):
    lmax = rcoeffs.shape[1] - 1
    alm = _make_alm(rcoeffs, lmax, 1, 1)
    alm = ducc0.sht.rotate_alm(
        alm, lmax, -angles[0], -angles[1], -angles[2], nthreads=nthreads
    )
    return _extract_alm(alm, lmax, 1, 1)


def SHRotateComplexCoef(ccoeffs, angles):
    lmax = ccoeffs.shape[1] - 1
    cilmR = _extractRealPart(ccoeffs)
    alm = _make_alm(cilmR, lmax, 1, 1)
    alm = ducc0.sht.rotate_alm(
        alm, lmax, -angles[0], -angles[1], -angles[2], nthreads=nthreads
    )
    cilmR = _extract_alm(alm, lmax, 1, 1)
    del alm
    cilmI = _extractImagPart(ccoeffs)
    alm = _make_alm(cilmI, lmax, 1, 1)
    alm = ducc0.sht.rotate_alm(
        alm, lmax, -angles[0], -angles[1], -angles[2], nthreads=nthreads
    )
    cilmI = _extract_alm(alm, lmax, 1, 1)
    del alm
    return _combineComplexCoef(cilmR, cilmI)


def MakeGridDH(
    cilm, lmax=None, norm=1, sampling=1, csphase=1, lmax_calc=None, extend=False
):
    if lmax is None:
        lmax = cilm.shape[1] - 1
    cilm = cilm[:, : lmax + 1, : lmax + 1]
    if lmax_calc is None:
        lmax_calc = cilm.shape[1] - 1
    alm = _make_alm(cilm, lmax_calc, norm, csphase)
    return _synthesize_DH(
        alm, lmax_calc, [2 * lmax + 2, sampling * (2 * lmax + 2)], extend
    )


def MakeGridDHC(
    cilm, lmax=None, norm=1, sampling=1, csphase=1, lmax_calc=None, extend=False
):
    if lmax is None:
        lmax = cilm.shape[1] - 1
    cilm = cilm[:, : lmax + 1, : lmax + 1]
    if lmax_calc is None:
        lmax_calc = cilm.shape[1] - 1
    cilmx = _extractImagPart(cilm)
    alm = _make_alm(cilmx, lmax_calc, norm, csphase)
    del cilmx
    res = 1j * _synthesize_DH(
        alm, lmax_calc, [2 * lmax + 2, sampling * (2 * lmax + 2)], extend
    )
    cilmx = _extractRealPart(cilm)
    alm = _make_alm(cilmx, lmax_calc, norm, csphase)
    del cilmx
    res += _synthesize_DH(
        alm, lmax_calc, [2 * lmax + 2, sampling * (2 * lmax + 2)], extend
    )
    return res


def SHExpandDH(grid, norm=1, sampling=1, csphase=1, lmax_calc=None):
    if grid.shape[1] != sampling * grid.shape[0]:
        raise RuntimeError("grid resolution mismatch")
    if lmax_calc is None:
        lmax_calc = grid.shape[0] // 2 - 1
    if lmax_calc > (grid.shape[0] // 2 - 1):
        raise RuntimeError("lmax_calc too high")
    lmax = grid.shape[0] // 2 - 1 if lmax_calc is None else lmax_calc
    alm = _analyze_DH(grid, lmax_calc)
    return _extract_alm(alm, lmax_calc, norm, csphase)


def SHExpandDHC(grid, norm=1, sampling=1, csphase=1, lmax_calc=None):
    if grid.shape[1] != sampling * grid.shape[0]:
        raise RuntimeError("grid resolution mismatch")
    if lmax_calc is None:
        lmax_calc = grid.shape[0] // 2 - 1
    if lmax_calc > (grid.shape[0] // 2 - 1):
        raise RuntimeError("lmax_calc too high")
    lmax = grid.shape[0] // 2 - 1 if lmax_calc is None else lmax_calc
    alm = _analyze_DH(grid.real, lmax_calc)
    cilmR = _extract_alm(alm, lmax_calc, norm, csphase)
    alm = _analyze_DH(grid.imag, lmax_calc)
    cilmI = _extract_alm(alm, lmax_calc, norm, csphase)
    del alm
    return _combineComplexCoef(cilmR, cilmI)


def MakeGridGLQ(
    cilm, zeros, lmax=None, norm=1, csphase=1, lmax_calc=None, extend=False
):
    if lmax is None:
        lmax = cilm.shape[1] - 1
    cilm = cilm[:, : lmax + 1, : lmax + 1]
    if lmax_calc is None:
        lmax_calc = cilm.shape[1] - 1
    alm = _make_alm(cilm, lmax_calc, norm, csphase)
    return _synthesize_GLQ(alm, lmax_calc, [lmax + 1, 2 * lmax + 1], extend)


def MakeGridGLQC(
    cilm, zeros, lmax=None, norm=1, csphase=1, lmax_calc=None, extend=False
):
    if lmax is None:
        lmax = cilm.shape[1] - 1
    cilm = cilm[:, : lmax + 1, : lmax + 1]
    if lmax_calc is None:
        lmax_calc = cilm.shape[1] - 1
    cilmx = _extractImagPart(cilm)
    alm = _make_alm(cilmx, lmax_calc, norm, csphase)
    del cilmx
    res = 1j * _synthesize_GLQ(alm, lmax_calc, [lmax + 1, 2 * lmax + 1], extend)
    cilmx = _extractRealPart(cilm)
    alm = _make_alm(cilmx, lmax_calc, norm, csphase)
    del cilmx
    res += _synthesize_GLQ(alm, lmax_calc, [lmax + 1, 2 * lmax + 1], extend)
    return res


def SHExpandGLQ(grid, weights, zeros, norm=1, csphase=1, lmax_calc=None):
    if lmax_calc is None:
        lmax_calc = grid.shape[0] - 1
    if lmax_calc > (grid.shape[0] - 1):
        raise RuntimeError("lmax_calc too high")
    alm = _analyze_GLQ(grid, lmax_calc)
    return _extract_alm(alm, lmax_calc, norm, csphase)


def SHExpandGLQC(grid, weights, zeros, norm=1, csphase=1, lmax_calc=None):
    if lmax_calc is None:
        lmax_calc = grid.shape[0] - 1
    if lmax_calc > (grid.shape[0] - 1):
        raise RuntimeError("lmax_calc too high")
    alm = _analyze_GLQ(grid.real, lmax_calc)
    cilmR = _extract_alm(alm, lmax_calc, norm, csphase)
    alm = _analyze_GLQ(grid.imag, lmax_calc)
    cilmI = _extract_alm(alm, lmax_calc, norm, csphase)
    del alm
    return _combineComplexCoef(cilmR, cilmI)


def MakeGradientDH(cilm, lmax=None, sampling=1, lmax_calc=None, extend=False):
    if lmax is None:
        lmax = cilm.shape[1] - 1
    cilm = cilm[:, : lmax + 1, : lmax + 1]
    if lmax_calc is None:
        lmax_calc = cilm.shape[1] - 1
    alm = _make_alm(cilm, lmax_calc, norm=1, csphase=1)
    return (1.0 / cilm[0, 0, 0]) * _synthesize_DH_deriv1(
        alm, lmax_calc, [2 * lmax + 2, sampling * (2 * lmax + 2)], extend
    )
