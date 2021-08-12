import numpy as np

try:
    import ducc0

    major, minor, patch = ducc0.__version__.split(".")
    if int(major) < 1 and int(minor) < 15:
        raise RuntimeError
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


def _rcilm2alm(cilm, lmax):
    alm = np.empty((_nalm(lmax, lmax),), dtype=np.complex128)
    alm[0:lmax + 1] = cilm[0, :, 0]
    ofs = lmax + 1
    for m in range(1, lmax + 1):
        alm[ofs:ofs + lmax + 1 - m].real = cilm[0, m:, m]
        alm[ofs:ofs + lmax + 1 - m].imag = cilm[1, m:, m]
        ofs += lmax + 1 - m
    return alm


def _ralm2cilm(alm, lmax):
    cilm = np.zeros((2, lmax + 1, lmax + 1), dtype=np.float64)
    cilm[0, :, 0] = alm[0:lmax + 1].real
    ofs = lmax + 1
    for m in range(1, lmax + 1):
        cilm[0, m:, m] = alm[ofs:ofs + lmax + 1 - m].real
        cilm[1, m:, m] = alm[ofs:ofs + lmax + 1 - m].imag
        ofs += lmax + 1 - m
    return cilm


def _apply_norm(alm, lmax, norm, csphase, reverse):
    lnorm = _get_norm(lmax, norm)
    if reverse:
        lnorm = 1.0 / lnorm
    alm[0:lmax + 1] *= lnorm[0:lmax + 1]
    lnorm *= np.sqrt(2.0) if reverse else (1.0 / np.sqrt(2.0))
    mlnorm = -lnorm
    ofs = lmax + 1
    for m in range(1, lmax + 1):
        if csphase == 1:
            if m & 1:
                alm[ofs:ofs + lmax + 1 - m].real *= mlnorm[m:]
                alm[ofs:ofs + lmax + 1 - m].imag *= lnorm[m:]
            else:
                alm[ofs:ofs + lmax + 1 - m].real *= lnorm[m:]
                alm[ofs:ofs + lmax + 1 - m].imag *= mlnorm[m:]
        else:
            alm[ofs:ofs + lmax + 1 - m].real *= lnorm[m:]
            alm[ofs:ofs + lmax + 1 - m].imag *= mlnorm[m:]
        ofs += lmax + 1 - m
    if norm == 3:  # special treatment for unnormalized a_lm
        r = np.arange(lmax + 1)
        fct = np.ones(lmax + 1)
        ofs = lmax + 1
        if reverse:
            alm[0:lmax + 1] /= np.sqrt(2)
            for m in range(1, lmax + 1):
                fct[m:] *= np.sqrt((r[m:] + m) * (r[m:] - m + 1))
                alm[ofs:ofs + lmax + 1 - m] /= fct[m:]
                ofs += lmax + 1 - m
        else:
            alm[0:lmax + 1] *= np.sqrt(2)
            for m in range(1, lmax + 1):
                fct[m:] *= np.sqrt((r[m:] + m) * (r[m:] - m + 1))
                alm[ofs:ofs + lmax + 1 - m] *= fct[m:]
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
    alm = np.empty((_nalm(lmax, lmax),), dtype=np.complex128)
    fct = (-1) ** np.arange(lmax + 1)
    alm[0:lmax + 1] = cilm[0, :, 0].real
    ofs = lmax + 1
    for m in range(1, lmax + 1):
        tmp = np.conj(cilm[1, m:, m])
        tmp *= fct[m]
        tmp += cilm[0, m:, m]
        tmp *= 1.0 / np.sqrt(2.0)
        alm[ofs:ofs + lmax + 1 - m] = np.conj(tmp)
        ofs += lmax + 1 - m
    return alm


def _ccilm2almi(cilm):
    lmax = cilm.shape[1] - 1
    alm = np.empty((_nalm(lmax, lmax),), dtype=np.complex128)
    fct = (-1) ** np.arange(lmax + 1)
    alm[0:lmax + 1] = cilm[0, :, 0].imag
    ofs = lmax + 1
    for m in range(1, lmax + 1):
        tmp = np.conj(cilm[1, m:, m])
        tmp *= -fct[m]
        tmp += cilm[0, m:, m]
        tmp *= 1.0 / np.sqrt(2.0)
        alm[ofs:ofs + lmax + 1 - m] = tmp.imag + 1j * tmp.real
        ofs += lmax + 1 - m
    return alm


def _addRealpart(cilm, alm):
    lmax = cilm.shape[1] - 1
    cilm[0, :, 0].real += alm[0:lmax + 1].real
    ofs = lmax + 1
    for m in range(1, lmax + 1):
        tmp = alm[ofs:ofs + lmax + 1 - m] / np.sqrt(2.0)
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
    cilm[0, :, 0].imag += alm[0:lmax + 1].real
    ofs = lmax + 1
    for m in range(1, lmax + 1):
        tmp = alm[ofs:ofs + lmax + 1 - m] / np.sqrt(2.0)
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


# dj_matrix is ignored
def SHRotateRealCoef(rcoeffs, angles, dj_matrix=None):
    lmax = rcoeffs.shape[1] - 1
    alm = _make_alm(rcoeffs, lmax, 1, 1)
    alm = ducc0.sht.rotate_alm(
        alm, lmax, -angles[0], -angles[1], -angles[2], nthreads=nthreads
    )
    return _extract_alm(alm, lmax, 1, 1)


# dj_matrix is ignored
def SHRotateComplexCoef(ccoeffs, angles, dj_matrix=None):
    lmax = ccoeffs.shape[1] - 1
    alm = _ccilm2almr(ccoeffs)
    alm = _apply_norm(alm, lmax, 1, 1, False)
    alm = ducc0.sht.rotate_alm(
        alm, lmax, -angles[0], -angles[1], -angles[2], nthreads=nthreads
    )
    alm = _apply_norm(alm, lmax, 1, 1, True)
    res = np.zeros((2, lmax + 1, lmax + 1), dtype=np.complex128)
    _addRealpart(res, alm)
    alm = _ccilm2almi(ccoeffs)
    alm = _apply_norm(alm, lmax, 1, 1, False)
    alm = ducc0.sht.rotate_alm(
        alm, lmax, -angles[0], -angles[1], -angles[2], nthreads=nthreads
    )
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
    if lmax is None:
        lmax = cilm.shape[1] - 1
    if lmax_calc is None:
        lmax_calc = cilm.shape[1] - 1
    cilm = cilm[:, : lmax_calc + 1, : lmax_calc + 1]
    alm = _make_alm(cilm, lmax_calc, norm, csphase)
    out = np.empty([2 * lmax + 2 + extend, sampling * (2 * lmax + 2) + extend])
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
    if lmax is None:
        lmax = cilm.shape[1] - 1
    if lmax_calc is None:
        lmax_calc = cilm.shape[1] - 1
    cilm = cilm[:, : lmax_calc + 1, : lmax_calc + 1]
    alm = _ccilm2almi(cilm)
    alm = _apply_norm(alm, lmax, norm, csphase, False)
    res = np.empty(
        [2 * lmax + 2 + extend, sampling * (2 * lmax + 2) + extend],
        dtype=np.complex128,
    )
    _synthesize_DH(alm, lmax_calc, extend, res.imag)
    alm = _ccilm2almr(cilm)
    alm = _apply_norm(alm, lmax, norm, csphase, False)
    _synthesize_DH(alm, lmax_calc, extend, res.real)
    return res


def SHExpandDH(grid, norm=1, sampling=1, csphase=1, lmax_calc=None):
    if grid.shape[1] != sampling * grid.shape[0]:
        raise RuntimeError("grid resolution mismatch")
    if lmax_calc is None:
        lmax_calc = grid.shape[0] // 2 - 1
    if lmax_calc > (grid.shape[0] // 2 - 1):
        raise RuntimeError("lmax_calc too high")
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
    res = np.zeros((2, lmax + 1, lmax + 1), dtype=np.complex128)
    alm = _analyze_DH(grid.real, lmax_calc)
    alm = _apply_norm(alm, lmax, norm, csphase, True)
    _addRealpart(res, alm)
    alm = _analyze_DH(grid.imag, lmax_calc)
    alm = _apply_norm(alm, lmax, norm, csphase, True)
    _addImagpart(res, alm)
    return res


# zero is ignored (they are computed internally)
def MakeGridGLQ(
    cilm, lmax=None, zero=None, norm=1, csphase=1, lmax_calc=None, extend=False
):
    if lmax is None:
        lmax = cilm.shape[1] - 1
    if lmax_calc is None:
        lmax_calc = cilm.shape[1] - 1
    cilm = cilm[:, : lmax_calc + 1, : lmax_calc + 1]
    alm = _make_alm(cilm, lmax_calc, norm, csphase)
    out = np.empty([lmax + 1, (2 * lmax + 1) + extend])
    return _synthesize_GLQ(alm, lmax_calc, extend, out)


# zero is ignored (they are computed internally)
def MakeGridGLQC(
    cilm, lmax=None, zero=None, norm=1, csphase=1, lmax_calc=None, extend=False
):
    if lmax is None:
        lmax = cilm.shape[1] - 1
    if lmax_calc is None:
        lmax_calc = cilm.shape[1] - 1
    cilm = cilm[:, : lmax_calc + 1, : lmax_calc + 1]
    alm = _ccilm2almi(cilm)
    alm = _apply_norm(alm, lmax, norm, csphase, False)
    res = np.empty([lmax + 1, 2 * lmax + 1 + extend], dtype=np.complex128)
    _synthesize_GLQ(alm, lmax_calc, extend, res.imag)
    alm = _ccilm2almr(cilm)
    alm = _apply_norm(alm, lmax, norm, csphase, False)
    _synthesize_GLQ(alm, lmax_calc, extend, res.real)
    return res


# weights and zeros are ignored (they are computed internally)
def SHExpandGLQ(
    grid, weights=None, zeros=None, norm=1, csphase=1, lmax_calc=None
):
    if lmax_calc is None:
        lmax_calc = grid.shape[0] - 1
    if lmax_calc > (grid.shape[0] - 1):
        raise RuntimeError("lmax_calc too high")
    alm = _analyze_GLQ(grid, lmax_calc)
    return _extract_alm(alm, lmax_calc, norm, csphase)


# weights and zeros are ignored (they are computed internally)
def SHExpandGLQC(
    grid, weights=None, zeros=None, norm=1, csphase=1, lmax_calc=None
):
    if lmax_calc is None:
        lmax_calc = grid.shape[0] - 1
    if lmax_calc > (grid.shape[0] - 1):
        raise RuntimeError("lmax_calc too high")
    res = np.zeros((2, lmax_calc + 1, lmax_calc + 1), dtype=np.complex128)
    alm = _analyze_GLQ(grid.real, lmax_calc)
    alm = _apply_norm(alm, lmax_calc, norm, csphase, True)
    _addRealpart(res, alm)
    alm = _analyze_GLQ(grid.imag, lmax_calc)
    alm = _apply_norm(alm, lmax_calc, norm, csphase, True)
    _addImagpart(res, alm)
    return res


def MakeGradientDH(
    cilm, lmax=None, sampling=1, lmax_calc=None, extend=False, radius=None
):
    if lmax is None:
        lmax = cilm.shape[1] - 1
    cilm = cilm[:, : lmax + 1, : lmax + 1]
    if lmax_calc is None:
        lmax_calc = cilm.shape[1] - 1
    alm = _make_alm(cilm, lmax_calc, norm=1, csphase=1)
    res = np.empty(
        (2, 2 * lmax + 2 + extend, sampling * (2 * lmax + 2) + extend)
    )
    res = _synthesize_DH_deriv1(alm, lmax_calc, extend, res)
    if radius is not None:
        res *= 1.0 / radius
    return res
