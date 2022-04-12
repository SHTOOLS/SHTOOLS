import numpy as np
import pyshtools as pysh
from time import time

nthreads = 1


def _l2error(a, b):
    return np.sqrt(np.sum(np.abs(a - b) ** 2) / np.sum(np.abs(a) ** 2))


# force SHTOOLS to deallocate temporary buffers
def flush_buffers(grd):
    degrees = np.arange(1, dtype=float)
    degrees[0] = np.inf
    power = degrees ** (-2)

    pysh.backends.select_preferred_backend("shtools")
    clm = pysh.SHCoeffs.from_random(power, seed=12345)
    grid2 = clm.expand(grid=grd)
    _ = grid2.expand()
    clm = pysh.SHCoeffs.from_random(power, seed=12345, kind="complex")
    grid2 = clm.expand(grid=grd)
    _ = grid2.expand()


def test_SHT(lmax, grd, csphase, normalization, extend):
    degrees = np.arange(lmax + 1, dtype=float)
    degrees[0] = np.inf
    power = degrees ** (-2)

    clm = pysh.SHCoeffs.from_random(power, seed=12345)

    clm = clm.convert(normalization=normalization, csphase=csphase, lmax=lmax)

    pysh.backends.select_preferred_backend("ducc", nthreads=nthreads)
    t0 = time()
    grid = clm.expand(grid=grd, extend=extend)
    cilm = grid.expand()
    tducc = time() - t0
    pysh.backends.select_preferred_backend("shtools")
    t0 = time()
    grid2 = clm.expand(grid=grd, extend=extend)
    cilm2 = grid2.expand()
    tshtools = time() - t0

    flush_buffers(grd)

    return (
        _l2error(grid.to_array(), grid2.to_array())
        + _l2error(cilm.to_array(), cilm2.to_array()),
        tshtools / tducc,
    )


def test_SHTducc(lmax, grd, nthreads):
    degrees = np.arange(lmax + 1, dtype=float)
    degrees[0] = np.inf
    power = degrees ** (-2)

    clm = pysh.SHCoeffs.from_random(power, seed=12345)

    clm = clm.convert(normalization="ortho", csphase=1, lmax=lmax)

    pysh.backends.select_preferred_backend("ducc", nthreads=nthreads)
    t0 = time()
    grid = clm.expand(grid=grd, extend=False)
    cilm = grid.expand(normalization="ortho", csphase=1)
    tducc = time() - t0

    return _l2error(clm.to_array(), cilm.to_array()), tducc


def test_SHTC(lmax, grd, csphase, normalization, extend):
    degrees = np.arange(lmax + 1, dtype=float)
    degrees[0] = np.inf
    power = degrees ** (-2)

    clm = pysh.SHCoeffs.from_random(power, seed=12345, kind="complex")

    clm = clm.convert(normalization=normalization, csphase=csphase, lmax=lmax)

    pysh.backends.select_preferred_backend("ducc", nthreads=nthreads)
    t0 = time()
    grid = clm.expand(grid=grd, extend=extend)
    cilm = grid.expand(normalization=normalization, csphase=csphase)
    tducc = time() - t0
    pysh.backends.select_preferred_backend("shtools")
    t0 = time()
    grid2 = clm.expand(grid=grd, extend=extend)
    cilm2 = grid2.expand(normalization=normalization, csphase=csphase)
    tshtools = time() - t0

    flush_buffers(grd)

    return (
        _l2error(grid.to_array(), grid2.to_array())
        + _l2error(cilm.to_array(), cilm2.to_array()),
        tshtools / tducc,
    )


def test_SHT_deriv(lmax, grd, csphase, extend):
    degrees = np.arange(lmax + 1, dtype=float)
    degrees[0] = 1.0
    power = degrees ** (-2)

    clm = pysh.SHCoeffs.from_random(power, seed=12345)

    clm = clm.convert(csphase=csphase, lmax=lmax)

    pysh.backends.select_preferred_backend("ducc", nthreads=nthreads)
    t0 = time()
    grad = clm.gradient(extend=extend, radius=3.4)
    tducc = time() - t0
    pysh.backends.select_preferred_backend("shtools")
    t0 = time()
    grad2 = clm.gradient(extend=extend, radius=1.0)
    tshtools = time() - t0

    flush_buffers(grd)

    return (
        _l2error(3.4 * grad.phi.to_array(), grad2.phi.to_array())
        + _l2error(3.4 * grad.theta.to_array(), grad2.theta.to_array()),
        tshtools / tducc,
    )


def test_rot(lmax, alpha, beta, gamma):
    degrees = np.arange(lmax + 1, dtype=float)
    degrees[0] = np.inf
    power = degrees ** (-2)

    clm = pysh.SHCoeffs.from_random(power, seed=12345)

    pysh.backends.select_preferred_backend("ducc", nthreads=nthreads)
    t0 = time()
    clm_rotated = clm.rotate(alpha, beta, gamma, degrees=True)
    tducc = time() - t0
    pysh.backends.select_preferred_backend("shtools")
    t0 = time()
    clm_rotated2 = clm.rotate(alpha, beta, gamma, degrees=True)
    tshtools = time() - t0

    return (
        _l2error(clm_rotated.to_array(), clm_rotated2.to_array()),
        tshtools / tducc,
    )


def test_rotc(lmax, alpha, beta, gamma):
    degrees = np.arange(lmax + 1, dtype=float)
    degrees[0] = np.inf
    power = degrees ** (-2)

    clm = pysh.SHCoeffs.from_random(power, seed=12345, kind="complex")

    pysh.backends.select_preferred_backend("ducc", nthreads=nthreads)
    t0 = time()
    clm_rotated = clm.rotate(alpha, beta, gamma, degrees=True)
    tducc = time() - t0
    pysh.backends.select_preferred_backend("shtools")
    t0 = time()
    clm_rotated2 = clm.rotate(alpha, beta, gamma, degrees=True)
    tshtools = time() - t0

    return (
        _l2error(clm_rotated.to_array(), clm_rotated2.to_array()),
        tshtools / tducc,
    )


def test_rot2(lmax, alpha, beta, gamma):
    degrees = np.arange(lmax + 1, dtype=float)
    degrees[0] = np.inf
    power = degrees ** (-2)

    clm = pysh.SHCoeffs.from_random(power, seed=12345)

    pysh.backends.select_preferred_backend("ducc", nthreads=nthreads)
    t0 = time()
    clm_rotated = clm.rotate(alpha, beta, gamma, degrees=True)
    clm_rotated = clm_rotated.rotate(-gamma, -beta, -alpha, degrees=True)
    tducc = time() - t0

    return _l2error(clm.to_array(), clm_rotated.to_array()), tducc


lmax_list = [127, 255, 511, 1023]

print("SHRealCoeff rotation tests:")
for lmax in lmax_list:
    for alpha in [47]:
        res = test_rot(lmax, alpha, 27, 59)
        print(
            "lmax={:4}: L2 error={:e}, speedup factor={:f}".format(
                lmax, res[0], res[1]
            )
        )

print("SHComplexCoeff rotation tests:")
for lmax in lmax_list:
    for alpha in [47]:
        res = test_rotc(lmax, alpha, 27, 59)
        print(
            "lmax={:4}: L2 error={:e}, speedup factor={:f}".format(
                lmax, res[0], res[1]
            )
        )

lmax_list = [80]

print("SHT tests unnorm:")
for grid in ["GLQ", "DH", "DH2"]:
    for csphase in [-1, 1]:
        for norm in ["unnorm"]:
            for extend in [True, False]:
                for lmax in [5, 10, 20, 85]:
                    res = test_SHT(lmax, grid, csphase, norm, extend)
                    print(
                        "{:3}, CS={:2}, norm={:7}, extend={:5}, lmax={:4}:  "
                        "L2 error={:e}, speedup factor={:f}".format(
                            grid, csphase, norm, extend, lmax, res[0], res[1]
                        )
                    )

lmax_list = [127, 255, 511, 1023, 2047]

print("SHT tests:")
for grid in ["GLQ", "DH", "DH2"]:
    for csphase in [-1, 1]:
        for norm in ["ortho", "4pi", "schmidt"]:
            for extend in [True, False]:
                for lmax in lmax_list:
                    res = test_SHT(lmax, grid, csphase, norm, extend)
                    print(
                        "{:3}, CS={:2}, norm={:7}, extend={:5}, lmax={:4}:  "
                        "L2 error={:e}, speedup factor={:f}".format(
                            grid, csphase, norm, extend, lmax, res[0], res[1]
                        )
                    )
print("SHTC tests:")
for grid in ["GLQ", "DH", "DH2"]:
    for csphase in [-1, 1]:
        for norm in ["ortho", "4pi", "schmidt"]:
            for extend in [True, False]:
                for lmax in lmax_list:
                    res = test_SHTC(lmax, grid, csphase, norm, extend)
                    print(
                        "{:3}, CS={:2}, norm={:7}, extend={:5}, lmax={:4}:  "
                        "L2 error={:e}, speedup factor={:f}".format(
                            grid, csphase, norm, extend, lmax, res[0], res[1]
                        )
                    )
print("SHT gradient tests:")
for grid in ["DH", "DH2"]:
    for csphase in [-1, 1]:
        for extend in [True, False]:
            for lmax in lmax_list:
                res = test_SHT_deriv(lmax, grid, csphase, extend)
                print(
                    "{:3}, CS={:2}, extend={:5}, lmax={:4}:  L2 error={:e}, "
                    "speedup factor={:f}".format(
                        grid, csphase, extend, lmax, res[0], res[1]
                    )
                )

print("DUCC: forward/backward rotation with high band limits:")
for lmax in [4095]:
    for alpha in [47]:
        res = test_rot2(lmax, alpha, 27, 59)
        print(
            "lmax={:4}: L2 error={:e}, time={:f}".format(lmax, res[0], res[1])
        )

print("DUCC: forward/backward SHT with high band limits:")
for lmax in [8191]:
    res = test_SHTducc(lmax, "GLQ", nthreads=nthreads)
    print("lmax={:4}: L2 error={:e}, time={:f}".format(lmax, res[0], res[1]))
