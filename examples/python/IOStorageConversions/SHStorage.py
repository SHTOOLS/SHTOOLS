"""
This script tests the conversions between real and complex spherical harmonics
coefficients
"""
import numpy as np
import pyshtools as pysh

pysh.utils.figstyle()


def main():
    test_SHStorage()


def test_SHStorage():
    # ---- input parameters ----
    lmax = 1
    mask = np.zeros((2, lmax + 1, lmax + 1), dtype=bool)
    for l in np.arange(lmax + 1):
        mask[:, l, :l + 1] = True
    mask[1, :, 0] = False
    # ---- creating random coefficients ----
    coeffs = np.random.normal(size=(2, lmax + 1, lmax + 1))
    coeffs[np.invert(mask)] = 0.

    print('\n---- testing SHCilmToCindex and SHCindexToCilm ----')
    coeffs_indexed = pysh.shio.SHCilmToCindex(coeffs)
    coeffs_recomp = pysh.shio.SHCindexToCilm(coeffs_indexed)
    print('input coeffs (l={:d}):'.format(lmax))
    print(coeffs)
    print('indexed coeffs:')
    print(coeffs_indexed)
    print('recomputed coeffs:')
    print(coeffs_recomp)

    print('\n---- testing SHCilmToVector and SHVectorToCilm ----')
    coeffs_indexed = pysh.shio.SHCilmToVector(coeffs)
    coeffs_recomp = pysh.shio.SHVectorToCilm(coeffs_indexed)
    print('\ninput coeffs (l={:d}):'.format(lmax))
    print(coeffs)
    print('\nindexed coeffs:')
    print(coeffs_indexed)
    print('\nrecomputed coeffs:')
    print(coeffs_recomp)


# ==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
