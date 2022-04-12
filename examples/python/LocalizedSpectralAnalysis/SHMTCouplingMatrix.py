"""
This script tests the SHMTCouplingMatrix routine
"""
import numpy as np
import pyshtools as pysh

pysh.utils.figstyle()


def main():
    test_CouplingMatrix()


def test_CouplingMatrix():
    print('\n---- testing SHMTCouplingMatrix ----')
    lmax = 30
    lwin = 10
    nwins = 1

    sqrt_taper_power = np.zeros((lwin+1, nwins))
    sqrt_taper_power[:, 0] = np.hanning(lwin+1)
    Mmt = pysh.spectralanalysis.SHMTCouplingMatrix(lmax, sqrt_taper_power)
    print(Mmt)


# ==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
