"""
This script tests the python class interface
"""
import numpy as np
import pyshtools as pysh

pysh.utils.figstyle()


# ==== MAIN FUNCTION ====


def main():
    example1()


def example1():
    # generate random spherical harmonics coefficients with a given
    # power spectrum and plot its bandspectrum
    degrees = np.arange(201)
    scale = 10
    power = 1. / (1. + (degrees / scale)**2)**2

    coeffs = pysh.SHCoeffs.from_random(power)
    coeffs.plot_spectrum(show=False, fname='power.png')

    # expand coefficients into a spatial grid and plot it
    grid1 = coeffs.expand(grid='DH1')
    grid1.plot(show=False, fname='DHGrid_unrotated.png')

    grid2 = coeffs.expand(grid='GLQ')
    grid2.plot(show=False, fname='GLQGrid.png')

    # rotate coefficients, expand to grid and plot again
    coeffs_rot = coeffs.rotate(40., 0., 0., degrees=True)
    grid3 = coeffs_rot.expand(grid='DH1')
    grid3.plot(show=False, fname='DHGrid_rotated.png')


# ==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
