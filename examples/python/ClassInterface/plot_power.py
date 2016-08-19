#!/usr/bin/env python
"""Test the different Power Spectra."""

import pyshtools
import numpy as np


def example():
    """Plot the 2d lm power spectrum of different data."""
    # Random model:
    # power = 1. / (1. + np.arange(100))
    # coeffs = pyshtools.SHCoeffs.from_random(power)

    # Earth topo:
    topo = np.loadtxt('../../ExampleDataFiles/topo.dat.gz')
    grid = pyshtools.SHGrid.from_array(topo, grid='DH')
    coeffs = grid.expand()

    # Mars topo:
    # path_coeffs = '../../ExampleDataFiles/MarsTopo719.shape'
    # coeffs = pyshtools.SHCoeffs.from_file(path_coeffs, 101)

    # Synthetic example:
    # coeffs = pyshtools.SHCoeffs.from_zeros(51)
    # coeffs.set_coeffs([1, 1, 1, 1], [0, 1, 2, 3], [0, 0, 0, 0])
    # coeffs = coeffs.rotate(50., 30., 0.)

    #coeffs.set_coeffs(0, 0, 0)
    # coeffs.plot_powerperlm(show=False, loglog=True)
    coeffs.plot_anisotropyspectrum(log=False)
    coeffs.plot_symmetries(lmax=40, with_grid=True)


def main():
    """Execute the example."""
    example()


if __name__ == "__main__":
    main()
