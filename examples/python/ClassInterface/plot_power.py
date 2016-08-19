#!/usr/bin/env python
"""Test the different Power Spectra."""

import pyshtools


def example():
    """Plot the 2d lm power spectrum of different data."""
    # topo = np.loadtxt('../../ExampleDataFiles/topo.dat.gz')
    # grid = pyshtools.SHGrid.from_array(topo, grid='DH')
    # coeffs = grid.expand()
    path_coeffs = '../../ExampleDataFiles/MarsTopo719.shape'
    coeffs = pyshtools.SHCoeffs.from_file(path_coeffs, 719)

    coeffs.set_coeffs(0, 0, 0)
    coeffs.plot_powerperlm(show=False, loglog=True)


def main():
    """Execute the example."""
    example()


if __name__ == "__main__":
    main()
