"""Test the phaseonly keyword."""

import matplotlib.pyplot as plt
import pyshtools as pysh
import numpy as np


def example():
    """Plot random phase and Gaussian random variable spectra."""
    ldata = 200
    degrees = np.arange(ldata+1, dtype=np.float64)
    degrees[0] = np.inf
    power = degrees**(-1)

    clm1 = pysh.SHCoeffs.from_random(power, exact_power=False)
    clm2 = pysh.SHCoeffs.from_random(power, exact_power=True)

    fig, ax = plt.subplots()
    ax.plot(clm1.spectrum(unit='per_l'), label='Normal distributed power')
    ax.plot(clm2.spectrum(unit='per_l'), label='Exact power')
    ax.set(xscale='log', yscale='log', xlabel='degree l',
           ylabel='power per degree l')
    ax.grid(which='both')
    ax.legend()

    plt.savefig('exact_power.png')


def main():
    example()


if __name__ == "__main__":
    main()
