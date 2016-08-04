## SHTOOLS - Tools for working with spherical harmonics ##

[![Join the chat at https://gitter.im/SHTOOLS/SHTOOLS](https://badges.gitter.im/SHTOOLS/SHTOOLS.svg)](https://gitter.im/SHTOOLS/SHTOOLS?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

### What is SHTOOLS? ###

SHTOOLS is an archive of Fortran 95 and Python software that can be used to perform spherical harmonic transforms and reconstructions, rotations of data expressed in spherical harmonics, and multitaper spectral analyses on the sphere.

### What makes SHTOOLS different? ###

SHTOOLS is extremely versatile:

* It can accommodate any standard normalization of the spherical harmonic functions ("geodesy" 4&pi; normalized,  Schmidt semi-normalized, orthonormalized, and unnormalized).

* Both real and complex spherical harmonics are supported.

* Spherical harmonic transforms are calculated by exact quadrature rules using either the sampling theorem of *Driscoll and Healy* (1994) where data are equally sampled (or spaced) in latitude and longitude, or Gauss-Legendre quadrature.

* One can choose to use or exclude the Condon-Shortley phase factor of (-1)<sup>m</sup> with the associated Legendre functions.

* The spherical harmonic transforms are accurate to approximately degree 2800, corresponding to a spatial resolution of better than 4 arc minutes.

* The fortran routines are OpenMP compatible and OpenMP thread-safe.

* Routines are included for performing localized multitaper spectral analyses.

* Routines are included for performing standard gravity and magnetic field calculations, such as computation of the geoid and the determination of the potential associated with finite-amplitude topography.

* The routines are fast. Spherical harmonic transforms and reconstructions take on the order of 1 second for bandwidths less than 600 and about 3 minutes for bandwidths close to 2800.

### How do I install SHTOOLS? ###

The most recent release of SHTOOLS can be downloaded from [GitHub](https://github.com/SHTOOLS/SHTOOLS/releases).  Installation of SHTOOLS can be as simple as executing the following command in a unix terminal

    make
    
To compile the OpenMP version of the fortran routines, enter<

	make fortran-mp
    
The Fortran 95 and Python test suites can be run by entering

    make fortran-tests
    make python2-tests
    make python3-tests

Full installation instructions can be found [here](www/install.html), and additional information can be found on the GitHub [wiki](https://github.com/SHTOOLS/SHTOOLS/wiki). You can keep up to date by following SHTOOLS on [Twitter](https://twitter.com/SH_tools).

### How do I use SHTOOLS? ###

SHTOOLS can be invoked in any Fortran 95, Python 2, or Python 3 program. The base SHTOOLS software is written in Fortran 95, and Python wrappers allow simple access to all fortran-compiled routines. SHTOOLS is open source (revised BSD license) and makes use of the freely available Fourier transform package [FFTW](http://www.fftw.org) and the linear algebra packages [LAPACK](http://www.netlib.org/lapack/) and [BLAS](http://www.netlib.org/blas/).

### Python tutorials

* [tutorial 1: Simple Spherical Harmonics Expansions](examples/notebooks/tutorial_1.ipynb)
* [tutorial 2: Localized Spectral Analysis on the Sphere](examples/notebooks/tutorial_2.ipynb)
* [tutorial 3: The SHTOOLS Class Interface](examples/notebooks/tutorial_3.ipynb)
* [tutorial 4: Spherical Harmonics Normalization and Parseval's theorem](examples/notebooks/tutorial_4.ipynb)
* [tutorial 5: Multitaper Windows - Class Interface](examples/notebooks/tutorial_5.ipynb)
