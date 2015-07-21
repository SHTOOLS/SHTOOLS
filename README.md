## SHTOOLS - Tools for working with spherical harmonics ##

### What is SHTOOLS? ###
		
SHTOOLS is an archive of Fortran 95 and Python software that can be used to perform spherical harmonic transforms and reconstructions, rotations of data expressed in spherical harmonics, and multitaper spectral analyses on the sphere.

### What makes SHTOOLS different? ###

SHTOOLS is extremely versatile:

* It can accommodate any standard normalization of the spherical harmonic functions ("geodesy" 4&pi; normalized,  Schmidt semi-normalized, orthonormalized, and unnormalized).
		
* Both real and complex spherical harmonics are supported.

* Spherical harmonic transforms are calculated by exact quadrature rules using either the sampling theorem of *Driscoll and Healy* (1994) where data are equally sampled (or spaced) in latitude and longitude, or Gauss-Legendre quadrature.

* One can choose to use or exclude the Condon-Shortley phase factor of (-1)<sup>m</sup> with the associated Legendre functions.

* The spherical harmonic transforms are proven to be accurate to approximately degree 2800, corresponding to a spatial resolution of better than 4 arc minutes.

* Routines are included for performing localized multitaper spectral analyses.

* Routines are included for performing standard gravity and magnetic field calculations, such as computation of the geoid and the determination of the potential associated with finite-amplitude topography.

* The routines are fast. Spherical harmonic transforms and reconstructions take on the order of 1 second for bandwidths less than 600 and about 3 minutes for bandwidths close to 2800.
		
### How do I use SHTOOLS? ###

SHTOOLS can be used in any Fortran 95 or Python program. The base SHTOOLS software is written in Fortran 95, and Python wrappers allow simple access to all fortran-compiled routines. SHTOOLS makes use of the freely available Fourier transform package [FFTW](http://www.fftw.org) and the linear algebra packages [LAPACK](http://www.netlib.org/lapack/) and [BLAS](http://www.netlib.org/blas/).

### How do I install SHTOOLS? ###

The most recent release of SHTOOLS can be downloaded from this link on [GitHub](https://github.com/SHTOOLS/SHTOOLS/releases). Installation of SHTOOLS can be as simple as executing the following command in a unix terminal:

    make
    
To run the Fortran 95 and Python tests, enter

    make fortran-tests
    make python-tests

More information can be found on the GitHub [wiki](https://github.com/SHTOOLS/SHTOOLS/wiki) and in the SHTOOLS [documentation](www/documentation.html). Keep up to date by following SHTOOLS on [Twitter](https://twitter.com/SH_tools).

### How much does it cost? ###

The SHTOOLS software package is entirely free and open source. Two free Fortran 95 compilers exist ([gfortran](http://gcc.gnu.org/) and [g95](http://www.g95.org/)), and [Python](https://www.python.org/) is readily available. SHTOOLS can be modified and distributed according to the revised BSD license.