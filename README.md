# SHTOOLS - Tools for working with spherical harmonics #

SHTOOLS is an archive of fortran 95 based software that can be used to perform spherical harmonic transforms and reconstructions, rotations of data expressed in spherical harmonics, and multitaper spectral analyses on the sphere.
 
## What makes SHTOOLS different? ##

While several collections of code exist for working with data expressed in spherical harmonics, this one is unique for several reasons:
	
* It can accommodate any standard normalization of the spherical harmonic functions ("geodesy" 4-pi; normalized,  Schmidt semi-normalized, orthonormalized, and unnormalized).
		
* Both real and complex spherical harmonics are supported.

* Spherical harmonic transforms are calculated by exact quadrature rules using either (1) the sampling theorem of *Driscoll and Healy* (1994) where data are equally sampled (or spaced) in latitude and longitude, or (2) Gauss-Legendre quadrature. A least squares inversion routine for irregularly sampled data is included.

* One can choose to use or exclude the Condon-Shortley phase factor of (-1)^m with the associated Legendre functions.

* The spherical harmonic transforms are proven to be accurate to approximately degree 2800, corresponding to a spatial resolution of better than 4 arc minutes.

* Routines are included for performing localized multitaper spectral analyses.

* Routines are included for performing standard gravity and magnetic field calculations, such as computation of the geoid and the determination of the potential associated with finite-amplitude topography.

* The routines are fast. Spherical harmonic transforms and reconstructions take on the order of 1 second for bandwidths less than 600 and about 3 minutes for bandwidths close to 2800.
		
## How do I use SHTOOLS? ##

SHTOOLS is invoked by standard subroutine and function calls in your fortran 90/95 program. For some routines, it will be necessary to have installed the freely available Fourier transform package [FFTW](http://www.fftw.org), and the linear algebra packages [LAPACK](http://www.netlib.org/lapack/) and [BLAS](http://www.netlib.org/blas/).

## How much does it cost? ##

The SHTOOLS software package is free. Two free Fortran 95 compilers exist (gfortran as part of [GCC](http://gcc.gnu.org/) and [g95](http://www.g95.org/)). SHTOOLS can be modified and distributed according to the revised BSD license.