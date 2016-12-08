# SHTOOLS

SHTOOLS is an archive of Fortran 95 and Python software that can be used to perform spherical harmonic transforms and reconstructions, rotations of data expressed in spherical harmonics, and multitaper spectral analyses on the sphere.

# Features

SHTOOLS is extremely versatile:

* It can accommodate any standard normalization of the spherical harmonic functions ("geodesy" 4-pi; normalized,  Schmidt semi-normalized, orthonormalized, and unnormalized).
* Both real and complex spherical harmonics are supported.
* Spherical harmonic transforms are calculated by exact quadrature rules using either the sampling theorem of *Driscoll and Healy* (1994) where data are equally sampled (or spaced) in latitude and longitude, or Gauss-Legendre quadrature.
* One can choose to use or exclude the Condon-Shortley phase factor of (-1)^m with the associated Legendre functions.
* The spherical harmonic transforms are proven to be accurate to approximately degree 2800, corresponding to a spatial resolution of better than 4 arc minutes.
* Routines are included for performing localized multitaper spectral analyses.
* Routines are included for performing standard gravity and magnetic field calculations, such as computation of the geoid and the determination of the potential associated with finite-amplitude topography.
* The routines are fast. Spherical harmonic transforms and reconstructions take on the order of 1 second for bandwidths less than 600 and about 3 minutes for bandwidths close to 2800.

# Usage

To call the SHTOOLS routines from within a Fortran 95 program, you will need to place the command

    use SHTOOLS

immediately after the program/subroutine/function name. To compile the program, use 

    gfortran -I/usr/local/include -m64 -fPIC -O3 -L/usr/local/lib -lSHTOOLS -lfftw3 -lm -llapack -lblas

To use SHTOOLS in Python, enter the following:

    import sys
    sys.path.append('/usr/local/lib/python2.7/site-packages')
    import pyshtools as shtools

# License

The SHTOOLS software package is entirely free and open source. It can be modified and distributed according to the revised BSD license.

# See also

[SHTOOLS](shtools.oca.eu)
