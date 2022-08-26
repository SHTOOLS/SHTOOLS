# SHTOOLS

SHTOOLS is an archive of Fortran 95 and Python software that can be used to perform spherical harmonic transforms and reconstructions, rotations of data expressed in spherical harmonics, and multitaper spectral analyses on the sphere.

# Features

* Supports all standard normalizations and phase conventions of the spherical harmonic functions.
* Use of both regularly sampled geographic grids and grids appropriate for Gauss-Legendre quadrature.
* Spherical harmonic transforms proven to be accurate up to about degree 2800.
* Perform localized multitaper spectral analyses, or expand functions in terms of localized Slepian basis.
* Perform basic operations on global gravity and magnetic field data.
* OpenMP compatible and OpenMP thread-safe versions of the Fortran routines.

# Usage

To call the SHTOOLS routines from within a Fortran 95 program, you will need to place the command

    use SHTOOLS

immediately after the program/subroutine/function name. To compile the program, it will be necessary to link to LAPACK, BLAS, and FFTW3 compatible libraries.

# License

The SHTOOLS software package is entirely free and open source. It can be modified and distributed according to the 3-clause BSD license.

# See also

[SHTOOLS](shtools.oca.eu)
