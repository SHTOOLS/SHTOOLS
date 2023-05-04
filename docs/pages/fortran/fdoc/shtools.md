---
title: SHTOOLS (Fortran)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: fortran_sidebar
permalink: shtools.html
summary:
tags: [fortran]
toc: false
editdoc: fdoc
---

SHTOOLS is a Fortran-95 library that can be used for spherical harmonic transforms, multitaper spectral analyses, expansions of gridded data into Slepian basis functions, and standard operations on global gravitational and magnetic field data.

## Features

* Supports all standard normalizations and phase conventions of the spherical harmonic functions.
* Use of both regularly sampled geographic grids and grids appropriate for Gauss-Legendre quadrature.
* Spherical harmonic transforms proven to be accurate up to about degree 2800.
* Perform localized multitaper spectral analyses, or expand gridded data in terms of localized Slepian basis functions.
* Perform basic operations on global gravity and magnetic field data.
* OpenMP compatible and OpenMP thread-safe versions of the Fortran routines.

## Usage

To call the SHTOOLS routines from within a Fortran 95 program, you will need to place the command

    use SHTOOLS

immediately after the program/subroutine/function name. To compile the program, it will be necessary to link to LAPACK, BLAS, and FFTW3 compatible libraries.

## License

The SHTOOLS software package is entirely free and open source. It can be modified and distributed according to the 3-clause BSD license.

## See also

[SHTOOLS](shtools.github.io/SHTOOLS/)
