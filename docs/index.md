---
title: "Spherical Harmonic Tools"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, Slepian functions, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: index.html
summary: SHTOOLS/pyshtools is an archive of Fortran 95 and Python software that can be used to perform spherical harmonic transforms, multitaper spectral analyses on the sphere, expansions of functions into Slepian bases, and standard operations on global gravitational and magnetic field data.
toc: false
---

<a href="https://twitter.com/pyshtools?ref_src=twsrc%5Etfw" class="twitter-follow-button" data-show-count="false">Follow @pyshtools</a><script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script>

## Features

SHTOOLS/pyshtools is extremely versatile:

* All standard normalizations of the spherical harmonic functions are supported: 4&pi; normalized, Schmidt semi-normalized, orthonormalized, and unnormalized.

* Both real and complex spherical harmonics are supported, and one can choose to either use or exclude the Condon-Shortley phase factor of (-1)<sup>m</sup>.

* Spherical harmonic transforms are calculated by exact quadrature rules using either the sampling theorem of *Driscoll and Healy* (1994) or Gauss-Legendre quadrature.

* The spherical harmonic transforms are fast and accurate to approximately degree 2800.

* Localized multitaper spectral analyses and expansions of functions in localized Slepian bases are easily performed.

* Standard operations on global gravitational and magnetic field data are supported.

* The Fortran routines are OpenMP compatible and OpenMP thread-safe.

* Standard data formats such as *xarray* and *netcdf* are supported.

## Installation

The Python components of SHTOOLS can be installed using the Python package manager `pip`. Binaries are pre-built for linux and macOS architectures, and you need only to execute one of the following commands in a unix terminal:

```bash
pip install pyshtools
pip install --upgrade pyshtools  # to upgrade a pre-existing installation
```

To install the Fortran 95 components for use in your Fortran programs, execute one or both of the following commands in the SHTOOLS directory

```bash
make fortran
make fortran-mp  # for OpenMP support
```

or alternatively install using the macOS package manager brew

```bash
brew tap shtools/shtools
brew install shtools
```

## Using

SHTOOLS/pyshtools can be called from any Fortran 95 or Python program. The core software is written in Fortran 95, and Python wrappers allow simple access to the fortran-compiled routines. A variety of Python notebooks and example files are included that demonstrate the major features of the library. When building from source, it will be necessary to link to LAPACK, BLAS, and FFTW compatible libraries. SHTOOLS is open source software (3-clause BSD license).

## Reference

Mark A. Wieczorek and Matthias Meschede (2018). SHTools --- Tools for working with spherical harmonics, *Geochemistry, Geophysics, Geosystems*, 19, 2574-2592, doi:[10.1029/2018GC007529](https://doi.org/10.1029/2018GC007529).
