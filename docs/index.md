---
title: "SHTOOLS - Tools for working with spherical harmonics"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: index.html
summary: SHTOOLS is an archive of Python and Fortran 95 software that can be used to perform spherical harmonic transforms and reconstructions, rotations of data expressed in spherical harmonics, and multitaper spectral analyses on the sphere.
toc: false
---

![Logo](images/company_logo.png)

<a href="https://twitter.com/pyshtools?ref_src=twsrc%5Etfw" class="twitter-follow-button" data-show-count="false">Follow @pyshtools</a><script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script>

## Features

SHTOOLS is extremely versatile:

* All standard normalizations of the spherical harmonic functions are supported (4&pi; normalized, Schmidt semi-normalized, orthonormalized, and unnormalized).

* Both real and complex spherical harmonics are supported.

* Spherical harmonic transforms are calculated by exact quadrature rules using either the sampling theorem of *Driscoll and Healy* (1994) or Gauss-Legendre quadrature.

* One can choose to use or exclude the Condon-Shortley phase factor of (-1)<sup>m</sup> with the associated Legendre functions.

* Localized multitaper spectral analyses are easily performed.

* Routines are included for performing standard gravity and magnetic field calculations.

* The Fortran routines are OpenMP compatible and OpenMP thread-safe.

* The spherical harmonic transforms are accurate to approximately degree 2800.

* The routines are fast. Spherical harmonic transforms and reconstructions take on the order of 1 second for bandwidths close to 800 and about 30 seconds for bandwidths close to 2600.

## Installation

The Python components of SHTOOLS can be installed using the Python package manager `pip`. Binaries are pre-built for linux, macOS, and windows architectures, and you need only to execute the following command in a unix terminal:

```bash
pip install pyshtools
```

To install the Fortran 95 components for use in your Fortran programs, installation can be as simple as executing the following command in the SHTOOLS directory

```bash
make fortran
make fortran-mp  # (for OpenMP)
```

or by using the brew package manager (macOS)

```bash
brew tap shtools/shtools
brew install shtools
```

## Using

SHTOOLS can be invoked in any Fortran 95 or Python program. The core software is written in Fortran 95, and Python wrappers allow simple access to the fortran-compiled routines. A variety of Python notebooks and example files are included that demonstrate the major features of the library.

## Acknowledgments

SHTOOLS is open source software (3-clause BSD license) and makes use of the freely available packages [FFTW](http://www.fftw.org), [LAPACK](http://www.netlib.org/lapack/) and [BLAS](http://www.netlib.org/blas/).
