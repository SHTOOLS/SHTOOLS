---
title: "Spherical Harmonic Tools"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, Slepian functions, fortran, Python, gravity, magnetic field
sidebar: fortran_sidebar
permalink: index-fortran.html
summary: SHTOOLS is an archive of Fortran 95 software that can be used to perform spherical harmonic transforms, multitaper spectral analyses, expansions of functions into Slepian bases, and standard operations on global gravitational and magnetic field data.
toc: false
folder: fortran
---

{% include note.html content="You are reading the Fortran 95 **SHTOOLS** documentation. Click [here](index.html) to access the **pyshtools** documentation." %}

## Features

* Supports all standard normalizations and phase conventions of the spherical harmonic functions.

* Use of both regularly sampled geographic grids and grids appropriate for Gauss-Legendre quadrature.

* Spherical harmonic transforms proven to be accurate up to about degree 2800.

* Perform localized multitaper spectral analyses, or expand functions in terms of localized Slepian bases.

* Perform basic operations on global gravity and magnetic field data.

* OpenMP compatible and OpenMP thread-safe versions of the Fortran routines.

## Easy installation

To install the Fortran 95 components for use in your Fortran programs, execute one or both of the following commands in the SHTOOLS directory

```bash
make fortran
make fortran-mp  # for OpenMP support
```

Alternatively, install using the [brew](http://brew.sh/) package manager (macOS)

```bash
brew tap shtools/shtools
brew install shtools
brew install shtools --with-openmp  # to install shtools with the OpenMP components
```

or the [macports](https://www.macports.org/) package manager (macOS)
```bash
sudo port install shtools
sudo port install shtools +openmp  # to install shtools with the OpenMP components
```

## Permissive licensing

SHTOOLS is open source software and is distributed under the 3-clause BSD license.

## Reference

Mark A. Wieczorek and Matthias Meschede (2018). SHTools --- Tools for working with spherical harmonics, *Geochemistry, Geophysics, Geosystems*, 19, 2574-2592, doi:[10.1029/2018GC007529](https://doi.org/10.1029/2018GC007529).
