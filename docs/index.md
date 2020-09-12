---
title: "Spherical Harmonic Tools"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, Slepian functions, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: index.html
summary: pyshtools is an archive of Python software that can be used to perform spherical harmonic transforms, multitaper spectral analyses, expansions of functions into Slepian bases, and standard operations on global gravitational and magnetic field data.
toc: false
---

{% include note.html content="You are reading the documentation for **pyshtools**. Click [here](index-fortran.html) to access the Fortran 95 **SHTOOLS** documentation." %}

## Features

* Supports all standard normalizations and phase conventions of the spherical harmonic functions.

* Effortless conversion between real and complex harmonics, between phase conventions, and between 4&pi; normalized, Schmidt semi-normalized, orthonormalized, and unnormalized harmonics.

* Use of both regularly sampled geographic grids and grids appropriate for Gauss-Legendre quadrature.

* Spherical harmonic transforms proven to be accurate up to about degree 2800.

* Perform localized multitaper spectral analyses, or expand functions in terms of localized Slepian bases.

* Support for standard data and file formats, including *xarray* and *netcdf*.

* Import research-grade gravity, topography, and magnetic field datasets with a single command.

* Creation of publication quality maps using [Cartopy](https://scitools.org.uk/cartopy) and [pygmt](https://www.pygmt.org/).

## Easy installation

The Python components of SHTOOLS can be installed using the Python package manager `pip` or `conda`. Binaries are pre-built for Linux, macOS and Windows architectures, and you need only to execute one of the following commands in a unix terminal:

```bash
pip install pyshtools
pip install --upgrade pyshtools  # to upgrade a pre-existing installation
conda install -c conda-forge pyshtools  # Linux and macOS only
```

## Permissive licensing

SHTOOLS is open source software and is distributed under the 3-clause BSD license.

## Reference

Mark A. Wieczorek and Matthias Meschede (2018). SHTools --- Tools for working with spherical harmonics, *Geochemistry, Geophysics, Geosystems*, 19, 2574-2592, doi:[10.1029/2018GC007529](https://doi.org/10.1029/2018GC007529).
