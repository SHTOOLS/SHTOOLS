---
title: "Frequently asked questions"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: fortran_sidebar
permalink: fortran-faq.html
summary: 
toc: true
folder: fortran
---

## I have a question

Before contacting us, first

* read the documentation on this web site,
* read the rest of this FAQ,
* ensure that you are using the current version of SHTOOLS,
* uninstall SHTOOLS using `make clean` and `pip uninstall pyshtools`, and then recompile/install SHTOOLS,
* consult the [SHTOOLS issues at GitHub](https://github.com/SHTOOLS/SHTOOLS/issues)

If at this point your problem is not resolved, you can

* open an [issue on GitHub](https://github.com/SHTOOLS/SHTOOLS/issues),
* ask your question on the [SHTOOLS gitter chat forum](https://gitter.im/SHTOOLS/SHTOOLS), or
* contact us by email.

## I have a suggestion

We often make improvements based on user suggestions. Nevertheless, please realize that the developers are very busy and that we are not paid to develop or maintain this archive. We suggest that you open an issue on GitHub describing your suggestion, and we will then flag the issue as a future enhancement.

## How do you pronounce "shtools"?

It is pronounced S-H-Tools.

## I would like to contribute

Please see [this page](how-to-contribute.html).

## Can I use the SHTOOLS library with C, F77, and Matlab?

A C-interface to the Fortran-95 SHTOOLS library does exist. However, it is currently a work in progress and is not well documented. A working C example program can be found in the directory `examples/cpp`, and the required header file can be found at `includes/shtools.h`.

As for Fortran 77 and Matlab, we would be happy to add instructions on how to interface with these (or other) languages if anyone provides them to us.

## Which FFT libraries work with SHTOOLS?

SHTOOLS was developed initially to use the [FFTW](http://www.fftw.org) library. Since then, Intel's [MKL](https://software.intel.com/en-us/mkl) has added wrapper functions to their FFT routines that use the same syntax as FFTW. SHTOOLS can be linked either to FFTW or MKL with no impact on performance.

## Does SHTOOLS work with Python 2.7?

The last version of pyshtools that supported Python 2.7 was version 4.5. Precompiled binaries of this release can be installed using the command
```bash
pip install pyshtools==4.5
```
Occassionally, critical updates may be made to the Python 2.7 code in the `python2.7` branch of the GitHub project. Though Python 2.7 compatible binaries will not be distributed for any of these updates, these can be installed from souce using the command
```bash
pip install git+https://github.com/SHTOOLS/SHTOOLS.git@python2.7
```

## Will you help me with my homework?

No.

## I don't understand Fortran and Python. Will you explain how to modify the example codes?

No.

## How do I make images with the output from SHTOOLS?

If you are using the Fortran version of SHTOOLS, the output is typically in the form of an ASCII or binary raster file. These can be read by any standard graphics package, such as the free unix-based command line software [GMT](https://www.generic-mapping-tools.org/).

If you are using the Python version of SHTOOLS, the classes for working with spherical harmonic coefficients and grids contain methods for making publication quality graphics that make use of the `matplotlib` and `pygmt` packages.

## Who maintains SHTOOLS?

This software package was created initially in 2004 by [Mark Wieczorek](https://www.oca.eu/fr/mark-wieczorek) and Matthias Meschede was responsible for the initial Python implementation of version 3. A list of all contributors can be found [here](contributors.html).

## How do I cite SHTOOLS in a publication?

SHTOOLS can be cited in two ways. First, one could cite the shtools paper that was published in *Geochemistry, Geophysics, Geosystems* (the reference is on the main web page). Secondly, one could cite a specific version of the code using the [Zenodo](https://zenodo.org/) DOI (digital object identifier) that is provided in the release notes.
