---
title: "Frequently asked questions"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: faq.html
summary: 
toc: true
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

## I would like to contribute

Please see [this page](how-to-contribute.html).

## Can I use the SHTOOLS library with C, F77, and Matlab?

Probably, but we have not yet implemented this. If you get this to work, let us know how you did it and we will add the instructions and source files to the distribution.

## How do I cite SHTOOLS in a publication?

Each SHTOOLS release has a DOI (digital object identifier) at [Zenodo](http://zenodo.org/). The suggested citation will always be provided in the release notes of each release.

## Where can I find more information about spherical harmonics?

Two online resources are:

* [Mathworld - Spherical Harmonic](http://mathworld.wolfram.com/SphericalHarmonic.html)
* [Wikipedia - Spherical harmonics](http://en.wikipedia.org/wiki/Spherical_harmonics)

## Will you help me with my homework?

No.

## I don't understand Fortran and Python. Will you explain how to modify the example codes?

No.

## How do I make images with the output from SHTOOLS?

If you are using the Fortran version of SHTOOLS, the output is typically in the form of an ASCII raster file. These can be read by any standard graphics package, such as the free unix-based command line software [GMT](http://gmt.soest.hawaii.edu/).

If you are using the Python version of SHTOOLS, the output can be visualized by use of the `matplotlib` package.

## Who maintains SHTOOLS?

This software package was created initially in 2004 by [Mark Wieczorek](https://www.oca.eu/fr/mark-wieczorek) who is the lead developer. Matthias Meschede is responsible for the initial Python implementation. A list of all contributors can be found [here](contributors.html).
