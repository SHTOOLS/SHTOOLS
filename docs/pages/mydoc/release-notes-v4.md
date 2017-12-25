---
title: "SHTOOLS release notes: Version 4"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: release-notes-v4.html
summary:
toc: true
---

## Version 4.1

This version adds improved functionality to SHTOOLS and fixes a couple of minor bugs. In addition, this release will be the first where pre-built wheels for unix/macOS/windows will be distributed via PYPI.

**Change log:**

* Added an optional argument `lmax` to `SHCoeffs.from_array()`.
* Coefficients are zero-padded when `lmax` is greater than the maximum degree in `SHCoeffs.to_array()`.
* The method `pad()` was added to the `SHCoeffs` class that zero pads or truncates the coefficients to a different `lmax`.
* Fixed the method `SHCoeffs.from_file()` such that the maximum spherical harmonic degree of the class is the maximum spherical harmonic degree of the coeffs (and not `lmaxin` as before).
* Fixed formatting issues with error messages in `SHCoeffs`.
* Removed print statements from the fortran code in `BAtoHilm` and `BAtoHilmRohH` that served no purpose.
* Fixed a bug in the argument order of the python wrappers of `CilmPlusRhoDH` and `BAtoHilmRhoDH`.
* Fixed the makefile to remove the `dist` directory during clean.
* Fixed a bug in the python routine `cross_spectrum()`, where the numpy arange function was incorrectly called.
* Fixed the`SHWindow` plotting methods to work when the number of rows is equal to 1.
* Conditional tests in the routine `Wigner3j` were reordered to avoid a division by zero.
* Numpy's auto-configuration is now used to detect the LAPACK libraries.
* Many minor updates to the python documentation and unix man pages.

**Citation:**

M. A. Wieczorek, M. Meschede, E. Sales de Andrade, I. Oshchepkov, B. Xu, and A. Walker (2017). SHTOOLS: Version 4.1, Zenodo, doi:[10.5281/zenodo.1067108](http://doi.org/10.5281/zenodo.1067108)


## Version 4.0
This is a major update that fixes bugs, adds new functionality, and improves Python error handing. All users are requested to upgrade to 4.0.

**Change log:**

* Instead of executing a Fortran STOP, which kills the Python kernel, the Fortran subroutines now return an `exitstatus` that allows Python to raise an exception. This technique does not work with the few Fortran functions that pyshtools calls, but these functions are relatively benign, and will soon be phased out for Python native functions.
* The Fortran `powerspectrum` routines have been removed from pyshtools, and have been replaced with Python native routines `spectrum` and `cross_spectrum`. The Python routines allow to specify the normalization, whether the output should be power, energy or l2norm, and whether the spectrum is per degree, per coefficient, or per log bandwidth.
* The method `plot_spectrum2d()` was added to the class `SHCoeffs` to plot the power as a function of degree and order.
* All pyshtools modules have been converted into proper Python subpackages. The subpackage `localizedpsectralanalysis` has been merged into `spectralanalysis`, and the subpackage `other` has been renamed `utils`.
* The Python class method `SHCoeffs.expand()` now can evaluate the function either on an SHGrid or for a list of latitude and longitude points. As part of this change, a new fortran function `MakeGridPointC` was created for complex coefficients.
* The majority of the methods for the classes `SHCoeffs`, `SHGrid` and `SHWindow` have been rename for consistency (see documentation!). Also, the classes now give the option of reading or saving to files as numpy arrays.
* Added new Python function `read_icgen_gfc` for reading ICGEM-format gravity coefficient files.
* The operator `pow` was added to the class `SHCoeffs`.
* All methods in the pyshtools classes now return copies by default, which can be modified by the optional argument `copy`.
* Added `pot` as a mandatory return argument for the Python routine `MakeGravGridDH`.
* Several minor modifications and bug fixes were made to the makefiles to improve compatibility and to allow the use of `make -j`.
* The routines `other.EigValSym`, `other.EigValVecSym`, `other.EigValVecSymTri`, `other.RandomGaussian`, `other.RandomN` and `other.PreGLQ` were removed from pyshtools, as these can be found in other scipy packages.
* The SHTOOLS routine `DHaj` was added to the pyshtools subpackage `utils`.
* Python docstrings have been streamlined and standardized.
* ...plus, many minor changes and optimizations...

**Citation:**

M. A. Wieczorek, M. Meschede, I. Oshchepkov, E. Sales de Andrade, and heroxbd  (2016). SHTOOLS: Version 4.0. Zenodo. doi:[10.5281/zenodo.206114](http://doi.org/10.5281/zenodo.206114)
