---
title: "SHTOOLS release notes: Version 3"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: fortran_sidebar
permalink: fortran-release-notes-v3.html
summary:
toc: true
folder: fortran
---

## Version 3.4

This release adds missing functionality to the SHGrids, SHCoeffs, and SHWindow classes, and adds support for PyPI.

**Change log:**

* Add pyshtools to PyPI repository. Can now be installed using `pip install pyshtools`.
* Add new function `SHBiasKMask` which is the arbitrary window counterpart to the spherical cap window `SHBiasK`.
* Add `get_biasedpowerspectrum()` method to `SHWindows` for arbitrary windows.
* Add `copy()` method to all classes, which returns a deep copy of the instance.
* Add `__sub__`, `__add__`, `__rsub__`, `__radd__`, `__mul__`, `__div__`, `__truediv__`, and `__pow__` operators for two sets of coefficient or grid classes, or one coefficient or grid class and a scalar.
* Add `nwinrot` option when rotating spherical cap windows in `SHWindow` that will rotate only the first `nwinrot` windows.
* Add degrees option to `get_lats()` and `get_lons()` methods.
* Add the constructor `from_file()` to initialize an `SHGrid` with a numpy formatted data file. Add option to read coeffs and grids from a numpy formatted binary file. Add`tofile()` methods to output raw grid and coefficient data as either text or binary formatted files.
* Update Intro 1 notebook and add example to Intro notebook 2 showing how to use arbitrary localization windows.
* Convert notebooks to html and add links to web documentation.
* Add option `fixed_power` to `SHCoeffs.from_random()` method to generate random coefficients that fit exactly the expected power spectrum.
* Add Earth topography coefficients referenced to mean sea level to the example files, expanded to degree 300: `srtmp300.msl`.

**Citation:**

M. A. Wieczorek, M. Meschede, I. Oshchepkov, E. Sales de Andrade  (2016). SHTOOLS: Version 3.4. Zenodo. doi:[10.5281/zenodo.61180](https://doi.org/10.5281/zenodo.61180)

## Version 3.3

This is a major upgrade to SHTOOLS. Full support for Python 3 has been added, a `setup.py` file has been added for easy installation, Python notebook tutorials have been created, and full support for three major Python classes has been provided (`SHCoeffs`, `SHGrid` and `SHWindow`). One bug in the Fortran code has been fixed, as well as several minor issues with the Python wrapper functions.

**Change log:**

* Added full support for Python 3.
* Fixed a critical, but rare, bug in `MakeGridDH`, `MakeGravGridDH`, `MakeMagGridDH`, and `MakeGravGradGridDH`. In these routines, the rows of the output grid are calculated by inverse Fourier transforming a vector that depends upon the spherical harmonic coefficients. This vector, when using the complex-to-real FFTW routines, includes one element that corresponds to the coefficients with `m=lmax+1`. This element of the array was not properly initialized to zero in the Fortran code, and if the Fortran compiler did not explicitly zero all new arrays, this could have resulted in incorrect output. Given that this term is 1 index above the Nyquist frequency (`m=lmax`), if this element were not zero, each column of the grid would contain a component that oscillates from -1 to +1 (scaled by the magnitude of the element). Even when this element was not initialized, a subsequent spherical harmonic expansion of the grid would usually give correct results.
* Added a `setup.py` file for easy installation.
* Makefiles were extensively modified to simplify the Fortran and Python builds: `make all2` and `make all3` were removed and replaced by flags that make use of the precompiler to resolve underscore problems; make install now places compiled module files in /usr/local/includes.
* Makefiles were improved to minimize problems when installing both fortran and fortran-mp components. When making the latter, all object and module files in the directory `src` are first deleted.
* The namespace of pyshtools was reorganized to list the routines by submodule name. A list of all routines is given in the submodule shtools.
* Add full support for the Python classes `SHCoeffs`, `SHGrid`, and `SHWindow`.
* Added an `info()` method to shtools constants that prints an info string.
* Changed the Python wrapper so that the output arrays in the SHExpand routines correspond to the optional input variable`lmax_calc`. Modified the dimensions by 1 for the output power spectrum in the wrapper functions for `SHBias` and `SHBiasK`.
* Changed the primary input parameter of `SHMTCouplingMatrix` to be a matrix of power spectra of the localization windows instead of a matrix of spherical-harmonic coefficients of spherical-cap localization windows. Furthermore, the output dimensions of the matrix have been switched.
* Added two fortran routines for performing multitaper spectral analyses when using windows generated from a mask, `SHMultiTaperMaskSE` and `SHMultiTaperCSE`.
* Minor modifications to the Python example scripts.
* Minor bug fix to the scripts that create the unix man documentation.
* Added Python notebook tutorials.
* Created an SHTOOLS development fund, funded by bitcoin donations.

**Citation:**

Mark Wieczorek et al. (2016). SHTOOLS: Version 3.3. Zenodo. doi:[10.5281/zenodo.60010](https://doi.org/10.5281/zenodo.60010)

## Version 3.2

**Change log:**

* Added the optional argument `centralmeridian` to `Curve2Mask` that accounts for cases where the curve makes a complete circle in longitude about the planet.
* Fixed in bug in the python implementation where the outputs `error` and `corr` of `SHAdmitCorr` were switched.
* Added OpenMP support. When compiling with `make fortran-mp`,`make fortran2-mp`, or `make fortran3-mp`, saved variables in the subroutines are defined as being `threadprivate`.
* Optimized performance of the routines `SHMultiTaperSE`, `SHMultiTaperCSE`, and `SHLocalizedAdmitCorr`.
* Minor documentation fixes.

**Citation:**

Mark Wieczorek et al.. (2016). SHTOOLS: Version 3.2. Zenodo. doi:[10.5281/zenodo.55790](https://doi.org/10.5281/zenodo.55790)

## Version 3.1

This release of SHTOOLS adds improved documentation for all Fortran 95 and Python routines, fixes several bugs, adds new functionalities, and adds additional example scripts.

**Change log:**
* Added OSX installation support via `brew`.
* Added `make install` that copies files to `/usr/local`.
* Reformatted all Fortran documentation files and rewrote the Python documentation and man pages.
* Added the routines `CilmMinus` and `CilmMinusRhoH`, which are the counterparts to `CilmPlus` and `CilmPlusRhoH`.
* Removed the following routines from the fortran documentation and Python wrappers: `DhAj`, `NGLQ`, `NGLQSH`, `NGLQSHN`.
* Removed the following redundant routines from SHTOOLS: `ComputeD0`, `SHMTVarOpt0`, `SHSjkPG0`.
* Renamed the routine `PreCompute` to `SHGLQ`.
* Renamed the routine `YilmIndex` to `YilmIndexVector`.
* Renamed the routine `Hilm` to `BAtoHilm`, and `HilmRhoH` to `BAtoHilmRhoH`.
* Renamed the routine `Wl` to `DownContFilterMA`, and `WlCurv` to`DownContFilterMC`.
* Added minimal support for three python classes: `SHGrid`, `SHCoeffs`, and `SHWindow`. These will be expanded upon in the next release.
* Added the routine `SHMTCouplingMatrix`.

**Citation:**

Mark A. Wieczorek, Matthias Meschede and Ilya Oshchepkov (2015). SHTOOLS - Tools for working with spherical harmonics (v3.1), ZENODO, doi:[10.5281/zenodo.20920](https://doi.org/10.5281/zenodo.20920).

## Version 3.0

This is a major release of SHTOOLS that adds full support for Python. In addition to Python support, this release contains minor bug fixes and improved documentation.

**Change log:**
* Added full python compatibility. This includes python wrappers for all SHTOOLS functions, python documentation for all functions, and a simple test suite to ensure that the SHTOOLS functions are working correctly.
* Moved the project from Sourceforge to GitHub. The Sourceforge project will no longer be maintained.
* Fixed bugs in `SHrtoc` and `SHctor` when the optional parameter `CONVENTION` was set to its default values, causing zeros to be returned.
* Modified `ComputeDMap` and `SJHReturnTapersMap` such that the parameter `SAMPLING` is optional.
* Fixed a bug in `Curve2Mask` that could give rise to vertical lines when the longitudes of the profile points were decreasing.
* Removed the routines `Import_Wisdom_From_File` and `Export_Wisdom_From_File`.

**Citation:**

Mark A. Wieczorek and Matthias Meschede (2015). SHTOOLS - Tools for working with spherical harmonics (v3.0), ZENODO, doi:[10.5281/zenodo.15967](https://doi.org/10.5281/zenodo.15967).


