---
title: "SHTOOLS release notes: Version 4"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: fortran_sidebar
permalink: fortran-release-notes-v4.html
summary:
toc: true
folder: fortran
---

## Version 4.12

**New Datasets!**

* Added new ultra-high degree shape models of the Moon in both principal axis and mean Earth/polar axis coordinate systems (LOLA_shape and LOLA_shape_pa) and moved MoonTopo2600p to historical.
* Added new ultra-high degree shape model of Mars (MOLA_shape) and moved MarsTopo2600 and MarsTopo719 to historical.
* Added a new ultra-high degree shape model of Mercury based on the USGS SPG DTM (USGS_SPG_shape).
* Added a new ultra-high degree shape model of Vesta based on the DLR SPG DTM (DLR_SPG_shape).
* Added a new ultra-high degree shape model of Ceres based on the DLR SPG DTM (DLR_SPG_shape) and a high-degree shape model based on the JPL SPC DTM (JPL_SPG_shape).
* Added historical degree-2 gravity models of Io, Europa, Ganymede, and Callisto derived from data collected by the Galileo mission.
* Added two shape models of Titan (Mitri2014_shape and Corlies2017_shape) and a degree 5 gravity model of Durante et al. 2019 (Durante2019_gravity).
* Added two gravity models of Enceladus (Iess2014_gravity and Park2024_gravity) and a high-degree shape model based on the JPL SPC DTM (JPL_SPG_shape).
* Added two shape models of Eros (NLR_shape and SPC_shape) and the JPL gravity model JGE15A01.
* Updated all constants to reflect the most recent datasets, and added new constant modules for Eros, Io, Europa, Ganymede, Callisto, Titan and Enceladus.
* Added the property `volume_equivalent_radius` to the constants modules for bodies that have a spherical harmonic shape model.
* Updated all datasets that download from zenodo to use the Pooch DOIDownloader.

**Bug fixes**

* When importing an xarray grid with `SHGrid.from_xarray()`, check if the first row is 90 or -90 and flip accordingly.
* Fixed a bug in `SHCoeffs.convert()` when coefficient errors are present.
* Fixed a bug in `spectralanalysis.cross_spectrum()` when using real unnormalized coefficients.

**Plotting improvements**

* Added the optional parameter `rectangle` to the `SHGrid.plotgmt()` method that allows the use of rectangular projections that are specified by the lower-left and upper-right coordinates.
* Added the optional parameter `cmap_background_foreground` to `SHGrid.plotgmt()` that controls how data are plotted when they exceed the limits of the colormap.
* Added the optional parameter `title_offset` to `SHGrid.plotgmt()` and `SHGrid.plot()` that specifies how much space to add between the plot and title.
* Cleaned up the `SHGrid.plotgmt()` shading options based on changes made in pygmt 0.7, including the use of xarrays directly instead of creating temporary files.
* Added the optional argument `cmap_rlimits` and `cmap_rlimits_complex` to `SHGrid.plot()`, `SHGrid.plotgmt()` and `SHGrid.plot3d()` to specify colorbar limits with respect to the maximum value of the data.
* Added the optional argument `cmap_scale` to `SHGrid.plot()` and `SHGrid.plotgmt()` to allow using either a linear or logarithmic color map.
* Add the optional argument `cb_power` to `SHGrid.plotgmt()` to allow plotting annotations as powers of 10.
* Added an option to `SHCoeffs.plot_spectrum()` to plot the error spectrum (if present) or not.

**Other changes**

* Added "ifx" to the list of possible Fortran compilers in the Makefile.
* Added the optional argument `r` to `SHMagCoeffs.expand()` and `SHGravCoeffs.expand()` that allows one to compute the field at an arbitrary list of (r, lat, lon) points.
* Allow the passing of `pathlib.Path` objects to all methods that require a filename.
* Added jupyter-core as a build dependency in `pyproject.toml`.
* Changed the behaviour of the fortan function `MakeCircleCoord` when the angular radius is zero.

M. A. Wieczorek, M. Meschede, A. Broquet, T. Brugere, A. Corbin, EricAtORS, A. Hattori, A. Kalinin, J. Kohler, D. Kutra, K. Leinweber, P. Lobo, J. Maia, D. Minton, I. Oshchepkov, P.-L. Phan, O. Poplawski, M. Reinecke, E. Sales de Andrade, E. Schnetter, S. Schröder, J. Sierra, A. Vasishta, A. Walker, xoviat, B. Xu (2024). SHTOOLS: Version 4.12, Zenodo, doi:[10.5281/zenodo.592762](https://doi.org/10.5281/zenodo.592762)


## Version 4.11

**Support for Python 3.12 using Meson**

This version no longer relies on `distutils` (which was deprecated in python 3.12) and instead makes use of [Meson](https://mesonbuild.com/) and [Meson-Python](https://meson-python.readthedocs.io) to build and test the pyshtools package. The package can be built from source using pip as before, however, if you need to create an editable install, it will be necessary to use the slightly modified command
```bash
pip install --no-build-isolation -e .
```
Please see the online documentation for instructions on how to run the test suites and benchmarks.

**Other changes**

* We no longer use `versioneer` to determine the package version, but instead set the version in the main `meson.build` file using `setuptools_scm`. At the present time, it is not possible to determine the version when using a source tarball, and for this case, the build will fail. Please ensure that when building from source that you are doing so from a git versioned repository.
* Fixed a problem with `SHCoeffs.volume()` when the coefficient normalization was `ortho`.
* Fixed a potential LAPACK underscore problem when compiling with `LAPACK_UNDERSCORE` specified.
* Minor changes were made to the python source files to ensure numpy v2 compatibility.

M. A. Wieczorek, M. Meschede, A. Broquet, T. Brugere, A. Corbin, EricAtORS, A. Hattori, A. Kalinin, J. Kohler, D. Kutra, K. Leinweber, P. Lobo, I. Oshchepkov, P.-L. Phan, O. Poplawski, M. Reinecke, E. Sales de Andrade, E. Schnetter, S. Schröder, J. Sierra, A. Vasishta, A. Walker, xoviat, B. Xu (2024). SHTOOLS: Version 4.11, Zenodo, doi:[10.5281/zenodo.592762](https://doi.org/10.5281/zenodo.592762)

## Version 4.10.4

**Bug fixes, minor enhancements, and deprecations**

* The module `pyshtools.shtools` has been removed, and is now accessible at `pyshtools.backends.shtools`.
* Fixed a bug in the Fortran source code of `MakeGravGradGridDH` and `MakeMagGradGridDH`, both of which are used in the pyshtools `tensor` method of the classes `SHGravCoeffs` and `SHMagCoeffs`. This bug only affected the southern hemisphere, and is most noticeable close to the south pole.
* Updated the urls for the Earth 2012/2014 datasets.
* Changed the order of the imports in `pyshtools/__init__.py` that led to a circular import problem on some systems.
* Updated how the `shtools` routines were wrapped using `functools.wraps`.

M. A. Wieczorek, M. Meschede, A. Broquet, T. Brugere, A. Corbin, EricAtORS, A. Hattori, A. Kalinin, J. Kohler, D. Kutra, K. Leinweber, P. Lobo, I. Oshchepkov, P.-L. Phan, O. Poplawski, M. Reinecke, E. Sales de Andrade, E. Schnetter, S. Schröder, J. Sierra, A. Vasishta, A. Walker, xoviat, B. Xu (2023). SHTOOLS: Version 4.10.4, Zenodo, doi:[10.5281/zenodo.592762](https://doi.org/10.5281/zenodo.592762)

## Version 4.10.3

**Minor packaging enhancement**

* The release modifies the way that the python doc strings are generated for the wrapped fortran functions. In previous releases the doc strings were generated at the time the package was built. Now, the doc strings are included directly in the repo, similar to the unix man pages. This minor change will help in making the macOS ARM conda builds, as well as in our transition from distutils to meson.

**Future deprecation**

The module `pyshtools.shtools` will be deprecated in the v4.11 release. This module represents 1 of 2 possible backends for pyshtools, and has been located at `pyshtools.backends.shtools` since version 4.9. Unless explicitly required, the user should avoid using the `backends` modules directly, and should instead call the routines that are located in the top level modules such as `pyshtools.expand` and `pyshtools.rotate`. Setting the backend by use of the routine `pyshtools.backends.selected_preferred_backend()` determines which backend to use when calling the routines in these top level modules.

M. A. Wieczorek, M. Meschede, T. Brugere, A. Corbin, A. Hattori, A. Kalinin, J. Kohler, D. Kutra, K. Leinweber, P. Lobo, I. Oshchepkov, P.-L. Phan, O. Poplawski, M. Reinecke, E. Sales de Andrade, E. Schnetter, S. Schröder, J. Sierra, A. Vasishta, A. Walker, xoviat, B. Xu (2023). SHTOOLS: Version 4.10.3, Zenodo, doi:[10.5281/zenodo.592762](https://doi.org/10.5281/zenodo.592762)

## Version 4.10.2

**Bug fixes and minor enhancements**

* Add the optional parameter `lmax` to `pysh.SHWindow.to_shgrid()`.
* Replace `libtool` by `ar` and `ranlib` in the fortran makefile.
* Replace `PWD` by `CURDIR` in the main makefile.
* Add method `.change_units` to the `SHMagCoeffs` class.
* Copy all routines from the `shtools` module to `backends.shtools` (the top-level shtools module will be deprecated in v4.11).
* Use `pkg_resources` instead of `setuptools.version.pkg_resources` in `setup.py`.
* Fix a couple bugs where `_np.int_` was mistakenly `_np.int`.
* Move unix man pages from section 1 to section 3.
* Modify fortran test programs to accept command line arguments, such as the location of the example data files and program input files.
* Remove `pypandoc` as a dependency in `setup.py` and don't convert the readme to reST format for pypi.
* Add a Ganymede gravity model to the `datasets` module.
* Minor updates to documentation and refactoring of the project `README`.

**Future deprecation**

The module `pyshtools.shtools` will be deprecated in the v4.11 release. This module represents 1 of 2 possible backends for pyshtools, and has been located at `pyshtools.backends.shtools` since version 4.9. Unless explicitly required, the user should avoid using the `backends` modules directly, and should instead call the routines that are located in the top level modules such as `pyshtools.expand` and `pyshtools.rotate`. Setting the backend by use of the routine `pyshtools.backends.selected_preferred_backend()` determines which backend to use when calling the routines in these top level modules.

M. A. Wieczorek, M. Meschede, T. Brugere, A. Corbin, A. Hattori, A. Kalinin, J. Kohler, D. Kutra, K. Leinweber, P. Lobo, I. Oshchepkov, P.-L. Phan, O. Poplawski, M. Reinecke, E. Sales de Andrade, E. Schnetter, S. Schröder, J. Sierra, A. Vasishta, A. Walker, xoviat, B. Xu (2023). SHTOOLS: Version 4.10.2, Zenodo, doi:[10.5281/zenodo.592762](https://doi.org/10.5281/zenodo.592762)

## Version 4.10.1

**Bug fixes and minor enhancements**

* Simplify the code for backend management and improve handing of the default backend when `ducc0` isn't installed.
* Changing the backend now changes which functions are referenced in the `expand` and `rotate` submodules.
* Fix `setup.py` to work with all versions of `setuptools`.
* Add the MarsTopo719 dataset for use in CI checks (using MarsTopo2600 often would timeout during download).
* Fix a bug in `SHGravRealCoeffs.expand` and `SHGravRealCoeffs.expand` that did not correctly compute the radius of the flattened ellipsoid when an array of latitudes was provided.
* Remove `-static` option from compiler options.
* Convert some strings to raw format when they contain latex backslashes.
* And other minor changes...

**Future deprecation**

The module `pyshtools.shtools` will be deprecated in the v4.11 release. This module represents 1 of 2 possible backends for pyshtools, and has been located at `pyshtools.backends.shtools` since version 4.9. Unless explicitly required, the user should avoid using the `backends` modules directly, and should instead call the routines that are located in the top level modules such as `pyshtools.expand` and `pyshtools.rotate`. Setting the backend by use of the routine `pyshtools.backends.selected_preferred_backend()` determines which backend to use when calling the routines in these top level modules.

M. A. Wieczorek, M. Meschede, T. Brugere, A. Corbin, A. Hattori, K. Leinweber, I. Oshchepkov, M. Reinecke, E. Sales de Andrade, E. Schnetter, S. Schröder, A. Vasishta, A. Walker, B. Xu, J. Sierra (2022). SHTOOLS: Version 4.10.1, Zenodo, doi:[10.5281/zenodo.592762](https://doi.org/10.5281/zenodo.592762)

## Version 4.10

**Enhancements**

* Change the preferred backend from 'shtools' to 'ducc' (when both are available).
* Link routines in top level modules `pyshtools.expand` and `pyshtools.rotate` to the corresponding backend routines.
* Add historical lunar topography dataset GLTM-2B.
* Add historical martian magnetic field models FSU50 and FSU90.
* Add new Mars gravity model MRO120F as well as several historical Mars gravity models.
* Add historical Venus topography datasets SHTJV360A01 and SHTJV360A02.
* Add Thebault2021 Earth magnetic field dataset.
* Add Mars topography dataset MarsTopo719, which is a truncated version of MarsTopo2600.
* Update urls for databases hosted at GSFC.
* Reorder optional arguments in docs for `makegravgradgriddh` and `makemaggravgradgrid` for consistency with code.
* Allow shtools and dov file formats to contain floats for degree and order.
* Minor changes and enhancements to the documentation.

**Bug fixes**

* Fix typo regarding `nthreads` in SHMagCoeffs.rotate() method.
* Fix bug with `SHGravCoeffs.admittance()` when using `function=geoid`.
* Fix bug in python wrapper of the routine `MakeGrid2D` concerning the mandatory variable `interval`.
* Add workaround to use pygmt with shading for versions >=0.4.
* Convert all grids to `float` before using the `ducc0` backend.
* `SHGeoid.to_netcdf()` now outputs double precision by default (consistent with the other grid classes).
* Fix bug with `SHWindow.multitaper_cross_spectrum()` when using arbitrary localization regions.
* Fix bug with the c-wrapper for `cMakeGradientDH` regarding the optional `radius` parameters.
* Minor changes to remove deprecation warnings.

**Future deprecation**

The module `pyshtools.shtools` will be deprecated in the v4.11 release. This module represents 1 of 2 possible backends for pyshtools, and has been located at `pyshtools.backends.shtools` since version 4.9. Unless explicitly required, the user should avoid using the `backends` modules directly, and should instead call the routines that are located in the top level modules such as `pyshtools.expand` and `pyshtools.rotate`. Setting the backend by use of the routine `pyshtools.backends.selected_preferred_backend()` determines which backend to use when calling the routines in these top level modules.

M. A. Wieczorek, M. Meschede, T. Brugere, A. Corbin, A. Hattori, K. Leinweber, I. Oshchepkov, M. Reinecke, E. Sales de Andrade, E. Schnetter, S. Schröder, A. Vasishta, A. Walker, B. Xu, J. Sierra (2022). SHTOOLS: Version 4.10, Zenodo, doi:[10.5281/zenodo.592762](https://doi.org/10.5281/zenodo.592762)

## Version 4.9

**Backends**

Implemented the option to use a different backend when performing certain operations requiring spherical harmonic transforms. At present, only 'shtools' (default) and the [Distinctly Useful Code Collection ('ducc')](https://gitlab.mpcdf.mpg.de/mtr/ducc) are supported.

* Introduced a new module `backends` that has functions allowing the user to control which backend is used. To set the backend for all subsequent operations, use `backends.select_preferred_backend()`.
* Added the optional parameters `backend` and `nthreads` ('ducc' only) to all methods of the pyshtools classes that allow multiple backends (such as `SHGrid` and `SHCoeffs`).
* Added a new `backends` web documentation page that describes the use of the new module.

**Plotting routines**

* Added the new methods `SHGrid.histogram()` and `SHGrid.plot_histogram()` for generating area-weighted histograms.
* Renamed a few instances of variable names used in the plotting routines from `title_labelsize` to `titlesize` for consistency. This affects the methods `Slepian.plot_coupling_matrix()`, `Slepian.plot_spectra()`, and `SHWindow.plot_coupling_matrix()`.
* Modified some plotting routines so that fontsizes can be specified using standard matplotlib strings.
* Replaced the deprecated matplotlib `get_geometry` with `get_subplotspec`.
* Replaced the pygmt option `I` with the alias `shading` for colorbars.

**IO and datasets**

* Updated the Venus rotation period using data from Margot et al. (2021).
* Renamed `constants` variables `r` to `mean_radius` (but kept `r` as an alias for backwards compatibility).
* Improved the reading of ICGEM formatted files of gravitational potential models. Fixed a bug where the time variable contribution was not computed when the coefficients had an associated `trnd` or `dot` term. When a line in the data section starts with an unknown key, a warning is now printed to the screen (which can be turned off by specifying `quiet=True`). For files formatted as `icgem2.0`, time variable terms were simply ignored if the specified epoch fell outside of the allowed range. Now, the routine will instead raise an error. Finally, the documentation was improved by describing the allowable keyword entries of the header and data section of the file.
* Added the option `encoding` for all routines and methods that read or write text-based spherical harmonic files.
* Hard coded all datasets to use `utf-8` in order to avoid problems with the XGM2019E dataset that has a character that can not be decoded by the GBK encoding that is the default in some Chinese installations.
* Added a few historical lunar gravity fields to the module `datasets.Moon.historical.gravity`.

**Other changes**

* Add a C wrapper for the function `MakeGradientDH()`.
* Changed the default behavior of the Fortran routine `MakeGradientDH()` and `SHCoeffs.gradient()`. The original behavior was to compute the gradient on a sphere of radius `r`, where `r` was the degree 0 coefficient of the function. The new behaviour is to compute the gradient on the unit sphere. This radius can be modified by supplying the optional argument `radius`.
* Added the option 'per_lm' for generating random spherical harmonic coefficients in the `.from_random()` methods of `SHCoeffs`, `SHGravCoeffs`, and `SHMagCoeffs`.
* Fixed two bugs related to complex spherical harmonic transforms. First, for complex grids, the last coefficient `coeffs[1, lmax, lmax]` was in error when `lmax` was odd and when using `DH` grids. Second, when using `SHGrid.expand()` with grid type `DH2` the parameter `sampling=2` was not passed to the Fortran routine.
* Updated `makefile install` to include example data files, and to place them in the correct directories when installing with homebrew.
* Updated a few dependencies, including `astropy>=4.0` and `pygmt>=0.3.0`.
* Fixed a floating point error-caused bug in `SHGrid` that could arise if the value input to arccos was greater than 1.
* Converted `np.float_` and `np.complex_` to `np.float64` and `np.complex128` to avoid numpy deprectation warning.
* Added `threadsafe` to numpy signature files.

M. A. Wieczorek, M. Meschede, T. Brugere, A. Corbin, A. Hattori, K. Leinweber, I. Oshchepkov, M. Reinecke, E. Sales de Andrade, E. Schnetter, S. Schröder, A. Vasishta, A. Walker, B. Xu (2021). SHTOOLS: Version 4.9, Zenodo, doi:[10.5281/zenodo.592762](https://doi.org/10.5281/zenodo.5254984)

## Version 4.8

* Several functions have been vectorized using `numpy.vectorize()`. These include: `spharm_lm()`, `legendre_lm()`, `MakeGridPoint()`, `MakeGridPointC()`, `DownContFilterMA()`, `DownContFilterMC()`, `NormalGravity()`, `SHConfidence()`, and `PlmIndex()`.
* A new Fortran routine `MakeGradientDH` was added to compute the horizontal gradient of a real scalar function. The method `.gradient()` was added to the `SHCoeffs` class, and a new class `SHGradient` was created to store and plot the two horizontal components of the gradient.
* Added new Fortran functions `MakeGravGridPoint` and `MakeMagGridPoint` to compute the gravity and magnetic field vector at a single point.
* Added the option to compute the gravity and magnetic field vectors at a single point using the python class methods `SHGravCoeffs.expand()` and `SHMagCoeffs.expand()`.
* The `plot_spectrum2d()` routines have been updated to include more plotting options, including placement of the origin, tick intervals, and colormaps. Most optional parameters are the same as in the `SHGrid.plot()` method.
* Added the option to including intensity shading in the `SHGrid.plotgmt()` routine. The shading can be derived from the gradient of the input grid (by setting `shading=True`) or from a different map by supplying an `SHGrid` class instance. Optional parameters include the azimuth of the shading (`shading_azimuth`), as well as the maximum amplitude of the intensity (`shading_amplitude`).
* Modified all the Fortran routines to use a slightly more efficient way to compute the radius of an ellipsoid as a function of geocentric latitude.
* Fixed a bug in `SHCoeffs.expand()` when `colat` was specified in radians.
* All declarations of integers in the Fortran code are now made using the types defined in the module `iso_fortran_env`. Furthermore, the python wrapper and signature files have been updated to be explicit when defining the Fortran variables.
* Fixed a bug where the old module name `constant` needed to be updated to `constants` in the method `SHCoeffs.centroid()`.
* Corrected the parameterization used when generating ellipsoids in `SHGrid.from_ellipsoid()`. Though this method was introduced in v4.7, it was not mentioned in the release notes.
* Changed the default behavior of `SHCoeffs.to_array()` so that the default value is not to return the errors by setting `errors=False`.
* Added the  optional attribute `name` to the coefficient classes `SHGrid`, `SHGravCoeffs`, `SHMagCoeffs` and `SlepianCoeffs`. All datasets now explicity set `name` to the function call of the dataset.
* Moved the file `shtools.h` from `src/` to `include/` and updated the Makefiles accordingly.

M. A. Wieczorek, M. Meschede, E. Sales de Andrade, A. Corbin, I. Oshchepkov, B. Xu, and A. Walker, A. Hattori, S. Schröder, K. Leinweber, A. Vasishta (2021). SHTOOLS: Version 4.8, Zenodo, doi:[10.5281/zenodo.592762](https://doi.org/10.5281/zenodo.4467730)

## Version 4.7.1

* Minor modifications were made to the Makefiles in order to submit shtools to the homebrew-core and macports package managers. Relative paths were removed in a few cases by explicitly passing variables such as `$(MODPATH)$` to all dependent sub-makefiles. Default variables are no longer set in the sub-makefiles, as these are not intended to be used independently: All variables are passed directly from the main Makefile. Renamed the directory `modules` to `include` to be consistent with macports and homebrew installations. The `F95FLAGS` are set by searching if the compiler name contains the "short" compiler name. This allows recognizing "gfortran-10" as being "gfortran".
* Added a `.github` folder with templates for issues and releases checklists.
* Converted matplotlib relative font sizes (such as 'large') to points when passing font sizes to the Cartopy and pygmt plotting routines.
* Minor changes to the Travis configuration file, the conda `environment.yml` file, and fortran documentation and man pages.
* Added initial experimental support for C-binded SHTOOLS wrapper functions. This includes replacing assumed-size arrays with fixed-size arrays with additional arguments for each dimension. Though this is not yet documented, a working example can be found in the folder `examples/cpp`.

M. A. Wieczorek, M. Meschede, E. Sales de Andrade, A. Corbin, I. Oshchepkov, B. Xu, and A. Walker, A. Hattori, S. Schröder, K. Leinweber, A. Vasishta (2020). SHTOOLS: Version 4.7.1, Zenodo, doi:[10.5281/zenodo.592762](https://doi.org/10.5281/zenodo.4048072)

## Version 4.7

**Datasets**

The new datasets module allows users to easily download spherical harmonic coefficient datasets and return them as `SHCoeffs`, `SHGravCoeffs` or `SHMagCoeffs` class instances.

To load a dataset, call the relevant method as in these examples:
```
    hlm = pysh.datasets.Venus.VenusTopo719()  # Venus shape
    clm = pysh.datasets.Earth.EGM2008()  # Earth gravity
    glm = pysh.datasetes.Earth.WDMAM2_800()  # Earth magnetic field
    clm = pysh.datasets.Moon.GRGM1200B()  # Gravity of the Moon
```

**Better IO routines**

* Added the functions (in the module `shio`) `shwrite()`,`read_dov()` `write_dov()`, `read_bshc()`, `write_bshc()` and `write_igcem_gfc()` to read and write 'shtools', 'dov', 'bshc', and 'icgem' files.
* Added the function `shio.read_igrf()` for reading IGRF formatted files, and returning coefficients for a specified year.
* The `SHCoeffs`, `SHMagCoeffs` and `SHGravCoeffs` methods `to_file()` and `from_file()` now accept all file formats.
* Added support for reading gzip and zip files in `shread`, `SHCoeffs.from_file()`, `SHGravCoeffs.from_file()`, `SHMagCoeffs.from_file()`, and `read_icgem_gfc()`
* Fixed a minor bug where netcdf files would not accept boolean attributes.

**Amittance and correlation methods**

* Added the methods `admittance()`, `correlation()` and `admitcorr()` for the classes `SHCoeffs`, `SHGravCoeffs`, and `SHMagCoeffs` to compute the admittance and/or correlation with another function.
* Added the methods `plot_admittance()`, `plot_correlation()` and `plot_admitcorr()` to easily plot these functions.

**Better plotting routines**

* Added the option `legend_loc` to most plotting routines to allow fine control over where the legend is placed.
* Minor bug fixes concerning colorbar parameters `cb_offset` and `cb_triangles`.

**Better treatmentment of uncertainties**

* Added the option to include error coefficients in the class `SHCoeffs`.
* Added the boolean option `errors` to the method `to_array()` in order to control whether the error coefficients are returned with the function spherical harmonic coefficients.
* Added the option `legend_error` to the `SHCoeffs`, `SHMagCoeffs` and `SHGravCoeffs` method `plot_spectrum()` to provide a customized legend entry for the error spectrum.

**New attributes for `SHCoeffs`, `SHGravCoeffs` and `SHMagCoeffs`**

* Added the attribute `error_kind` to specify the type of errors.
* Added the attribute `units` to all grid and coefficient classes.
* Added the attribute `epoch` to `SHGravCoeffs`, `SHGravGrid` , `SHGeoid` and `SHTensor`.
* Added the attribute `year` to `SHMagCoeffs`, `SHMagGrid` , and `SHTensor`.

**Improved Documentation**

* The web documentation has been broken into two separate components: pyshtools (python) and SHTOOLS (Fortran 95).
* Reorganized the web documentation for clarity (re-organization of tutorials and guides, creation of a separate page for shtools grid formats, creation of separate pages for datasets, constants, and spherical harmonic coefficient file coeeficients).
* The python tutorial notebooks are now rendered by the jupyter nbviewer web page. From this viewer, the user can easily download the notebook, or run it in a binder session.
* Updated the documentation for installing pyshtools using Conda.

**Initial support for fpm**

Initial experimental support is added for use with fpm: the [fortran package manager](https://github.com/fortran-lang/fpm).

To install as a stand-alone project, it is only necesssary to use the command
```
fpm build
```
This will place the necessary .mod and .a files in a subdirectory of `build`.

To include shtools as a dependency in a project that compiles with fpm, you only need to add the following to the `fpm.toml` file:
```
[dependencies]
SHTOOLS = {git="https://github.com/SHTOOLS/SHTOOLS.git"}
```

In the current state of fpm (which is undergoing active development), it is not possible to link to system wide libraries, such as fftw and lapack, which are required by shtools.

**Other changes**

* Added error checks to the pyshtools function `YilmIndexVector`.
* Renamed the `constant` module to `constants`, and reogranized the constants in a more logical way (i.e., `constants.Mars.r` instead of `constants.Mars.r_mars`).
* Added a gmt xarray accessor for use with pygmt.
* Fixed a bug in `Curve2Mask` python wrapper when using extended grids, and fixed a bug in the fortran code when the input file contained points at exactly 0 or 360 degree.
* pyshtools versioning is now done using `versioneer`, instead of the homemade system that was in the setup.py (which was somewhat complicated and needed to set `ISREALESED` to True or False). Versioneer gets the version number automatically from git tags.

M. A. Wieczorek, M. Meschede, E. Sales de Andrade, I. Oshchepkov, B. Xu, and A. Walker, A. Hattori, S. Schröder, K. Leinweber, A. Vasishta (2020). SHTOOLS: Version 4.7, Zenodo, doi:[10.5281/zenodo.592762](https://doi.org/10.5281/zenodo.4028484)

## Version 4.6

**New extended grids**

All grid formats now allow to compute the redundant values at 360 E longitude (GLQ and DH), as well as at 90 S (DH only). These *extended* grids are now the default in pyshtools, but remain optional in the Fortran 95 routines. The use of extended grids is controlled by the optional argument `extend`. The purpose of these extended grids is to better integrate with the plotting routines that require these points (i.e.., Cartopy and pygmt).

**Improved plotting and map projections**

The plotting routine `SHGrid.plot()` has been refactored to allow support for projections using `Cartopy` and `pygmt`.

* An incorrect 0.5 pixel offset was fixed when plotting grids via matplotlib, and grids now correctly plot both 0 and 360 degrees using the new "extended" grids of `SHGrid`.
* Support was added for Cartopy projections, by specifying: `SHGrid.plot(projection=ccrs.ProjectionName())`.
* The argument `colorbar` now takes the options 'top', 'bottom', 'left' or 'right'.
* Improved plotting and placement of colorbars. New optional arguments include `cb_label` for labels, `cb_ylabel` for a label on the y axis of the colorbar, `cb_tick_interval` for specifying the major tick interval, `cb_minor_tick_interval` for specifying minor tick intervals, `cb_triangles` for plotting upper/lower limit triangles at the ends of the colorbar, `cb_width` to specify the colorbar width, and `cb_offset` to override the default spacing between the map and colorbar.
* Improved colormap handling: New optional arguments include `cmap_limits` to specify the lower and upper bounds of the data, as well as an interval for constant color intervals, and `cmap_reverse` to reverse the colormap.
* Improved handling of ticks and annotations: The optional argument `ticks` specifies which ticks and annotations to show, using a syntax from the generic mapping tools (i.e., `'WSen'`).
* Experimental support for pygmt using the routine `SHGrid.plotgmt()`. This function takes nearly the same arguments as `plot()`. As soon as pygmt implements projection classes (https://github.com/GenericMappingTools/pygmt/pull/379), this will be incorporated into the `plot` function in the same manner as Cartopy was.
* All gravity, magnetics, tensor, localization windows and slepian function plotting routines incorporate these changes.
* Added a new introductory notebook that shows how to use all features of the `plot()` function.

**Improved integration with xarray DataArrays, xarray DataSets, and netcdf files**

* Added the methods `to_netcdf()` and `from_netcdf()` to the `SHCoeffs`, `SHGravCoeffs` and `SHMagCoeffs` classes.
* Added the method `SHGrid.from_xarray()` to initialize a grid from an xarray DataArray.
* Added improved descriptive attributes for netcdf files that mirror these [conventions](http://cfconventions.org/cf-conventions/cf-conventions.html).
* Added the method `SHGeoid.to_xarray()` to export an xarray DataArray and `to_netcdf()` to export a netcdf object readable by the generic mapping tools.
* Added the methods `SHGravGrid.to_xarray()` and `SHMagGrid.to_xarray()` to export all gridded data (radial, theta, phi, total, and potential) as an xarray DataSet.
* Added the methods `SHGravTensor.to_xarray()` and `SHMagTensor.to_xarray()`to export all gridded data (Vxx, invariants, eigenvalues) as an xarray DataSet.

**Gravity routine improvements**

* Added the method `SHGravCoeffs.center_of_mass` to calculate the center of mass of a body.
* Added the method `SHGravCoeffs.inertia_tensor()` to calculate the moment of inertia tensor.
* Added the Earth dynamical flattening constant H (IERS Conventions 2010) to the `constant` module.
* The `read_icgem_gfc()` function was extended with the option `encoding` as some models in ICGEM are not in UTF-8.
* Addded the method `centroid()` to the class `SHCoeffs`. The centroid is computed as the center of mass of a homogeneous object.

**Other changes**

* New methods `SHGrid.to_real()` and `SHGrid.to_imag()` return the real and imaginary components of a complex `SHGrid` instance.
* Added an optional argument `copy` to `SHCoeffs.pad()`.
* Fixed bugs in the Fortran code of `PlBar_d1` and `PlON_d1` when calculating the Legendre polynomials at the north and south pole.
* Spherical harmonic coefficients can be read remotely by specifying a URL as the filename. This functionality uses `requests.get()`, and has been implemented in the function `shread()` and the `SHCoeffs` method `from_file()`.
* Fixed a bug in the fortran code of `Curve2Mask`. As part of this fix, the optional parameter `centralmeridian` has been removed as it is no longer required. The longitudes of the curve can possess values from -360 to 720 degrees, and the routine searches for discontinuities that may occur between two successive points as the longitudes pass from 360 to 0, or -180 to 180 degrees.

M. A. Wieczorek, M. Meschede, E. Sales de Andrade, I. Oshchepkov, B. Xu, and A. Walker, A. Hattori, S. Schröder, K. Leinweber, A. Vasishta (2020). SHTOOLS: Version 4.6, Zenodo, doi:[10.5281/zenodo.3698050](https://doi.org/10.5281/zenodo.3698050)

## Version 4.5

**SHCoeffs**
* Added `cross_spectrum()`, `plot_cross_spectrum()` and `plot_cross_spectrum2d()` methods to the `SHCoeffs` class.
* Added `from_cap()` constructor to create coefficients of a spherical cap.
* Fixed a small bug when rotating coefficients in the `SHCoeffs`, `SHMagCoeffs` and `SHGravCoeffs` classes when the input `csphase` is different from the default value.

**SHGrid**
* Added the method `from_zeros()` to initialize a grid with zeros.
* Added the method `from_cap()` to initial a grid with a spherical cap.
* Added support for saving gridded data to netcdf and xarray formats. `to_netcdf()` exports data to netcdf format, and when saved to file can be used directly with GMT (generic mapping tools) where they are known as 'grd' files. `to_xarray()` exports data to an xarray DataArray.
* Fixed a type error when expanding complex coefficients at arbitrary points.

**Slepian functions**
* Added the fortran function `SHSCouplingMatrix` and `Slepain` class method `coupling matrix()` for computing the coupling matrix that relates the Slepian expansion power spectrum to the global power spectrum.
* Added the fortran function `SHSCouplingMatrixCap` which is optimized for working with spherical cap Slepian functions.
* Improved the plotting capabilities of the method `plot_coupling_matrix`, such as the addition of colorbars, and the option to normalize to maximum value to unity.
* Added the option `taper_degrees()` and `slepian_degrees()` to the `SHWindow` and `Slepian` classes, respectively, to allow for the construction of Slepian functions that exclude certain spherical harmonic degrees.
* Added the fortran function `SHMTVar` and the `SHWindow` class method `variance()` to compute the variance of a multitaper spectral estimate (based on `SHMTVarOpt`).
* Added the fortran function `SHSlepianVar` and the `Slepian` class method `variance()` to compute the variance of a Slepian expansion spectral estimate.
* Improved the handing of the optional parameter `weights` in all methods of `SHWindow`.
* **Changed the name of the optional parameter `nwin` of `SHWindow.coupling_matrix` to `k` for consistency with the localized spectral analysis routines.**

**Fortran 95**
* The Fortran code was modified to be strictly compliant with the f95 standard (`-std=f95` in the `gfortran` compiler). Double precision, double complex, and long integers are defined as `real(dp)`, `complex(dp)`, and `integer(int8)`, and the types are defined in a new module `ftypes.f95`.
* A few intrinsic function calls have been renamed to conform with the standards, and all double precision constants are now defined by appending `_dp` to them.

**FFTW**
* The syntax of the `fftw` routines has been updated. In particular, the old `call dfftw_execute(plan)` statements now include all their dependent variables using the new syntax `call fftw_execute_dft(plan, grid, coef)`. Importantly, the old syntax caused the GCC9 optimizer to break the spherical transform routines, generating meaningless output for large parts of the grids or coefficients.
* The fortran routines now access the FFTW library using Fortran 2003 standards. Though most compilers do not yet support all features of F2003, the few features used here (such as the module `iso_c_binding`) are claimed to be supported by the majority of modern compilers. This allows us to access the FFTW routines without use of fortran bindings, which are not always included in compiled versions of FFTW.
* All FFTW routines are explicitly defined in an interface block in `FFTW3.f95`. For simplicity, the use of the optional parameter `FFTW_UNDERSCORE` used in the Makefile has been deprecated.

M. A. Wieczorek, M. Meschede, E. Sales de Andrade, I. Oshchepkov, B. Xu, and A. Walker (2019). SHTOOLS: Version 4.5, Zenodo, doi:[10.5281/zenodo.2350781](https://doi.org/10.5281/zenodo.2350781)

## Version 4.4

**New Slepian expansion routines**

* Added the new Fortran 95 functions `SHRotateTapers` for rotating spherical cap Slepian functions, `SlepianCoeffs` for expanding a function in a Slepian basis, and `SlepianCoeffsToSH` to convert the Slepian coefficients into spherical harmonic coefficients.
* Added two python classes `Slepian` and `SlepianCoeffs` for managing Slepian basis functions and the Slepian expansion coefficients of a function.

**Legendre and spherical harmonic convenience functions**

* Added 4 python convenience functions for computing Legendre and spherical harmonic functions: `legendre()` to compute all the legendre functions, `legendre_lm()` to compute the legendre function for a specified degree and order, `spharm()` to compute all the spherical harmonic functions, and `spharm_lm()` to compute the spherical harmonic functions for a specific degree and order.

**Other improvements**

* Added the option to perform a weighted least squares inversion in the fortran routine `SHExpandLSQ`, and added a python wrapped function `SHExpandWLSQ()`.
* Added `min()` and `max()` methods to the `SHGrid` class, to return the minimum and maximum value of the gridded data.
* Added the property `mass` to the `SHGravCoeffs` class, which is computed by the input GM and the codata value for `G`.
* Add the optional parameter `omega` to `SHGravCoeffs.geoid()` to override the value provided in the class instance.
* Fixed some minor typos and usability issues, and cleaned up the python wrapper and signature files.

M. A. Wieczorek, M. Meschede, E. Sales de Andrade, I. Oshchepkov, B. Xu, and A. Walker (2018). SHTOOLS: Version 4.4, Zenodo, doi:[10.5281/zenodo.592762](https://doi.org/10.5281/zenodo.592762)

## Version 4.3

**New Gravity and Magnetic field classes**

* Added gravity classes `SHGravCoeffs`, `SHGravGrid`, `SHGravTensor` and `SHGeoid`.
* Added magnetic field classes `SHMagCoeffs`, `SHMagGrid`, and `SHMagTensor`.
* Added new fortran subroutine `MakeMagGradGridDH`, which is analogous to `MakeGravGradGridDH`.
* Reorded the arguments of `CilmPlusRhoHDH` to be consistent with `CilmPlusDH`.
* The python routine `MakeMagGridDH` now also outputs the magnetic potential as a grid.

**Better figures**

* Addition of the function `pyshtools.utils.figstyle()`, that sets several matplotlib parameters for better figures. This function takes as optional parmeters the maximum useable width of a journal page, the relative width of the figure with respect to this value, and the screen resolution in dpi.
* Most plotting routines have optional parameters to set minor tick intervals, grids, label font size, and tick font size.
* Degree symbols are plotted on tick labels for maps.
* `examples/python/Common/FigStyle.py` was removed from the examples.
* Added the options `vmin` and `vmax` to the plotting methods `SHCoeffs.plot_spectrum2d()` in order to specify the limits of the color scale.
* Added the option to plot colorbars on `SHGrid` plots, along with the option to specify their orientation and a text label.
* All notebooks have been updated.

**New constant subpackage**

* The `constant` subpackage has been completely rewritten and now makes use of the [astropy](https://docs.astropy.org/en/latest/constants/index.html) `Constant` class. This class has attributes `name`, `value`, `uncertainty`, `unit`, and `reference`. The naming of the constants has changed in some cases for consistency. A few constants that are not necessary were removed. Many of the constants were updated with more recent values. `Constants` can be used in arithmetic operations with either other `Constants` or with objects of the astropy class `Quantity`.
* Constants are organized into modules for each of the planets (`Mercury`, `Venus`, `Earth`, `Moon`, and `Mars`), and for convenience, these are all added to the main namespace. The fundamental constants `G` and `mu0` from the astropy constants package were added (as taken from CODATA 2014).

**Other changes**

* Fixed a bug in how the random coeffcients were determined for unnormalized coeffcients in SHCoeffs.
* Optional parameter `seed` added to `SHCoeffs.from_random()` to allow for reproducibility.
* One can now specify colat instead of lat for the method SHCoeffs.expand().
* Added `__repr__` methods to all pyshtools classes.
* Changed the mathematical operators of `SHCoeffs` such that addition and subtraction of a constant only affects the degree 0 term.
* Added the optional parameter `lmax` to `SHCoeffs.plot_spectrum()` and `SHCoeffs.plot_spectrum2d()`.
* Fixed a bug in `SHCoeffs.pad()` where the attribute `mask` was not similarly padded.
* For mathematical operations with `SHCoeffs` grids, it is now required that the
two class instances have the same `lmax`.
* Clarified the documentation of `SHRotateCoef` to point out that this is only valid for intrinsically real functions that are expressed in complex harmonics.
* Added the method `volume()` to the class SHCoeffs, that calculates the volume of the object.
* Added the attributes `area` and `shannon` to `SHWindow`, which provides the area of the concentration domain and the shannon number.
* Removed python installation support from Makefile: use pip instead.

**Citation:**

M. A. Wieczorek, M. Meschede, E. Sales de Andrade, I. Oshchepkov, B. Xu, and A. Walker (2018). SHTOOLS: Version 4.3, Zenodo, doi:[10.5281/zenodo.1346663](https://doi.org/10.5281/zenodo.1346663)


## Version 4.2

**Change log:**

* Full support added for the use of unnormalized harmonics in the classes `SHCoeffs` and `SHGrid`. To make use of this normalization, just specify `normalization='unnorm'`.
* Added a new python native routine `mag_spectrum()` in the subpackage `gravmag` that replaces the original fortran wrapped routines. The old python wrapped functions have been removed from pyshtools. This new routine is nearly the same as `spectrum()`, and further allows one to compute the spectrum of the potential or magnetic intensity.
* The old fortran based `SHRead` function has been replaced with a python native version `shread()`. The functionality is nearly identical as before, and combines the previous routines `SHRead`, `SHReadError`, `SHReadErrorH` and `SHReadH` into one. Differences include: (1) It is no longer necessary to specify the lmax of the file: This is determined automatically by reading the file from the end, (2) both real and complex coefficients are supported, (3) a header line can be output, but it is a simple list of type str that will need to be converted to the correct format by the user, and (4) "comment" lines are read and ignored: A valid line is one where there are 4 or more words, and where the first two words are integers.
* A new python native function `convert()` was added in the subpackage `shio` that converts between arrays of spherical harmonic coefficients with different normalizations. The class `SHCoeffs` was then simplified by using this external function for all conversions involving `SHCoeffs` class instances.
* The optional parameter lmax was added to `SHCoeffs.spectrum()`.
* When plotting grid from the class `SHGrid`, one can now specify the label to use for the x and y axes with `xlabel` and `ylabel`, as well as the interval to use when plotting ticks on both axes using `tick_interval`.
* The pyshtools rotation routines now allow you to specify the optional parameter `convention` to treat Euler angles in either the `x` or `y` conventions (i.e., which axes to use for the second rotation). Furthermore, the optional argument `body` allows you to specify if you want to rotate the body (True), or coordinate system (False, default). The tutorial number 3 was updated to clear up some inconsistencies in how the angles were defined.
* New optional parameters added to `SHWindow.plot_windows()` and `SHWindow.plot_spectra()` that include `xlim` and `ylim` for the limits when plotting spectra, `maxcolumns` for the number of columns to use when plotting several windows, and `lmax` which controls the grid spacing when plotting the windows.
* Added the optional argument `lmax` to `SHCoeffs.from_random()` that allows you to create coefficients with maximum bandwidths that are either greater or less than the bandwidth of the input power spectrum.
* Added a warning message when using `SHCoeffs.rotate()` with degrees greater than 1200, as the routine is not accurate beyond this value.
* Added an optional argument `legend` to `SHWindow.plot_windows()` to control whether the legend is plotted or not.
* Fixed a minor bug in ClassExample.py file concerning the use of `SHCoeffs.rotate()`.
* Added support for plotting to an already existing figure by allowing the user to specify an existing matplotlib axes.
* Removed some non-standard ascii dashes in the documentation files, and forced all doc files to be opened as `utf-8`.
* HTML documentation was completely redone using Jekyll. The markdown source files are now located in `pages/mydoc`. A static html web site is built using `jekyll`, whose files are located in `doc`. Github will automatically create the static pages and serve them on [shtools.github.io/SHTOOLS](https://shtools.github.io/SHTOOLS/). To build the static pages yourself, it is only necessary to execute `bundle exec jekyll build` in the directory `doc`, which will build the site into `_site` in the same directory. Alternatively, `make www` in the main directory will create a static site in the top-level directory `www` that could be used to deploy on a different web server. The site is based on the template [Jekyll documentation theme](https://github.com/tomjoht/documentation-theme-jekyll) by @tomjoht.

**Citation:**

M. A. Wieczorek, M. Meschede, E. Sales de Andrade, I. Oshchepkov, B. Xu, and A. Walker (2018). SHTOOLS: Version 4.2, Zenodo, doi:[10.5281/zenodo.1250054](https://doi.org/10.5281/zenodo.1250054)


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

M. A. Wieczorek, M. Meschede, E. Sales de Andrade, I. Oshchepkov, B. Xu, and A. Walker (2017). SHTOOLS: Version 4.1, Zenodo, doi:[10.5281/zenodo.1067108](https://doi.org/10.5281/zenodo.1067108)


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

M. A. Wieczorek, M. Meschede, I. Oshchepkov, E. Sales de Andrade, and heroxbd  (2016). SHTOOLS: Version 4.0. Zenodo. doi:[10.5281/zenodo.206114](https://doi.org/10.5281/zenodo.206114)
