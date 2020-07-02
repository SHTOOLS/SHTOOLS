![LOGO](misc/logo.png)

[![Documentation](https://img.shields.io/badge/documentation-shtools.github.io%2FSHTOOLS%2F-yellow.svg)](https://shtools.github.io/SHTOOLS/)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.592762.svg)](https://doi.org/10.5281/zenodo.592762)
[![Paper](https://img.shields.io/badge/paper-10.1029/2018GC007529-orange.svg)](https://doi.org/10.1029/2018GC007529)
[![Chat on matrix](https://img.shields.io/badge/chat-on_[matrix]-4bb596.svg)](https://matrix.to/#/!SrkiFczPSWmYrlSNYF:matrix.org?via=matrix.org)
[![Chat at gitter](https://badges.gitter.im/SHTOOLS/SHTOOLS.svg)](https://gitter.im/SHTOOLS/SHTOOLS?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Twitter](https://img.shields.io/twitter/follow/pyshtools.svg?style=social&label=Follow)](https://twitter.com/intent/follow?screen_name=pyshtools)

SHTOOLS/pysthools is a Fortran-95/Python library that can be used to perform
spherical harmonic transforms, multitaper spectral analyses, expansions of functions into Slepian bases, and standard operations on global gravitational and magnetic field data.

### FEATURES ###

* Supports all standard normalizations and phase conventions of the spherical harmonic functions.

* Effortless conversion between real and complex harmonics, between phase conventions, and between 4&pi; normalized, Schmidt semi-normalized, orthonormalized, and unnormalized harmonics.

* Use of both regularly sampled geographic grids and grids appropriate for Gauss-Legendre quadrature.

* Spherical harmonic transforms proven to be accurate up to about degree 2800.

* Perform localized multitaper spectral analyses, or expand functions in terms of localized Slepian bases.

* Support for standard data and file formats, including *xarray* and *netcdf*.

* Import research-grade gravity, topography, and magnetic field datasets with a single command.

* Creation of publication quality maps using [Cartopy](https://scitools.org.uk/cartopy) and [pygmt](https://www.pygmt.org/).

* OpenMP compatible and OpenMP thread-safe versions of the Fortran routines.

### INSTALLATION ###
#### pyshtools for Python ####

Binary install using pip or conda (Linux, macOS and Windows):
```bash
pip install pyshtools
pip install --upgrade pyshtools  # to upgrade a pre-existing installation
conda install -c conda-forge pyshtools  # Linux and macOS only
```

Build from source:
```bash
pip install pyshtools --no-binary pyshtools
```

To install the develop branch from source:
```bash
pip install git+https://github.com/SHTOOLS/SHTOOLS@develop
```

#### pyshtools for developers ####
Linux requirements:
```bash
sudo apt-get install libblas-dev liblapack-dev g++ gfortran libfftw3-dev tcsh
```
macOS requirements:
```bash
brew install fftw
# for lapack and blas, link to the system '-framework Accelerate'
```

Clone the shtools repo and then install:
```bash
git clone https://github.com/SHTOOLS/SHTOOLS.git
cd shtools
git checkout develop
pip install .  # install into the active python environment lib folder, or
pip install -e .  # install into the SHTOOLS/pyshtools folder and link to the active python environment
```

#### Fortran Library ####
Clone the shtools repo, and then execute one of the following commands in the shtools directory:
```bash
make fortran
make fortran-mp  # for OpenMP Fortran routines
```
Or use the [brew](http://brew.sh/) package manager (macOS):
```bash
brew tap shtools/shtools
brew install shtools
brew install shtools --with-openmp  # to install shtools with the OpenMP components.
```

Further installation instructions and options can be found in the [web documentation](https://shtools.github.io/SHTOOLS/).

### HOW TO USE ###

SHTOOLS can be invoked in any Fortran 95 or Python program. The core software is written in Fortran 95, and Python wrappers and dedicated classes allow simple access to the fortran-compiled routines. A variety of Python tutorials and guides are included that demonstrate the major features of the library.

To get started, check out the following Python tutorial notebooks:

* [Spherical harmonic coefficients and grids](https://shtools.github.io/SHTOOLS/pages/mydoc/notebooks/grids-and-coefficients.html)
* [Localization windows and spectral analysis](https://shtools.github.io/SHTOOLS/pages/mydoc/notebooks/localized-spectral-analysis.html)
* [Gravity and magnetic fields](https://shtools.github.io/SHTOOLS/pages/mydoc/notebooks/gravity-and-magnetic-fields.html)
* [Plotting maps](https://shtools.github.io/SHTOOLS/pages/mydoc/notebooks/plotting-maps.html)

### DEVELOPERS ###

We work on the `develop` branch and only push releases to `master`. Please base all pull requests on `develop`.

### REFERENCE ###

Mark A. Wieczorek and Matthias Meschede (2018). SHTools --- Tools for working with spherical harmonics, *Geochemistry, Geophysics, Geosystems*, 19, 2574-2592, doi:[10.1029/2018GC007529](https://doi.org/10.1029/2018GC007529).
