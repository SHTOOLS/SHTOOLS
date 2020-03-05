![LOGO](misc/logo.png)

[![Documentation](https://img.shields.io/badge/documentation-shtools.github.io%2FSHTOOLS%2F-yellow.svg)](https://shtools.github.io/SHTOOLS/)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.592762.svg)](https://doi.org/10.5281/zenodo.592762)
[![Paper](https://img.shields.io/badge/paper-10.1029/2018GC007529-orange.svg)](https://doi.org/10.1029/2018GC007529)
[![Join the chat at https://gitter.im/SHTOOLS/SHTOOLS](https://badges.gitter.im/SHTOOLS/SHTOOLS.svg)](https://gitter.im/SHTOOLS/SHTOOLS?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Twitter](https://img.shields.io/twitter/follow/pyshtools.svg?style=social&label=Follow)](https://twitter.com/intent/follow?screen_name=pyshtools)

SHTOOLS/pysthools is a Fortran-95/Python library that can be used to perform
spherical harmonic transforms, multitaper spectral analyses on the sphere, expansions of functions into Slepian bases, and standard operations on global gravitational and magnetic field data.

### FEATURES ###

* A wide range of supported spherical harmonic functions:
   * real and complex,
   * all standard normalizations (Geodesy 4&pi;, Schmidt semi-normalized, orthonormalized, unnormalized),
   * Condon-Shortley phase factor of (-1)<sup>m</sup>.

* Clean implementation of the spherical harmonic transforms:
  * exact quadrature rules using the sampling theorem of *Driscoll and Healy* (1994) or Gauss-Legendre quadrature,
  * accurate and fast to approximately degree 2800,
  * OpenMP compatible and OpenMP thread-safe.

* Selected applications and routines:
   * global spectral analysis, spherical harmonic rotations, and Wigner 3j symbols,
   * localized multitaper spectral analyses, and expansions in Slepian basis functions,
   * standard gravity and magnetic field calculations, computation of the geoid, and finite-amplitude potential from topography.

* Integrated support for working with xarray and netcdf data.

* Publication quality maps using [Cartopy](https://scitools.org.uk/cartopy) or [pygmt](https://www.pygmt.org/dev/).

* SHTOOLS is open source software (3-clause BSD license).

### INSTALLATION ###
#### pyshtools for Python ####

Binary install for linux and macOS:
```bash
pip install pyshtools
pip install --upgrade pyshtools  # upgrade to a new version
```

Build from source:
```bash
pip install pyshtools --no-binary pyshtools
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
To install the develop branch:
```bash
pip install git+https://github.com/SHTOOLS/SHTOOLS@develop
```
Alternatively, clone the shtools repo and then install:
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
brew install shtools --with-openmp # to install shtools with the OpenMP components.
```

More installation instructions and options can be found in the [web documentation](https://shtools.github.io/SHTOOLS/).

### HOW TO USE ###

SHTOOLS can be invoked in any Fortran 95 or Python program. The core software is written in Fortran 95, and Python wrappers allow simple access to the fortran-compiled routines. A variety of Python notebooks and example files are included that demonstrate the major features of the library.

To get started, check out the following Python tutorial notebooks:

* [Introduction to Grids and Spherical Harmonic Coefficients](https://shtools.github.io/SHTOOLS/pages/mydoc/notebooks/Introduction-1.html)
* [Introduction to Localization Windows and Spectral Analysis](https://shtools.github.io/SHTOOLS/pages/mydoc/notebooks/Introduction-2.html)
* [Introduction to Gravity and Magnetic Fields](https://shtools.github.io/SHTOOLS/pages/mydoc/notebooks/Introduction-3.html)
* [Introduction to Plotting Maps](https://shtools.github.io/SHTOOLS/pages/mydoc/notebooks/Introduction-4.html)

### DEVELOPERS ###

We work on the `develop` branch and only push releases to `master`. Please base all pull requests on `develop`.

### CITATION ###
Mark A. Wieczorek and Matthias Meschede (2018). SHTools --- Tools for working with spherical harmonics, *Geochemistry, Geophysics, Geosystems*, 19, 2574-2592, doi:[10.1029/2018GC007529](https://doi.org/10.1029/2018GC007529).
