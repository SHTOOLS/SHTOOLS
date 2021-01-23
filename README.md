![LOGO](misc/logo.png)

[![Documentation](https://img.shields.io/badge/documentation-shtools.github.io%2FSHTOOLS%2F-yellow.svg)](https://shtools.github.io/SHTOOLS/)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.592762.svg)](https://doi.org/10.5281/zenodo.592762)
[![Paper](https://img.shields.io/badge/paper-10.1029/2018GC007529-orange.svg)](https://doi.org/10.1029/2018GC007529)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/SHTOOLS/SHTOOLS/master?filepath=examples%2Fnotebooks%2F)
[![Chat on matrix](https://img.shields.io/badge/chat-on_[matrix]-4bb596.svg)](https://matrix.to/#/#pyshtools:matrix.org?via=matrix.org&via=gitter.im)
[![Chat at gitter](https://badges.gitter.im/SHTOOLS/SHTOOLS.svg)](https://gitter.im/SHTOOLS/SHTOOLS?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Twitter](https://img.shields.io/twitter/follow/pyshtools.svg?style=social&label=Follow)](https://twitter.com/intent/follow?screen_name=pyshtools)

SHTOOLS/pyshtools is a Fortran-95/Python library that can be used to perform
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
#### pyshtools (for Python) ####

Install using `conda`:
```bash
conda install -c conda-forge pyshtools  # Linux and macOS only
conda update -c conda-forge pyshtools  # to upgrade a pre-existing installation
```

Install using `pip`:
```bash
pip install pyshtools
pip install --upgrade pyshtools  # to upgrade a pre-existing installation
pip install pyshtools --no-binary pyshtools  # build from source
pip install git+https://github.com/SHTOOLS/SHTOOLS@develop  # install the develop branch from source
```

For developers, install the requirements
```bash
# Linux: install gfortran, fftw3, blas, and lapack
sudo apt-get install g++ gfortran libfftw3-dev libblas-dev liblapack-dev
# macOS: install fftw using brew or macports
brew install fftw
sudo port install fftw-3
# macOS: for LAPACK, link to the system '-framework Accelerate' or install openblas
```

then clone the shtools repo and install manually:
```bash
git clone https://github.com/SHTOOLS/SHTOOLS.git
cd shtools
git checkout develop
pip install -e .  # install into the shtools folder and link to the active python environment
```

#### SHTOOLS (for Fortran 95) ####

Install using the [brew](http://brew.sh/) package manager (macOS, linux, windows):
```bash
brew install shtools
```

Install using the [macports](https://www.macports.org/) package manager (macOS)
```bash
sudo port install shtools
```

Install from source. Clone or download the shtools repo, and then execute one (or both) of the following commands in the shtools directory:
```bash
make fortran
make fortran-mp  # for OpenMP Fortran routines
```

Further installation instructions and options can be found in the [web documentation](https://shtools.github.io/SHTOOLS/).

### HOW TO USE ###

SHTOOLS can be invoked in any Fortran 95 or Python program. The core software is written in Fortran 95, and Python wrappers and dedicated classes allow simple access to the fortran-compiled routines. A variety of Python tutorials and guides are included that demonstrate the major features of the library.

To get started, click on the following Python tutorials and then run them interactively in Binder:

* [Spherical harmonic coefficients and grids](https://nbviewer.jupyter.org/github/SHTOOLS/SHTOOLS/blob/master/examples/notebooks/grids-and-coefficients.ipynb)
* [Localization windows and spectral analysis](https://nbviewer.jupyter.org/github/SHTOOLS/SHTOOLS/blob/master/examples/notebooks/localized-spectral-analysis.ipynb)
* [Gravity and magnetic fields](https://nbviewer.jupyter.org/github/SHTOOLS/SHTOOLS/blob/master/examples/notebooks/gravity-and-magnetic-fields.ipynb)
* [Plotting maps](https://nbviewer.jupyter.org/github/SHTOOLS/SHTOOLS/blob/master/examples/notebooks/plotting-maps.ipynb)

### DEVELOPERS ###

We work on the `develop` branch and only push releases to `master`. Please base all pull requests on `develop`.

### REFERENCE ###

Mark A. Wieczorek and Matthias Meschede (2018). SHTools --- Tools for working with spherical harmonics, *Geochemistry, Geophysics, Geosystems*, 19, 2574-2592, doi:[10.1029/2018GC007529](https://doi.org/10.1029/2018GC007529).
