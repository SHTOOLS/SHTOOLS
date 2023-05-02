<p align="center">
    <img alt="SHTOOLS LOGO" src="misc/logo.png" width="100%">
    <a href="https://shtools.github.io/SHTOOLS/"><img alt="Documentation" src="https://img.shields.io/badge/documentation-shtools.github.io%2FSHTOOLS%2F-yellow.svg"></a>
    <a href="https://doi.org/10.5281/zenodo.592762"><img alt="DOI" src="https://zenodo.org/badge/doi/10.5281/zenodo.592762.svg"></a>
    <a href="https://doi.org/10.1029/2018GC007529"><img alt="Paper" src="https://img.shields.io/badge/paper-10.1029/2018GC007529-orange.svg"></a>
    <a href="https://mybinder.org/v2/gh/SHTOOLS/SHTOOLS/master?filepath=examples%2Fnotebooks%2F"><img alt="Binder" src="https://mybinder.org/badge_logo.svg"></a>
    <a href="https://matrix.to/#/%23pyshtools:matrix.org"><img alt="Chat on matrix" src="https://img.shields.io/badge/chat-on_[matrix]-4bb596.svg"></a>
    <a href="https://app.gitter.im/#/room/#pyshtools:matrix.org"><img alt="Chat at gitter" src="https://badges.gitter.im/SHTOOLS/SHTOOLS.svg"></a>
    <img alt="License" src="https://img.shields.io/badge/License-BSD_3--Clause-brightgreen.svg">
    <a href="https://fosstodon.org/@shtools"><img alt="Mastodon Follow" src="https://img.shields.io/mastodon/follow/108112567255227248?domain=https%3A%2F%2Ffosstodon.org&style=social"></a>
</p>

SHTOOLS/pyshtools is a Fortran-95/Python library that can be used for spherical harmonic transforms, multitaper spectral analyses, expansions of gridded data into Slepian basis functions, standard operations on global gravitational and magnetic field data.

### TABLE OF CONTENTS
<!-- toc -->

- [FEATURES](#features)
- [HOW TO USE](#how-to-use)
- [INSTALLATION](#installation)
  - [pyshtools (for Python)](#pyshtools-for-python)
    - [Install using `conda`](#install-using-conda)
    - [Install using `pip`](#install-using-pip)
    - [For Developers](#for-developers)
  - [SHTOOLS (for Fortran 95)](#shtools-for-fortran-95)
    - [Install using the `brew` package manager (MacOS, Linux, Windows)](#install-using-the-brew-package-manager-macos-linux-windows)
    - [Install using the `macports` package manager (MacOS)](#install-using-the-macports-package-manager-macos)
    - [Install from source](#install-from-source)
- [CONTRIBUTING](#contributing)
- [MAINTAINERS AND CONTRIBUTORS](#maintainers-and-contributors)
- [LICENSE](#license)
- [REFERENCES](#references)

<!-- tocstop -->

### FEATURES

* Supports all standard normalizations and phase conventions of the spherical harmonic functions.

* Effortless conversion between real and complex harmonics, and between different normalization and phase conventions.

* Use of both regularly sampled geographic grids and grids appropriate for Gauss-Legendre quadrature.

* Spherical harmonic transforms proven to be accurate up to about degree 2800 for the native Fortran 95 backend and beyond using the [DUCC0](https://gitlab.mpcdf.mpg.de/mtr/ducc) backend.

* Perform localized multitaper spectral analyses, or expand gridded data in terms of localized Slepian basis functions.

* Support for standard data and file formats, including *xarray* and *netcdf*.

* Import research-grade gravity, topography, and magnetic field datasets with a single command.

* Creation of publication quality maps using [Cartopy](https://scitools.org.uk/cartopy) and [PyGMT](https://scitools.org.uk/cartopy/docs/latest).

* Support multithreaded programming using the [OpenMP](https://www.openmp.org) API.

### HOW TO USE

A variety of Python tutorials and guides are available to explain the main library features. To get started, click on the following Python tutorials and run them interactively in Binder:

* [Spherical harmonic coefficients and grids](https://nbviewer.jupyter.org/github/SHTOOLS/SHTOOLS/blob/master/examples/notebooks/grids-and-coefficients.ipynb)
* [Localization windows and spectral analysis](https://nbviewer.jupyter.org/github/SHTOOLS/SHTOOLS/blob/master/examples/notebooks/localized-spectral-analysis.ipynb)
* [Gravity and magnetic fields](https://nbviewer.jupyter.org/github/SHTOOLS/SHTOOLS/blob/master/examples/notebooks/gravity-and-magnetic-fields.ipynb)
* [Plotting maps](https://nbviewer.jupyter.org/github/SHTOOLS/SHTOOLS/blob/master/examples/notebooks/plotting-maps.ipynb)

### INSTALLATION

SHTOOLS can be invoked in any Fortran 95 or Python program. The core software is written in Fortran 95, and Python wrappers and dedicated classes allow simple access to the fortran-compiled routines. To install it, run these commands below:

#### pyshtools (for Python)

##### Install using `conda`:
```bash
conda install -c conda-forge pyshtools  # Linux and macOS only
conda update -c conda-forge pyshtools  # to upgrade a pre-existing installation
```

##### Install using `pip`:
```bash
pip install pyshtools
pip install --upgrade pyshtools  # to upgrade a pre-existing installation
pip install pyshtools --no-binary pyshtools  # build from source
pip install git+https://github.com/SHTOOLS/SHTOOLS@develop  # install the develop branch from source
```

##### For developers:
Install the requirements:
```bash
# Linux: install gfortran, fftw3, blas, and lapack
sudo apt-get install g++ gfortran libfftw3-dev libblas-dev liblapack-dev
# macOS: install fftw using brew or macports
brew install fftw
sudo port install fftw-3
# macOS: for LAPACK, link to the system '-framework Accelerate' or install openblas
```

Then clone the shtools repo and install manually:
```bash
git clone https://github.com/SHTOOLS/SHTOOLS.git
cd shtools
git checkout develop
pip install -e .  # install into the shtools folder and link to the active python environment
```

#### SHTOOLS (for Fortran 95)

##### Install using the [brew](http://brew.sh/) package manager (MacOS, Linux, Windows):
```bash
brew install shtools
```

##### Install using the [macports](https://www.macports.org/) package manager (MacOS):
```bash
sudo port install shtools
```

##### Install from source:
Clone or download the this repo, and then execute one (or both) of the following commands in the `shtools` directory:
```bash
make fortran
make fortran-mp  # for OpenMP Fortran routines
```

Further installation instructions and options can be found in the [web documentation](https://shtools.github.io/SHTOOLS/).

### CONTRIBUTING

We work on the `develop` branch and only push releases to `master`. Please base all pull requests on `develop`.

### CONTRIBUTORS

For the full list of contributors, see the [AUTHORS](./AUTHORS.md) file.

### LICENSE
This project uses the BSD 3-Clause license, as found in the [LICENSE](./LICENSE.txt) file.

### REFERENCES
- Mark A. Wieczorek and Matthias Meschede (2018). SHTools --- Tools for working with spherical harmonics, *Geochemistry, Geophysics, Geosystems*, 19, 2574-2592, doi:[10.1029/2018GC007529](https://doi.org/10.1029/2018GC007529).
- Mark Wieczorek, *et al*. (2019). SHTOOLS/SHTOOLS. Zenodo, doi:[10.5281/zenodo.3457861](https://doi.org/10.5281/zenodo.3457861)
