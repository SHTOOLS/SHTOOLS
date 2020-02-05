---
title: "Installing pyshtools"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: installing-python.html
summary: The python components of SHTOOLS can be installed using the pip package manager.
toc: true
---

## Python package installer (pip)

The easiest way to install the Python components of SHTOOLS (pyshtools) is to use `pip`. On Linux and macOS architectures, the binary wheels can be installed by executing one of the following commands
```bash
pip install pyshtools
pip install --upgrade pyshtools  # to upgrade a pre-existing installation
```
If you wish to instead compile the archive yourself, first make sure that you have the required dependencies installed for the Fortran-95 components. On most linux distributions, this can be accomplished using
```bash
sudo apt-get install libblas-dev liblapack-dev g++ gfortran libfftw3-dev tcsh
```
or on macOS using [brew](https://brew.sh/)
```bash
brew install fftw
```
(Note that SHTOOLS supports the use of any FFTW3-compatible library, such as Intel's [MKL](https://software.intel.com/en-us/mkl)). Then build from source using the command
```bash
pip install pyshtools --no-binary pyshtools
```

If you would like to modify the source code, you should clone the SHTOOLS repo:
```bash
git clone https://github.com/SHTOOLS/SHTOOLS.git
```
Once the repo is cloned, enter the directory, and use either the command
```bash
pip install .
```
to install pyshtools in the active Python environment lib folder, or
```bash
pip install -e .
```
to install the files in the current working directory and link them to the system Python directory.

To uninstall pyshtools from your system directory, use the command
```bash
pip uninstall pyshtools
```
Note that these commands will install only the Python version that corresponds to the version of `pip` being used. On some systems, it may be necessary to specify explicitly `pip3.x`.

### Python dependencies

In order to use the most basic aspects of pyshtools, it will be necessary to install the python packages [numpy](https://numpy.org/), [scipy](https://www.scipy.org/), and [matplotlib](https://matplotlib.org/). Furthermore, [astropy(https://www.astropy.org/) is required for the planetary constants module, [xarray](https://xarray.pydata.org/en/stable/#) is required for netcdf file support, and [requests](https://2.python-requests.org/en/master/#) is required when reading files from urls. All of these packages should be installed automatically when installing pyshtools.

In addition to these packages, it will be necessary to install manually [cartopy](https://scitools.org.uk/cartopy/docs/latest/) and/or [pygmt](https://www.pygmt.org) in order to access the geographic projections of the plotting functions. Finally, the package [palettable](https://jiffyclub.github.io/palettable/) is required by one of the notebooks, and this is useful for providing access to a suite of scientific color maps.
