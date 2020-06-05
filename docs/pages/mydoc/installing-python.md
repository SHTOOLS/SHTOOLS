---
title: "Installing pyshtools"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: installing-python.html
summary: The python components of SHTOOLS can be installed using the pip package manager or conda.
toc: true
---

## Python package installer (pip)

On Linux, macOS and Windows architectures, the binary wheels can be installed using `pip` by executing one of the following commands:
```bash
pip install pyshtools
pip install --upgrade pyshtools  # to upgrade a pre-existing installation
```

## Conda package installer

The binary pre-compiled pyshtools library can also be installed with the conda package manager using the command:
```bash
conda install -c conda-forge pyshtools  # Linux and macOS only
```
The conda packages do not support Windows architectures at the present time.

## Build from source (pip)

If you wish to compile the archive yourself, first make sure that you have the required dependencies installed for the Fortran-95 components. On most Linux distributions, this can be accomplished using
```bash
sudo apt-get install libblas-dev liblapack-dev g++ gfortran libfftw3-dev tcsh
```
or on macOS using [brew](https://brew.sh/)
```bash
brew install fftw
```
Then build from source using the command
```bash
pip install pyshtools --no-binary pyshtools
```
Note that SHTOOLS supports the use of any FFTW3-compatible library, such as Intel's [MKL](https://software.intel.com/en-us/mkl).

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
Note that these commands will install only the Python version that corresponds to the version of `pip` being used. On some systems, it may be necessary to specify explicitly `pip3` or `pip3.x`.

## Python dependencies

In order to use the most basic aspects of pyshtools, it will be necessary to install the following python packages:

* [numpy](https://numpy.org/)
* [scipy](https://www.scipy.org/)
* [matplotlib](https://matplotlib.org/)
* [astropy](https://www.astropy.org/): required for the planetary constants module.
* [xarray](https://xarray.pydata.org/en/stable/#): required for netcdf file support.
* [requests](https://2.python-requests.org/en/master/#): required when reading files from urls.
* [pooch](https://www.fatiando.org/pooch/latest/index.html): required for reading datasets.
* [tqdm](https://tqdm.github.io/): required for showing progress bars when downloading datasets.

All of the above packages should be installed automatically when installing pyshtools with pip or conda. For more specialized operations, it will be necessary to install manually the following packages:

* [Cartopy](https://scitools.org.uk/cartopy/docs/latest/): required for Cartopy map projections.
* [pygmt](https://www.pygmt.org): required for pygmt map projections.
* [palettable](https://jiffyclub.github.io/palettable/): scientific color maps required by one of the tutorials.
