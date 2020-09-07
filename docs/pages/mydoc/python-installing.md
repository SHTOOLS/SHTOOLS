---
title: "Installing pyshtools"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: python-installing.html
summary: pyshtools can be installed using the conda or pip package manager.
toc: true
folder: mydoc
---

## Conda package installer

The binary pre-compiled pyshtools library, with all required dependencies, can be installed with the conda package manager using the command:
```bash
conda install -c conda-forge pyshtools  # Linux and macOS only
```
The conda packages do not support Windows architectures at the present time.

## Python package installer (pip)

On Linux, macOS and Windows architectures, the binary wheels can be installed using `pip` by executing one of the following commands:
```bash
pip install pyshtools
pip install --upgrade pyshtools  # to upgrade a pre-existing installation
```
In order to use the map projection routines, it will be necessary to install either Cartopy and/or pygmt separately, as described in the section [Python dependencies](#dependencies).

## Build from source (pip)

If you wish to compile the archive yourself, first make sure that you have the required dependencies installed for the Fortran-95 components, which inlcude FFTW3 and LAPACK compatible libraries (see [these instructions](fortran-installing.html)). This can be accomplished on most Linux distributions using
```bash
sudo apt-get install gcc gfortran libfftw3-dev libblas-dev liblapack-dev
```
or on macOS by using one of the following
```bash
brew install fftw  # using brew
sudo port install fftw-3  # using macports
conda install fftw  # using conda
```
Then build from source using the command
```bash
pip install pyshtools --no-binary pyshtools
```
Note that pyshtools supports the use of any FFTW3-compatible library, such as Intel's [MKL](https://software.intel.com/en-us/mkl).

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

## Python dependencies {#dependencies}

When installing pyshtools using `pip` or `conda`, the following packages should be installed automatically:

* [numpy](https://numpy.org/): required for all numerical calculations
* [scipy](https://www.scipy.org/): required for a few specialized functions.
* [matplotlib](https://matplotlib.org/): required for most plotting functions.
* [astropy](https://www.astropy.org/): required for the constants module.
* [xarray](https://xarray.pydata.org/en/stable/#): required for netcdf file support.
* [requests](https://2.python-requests.org/en/master/#): required when reading files from urls.
* [pooch](https://www.fatiando.org/pooch/latest/index.html): required for reading datasets.
* [tqdm](https://tqdm.github.io/): required for showing progress bars when downloading datasets.

When installing pyshtools using `conda`, the following will also be installed automatically:

* [cartopy](https://scitools.org.uk/cartopy/docs/latest/): required for Cartopy map projections. Cartopy requires (see below) *proj*, *geos*, *cython*, *pyshp*, *six*, and *shapely*.
* [pygmt](https://www.pygmt.org) (>=0.2): required for pygmt map projections. pygmt requires (see below) *gmt (>=6.1.1)*.
* [palettable](https://jiffyclub.github.io/palettable/): scientific color maps required by one of the tutorials.

The above three packages will need to be installed separately when installing pyshtools with `pip`, as described in the following subsections.

### How to install Cartopy

The easiest way to install Cartopy is with `conda`:
```bash
conda install -c conda-forge cartopy
```
Cartopy can also be installed using `pip`, but there are several dependencies that need to be installed first. On macOS, this can be accomplished using
```bash
brew install proj geos
pip install --upgrade cython numpy pyshp six
# shapely needs to be built from source to link to geos. If it is already
# installed, uninstall it by: pip3 uninstall shapely
pip install shapely --no-binary shapely
pip install cartopy
```
See [these instructions](https://scitools.org.uk/cartopy/docs/latest/installing.html#installing) for further details.

### How to install pygmt
In order to use the *pygmt* plotting routines, it will be necessary to install both *pygmt (>=0.2)* and the *gmt (>=6.1.1)* library. This is most easily achieved using conda with
```bash
conda install -c conda-forge pygmt gmt
```
Alternatively, *pygmt* can be installed using `pip`
```bash
pip install pygmt
```
The *gmt* library will then need to be installed using other means, such as with brew, macports or apt-get:
```bash
brew install gmt  # using brew on macOS
sudo port install gmt6  # using macports on macOS
sudo apt-get install gmt  # using apt-get on linux
```