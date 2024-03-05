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
conda install -c conda-forge pyshtools
conda update -c conda-forge pyshtools  # to upgrade a pre-existing installation
```

## Python package installer (pip)

The binary wheels can be installed using `pip` by executing one of the following commands:
```bash
pip install pyshtools
pip install --upgrade pyshtools  # to upgrade a pre-existing installation
```
In order to use the map projection routines, it will be necessary to install either Cartopy and/or pygmt separately, as described in the section [Python dependencies](#dependencies).

## Build from source

If you wish to compile the archive yourself, you will need to make sure that you have the required dependencies installed, which include a fortran and C compiler, and FFTW3, BLAS and LAPACK compatible libraries. First, ensure that the relevant build utilities and compilers are installed:
```bash
# Debian, Ubuntu and derivatives
sudo apt-get install build-essential cmake gfortran
# Fedora, Centos, RHEL and derivatives
sudo dnf group install "C Development Tools and Libraries" "Development Tools"
sudo dnf install cmake gcc-fortran
# macOS
xcode-select --install
```
Next, clone the SHTOOLS repo and enter the directory
```bash
git clone https://github.com/SHTOOLS/SHTOOLS.git
cd shtools
```
All dependencies can be then be installed using the conda `environment.yml` file that is found in the top level directory:
```bash
conda create -n your_env_name python=3.xx  # create a new conda environment, if desired
conda env update -n your_env_name -f environment.yml
conda activate your_env_name  # activate the new environment
```
To build pyshtools from source and install in the active Python environment lib folder, use the command
```bash
pip install .
```
To instead install the files in the current working directory in editable mode and link them to the system Python directory, use the command
```bash
pip install --no-build-isolation -e .
```

When pyshtools is installed in editable mode, a number of tests and benchmarks can be run from within the repo. First, you will need to find out the name of the build directory by typing
```bash
ls build
```
which should return something like `cp312`. Then, the tests and benchmarks can be run by executing the following commands in the top level directory
```bash
meson test -C build/cp312 --suite python
meson test -C build/cp312 --suite python-notebooks
meson test -C build/cp312 --suite fortran
meson test -C build/cp312 --benchmark --suite python
meson test -C build/cp312 --benchmark --suite fortran
```

To uninstall pyshtools from your system directory, use the command
```bash
pip uninstall pyshtools
```
Note that on some systems, it may be necessary to specify explicitly `pip3` or `pip3.x` for the `pip` command.

## Python dependencies {#dependencies}

When installing pyshtools using `pip` or `conda`, the following packages should be installed automatically:

* [numpy](https://numpy.org/): required for all numerical calculations.
* [meson-python](https://meson-python.readthedocs.io/en/latest/#): required to build pyshtools.
* [setuptools_scm](https://setuptools-scm.readthedocs.io/en/latest/): required to obtain the pyshtools version number.
* [scipy](https://www.scipy.org/): required for a few specialized functions.
* [matplotlib](https://matplotlib.org/): required for most plotting functions.
* [astropy](https://www.astropy.org/): required for the constants module.
* [xarray](https://xarray.pydata.org/en/stable/#): required for netcdf file support.
* [requests](https://2.python-requests.org/en/master/#): required when reading files from urls.
* [pooch](https://www.fatiando.org/pooch/latest/index.html): required for reading datasets.
* [tqdm](https://tqdm.github.io/): required for showing progress bars when downloading datasets.

When installing pyshtools using `conda`, the following will also be installed automatically:

* [fftw](https://www.fftw.org/): required for the fortran library.
* [blas-devel](https://anaconda.org/conda-forge/blas-devel): required for the fortran components.
* [cartopy](https://scitools.org.uk/cartopy/docs/latest/): required for Cartopy map projections.
* [pygmt](https://www.pygmt.org) (>=0.7): required for pygmt map projections. pygmt requires (see below) *gmt (>=6.3.0)*.
* [ducc0](https://gitlab.mpcdf.mpg.de/mtr/ducc) (>=0.15): required for using the 'ducc' backend for spherical harmonic transforms.
* [palettable](https://jiffyclub.github.io/palettable/): scientific color maps required by one of the tutorials.

Thus, when installing pyshtools with `pip`, it will be necessary to install compatible versions of BLAS, LAPACK and FFTW3, as well as Cartopy, pygmt, and ducc0.

### How to install BLAS, LAPACK and FFTW3

The easiest way to install BLAS, LAPACK and FFTW is with `conda`:
```bash
conda install -c conda-forge blas-devel>=3.8 fftw>=3.3.8
```
If it is not possible to use conda, these can instead be installed using the system package manager:
```bash
sudo apt-get install libfftw3-dev libblas-dev liblapack-dev  # Debian, Ubuntu and derivatives
sudo dnf install blas-devel lapack-devel fftw-devel  # Fedora, Centos, RHEL and derivatives
brew install fftw  # macOS
```
Note that pyshtools supports the use of any FFTW3-compatible library, such as Intel's [MKL](https://software.intel.com/en-us/mkl).

### How to install pygmt

In order to use the *pygmt* plotting routines, it will be necessary to install both *pygmt (>=0.7)* and the *gmt (>=6.3.0)* library. This is most easily achieved using conda with
```bash
conda install -c conda-forge pygmt gmt
```
Alternatively, *pygmt* can be installed using `pip`
```bash
pip install pyshtools[pygmt]  # installs pygmt at the same time as pyshtools
pip install pygmt  # install pygmt with pip
```
For this case, however, the *gmt* library will need to be installed using other means, such as with your system package manager:
```bash
sudo apt-get install gmt  # Debian, Ubuntu and derivatives
sudo dnf install gmt-devel  # Fedora, Centos, RHEL and derivatives
brew install gmt  # macOS using brew
sudo port install gmt6  # macOS using macports
```

### How to install ducc

To make use of the 'ducc' backend for the spherical harmonic transforms, it will be necessary to install the *ducc0 (>=0.15)* package using either `pip` or `conda`:
```bash
conda install -c conda-forge ducc0>=0.15  # install using conda
pip install ducc0>=0.15  # install using pip
pip install pyshtools[ducc]  # installs ducc at the same time as pyshtools
pip install ducc0>=0.15 -no-binary ducc0  # install ducc from source
```
Note: By installing *ducc0* from source, it might be possible to benefit from the use of AVX instructions that can improve execution speeds by a factor of about 2.
