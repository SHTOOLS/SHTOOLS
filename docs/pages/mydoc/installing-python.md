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
If you wish to instead compile the archive yourself, first make sure that you have the necessary dependencies installed. On most linux distributions, this can be accomplished using
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
