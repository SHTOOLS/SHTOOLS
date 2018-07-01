![LOGO](misc/logo.png)

[![Documentation](https://img.shields.io/badge/documentation-shtools.github.io%2FSHTOOLS%2F-yellow.svg)](https://shtools.github.io/SHTOOLS/)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.592762.svg)](http://dx.doi.org/10.5281/zenodo.592762)
[![Paper](https://img.shields.io/badge/paper-10.1029/2018GC007529-orange.svg)](http://doi.org/10.1029/2018GC007529)
[![Join the chat at https://gitter.im/SHTOOLS/SHTOOLS](https://badges.gitter.im/SHTOOLS/SHTOOLS.svg)](https://gitter.im/SHTOOLS/SHTOOLS?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Twitter](https://img.shields.io/twitter/follow/pyshtools.svg?style=social&label=Follow)](https://twitter.com/intent/follow?screen_name=pyshtools)

SHTOOLS/pysthools is a Fortran-95/Python library that can be used to perform
spherical harmonic transforms and reconstructions, rotations of data expressed
in spherical harmonics, and multitaper spectral analyses on the sphere.

### FEATURES ###

* A wide range of supported spherical harmonic functions:
   * real and complex,
   * all standard normalizations (Geodesy 4&pi;, Schmidt semi-normalized, orthonormalized, unnormalized),
   * Condon-Shortley phase factor of (-1)<sup>m</sup>.

* Selected applications and routines:
   * global spectral analysis, spherical harmonic rotations, Wigner 3j symbols,
   * localized multitaper spectral analyses, optimal window generation, spherical harmonic coupling matrices,
   * standard gravity and magnetic field calculations, computation of the geoid, finite-amplitude potential from topography.

* Clean implementation of the spherical harmonic transforms:
  * exact quadrature rules using the sampling theorem of *Driscoll and Healy* (1994) or Gauss-Legendre quadrature,
  * accurate and fast to approximately degree 2800 (spatial resolution higher than 4 arc minutes),
  * Fortran 95 routines are OpenMP compatible and OpenMP thread-safe.

### INSTALLATION ###
#### pyshtools for Python ####

Binary install for linux/macOS/windows:
```bash
pip install pyshtools
```
Build from source:
```bash
pip install pyshtools --no-binary pyshtools
```

#### pyshtools for Anaconda Python ####

Install fftw3 with fortran bindings and then install pyshtools using `pip`:

```bash
conda install -c eumetsat fftw3
pip install pyshtools
```

#### pyshtools for Python (developer install) ####
Linux requirements:
```bash
sudo apt-get install libblas-dev liblapack-dev g++ gfortran libfftw3-dev tcsh
```
macOS requirements:
```bash
brew install fftw --with-fortran
```
To install the develop branch use:
```bash
pip install git+https://github.com/SHTOOLS/SHTOOLS@develop
```
Alternatively, clone the shtools repo
```bash
git clone https://github.com/SHTOOLS/SHTOOLS.git
```
and then execute one of the following commands in the shtools directory:
```bash
pip install .  # installs into the active python environment lib folder
pip install -v -e .  # installs into the SHTOOLS/pyshtools folder and links to the active python environment
```

#### Fortran Library ####
Clone the shtools repo, and then execute one of the following commands in the shtools directory:
```bash
make fortran
make fortran-mp  # OpenMP Fortran routines
```
Or use the [brew](http://brew.sh/) package manager (macOS):
```bash
brew tap shtools/shtools
brew install shtools
brew install shtools --with-openmp # to install shtools with the OpenMP components.
```

More installation instructions and options can be found in the [web documentation](https://shtools.github.io/SHTOOLS/) and GitHub 
[wiki](https://github.com/SHTOOLS/SHTOOLS/wiki).


### HOW TO USE ###

SHTOOLS can be invoked from Fortran 95, Python 2 or Python 3. The
base SHTOOLS software is written in Fortran 95, and the Python library allows
simple access to all fortran-compiled routines and offers helper routines as
well as simple interfaces.

To get started, check out the following Python tutorial notebooks:

* [Introduction 1: Grids and Spherical Harmonic Coefficients](https://shtools.github.io/SHTOOLS/pages/mydoc/notebooks/Introduction-1.html)
* [Introduction 2: Localization Windows and Spectral Analysis](https://shtools.github.io/SHTOOLS/pages/mydoc/notebooks/Introduction-2.html)

### ACKNOWLEDGMENTS ###
SHTOOLS is open source (revised BSD license) and makes use of the freely
available Fourier transform package
[FFTW](http://www.fftw.org) and the linear algebra packages
[LAPACK](http://www.netlib.org/lapack/) and
[BLAS](http://www.netlib.org/blas/).

### CITATION ###
Mark A. Wieczorek and Matthias Meschede (2018). SHTools --- Tools for working with spherical harmonics, *Geochemistry, Geophysics, Geosystems*, doi:[10.1029/2018GC007529](http://doi.org/10.1029/2018GC007529).

M. A. Wieczorek, M. Meschede, E. Sales de Andrade, I. Oshchepkov, B. Xu, and A. Walker (2018). SHTOOLS, *Zenodo*, doi:[10.5281/zenodo.592762](http://doi.org/10.5281/zenodo.592762).
