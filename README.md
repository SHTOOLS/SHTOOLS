![LOGO](misc/logo.png)

[![Join the chat at https://gitter.im/SHTOOLS/SHTOOLS](https://badges.gitter.im/SHTOOLS/SHTOOLS.svg)](https://gitter.im/SHTOOLS/SHTOOLS?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.592762.svg)](http://dx.doi.org/10.5281/zenodo.592762)

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

#### pyshtools for Anaconda Python ####

fftw3 with fortran bindings and pyshtools can be installed on anaconda like 
this

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
Clone the shtools repo
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

More installation instructions and options can be found in the [web documentation](https://shtools.oca.eu) and GitHub 
[wiki](https://github.com/SHTOOLS/SHTOOLS/wiki).


### HOW TO USE ###

SHTOOLS can be invoked from Fortran 95, Python 2 or Python 3. The
base SHTOOLS software is written in Fortran 95, and the Python library allows
simple access to all fortran-compiled routines and offers helper routines as
well as simple interfaces.

To get started, check out the following Python tutorial notebooks:

* Introduction 1: Grids and Spherical Harmonic Coefficients [\[ipynb\]](examples/notebooks/Introduction-1.ipynb)
* Introduction 2: Localization Windows and Spectral Analysis [\[ipynb\]](examples/notebooks/Introduction-2.ipynb)

You can keep up to date by following SHTOOLS on [Twitter](https://twitter.com/pyshtools).

### ACKNOWLEDGMENTS ###
SHTOOLS is open source (revised BSD license) and makes use of the freely
available Fourier transform package
[FFTW](http://www.fftw.org) and the linear algebra packages
[LAPACK](http://www.netlib.org/lapack/) and
[BLAS](http://www.netlib.org/blas/).

### CITATION ###
M. A. Wieczorek, M. Meschede, E. Sales de Andrade, I. Oshchepkov, B. Xu, and A. Walker (2017). SHTOOLS, Zenodo, doi:[10.5281/zenodo.592762](http://doi.org/10.5281/zenodo.592762)
