![LOGO](misc/logo.png)

[![Join the chat at https://gitter.im/SHTOOLS/SHTOOLS](https://badges.gitter.im/SHTOOLS/SHTOOLS.svg)](https://gitter.im/SHTOOLS/SHTOOLS?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.61180.svg)](http://dx.doi.org/10.5281/zenodo.61180)

SHTOOLS is a Fortran 95 / Python library that can be used to perform
spherical harmonic transforms and reconstructions, rotations of data expressed
in spherical harmonics, and multitaper spectral analyses on the sphere.

### FEATURES ###

* A wide range of supported spherical harmonic functions:
   * real and complex,
   * different normalizations (Geodesy 4&pi;, Schmidt semi-normalized, orthonormalized, unnormalized),
   * Condon-Shortley phase factor of (-1)<sup>m</sup>.

* Selected applications and routines:
   * global spectral analysis, spherical harmonic rotations, Wigner 3j symbols,
   * localized multitaper spectral analyses, optimal window generation, spherical harmonic coupling matrices,
   * standard gravity and magnetic field calculations, computation of the geoid, finite-amplitude potential from topography.

* Clean implementation of the spherical harmonic transforms:
  * Exact quadrature rules using either the sampling theorem of *Driscoll and Healy* (1994) where data are equally sampled (or spaced) in latitude and longitude, or Gauss-Legendre quadrature.

  * Accurate and fast to approximately degree 2800, corresponding to a spatial
    resolution higher than 4 arc minutes. Transforms and reconstructions take
    on the order of 1 second for bandwidths less than 600 and about 3 minutes
    for bandwidths close to 2800 on standard machines. The Fortran 95 routines are
    OpenMP compatible and OpenMP thread-safe.

### INSTALLATION ###
#### Requirements ####
Linux:
```bash
sudo apt-get install libblas-dev liblapack-dev g++ gfortran libfftw3-dev tcsh
```
OSX:
```bash
brew install fftw --with-fortran
```

#### Python Library ####
```bash
pip install pyshtools
```
Or, to install a developer version, [download](https://github.com/SHTOOLS/SHTOOLS/zipball/master) or clone the SHTOOLS repository, enter the SHTOOLS folder and then execute one of the following commands:
```bash
pip install .  # installs into the active python environment lib folder
pip install -v -e .  # installs into the SHTOOLS/pyshtools folder and links to the active python environment
```

#### Fortran Library ####
To install the Fortran 95 library, enter one of the following
```bash
make fortran
make fortran-mp  # Open-MP Fortran routines
```

Or, with OSX, use the [brew](http://brew.sh/) package manager:
```bash
brew tap shtools/shtools
brew install shtools
```
To also install the OpenMP components, add ```--with-openmp``` to the last command.

More installation instructions and options can be found in the [web documentation](https://shtools.ipgp.fr) and GitHub 
[wiki](https://github.com/SHTOOLS/SHTOOLS/wiki).


### HOW TO USE ###

SHTOOLS can be invoked from Fortran 95, Python 2, or Python 3. The
base SHTOOLS software is written in Fortran 95, and the Python library allows
simple access to all fortran-compiled routines and offers helper routines as
well as simple interfaces.

To get started, check out the following Python tutorial notebooks:

* Introduction 1: Grids and Spherical Harmonic Coefficients [\[ipynb\]](examples/notebooks/Introduction-1.ipynb)
* Introduction 2: Localization Windows and Spectral Analysis [\[ipynb\]](examples/notebooks/Introduction-2.ipynb)
* Tutorial 1: Simple Spherical Harmonic Analyses [\[ipynb\]](examples/notebooks/tutorial_1.ipynb)
* Tutorial 2: Localized Spectral Analysis on the Sphere [\[ipynb\]](examples/notebooks/tutorial_2.ipynb)
* Tutorial 3: The SHTOOLS Class Interface [\[ipynb\]](examples/notebooks/tutorial_3.ipynb)
* Tutorial 4: Spherical Harmonic Normalizations and Parseval's theorem [\[ipynb\]](examples/notebooks/tutorial_4.ipynb)
* Tutorial 5: Multitaper Spectral Analysis Class Interface [\[ipynb\]](examples/notebooks/tutorial_5.ipynb)
* Tutorial 6: 3D Spherical Harmonic Plots [\[ipynb\]](examples/notebooks/tutorial_6.ipynb)

You can keep up to date by following SHTOOLS on [Twitter](https://twitter.com/SH_tools).

### ACKNOWLEDGMENTS ###
SHTOOLS is open source (revised BSD license) and makes use of the freely
available Fourier transform package
[FFTW](http://www.fftw.org) and the linear algebra packages
[LAPACK](http://www.netlib.org/lapack/) and
[BLAS](http://www.netlib.org/blas/).

### CITATION ###
Wieczorek, M. A., M. Meschede, I. Oshchepkov, E. Sales de Andrade (2016). SHTOOLS: Version 3.4. Zenodo. doi:[10.5281/zenodo.61180](http://dx.doi.org/10.5281/zenodo.61180).
