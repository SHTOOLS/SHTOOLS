---
title: "Installing SHTOOLS"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: installing.html
summary: SHTOOLS can be installed in several ways. If you will be using only the Python components, you should use the pip package manager. If you will be writing and compiling Fortran 95 code, you should use either the brew package manager (on macOS) or compile manually using the Makefile.
toc: true
---

## Python package installer (pip)

The easiest way to install the Python components of SHTOOLS is to use `pip`. On Linux, macOS and windows machines, the binary wheels can be installed by executing the following command:
```bash
pip install pyshtools
```
If you wish to instead compile the archive yourself, first make sure that you have the necessary dependencies installed. On most linux distributions, this can be accomplished using
```bash
sudo apt-get install libblas-dev liblapack-dev g++ gfortran libfftw3-dev tcsh
```
or on macOS using [brew](https://brew.sh/)
```bash
brew install fftw --with-fortran
```
Then clone the SHTOOLS repo
```bash
git clone https://github.com/SHTOOLS/SHTOOLS.git
```
To install pyshtools in the active Python environment lib folder, execute this command in the main directory:
```bash
pip install .
```
To instead install the files in the current working directory and link them to the system Python directory, use
```bash
pip install -e .
```
To uninstall pyshtools from your system directory, use the command
```bash
pip uninstall pyshtools
```
Note that these commands will install only the Python version that corresponds to the version of `pip` being used. On some systems, it may be necessary to specify `pip2.7` or `pip3.x`.

## Fortran 95 library using brew (macOS)

If [brew](https://brew.sh/) is already installed, it is only necessary to enter the following commands in the terminal:
```bash
brew tap shtools/shtools
brew install shtools
```
To also install both the standard Fortran 95 library along with the OpenMP components, add the option `--with-openmp` to the last command:
```bash
brew install shtools --with-openmp
```

## Fortran 95 library using the Makefile

Before installing the Fortran 95 components of SHTOOLS, it will be necessary to have a Fortran 95 compiler and the [FFTW](http://www.fftw.org), [LAPACK](http://www.netlib.org/lapack/), and [BLAS](http://www.netlib.org/blas/) libraries installed on your computer. After this is done, the Fortran 95 components of SHTOOLS can be compiled in most cases by executing the following command in a unix shell in the main directory:
```bash
make
```
To compile the Fortran 95 components with OpenMP use
```bash
make fortran-mp
```
To compile and run the Fortran 95 test suites, use
```bash
make fortran-tests
```
To delete the compiled archive, module files, object files, and test suite files, use
```bash
make clean
```
By default, the `Makefile` will use the `gfortran` compiler. Different compilers and options can be specified by passing optional arguments to the `Makefile` using the syntax
```bash
make F95 = "MyCompiler" F95FLAGS = "MyCompilerFlags"
```
where "`MyCompiler`" and "`MyCompilerFlags`" are to be replaced by the path of the compiler name and options, respectively. Supported options include:
```bash
F95 = "Name (including path) of the f95 compiler"
F95FLAGS = "F95 compiler options"
FFTW = Name and path of the FFTW3 library of the form "-Lpath -lfftw3"
LAPACK = Name and path of the LAPACK library of the form "-Lpath -llapack"
BLAS = Name and path of the BLAS library of the form "-Lpath -lblas"
```
Successful compilation will create the library file `libSHTOOLS.a` (and `libSHTOOLS-mp.a` when compiling with OpenMP) in the directory `lib`, and will place a few compiled module files in the directory `modules`. If you need to recompile SHTOOLS a second time using a different set of compiler flags, it will be necessary to first remove all the previously compiled object files by using `make clean`.

To make all files available at a system level, execute
```bash
make install
```
This will move the compiled SHTOOLS files and documentation to
```
SYSLIBPATH (libSHTOOLS.a, libSHTOOLS-mp.a)
SYSMODPATH (fftw3.mod, planetsconstants.mod, shtools.mod)
SYSSHAREPATH/shtools/examples (example files)
SYSSHAREPATH/man/man1 (man pages)
SYSDOCPATH/shtools (index.html, web documentation)
```
The locations of the above directories can be set as optional arguments passed to the `Makefile`, and the default locations are
```bash
SYSLIBPATH = /usr/local/lib
SYSMODPATH = /usr/local/include
SYSSHAREPATH = /usr/local/share
SYSDOCPATH = /usr/local/share/doc
```
To remove all installed SHTOOLS files, use
```bash
make uninstall
```
To access the unix man pages, it will be necessary to add `SYSSHAREPATH/man` to your man path.

## Fortran 95 compiler specific flags and optimizations

Default compiler options are specified in the main `Makefile` for a few common compilers (`gfortran`, Absoft `f95`, `g95`, and `ifort`). If it is necessary to change these, consider the following guidelines.

One should always use some form of optimization when compiling SHTOOLS, such as by specifying the option
```bash
-O3
```
Performance will likely be increased by 10s of percent by specifying the compiler to target the host architecture
```bash
-march=host  (f95)
-march=native  (gfortran)
```
and to use fast math
```bash
-speed_math=11  (f95)
-ffast-math  (gfortran)
```
The biggest difficulty in compiling SHTOOLS is setting the compiler flags so that the external subroutine names are in a format that is compatible with the already compiled FFTW and LAPACK libraries. In general, it is necessary to ensure that the SHTOOLS subroutine names are in lower case and have the right number of underscores appended to them.

For Absoft ProFortran, this is achieved by setting
```bash
-YEXT_NAMES=LCS -YEXT_SFX=_
```
For `g95`, it will be necessary to use one of the following:
```bash
-fno-second-underscore (most likely)
-fno-underscoring
```
For `gfortran`, it is generally not necessary to use any special flags, though it could arise that one of the following might be necessary:
```bash
-fno-underscoring
-fsecond-underscore
```
For the Intel Fortran compiler `ifort`, it will be necessary to use
```bash
-free -Tf
```
in order that the compiler recognizes files with the extension `.f95` as Fortran 95 files. In this case, the f95 file should come directly after the option `-Tf`.

Setting the right compiler flags is more complicated when the FFTW and LAPACK libraries have different naming and underscoring conventions. In order to accommodate this case, underscores can be added explicitly to either the LAPACK or FFTW subroutine names in the SHTOOLS source code by specifying the optional `make` arguments when building the archive:
```bash
FFTW_UNDERSCORE = 1 (to add an extra underscore to the FFTW routine names)
LAPACK_UNDERSCORE = 1 (to add an extra underscore to the LAPACK routine names)
```
For both cases, compiler flags should probably be set so that underscores are not appended to routine names. See the [FAQ](faq.html) for further information.

To generate 64 bit code, use the compiler option
```bash
-m64
```
For this case, it will be necessary to use 64-bit compiled FFTW and LAPACK libraries.


