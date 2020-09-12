---
title: "Installing SHTOOLS"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: fortran_sidebar
permalink: fortran-installing.html
summary: SHTOOLS can be installed using the brew and macports package managers, or manually using the Makefile.
toc: true
folder: fortran
---

## brew (macOS)

If the [brew](https://brew.sh/) package manager is already installed, it is only necessary to enter the following commands in the terminal:
```bash
brew tap shtools/shtools
brew install shtools
```
To install the OpenMP components along with the standard Fortran 95 library, add the option `--with-openmp` to the last command:
```bash
brew install shtools --with-openmp
```
The shtools library `libSHTOOLS.a` will be installed in the directory `/usr/local/lib`, and the compiled module files will be installed in `/usr/local/include`.

To install the example data files and test programs in the directory `/usr/local/share/shtools/` use:
```bash
brew install shtools --with-examples
```
To run the test suite, use
```bash
brew test shtools
```
The output of the tests can be inspected in the user's `Library/Logs/Homebrew/shtools` directory.

## macports (macOS)

If the [macports](https://www.macports.org/) package manager is already installed, it is only necessary to enter the following command in the terminal:
```bash
sudo port install shtools
```
To install the OpenMP components along with the standard Fortran 95 library, add the option `+openmp` to the last command:
```bash
sudo port install shtools +openmp
```
The shtools library `libSHTOOLS.a` will be installed in the directory `/opt/local/lib`, and the compiled module files will be installed in `/opt/local/include`. To run the test suite, which is located in `/opt/local/share/examples/shtools`, use the command
```bash
sudo port test shtools
```

## Using the Makefile

Before trying to install the Fortran 95 components of SHTOOLS, it will be necessary to have a Fortran 95 compiler and [LAPACK](https://www.netlib.org/lapack/), [BLAS](https://www.netlib.org/blas/) and [FFTW3](http://www.fftw.org)-compatible libraries. On most linux distributions, this can be accomplished using
```bash
sudo apt-get install gcc gfortran libfftw3-dev libblas-dev liblapack-dev
```
or on macOS using
```bash
brew install fftw  # using brew
sudo port install fftw-3  # using macports
conda install fftw  # using conda
# lapack and blas can be accessed by linking to the system '-framework Accelerate'
```

The Fortran 95 components of SHTOOLS can then be compiled in most cases by executing the following command in a unix shell in the main directory:
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
```bash
SYSLIBPATH  # libSHTOOLS.a, libSHTOOLS-mp.a
SYSMODPATH  # ftypes.mod fftw3.mod, planetsconstants.mod, shtools.mod
SYSSHAREPATH/shtools/examples  # example files
SYSSHAREPATH/man/man1  # man pages
SYSDOCPATH/shtools  # index.html, web documentation 
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
-march=host  # f95
-march=native  # gfortran
```
and to use fast math
```bash
-speed_math=11  # f95
-ffast-math  # gfortran
```
The biggest difficulty in compiling SHTOOLS is setting the compiler flags so that the external subroutine names are in a format that is compatible with the already compiled FFTW and LAPACK libraries. In general, it is necessary to ensure that the SHTOOLS subroutine names are in lower case and have the right number of underscores appended to them.

For Absoft ProFortran, this is achieved by setting
```bash
-YEXT_NAMES=LCS -YEXT_SFX=_
```
For `g95`, it will be necessary to use one of the following:
```bash
-fno-second-underscore # most likely
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

Setting the right compiler flags is more complicated when the FFTW and LAPACK libraries have different naming and underscoring conventions. In order to accommodate this case, underscores can be added explicitly to the LAPACK subroutine names in the SHTOOLS source code by specifying the optional `make` arguments when building the archive:
```bash
LAPACK_UNDERSCORE=1  # add an extra underscore to the LAPACK routine names
```
For this case, compiler flags should probably be set so that underscores are not appended to routine names. See [Fortran 95 problems](fortran-95-problems.html) for further information.

To generate 64 bit code, use the compiler option
```bash
-m64
```
For this case, it will be necessary to use 64-bit compiled FFTW and LAPACK libraries.
