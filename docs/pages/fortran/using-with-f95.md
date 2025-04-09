---
title: "Using SHTOOLS"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: fortran_sidebar
permalink: using-with-f95.html
summary: SHTOOLS subroutines and functions can be accessed easily from any Fortran 95 program. It is only necessary to use the SHTOOLS module and link to the compiled archive.
toc: true
folder: fortran
---

## Calling SHTOOLS routines

To access the SHTOOLS functions and subroutines in your code, it is necessary to place the statement
```fortran
use SHTOOLS
```
directly after the program, subroutine, or function declaration (i.e., before an `implicit none` statement). The SHTOOLS module contains an interface block that declares the subroutines and functions used in this archive and allows for the use of implicitly shaped arrays. It should be noted that all arrays passed to a subroutine or function can be *larger* than that specified in the documentation.

## Linking to the SHTOOLS library

When compiling a code that references SHTOOLS, it will be necessary to link to the SHTOOLS library and module files. For this purpose, it may be useful to define the environment variables `SHTOOLSMODPATH` and `SHTOOLSLIBPATH` using the default file locations. In a C-shell, this is accomplished using
```bash
setenv SHTOOLSLIBPATH = "/usr/local/lib"
setenv SHTOOLSMODPATH = "/usr/local/include"
```
whereas for a bash shell the syntax is
```bash
export SHTOOLSLIBPATH = "/usr/local/lib"
export SHTOOLSMODPATH = "/usr/local/include"
```
In addition, most routines require linking to libraries that are compatible with the fast Fourier transform package [FFTW](http://www.fftw.org), and the linear algebra packages [LAPACK](https://www.netlib.org/lapack/) and [BLAS](https://www.netlib.org/blas/). SHTOOLS is compatible with the FFT routines in Intel's [MKL](https://software.intel.com/en-us/mkl) library. Furthermore, it is noted that there are many different system optimized versions of LAPACK and BLAS, some of which may be pre-installed on your machine.

Typical examples of compiling and linking a program `MyProgram.f95` to the necessary library and module files are given below for several common compilers.

### gfortran
```bash
gfortran MyProgram.f95 -I$SHTOOLSMODPATH -L$SHTOOLSLIBPATH -lSHTOOLS -lfftw3 -lm -llapack -lblas -O3 -o MyProgram
```

### Absoft Pro Fortran (f95)
```bash
f95 MyProgram.f95 -p $SHTOOLSMODPATH -L$SHTOOLSLIBPATH -YEXT_NAMES=LCS -YEXT_SFX=_ -lSHTOOLS -lfftw3 -lm -llapack -lblas -O3 -o MyProgram
```

### g95
```bash
g95 MyProgram.f95 -I$SHTOOLSMODPATH -L$SHTOOLSLIBPATH -fno-second-underscore -lSHTOOLS -lfftw3 -lm -llapack -lblas -O3 -o MyProgram
```

### Intel Fortran (ifort)
```bash
ifort -fpp -free -I$SHTOOLSMODPATH -L$SHTOOLSLIBPATH -lSHTOOLS -lfftw3 -lm -llapack -lblas -O3 -Tf MyProgram.f95 -o MyProgram
```
Note that the position of the source file in the above examples might be important for some compilers. Note also that on macOS, linking to the preinstalled LAPACK and BLAS libraries can be accomplished using `-framework Accelerate`. (See [installing SHTOOLS](installing-fortran.html) for more details).

## OpenMP
The Fortran 95 routines in the OpenMP version of SHTOOLS are OpenMP compatible and OpenMP thread-safe. To use these routines in a program that uses OpenMP, it is necessary to link to the library `libSHTOOLS-mp.a` and to add one of the following compiler options:
```bash
-fopenmp  # gfortran
-openmp  # Absoft Pro Fortran
-qopenmp  # ifort
```

## Fortran array indices
Fortran arrays usually start with an index of 1. However, with spherical harmonic functions and coefficients a degree-0 term exists. In order to deal with this, all arrays that are a function of spherical harmonic degree `l` and order `m` have 1 added to each of their indices. For instance, the real Fortran array `clm` is related to the spherical harmonic coefficients $$C_{lm}$$ by

$$ C_{lm} = \left \lbrace
\begin{array}{ll}
\mbox{clm}(1,l+1,m+1)  & \mbox{if } m\ge 0 \\
\mbox{clm}(2,l+1,m+1)  & \mbox{if } m < 0.\\
\end{array}
\right. $$

In this notation, the positive and negative angular orders correspond to the cosine and sine coefficients, respectively (see the page [real spherical harmonics](real-spherical-harmonics.html) for more information).

## Using optional parameters
Many of the subroutines and functions in this archive can accept one or more optional parameters. To specify these parameters, it is only necessary to use the syntax
```fortran
call SHRead (FILENAME, CILM, LMAX, SKIP=1, ERROR=errorcoef)
```
where `SKIP` and `ERROR` are the names of two optional parameters, and the constant `1` and variable `errorcoef` are their respective arguments.

## Documentation

Documentation for the Fortran 95 subroutines and functions in SHTOOLS can be accessed by their unix man pages, using all lower case letters. As an example, to access the `MakeGridDH` man page, use
```bash
man makegriddh
```
Alternatively, the man pages can be accessed from the *Fortran 95 Reference* menu on this web site.
