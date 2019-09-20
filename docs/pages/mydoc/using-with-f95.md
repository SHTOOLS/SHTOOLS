---
title: "Using SHTOOLS with Fortran 95"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: using-with-f95.html
summary: SHTOOLS subroutines and functions can be accessed easily from any Fortran 95 program. It is only necessary to use the SHTOOLS module and link to the compiled archive.
toc: true
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
In addition, most routines require linking to the fast Fourier transform package [FFTW](http://www.fftw.org), and the linear algebra packages [LAPACK](https://www.netlib.org/lapack/) and [BLAS](https://www.netlib.org/blas/). Typical examples of compiling and linking a program `MyProgram.f95` to the necessary library and module files are given below for several common compilers.

### gfortran
```bash
gfortran MyProgram.f95 -I$SHTOOLSMODPATH -L$SHTOOLSLIBPATH -lSHTOOLS -lfftw3 -lm -llapack -lblas -O3 -m64 -o MyProgram
```

### Absoft Pro Fortran (f95)
```bash
f95 MyProgram.f95 -p $SHTOOLSMODPATH -L$SHTOOLSLIBPATH -YEXT_NAMES=LCS -YEXT_SFX=_ -lSHTOOLS -lfftw3 -lm -llapack -lblas -O3 -m64 -o MyProgram
```

### g95
```bash
g95 MyProgram.f95 -I$SHTOOLSMODPATH -L$SHTOOLSLIBPATH -fno-second-underscore -lSHTOOLS -lfftw3 -lm -llapack -lblas -O3 -m64 -o MyProgram
```

### Intel Fortran (ifort)
```bash
ifort -fpp -free -I$SHTOOLSMODPATH -L$SHTOOLSLIBPATH -lSHTOOLS -lfftw3 -lm -llapack -lblas -O3 -m64 -Tf MyProgram.f95 -o MyProgram
```
Note that the position of the source file in the above examples might be important for some compilers. It may be necessary to modify some options in order to properly link to both the LAPACK and FFTW libraries (see [installing SHTOOLS](installing.html) for more details). If the library ATLAS exists on your system, this could be used instead of BLAS.

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
Alternatively, the man pages can be accessed from the *Fortran 95 Routines* menu on this web site.
