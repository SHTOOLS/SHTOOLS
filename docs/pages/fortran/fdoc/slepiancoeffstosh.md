---
title: SlepianCoeffsToSH (Fortran)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: fortran_sidebar
permalink: slepiancoeffstosh.html
summary:
tags: [fortran]
toc: false
editdoc: fdoc
---

Convert a function expressed in Slepian coefficients to spherical harmonic coefficients.

## Usage

call SlepianCoeffsToSH(`film`, `falpha`, `galpha`, `lmax`, `nmax`, `exitstatus`)

## Parameters

`film` : output, real(dp), dimension (2, `lmax`+1, `lmax`+1)
:   The spherical harmonic coefficients of the global function.

`falpha` : input, real(dp), dimension (`nmax`)
:   A vector containing the Slepian coefficients of the function.

`galpha` : input, real(dp), dimension ((`lmax`+1)**2, `nmax`)
:   An array containing the spherical harmonic coefficients of the Slepian functions. Each column corresponds to a single function of which the spherical harmonic coefficients can be unpacked with `SHVectorToCilm`.

`lmax` : input, integer
:   The spherical harmonic bandwidth of the Slepian functions.

`nmax` : input, integer
:   The number of expansion coefficients to compute. This must be less than or equal to (`lmax`+1)\*\*2.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

## Description

`SlepianCoeffsToSH` will compute the spherical harmonic coefficients of a global function `film` given the Slepian functions `galpha` and the corresponding Slepian coefficients `falpha`. The Slepian functions are determined by a call to either (1) `SHReturnTapers` and then `SHRotateTapers`, or (2) `SHReturnTapersMap`. Each row of `galpha` contains the (`lmax`+1)**2 spherical harmonic coefficients of a Slepian function that can be unpacked using `SHVectorToCilm`. The Slepian functions must be normalized to have unit power (that is the sum of the coefficients squared is 1), and the spherical harmonic coefficients are calculated as

`f_ilm = sum_{alpha}^{nmax} f_alpha g(alpha)_lm`  

## See also

[slepiancoeffs](slepiancoeffs.html), [shslepianvar](shslepianvar.html), [shreturntapers](shreturntapers.html), [shreturntapersmap](shreturntapersmap.html), [shrotatetapers](shrotatetapers.html), [shvectortocilm](shvectortocilm.html), [shscouplingmatrix](shscouplingmatrix.html), [shscouplingmatrixcap](shscouplingmatrixcap.html), [shmtcouplingmatrix](shmtcouplingmatrix.html),
