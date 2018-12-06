---
title: SlepianCoeffs (Fortran)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: slepiancoeffs.html
summary:
tags: [fortran]
toc: false
editdoc: fdoc
---

Determine the expansion coefficients of a function for a given set of input Slepian functions.

## Usage

call SlepianCoeffs(`salpha`, `galpha`, `flm`, `lmax`, `nalpha`, `exitstatus`)

## Parameters

`salpha` : output, real\*8, dimension (`nalpha`)
:   A vector containing the Slepian coefficients of the input function `flm`.

`galpha` : input, real\*8, dimension ((`lmax`+1)**2, `nalpha`)
:   An array containing the spherical harmonic coefficients of the Slepian functions. Each column corresponds to a single function of which the spherical harmonic coefficients can be unpacked with `SHVectorToCilm`.

`flm` : input, real\*8, dimension (2, `lmax`+1, `lmax`+1)
:   The spherical harmonic coefficients of the global function to be expanded in Slepian functions.

`lmax` : input, integer
:   The spherical harmonic bandwidth of the Slepian functions.

`nalpha` : input, integer
:   The number of expansion coefficients to compute.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

## Description

`SlepianCoeffs` will compute the Slepian coefficients of a global input function `flm` given the Slepian functions `galpha`. The Slepian functions are determined by a call to either (1) `SHReturnTapers` and then `SHRotateTapers`, or (2) `SHReturnTapersMap`. Each row of `galpha` contains the (`lmax`+1)**2 spherical harmonic coefficients of a Slepian function that can be unpacked using `SHVectorToCilm`. The Slepian functions must be normalized to have unit power (that is the sum of the coefficients squared is 1), and the Slepian coefficients are calculated as

`s(alpha) = sum_{lm}^{lmax} f_lm g_lm(alpha)`  

## See also

[shreturntapers](shreturntapers.html), [shreturntapersmap](shreturntapersmap.html), [shrotatetapers](shrotatetapers.html), [shvectortocilm](shvectortocilm.html)
