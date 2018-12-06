---
title: SlepianCoeffs (Python)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: pyslepiancoeffs.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Determine the expansion coefficients of a function for a given set of input Slepian functions.

## Usage

`salpha` = SlepianCoeffs(`galpha`, `flm`, `nalpha`)

## Returns

`salpha` : float, dimension (`nalpha`)
:   A vector containing the Slepian coefficients of the input function `flm`.

## Parameters

`galpha` : float, dimension ((`lmax`+1)\*\*2, `nalpha`)
:   An array containing the spherical harmonic coefficients of the Slepian functions, where `lmax` is the spherical harmonic bandwidth of the functions. Each column corresponds to a single function of which the spherical harmonic coefficients can be unpacked with `SHVectorToCilm`.

`flm` : float, dimension (2, `lmax`+1, `lmax`+1)
:   The spherical harmonic coefficients of the global function to be expanded in Slepian functions.

`nalpha` : integer
:   The number of expansion coefficients to compute.

## Description

`SlepianCoeffs` will compute the Slepian coefficients of a global input function `flm` given the Slepian functions `galpha`. The Slepian functions are determined by a call to either (1) `SHReturnTapers` and then `SHRotateTapers`, or (2) `SHReturnTapersMap`. Each row of `galpha` contains the (`lmax`+1)\*\*2 spherical harmonic coefficients of a Slepian function that can be unpacked using `SHVectorToCilm`. The Slepian functions must be normalized to have unit power (that is the sum of the coefficients squared is 1), and the Slepian coefficients are calculated as

`s(alpha) = sum_{lm}^{lmax} f_lm g_lm(alpha)`  

## See also

[shreturntapers](pyshreturntapers.html), [shreturntapersmap](pyshreturntapersmap.html), [shrotatetapers](pyshrotatetapers.html), [shvectortocilm](pyshvectortocilm.html)
