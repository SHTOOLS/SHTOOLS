---
title: SlepianCoeffsToSH (Python)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: pyslepiancoeffstosh.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Convert a function expressed in Slepian coefficients to spherical harmonic coefficients.

## Usage

`flm` = SlepianCoeffsToSH(`salpha`, `galpha`, `nalpha`)

## Returns

`flm` :float, dimension (2, `lmax`+1, `lmax`+1)
:   The spherical harmonic coefficients of the global function.

## Parameters

`salpha` :float, dimension (`nalpha`)
:   A vector containing the Slepian coefficients of the function.

`galpha` : float, dimension ((`lmax`+1)\*\*2, `nalpha`)
:   An array containing the spherical harmonic coefficients of the Slepian functions, where `lmax` is the spherical harmonic bandwidth of the functions. Each column corresponds to a single function of which the spherical harmonic coefficients can be unpacked with `SHVectorToCilm`.

`nalpha` : input, integer
:   The number of expansion coefficients to compute.

## Description

`SlepianCoeffsToSH` will compute the spherical harmonic coefficients of a global function `flm` given the Slepian functions `galpha` and the corresponding Slepian coefficients `salpha`. The Slepian functions are determined by a call to either (1) `SHReturnTapers` and then `SHRotateTapers`, or (2) `SHReturnTapersMap`. Each row of `galpha` contains the (`lmax`+1)\*\*2 spherical harmonic coefficients of a Slepian function that can be unpacked using `SHVectorToCilm`. The Slepian functions must be normalized to have unit power (that is the sum of the coefficients squared is 1), and the spherical harmonic coefficients are calculated as

`f_lm = sum_{i}^{nalpha} s(alpha) g_lm(alpha)`  

## See also

[shreturntapers](pyshreturntapers.html), [shreturntapersmap](pyshreturntapersmap.html), [shrotatetapers](pyshrotatetapers.html), [shvectortocilm](pyshvectortocilm.html)
