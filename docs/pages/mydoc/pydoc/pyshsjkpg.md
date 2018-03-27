---
title: SHSjkPG (Python)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: pyshsjkpg.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Calculate the expectation of the product of two functions, each multiplied by a different data taper, for a given spherical harmonic degree and two different angular orders.

## Usage

`value` = SHSjkPG (`incspectra`, `l`, `m`, `mprime`, `hj_real`, `hk_real`, `mj`, `mk`, `lwin`, `hkcc`)

## Returns

`value` : complex
:   The expectation of the product of two functions, each multiplied by a different data taper, for a given spherical harmonic degree and two different angular orders.

## Parameters

`incspectra` : float, dimension (`l`+`lwin`+1)
:   The global cross-power spectrum of `f` and `g`.

`l` : integer
:   The spherical harmonic degree for which to calculate the expectation.

`m` :  integer
:   The angular order of the first localized function, `Phi`.

`mprime` : integer
:   The angular order of the second localized function, `Gamma`.

`hj_real` : float, dimension (`lwin`+1)
:   The real spherical harmonic coefficients of angular order `mj` used to localize the first function `f`. These are obtained by a call to `SHReturnTapers`.

`hk_real` : float, dimension (`lwin`+1)
:   The real spherical harmonic coefficients of angular order `mk` used to localize the second function `g`. These are obtained by a call to `SHReturnTapers`.

`mj` : integer
:   The angular order of the window coefficients `hj_real`.

`mk` : integer
:   The angular order of the window coefficients `hk_real`.

`lwin` : integer
:   the spherical harmonic bandwidth of the localizing windows `hj_real` and `hk_real`.

`hkcc` : integer
:   If 1, the function described in the `description` will be calculated as is. If 2, the second localized function `Gamma` will not have its complex conjugate taken.

## Description

`SHSjkPG` will calculate the expectation of two functions (`f` and `g`), each localized by a different data taper that is a solution of the spherical cap concentration problem, for a given spherical harmonic degree and two different angular orders. As described in Wieczorek and Simons (2007), this is the function

      /    m(j)       mprime(k)* \
     |  Phi      Gamma            |
      \    l          l          /

The global cross-power spectrum of `f` and `g` is input as `incspectra`, and the real coefficients of the two data tapers of angular order `mj` and `mk` (obtained by a call to `SHReturnTapers`) are specified by `hj_real` and `hk_real`. If `hkcc` is set to 1, then the above function is calculated as is. However, if this is set to 2, then the complex conjugate of the second localized function is not taken.

## References

Wieczorek, M. A. and F. J. Simons, Minimum-variance multitaper spectral estimation on the sphere, J. Fourier Anal. Appl., 13, doi:10.1007/s00041-006-6904-1, 665-692, 2007.

## See also

[shreturntapers](pyshreturntapers.html), [shmtvaropt](pyshmtvaropt.html)
