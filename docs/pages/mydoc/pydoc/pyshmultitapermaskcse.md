---
title: SHMultiTaperMaskCSE (Python)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: pyshmultitapermaskcse.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Perform a localized multitaper cross-spectral analysis using arbitrary windows derived from a mask.

## Usage

`mtse`, `sd` = SHMultiTaperMaskCSE (`sh1`, `sh2`, `tapers`, [`lmax1`, `lmax2`, `lmaxt`, `k`, `taper_wt`, `norm`, `csphase`])

## Returns

`mtse` : float, dimension (`lmax`-`lmaxt`+1)
:   The localized multitaper cross-power spectrum estimate. `lmax` is the smaller of `lmax1` and `lmax2`.

`sd` : float, dimension (`lmax`-`lmaxt`+1)
:   The standard error of the localized multitaper cross-power spectral estimates. `lmax` is the smaller of `lmax1` and `lmax2`.

## Parameters

`sh1` : float, dimension (2, `lmax1in`+1, `lmax1in`+1)
:   The spherical harmonic coefficients of the first function.

`sh2` : float, dimension (2, `lmax2in`+1, `lmax2in`+1)
:   The spherical harmonic coefficients of the second function.

`tapers` : float, dimension ((`lmaxtin`+1)**2, `kin`)
:   An array of the `k` windowing functions, arranged in columns, obtained from a call to `SHReturnTapersMap`. The spherical harmonic coefficients are packed according to the conventions in `SHCilmToVector`.

`lmax1` : optional, integer, default = `lmax1in`
:   The spherical harmonic bandwidth of `sh1`. This must be less than or equal to `lmax1in`.

`lmax2` : optional, integer, default = `lmax2in`
:   The spherical harmonic bandwidth of `sh2`. This must be less than or equal to `lmax2in`.

`lmaxt` : optional, integer, default = `lmaxtin`
:   The spherical harmonic bandwidth of the windowing functions in the array `tapers`.

`k` : optional, integer, default = `kin`
:   The number of tapers to be utilized in performing the multitaper spectral analysis.

`taper_wt` : optional, float, dimension (`kin`), default = -1
:   The weights used in calculating the multitaper spectral estimates and standard error. Optimal values of the weights (for a known global power spectrum) can be obtained from the routine `SHMTVarOpt`. The default value specifies not to use `taper_wt`.

`norm` : optional, intger, default = 1
:   1 (default) = 4-pi (geodesy) normalized harmonics; 2 = Schmidt semi-normalized harmonics; 3 = unnormalized harmonics; 4 = orthonormal harmonics.

`csphase` : optional, integer, default = 1
:   1 (default) = do not apply the Condon-Shortley phase factor to the associated Legendre functions; -1 = append the Condon-Shortley phase factor of (-1)^m to the associated Legendre functions.

## Description

`SHMultiTaperMaskCSE` will perform a localized multitaper cross-spectral analysis of two input functions expressed in spherical harmonics, `SH1` and `SH2`, using an arbitrary set of windows derived from a mask. The maximum degree of the localized multitaper power spectrum estimate is `lmax-lmaxt`, where `lmax` is the smaller of `lmax1` and `lmax2`. The matrix `tapers` contains the spherical harmonic coefficients of the windows and can be obtained by a call to `SHReturnTapersMap`. The coefficients of each window are stored in a single column, ordered according to the conventions used in `SHCilmToVector`.

If the optional array `taper_wt` is specified, then these weights will be used in calculating a weighted average of the individual `k` tapered estimates (`mtse`) and the corresponding standard error of the estimates (`sd`). If not present, the weights will all be assumed to be equal. When `taper_wt` is not specified, the mutltitaper spectral estimate for a given degree is calculated as the average obtained from the `k` individual tapered estimates. The standard error of the multitaper estimate at degree l is simply the population standard deviation, `S = sqrt(sum (Si - mtse)^2 / (k-1))`, divided by sqrt(`k`). See Wieczorek and Simons (2007) for the relevant expressions when weighted estimates are used.

The employed spherical harmonic normalization and Condon-Shortley phase convention can be set by the optional arguments `norm` and `csphase`; if not set, the default is to use geodesy 4-pi normalized harmonics that exclude the Condon-Shortley phase of (-1)^m.

## References

Wieczorek, M. A. and F. J. Simons, Minimum-variance multitaper spectral estimation on the sphere, J. Fourier Anal. Appl., 13, doi:10.1007/s00041-006-6904-1, 665-692, 2007.

## See also

[shmultitapermaskse](pyshmultitapermaskse.html), [shreturntapersmap](pyshreturntapersmap.html), [shcilmtovector](pyshcilmtovector.html)
