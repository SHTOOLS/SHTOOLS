---
title: SHMultiTaperSE (Python)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: pyshmultitaperse.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Perform a localized multitaper spectral analysis using spherical cap windows.

## Usage

`mtse`, `sd` = SHMultiTaperSE (`sh`, `tapers`, `taper_order`, [`lmax`,  `lmaxt`, `k`, `lat`, `lon`, `taper_wt`, `norm`, `csphase`])

## Returns

`mtse` : float dimension (`lmax`-`lmaxt`+1)
:   The localized multitaper power spectrum estimate. 

`sd` : float, dimension (`lmax`-`lmaxt`+1)
:   The standard error of the localized multitaper power spectral estimates.

## Parameters

`sh` : float, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The spherical harmonic coefficients of the function to be localized.

`tapers` : float, dimension (`lmaxtin`+1, `kin`)
:   An array of the `k` windowing functions, arranged in columns, obtained from a call to `SHReturnTapers`. Each window has non-zero coefficients for a single angular order that is specified in the array `taper_order`.

`taper_order` : integer, dimension (`kin`)
:   An array containing the angular orders of the spherical harmonic coefficients in each column of the array `tapers`.

`lmax` : optional, integer, default = `lmaxin`
:   The spherical harmonic bandwidth of `sh`. This must be less than or equal to `lmaxin`.

`lmaxt` : optional, integer, default = `lmaxtin`
:   The spherical harmonic bandwidth of the windowing functions in the array `tapers`.

`k` : optional, integer, default = `kin`
:   The number of tapers to be utilized in performing the multitaper spectral analysis.

`lat` : optional, float, default = 90
:   The latitude in degrees of the localized analysis. The default is to perform the spectral analysis at the north pole.

`lon` : optional, float, default = 0
:   The longitude in degrees of the localized analysis.

`taper_wt` : optional, float, dimension (`kin`), default = -1
:   The weights used in calculating the multitaper spectral estimates and standard error. Optimal values of the weights (for a known global power spectrum) can be obtained from the routine `SHMTVarOpt`. The default value specifies not to use `taper_wt`.

`norm` : optional, integer, default = 1
:   1 (default) = 4-pi (geodesy) normalized harmonics; 2 = Schmidt semi-normalized harmonics; 3 = unnormalized harmonics; 4 = orthonormal harmonics.

`csphase` : optional, integer, default = 1
:   1 (default) = do not apply the Condon-Shortley phase factor to the associated Legendre functions; -1 = append the Condon-Shortley phase factor of (-1)^m to the associated Legendre functions.

## Description

`SHMultiTaperSE` will perform a localized multitaper spectral analysis of an input function expressed in spherical harmonics. The maximum degree of the localized multitaper cross-power spectrum estimate is `lmax-lmaxt`. The coefficients and angular orders of the windowing coefficients (`tapers` and `taper_order`) are obtained by a call to `SHReturnTapers`. If `lat` and `lon` are specified, the symmetry axis of the localizing windows will be rotated to these coordinates. Otherwise, the localized spectral analysis will be centered over the north pole.

If the optional array `taper_wt` is specified, these weights will be used in calculating a weighted average of the individual `k` tapered estimates `mtse` and the corresponding standard error of the estimates `sd`. If not present, the weights will all be assumed to be equal. When `taper_wt` is not specified, the mutltitaper spectral estimate for a given degree is calculated as the average obtained from the `k` individual tapered estimates. The standard error of the multitaper estimate at degree `l` is simply the population standard deviation, `S = sqrt(sum (Si - mtse)^2 / (k-1))`, divided by `sqrt(k)`. See Wieczorek and Simons (2007) for the relevant expressions when weighted estimates are used.

The employed spherical harmonic normalization and Condon-Shortley phase convention can be set by the optional arguments `norm` and `csphase`; if not set, the default is to use geodesy 4-pi normalized harmonics that exclude the Condon-Shortley phase of (-1)^m.

## References

Wieczorek, M. A. and F. J. Simons, Minimum-variance multitaper spectral estimation on the sphere, J. Fourier Anal. Appl., 13, doi:10.1007/s00041-006-6904-1, 665-692, 2007.

## See also

[shmultitapercse](pyshmultitapercse.html), [shreturntapers](pyshreturntapers.html), [shreturntapersm](pyshreturntapersm.html), [shmtvaropt](pyshmtvaropt.html)
