---
title: SHMultiTaperCSE (Fortran)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: shmultitapercse.html
summary:
tags: [fortran]
toc: false
editdoc: fdoc
---

Perform a localized multitaper cross-spectral analysis using spherical cap windows.

## Usage

call SHMultiTaperCSE (`mtse`, `sd`, `sh1`, `lmax1`, `sh2`, `lmax2`, `tapers`, `taper_order`, `lmaxt`, `k`, `alpha`, `lat`, `lon`, `taper_wt`, `norm`, `csphase`, `exitstatus`)

## Parameters

`mtse` : output, real\*8, dimension (`lmax`-`lmaxt`+1)
:   The localized multitaper cross-power spectrum estimate. `lmax` is the smaller of `lmax1` and `lmax2`.

`sd` : output, real\*8, dimension (`lmax`-`lmaxt`+1)
:   The standard error of the localized multitaper cross-power spectral estimates. `lmax` is the smaller of `lmax1` and `lmax2`.

`sh1` : input, real\*8, dimension (2, `lmax1`+1, `lmax1`+1)
:   The spherical harmonic coefficients of the first function.

`lmax1` : input, integer
:   The spherical harmonic bandwidth of `sh1`.

`sh2` : input, real\*8, dimension (2, `lmax2`+1, `lmax2`+1)
:   The spherical harmonic coefficients of the second function.

`lmax2` : input, integer
:   The spherical harmonic bandwidth of `sh2`.

`tapers` : input, real\*8, dimension (`lmaxt`+1, `k`)
:   An array of the `k` windowing functions, arranged in columns, obtained from a call to `SHReturnTapers`. Each window has non-zero coefficients for a single angular order that is specified in the array `taper_order`.

`taper_order` : input, integer, dimension (`k`)
:   An array containing the angular orders of the spherical harmonic coefficients in each column of the array `tapers`.

`lmaxt` : input, integer
:   The spherical harmonic bandwidth of the windowing functions in the array `tapers`.

`k` : input, integer
:   The number of tapers to be utilized in performing the multitaper spectral analysis.

`alpha` : input, optional, real\*8, dimension (3)
:   The Euler rotation angles used in rotating the windowing functions. `alpha(1)=0`, `alpha(2)=-(90-lat)*pi/180`, `alpha(3)=-lon*pi/180`. Either `alpha` or `lat` and `lon` can be specified, but not both. If none of these are specified, the spectral analysis will be centered at the north pole.

`lat` : input, optional, real\*8
:   The latitude in degrees of the localized analysis. Either `alpha` or `lat` and `lon` can be specified but not both. If none of these are specified, the spectral analysis will be centered at the north pole.

`lon` : input, optional, real\*8
:   The longitude in degrees of the localized analysis. Either `alpha` or `lat` and `lon` can be specified, but not both. If none of these are specified, the spectral analysis will be centered at the north pole.

`taper_wt` : input, optional, real\*8, dimension (`k`)
:   The weights used in calculating the multitaper spectral estimates and standard error. Optimal values of the weights (for a known global power spectrum) can be obtained from the routine `SHMTVarOpt`.
	
`norm` : input, optional, integer, default = 1
:   1 (default) = 4-pi (geodesy) normalized harmonics; 2 = Schmidt semi-normalized harmonics; 3 = unnormalized harmonics; 4 = orthonormal harmonics.

`csphase` : input, optional, integer, default = 1
:   1 (default) = do not apply the Condon-Shortley phase factor to the associated Legendre functions; -1 = append the Condon-Shortley phase factor of (-1)^m to the associated Legendre functions.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

## Description

`SHMultiTaperCSE` will perform a localized multitaper cross-spectral analysis of two input functions expressed in spherical harmonics, `SH1` and `SH2`. The maximum degree of the localized multitaper power spectrum estimate is `lmax-lmaxt`, where `lmax` is the smaller of `lmax1` and `lmax2`. The coefficients and angular orders of the windowing coefficients (`tapers` and `taper_order`) are obtained by a call to `SHReturnTapers`. If `lat` and `lon` or `alpha` is specified, then the symmetry axis of the localizing windows will be rotated to these coordinates. Otherwise, the localized spectral analysis will be centered over the north pole.

If the optional array `taper_wt` is specified, then these weights will be used in calculating a weighted average of the individual `k` tapered estimates (`mtse`) and the corresponding standard error of the estimates (`sd`). If not present, the weights will all be assumed to be equal. When `taper_wt` is not specified, the mutltitaper spectral estimate for a given degree will be calculated as the average obtained from the `k` individual tapered estimates. The standard error of the multitaper estimate at degree l is simply the population standard deviation, `S = sqrt(sum (Si - mtse)^2 / (k-1))`, divided by sqrt(`k`). See Wieczorek and Simons (2007) for the relevant expressions when weighted estimates are used.

The employed spherical harmonic normalization and Condon-Shortley phase convention can be set by the optional arguments `norm` and `csphase`; if not set, the default is to use geodesy 4-pi normalized harmonics that exclude the Condon-Shortley phase of (-1)^m.

## References

Wieczorek, M. A. and F. J. Simons, Minimum-variance multitaper spectral estimation on the sphere, J. Fourier Anal. Appl., 13, doi:10.1007/s00041-006-6904-1, 665-692, 2007.

## See also

[shmultitaperse](shmultitaperse.html), [shreturntapers](shreturntapers.html), [shreturntapersm](shreturntapersm.html), [shmtvaropt](shmtvaropt.html)
