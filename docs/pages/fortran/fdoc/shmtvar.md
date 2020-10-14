---
title: SHMTVar (Fortran)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: fortran_sidebar
permalink: shmtvar.html
summary:
tags: [fortran]
toc: false
editdoc: fdoc
---

Calculate the theoretical variance of a multitaper spectral estimate for a given input power spectrum.

## Usage

call SHMTVar (`l`, `tapers`, `taper_order`, `lwin`, `kmax`, `sff`, `variance`, `taper_wt`, `unweighted_covar`, `nocross`, `exitstatus`)

## Parameters

`l` : input, integer(int32)
:   The spherical harmonic degree used to calculate the theoretical variance.

`tapers` : input, real(dp), dimension (`lwin`+1, `kmax`)
:   A matrix of localization functions obtained from `SHReturnTapers` or `SHReturnTapersM`.

`taper_order` : input, integer(int32), dimension (`kmax`)
:   The angular order of the windowing coefficients in `tapers`.

`lwin` : input, integer(int32)
:   The spherical harmonic bandwidth of the localizing windows.

`kmax` : input, integer(int32)
:   The maximum number of tapers to be used when calculating the variance.

`sff` : input, real(dp), dimension (`l`+`lwin`+1)
:   The global unwindowed power spectrum of the function to be localized.

`variance` : output, real(dp)
:   The theoretical variance of the multitaper spectral estimate for degree `l`.

`taper_wt` : optional, input, real(dp), dimension (`kmax`)
:   The weights to be applied to the multitaper spectral estimates.

`unweighted_covar` : optional, output, real(dp), dimension (`kmax`, `kmax`)
:   The unweighted covariance matrix of the `kmax` tapers (i.e., Fij in Wieczorek and Simons 2007).

`nocross` : optional, input, integer(int32), default = 0
:   If 1, only the diagonal terms of the covariance matrix Fij will be computed. If 0, all terms will be computed.

`exitstatus` : output, optional, integer(int32)
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

## Description

`SHMTVar` will determine the theoretical variance of a multitaper spectral estimate for a given input power spectrum, degree l, and optionally, a set of taper weights (see eq. C21 of Wieczorek and Simons 2007). The windowing functions are assumed to be solutions to the spherical-cap concentration problem, as determined by a call to `SHReturnTapers` or `SHReturnTapersM`. If `unweighted_covar` is specified, then the unweighted covariance matrix of the `kmax` tapers (i.e., Fij) will be output. If the optional argument `nocross` is set to 1, then only the diagnonal terms of `Fij` will be computed.

## References

Wieczorek, M. A. and F. J. Simons, Minimum-variance multitaper spectral estimation on the sphere, J. Fourier Anal. Appl., 13, doi:10.1007/s00041-006-6904-1, 665-692, 2007.

## See also

[shmtvaropt](shmtvaropt.html), [shreturntapers](shreturntapers.html), [shreturntapersm](shreturntapersm.html), [shmultitaperse](shmultitaperse.html), [shmultitapercse](shmultitapercse.html), [shlocalizedadmitcorr](shlocalizedadmitcorr.html), [shbiasadmitcorr](shbiasadmitcorr.html), [shbiask](shbiask.html), [shmtdebias](shmtdebias.html)
