---
title: SHMTVarOpt (Python)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: pyshmtvaropt.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Calculate the minimum variance and corresponding optimal weights of a localized multitaper spectral estimate.

## Usage

`var_opt`, `var_unit`, `weight_opt` = SHMTVarOpt (`l`, `tapers`, `taper_order`, `sff`, [`lwin`, `kmax`, `nocross`])

## Returns

`var_opt` : float, dimension (`kmax`)
:   The minimum variance of the multitaper spectral estimate for degree `l` using 1 through `kmax` tapers.

`var_unit` : float, dimension (`kmax`)
:   The variance of the multitaper spectral estimate using equal weights for degree `l` using 1 through `kmax` tapers.

`weight_opt` : float, dimension (`kmax`, `kmax`)
:   The optimal weights (in columns) that minimize the multitaper spectral estimate's variance using 1 through `kmax` tapers.

## Parameters

`l` : integer
:   The angular degree to determine the minimum variance and optimal weights.

`tapers` : float, dimension (`lwinin`+1, `kmaxin`)
:   A matrix of localization functions obtained from `SHReturnTapers` or `SHReturnTapersM`.

`taper_order` : integer, dimension (`kmaxin`)
:   The angular order of the windowing coefficients in TAPERS. If this matrix was created using `SHReturnTapersM`, then this array must be composed of zeros.

`sff` : float, dimension (`l`+`lwinin`+1)
:   The global unwindowed power spectrum of the function to be localized.

`lwin` : optional, integer, default = `lwinin`
:   The spherical harmonic bandwidth of the localizing windows.

`kmax` : optional, integer, default = `kmaxin`
:   The maximum number of tapers to be used when calculating the minimum variance and optimal weights.

`nocross` : optional, integer, default = 0
:   If 1, only the diagonal terms of the covariance matrix Fij will be computed. If 0, all terms will be computed.

## Description

`SHMTVarOpt` will determine the minimum variance that can be achieved by a weighted multitaper spectral analysis, as is described by Wieczorek and Simons (2007). The minimum variance is output as a function of the number of tapers utilized, from 1 to a maximum of `kmax`, and the corresponding variance using equal weights is output for comparison. The windowing functions are assumed to be solutions to the spherical-cap concentration problem, as determined by a call to `SHReturnTapers` or `SHReturnTapersM`. The minimum variance and weights are dependent upon the form of the global unwindowed power spectrum, `Sff`.

If the optional argument `nocross` is set to 1, then only the diagnonal terms of `Fij` will be computed.

## References

Wieczorek, M. A. and F. J. Simons, Minimum-variance multitaper spectral estimation on the sphere, J. Fourier Anal. Appl., 13, doi:10.1007/s00041-006-6904-1, 665-692, 2007.

## See also

[shreturntapers](pyshreturntapers.html), [shreturntapersm](pyshreturntapersm.html), [shmultitaperse](pyshmultitaperse.html), [shmultitapercse](pyshmultitapercse.html); [shlocalizedadmitcorr](pyshlocalizedadmitcorr.html), [shbiasadmitcorr](pyshbiasadmitcorr.html), [shbiask](pyshbiask.html), [shmtdebias](pyshmtdebias.html)
