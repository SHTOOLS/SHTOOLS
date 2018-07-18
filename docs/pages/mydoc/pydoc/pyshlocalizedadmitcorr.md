---
title: SHLocalizedAdmitCorr (Python)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: pyshlocalizedadmitcorr.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Calculate the localized admittance and correlation spectra of two functions at a given location using spherical cap localization windows.

## Usage

`admit`, `corr`, `admit_error`, `corr_error` = SHLocalizedAdmitCorr (`gilm`, `tilm`, `tapers`, `taper_order`, `lat`, `lon`, [`k`, `lwin`, `lmax`, `taper_wt`, `mtdef`, `k1linsig`])

## Returns

`admit` : float, dimension (`lmax`-`lwin`+1)
:   The admittance function, which is equal to `Sgt/Stt`.

`corr` : float, dimension (`lmax`-`lwin`+1)
:   The degree correlation function, which is equal to `Sgt/sqrt(Sgg Stt)`.

`admit_error` : float, dimension (`lmax`-`lwin`+1)
:   The standard error of the admittance function.

`corr_error` : float, dimension (`lmax`-`lwin`+1)
:   The standard error of the degree correlation function.

## Parameters

`gilm` : float, dimension (2, `lmaxgin`+1, `lmaxgin`+1)
:   The spherical harmonic coefficients of the function G.

`tilm` : float, dimension (2, `lmaxtin`+1, `lmaxtin`+1)
:   The spherical harmonic coefficients of the function T.

`tapers` : float, dimension (`lwinin`+1, `kin`)
:   A matrix of spherical cap localization functions obtained from `SHReturnTapers` or `SHReturnTapersM`.

`taper_order` : integer, dimension (`kin`)
:   The angular order of the windowing coefficients in `tapers`.

`lat` : float
:   The latitude of the localized analysis in degrees.

`lon` : float
:   The longitude of the localized analysis in degrees.

`k` : optional, integer, default = `kin`
:   The number of tapers to be used in the multitaper spectral analysis.

`lwin` : optional, integer, default = `lwinin`
:   The spherical harmonic bandwidth of the localizing windows.

`lmax` : optional, integer, default = min(`lmaxgin`, `lmaxtin`)
:   The maximum spherical harmonic degree of the input functions corresponding to `gilm` and `tilm`.

`taper_wt` : optional, float, dimension (`k`), default = -1 (not used)
:   The weights to be applied to the spectral estimates when calculating the admittance, correlation, and their associated errors. This must sum to unity. The default value specifies that taper weights will not be used.

`mtdef` : optional, integer, default = 1
:   1 (default): Calculate the multitaper spectral estimates Sgt, Sgg and Stt first, and then use these to calculate the admittance and correlation functions. 2: Calculate admittance and correlation spectra using each individual taper, and then average these to obtain the multitaper admittance and correlation functions.

`k1linsig` : optional, integer, default = -1
:   If equal to one, and only a single taper is being used, the errors in the admittance function will be calculated by assuming that the coefficients of `gilm` and `tilm` are related by a linear degree-dependent transfer function and that the lack of correlation is a result of uncorrelated noise. This is the square root of eq. 33 of Simons et al. 1997. The default value specifies that this is not set.

## Description

`SHLocalizedAdmitCorr` will calculate the localized admittance and degree correlation spectra of two functions at a given location. The windowing functions are solutions to the spherical-cap concentration problem (as calculated by `SHReturnTapers` or `SHReturnTapersM`), of which the best `k` concentrated tapers are utilized. If `k` is greater than 1, then estimates of the standard error for the admittance and correlation will be returned in the optional arrays `admit_error` and `corr_error`. The symmetry axis of the localizing windows are rotated to the coordinates (`lat`, `lon`) before performing the windowing operation.

The admittance is defined as `Sgt/Stt`, where `Sgt` is the localized cross-power spectrum of two functions `G` and `T` expressed in spherical harmonics. The localized degree-correlation spectrum is defined as `Sgt/sqrt(Sgg Stt)`, which can possess values between -1 and 1. Two methods are available for calculating the multitaper admittance and correlation functions. When `mtdef` is 1 (default), the multitaper estimates and errors of Sgt, Stt, and Sgg are calculated by calls to `SHMultiTaperSE` and `SHMultiTaperCSE`, and these results are then used to calculate the final admittance and correlation functions. When `mtdef` is 2, the admitance and correlation are calculated invidivually for each individual taper, and these results are then averaged.

If the optional parameter `k1linsig` is specified, and only a single taper is being used, the uncertainty in the admittance function will be calculated by assuming the two sets of coefficients are related by a linear degree-dependent transfer function and that the lack of correlation is a result of uncorrelated noise. 

When `mtdef` is 1, by default, the multitaper spectral estimates are calculated as an unweighted average of the individual tapered estimates. However, if the optional argument `taper_wt` is specified, a weighted average will be employed using the weights in this array. Minimum variance optimal weights can be obtained from the routines `SHMTVarOpt` if the form of the underlying global power spectrum is known. Taper weights can not be used when `mtdef` is 2

This routine assumes that the input functions and tapers are expressed using geodesy 4-pi normalized spherical harmonic functions that exclude the  Condon-Shortley phase factor of (-1)^m.

## References

Wieczorek, M. A. and F. J. Simons, Minimum-variance multitaper spectral estimation on the sphere, J. Fourier Anal. Appl., 13, doi:10.1007/s00041-006-6904-1, 665-692, 2007.

Simons, F. J., F. A. Dahlen and M. A. Wieczorek, Spatiospectral concentration on the sphere, SIAM Review, 48, 504-536, doi:10.1137/S0036144504445765, 2006. 

Wieczorek, M. A. and F. J. Simons, Localized spectral analysis on the sphere, 
Geophys. J. Int., 162, 655-675, 2005.

Simons, M., S. C. Solomon and B. H. Hager, Localization of gravity and topography: constrains on the tectonics and mantle dynamics of Venus, Geophys. J. Int., 131, 24-44, 1997.

## See also

[shreturntapers](pyshreturntapers.html), [shreturntapersm](pyshreturntapersm.html), [shmultitaperse](pyshmultitaperse.html), [shmultitapercse](pyshmultitapercse.html)

