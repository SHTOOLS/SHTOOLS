# SHLocalizedAdmitCorr

Calculate the localized admittance and correlation spectra of two functions at a given location using spherical cap localization windows.

# Usage

call SHLocalizedAdmitCorr (`tapers`, `taper_order`, `lwin`, `lat`, `lon`, `gilm`, `tilm`, `lmax`, `admit`, `corr`, `k`, `admit_error`, `corr_error`, `taper_wt`, `mtdef`, `k1linsig`, `exitstatus`)

# Parameters

`tapers` : input, real\*8, dimension (`lwin`+1, `k`)
:   A matrix of spherical cap localization functions obtained from `SHReturnTapers` or `SHReturnTapersM`.

`taper_order` : input, integer, dimension (`k`)
:   The angular order of the windowing coefficients in `tapers`.

`lwin` : input, integer
:   The spherical harmonic bandwidth of the localizing windows.

`lat` : input, real\*8
:   The latitude of the localized analysis in degrees.

`lon` : input, real\*8
:   The longitude of the localized analysis in degrees.

`gilm` : input, real\*8, dimension (2, `lmax`+1, `lmax`+1)
:   The spherical harmonic coefficients of the function G.

`tilm` : input, real\*8, dimension (2, `lmax`+1, `lmax`+1)
:   The spherical harmonic coefficients of the function T.

`lmax` : input, integer
:   The maximum spherical harmonic degree of the input functions corresponding to `gilm` and `tilm`.

`admit` : output, real\*8, dimension (`lmax`-`lwin`+1)
:   The admittance function, which is equal to `Sgt/Stt`.

`corr` : output, real\*8, dimension (`lmax`-`lwin`+1)
:   The degree correlation function, which is equal to `Sgt/sqrt(Sgg Stt)`.

`k` : input, integer
:   The number of tapers to be used in the multitaper spectral analysis.

`admit_error` : output, optional, real\*8, dimension (`lmax`-`lwin`+1)
:   The standard error of the admittance function.

`corr_error` : output, optional, real\*8, dimension (`lmax`-`lwin`+1)
:   The standard error of the degree correlation function.

`taper_wt` : input, optional, real\*8, dimension (`k`)
:   The weights to be applied to the spectral estimates when calculating the admittance, correlation, and their associated errors. This must sum to unity.

`mtdef` : input, optional, integer, default = 1
:   1 (default): Calculate the multitaper spectral estimates Sgt, Sgg and Stt first, and then use these to calculate the admittance and correlation functions. 2: Calculate admittance and correlation spectra using each individual taper, and then average these to obtain the multitaper admittance and correlation functions.

`k1linsig` : input, optional, integer
:   If equal to one, and only a single taper is being used, the errors in the admittance function will be calculated by assuming that the coefficients of `gilm` and `tilm` are related by a linear degree-dependent transfer function and that the lack of correlation is a result of uncorrelated noise. This is the square root of eq. 33 of Simons et al. 1997.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`SHLocalizedAdmitCorr` will calculate the localized admittance and degree correlation spectra of two functions at a given location. The windowing functions are solutions to the spherical-cap concentration problem (as calculated by `SHReturnTapers` or `SHReturnTapersM`), of which the best `k` concentrated tapers are utilized. If `k` is greater than 1, then estimates of the standard error for the admittance and correlation will be returned in the optional arrays `admit_error` and `corr_error`. The symmetry axis of the localizing windows are rotated to the coordinates (`lat`, `lon`) before performing the windowing operation.

The admittance is defined as `Sgt/Stt`, where `Sgt` is the localized cross-power spectrum of two functions `G` and `T` expressed in spherical harmonics. The localized degree-correlation spectrum is defined as `Sgt/sqrt(Sgg Stt)`, which can possess values between -1 and 1. Two methods are available for calculating the multitaper admittance and correlation functions. When `mtdef` is 1 (default), the multitaper estimates and errors of Sgt, Stt, and Sgg are calculated by calls to `SHMultiTaperSE` and `SHMultiTaperCSE`, and these results are then used to calculate the final admittance and correlation functions. When `mtdef` is 2, the admitance and correlation are calculated invidivually for each individual taper, and these results are then averaged.

If the optional parameter `k1linsig` is specified, and only a single taper is being used, the uncertainty in the admittance function will be calculated by assuming the two sets of coefficients are related by a linear degree-dependent transfer function and that the lack of correlation is a result of uncorrelated noise.

When `mtdef` is 1, by default, the multitaper spectral estimates are calculated as an unweighted average of the individual tapered estimates. However, if the optional argument `taper_wt` is specified, a weighted average will be employed using the weights in this array. Minimum variance optimal weights can be obtained from the routines `SHMTVarOpt` if the form of the underlying global power spectrum is known. Taper weights can not be used when `mtdef` is 2

This routine assumes that the input functions and tapers are expressed using geodesy 4-pi normalized spherical harmonic functions that exclude the  Condon-Shortley phase factor of (-1)^m.

# See also

[shreturntapers](shreturntapers.html), [shreturntapersm](shreturntapersm.html), [shmultitaperse](shmultitaperse.html), [shmultitapercse](shmultitapercse.html)

# References

Wieczorek, M. A. and F. J. Simons, Minimum-variance multitaper spectral estimation on the sphere, J. Fourier Anal. Appl., 13, doi:10.1007/s00041-006-6904-1, 665-692, 2007.

Simons, F. J., F. A. Dahlen and M. A. Wieczorek, Spatiospectral concentration on the sphere, SIAM Review, 48, 504-536, doi:10.1137/S0036144504445765, 2006. 

Wieczorek, M. A. and F. J. Simons, Localized spectral analysis on the sphere, 
Geophys. J. Int., 162, 655-675, 2005.

Simons, M., S. C. Solomon and B. H. Hager, Localization of gravity and topography: constrains on the tectonics and mantle dynamics of Venus, Geophys. J. Int., 131, 24-44, 1997.
