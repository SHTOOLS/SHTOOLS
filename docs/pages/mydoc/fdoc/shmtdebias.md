---
title: SHMTDebias (Fortran)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: shmtdebias.html
summary:
tags: [fortran]
toc: false
editdoc: fdoc
---

Invert for the global power spectrum given a multitaper spectrum estimate formed with spherical cap localization windows.

## Usage

call SHMTDebias (`mtdebias`, `mtspectra`, `lmax`, `tapers`, `lwin`, `k`, `nl`, `lmid`, `n`, `taper_wt`, `exitstatus`)

## Parameters

`mtdebias` : output, real\*8, dimension (2, ceiling((`lmax`+1)/`nl`))
:   The global power spectrum (column 1) and uncertainty (column 2). The midpoints of the N spherical harmonic bins are given in `lmid`.

`mtspectra` : input, real\*8, dimension (2, `lmax`+1)
:   The localized multitaper spectrum estimate and uncertainty, obtained from a routine such as `SHMultitaperCSE` or `SHMultitaperSE`.

`lmax` : input, integer
:   The spherical harmonic bandwidth of the localized multitaper spectrum estimates.

`tapers` : input, real\*8, dimension (`lmaxt`+1, `k`)
:   An array of the K windowing functions, arranged in columns, obtained from a call to `SHReturnTapers`.

`lwin` : input, integer
:   The spherical harmonic bandwidth of the windowing functions in the array `tapers`.

`k` : input, integer
:   The number of tapers utilized in the multitaper spectral analysis.

`nl` : input, integer
:   The global power spectrum is assumed to be constant within bins of spherical harmonic wdith `nl`. In addition, the global power spectrum will be assumed to be constant beyond `lmax`.

`lmid` : output, real\*8, dimension ((`lmax`+1)/`nl`+1)
:   The midpoints of the spherical harmonic bins for which the global power spectrum is constant.

`n` : output, integer
:   The number of bins used in `mtdebias` and `lmid`.

`taper_wt` : input, optional, real\*8, dimension (`k`)
:   The weights used in calculating the multitaper spectral estimates. Optimal values of the weights (for a known global power spectrum) can be obtained from the routine `SHMTVarOpt`.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

## Description

`SHMTDebias` will invert for the global power spectrum given a localized multitaper spectrum estimate formed from spherical cap localization windows. This linear inverse problem is inherently underdetermined, and in order to achive a unique solution it is assumed that the global spectrum is constant in bins of width `nl`, and that the global power spectrum is constant for degrees greater than `lmax`. In practice `nl` should be increased until the global power spectrum is everywhere positive (negative values would be unphysical) and the variances are reasonable. Further details can be found in Wieczorek and Simons (2007).

This set of linear equations is solved using the method of singular value decomposition as outlined in Press et al. (1992, pp. 670-672). Each value of the multitaper spectrum estimate `mtspectra(1,:)`, as well as the corresponding rows of the transformation matrix, is divided by the uncertainties of the estimate `mtspectra(2,:)`. The solution and uncertainty are given by eqs 15.4.17 and 15.4.19 of Press et al. (1992, p. 671), respectively.

If `taper_wt` is not specified, the weights will all be assumed to be equal to `1/K`.

## References

Press, W.H., S.A. Teukolsky, W.T. Vetterling, and B.P. Flannery, Numerical Recipes in FORTRAN: The Art of Scientific Computing, 2nd ed., Cambridge Univ. Press, Cambridge, UK, 1992.

Wieczorek, M. A. and F. J. Simons, Minimum-variance multitaper spectral estimation on the sphere, J. Fourier Anal. Appl., 13, 665-692, doi:10.1007/s00041-006-6904-1, 2007.

## See also

[shmultitaperse](shmultitaperse.html), [shmultitapercse](shmultitapercse.html), [shreturntapers](shreturntapers.html), [shmtvaropt](shmtvaropt.html),  [shmtcouplingmatrix](shmtcouplingmatrix.html)
