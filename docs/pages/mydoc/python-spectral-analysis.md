---
title: "Global and localized spectral analysis"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: python-spectral-analysis.html
summary: 
toc: true
folder: mydoc
---

<style>
table:nth-of-type(n) {
    display:table;
    width:100%;
}
table:nth-of-type(n) th:nth-of-type(2) {
    width:75%;
}
</style>

## Global spectral analysis

| Function name | Description |
| ------------- | ----------- |
| [spectrum](spectrum.html) | Calculate the spectrum of a real or complex function. |
| [cross_spectrum](cross_spectrum.html) | Calculate the cross-spectrum of a real or complex function. |
| [SHAdmitCorr](pyshadmitcorr.html) | Calculate the admittance and correlation spectra of two functions. |
| [SHConfidence](pyshconfidence.html) | Compute the probability that two sets of spherical harmonic coefficients are correlated at a given degree and for a given correlation coefficient. |

## Multitaper spectral estimation (spherical cap domain)

| Function name | Description |
| ------------- | ----------- |
| [SHMultiTaperSE](pyshmultitaperse.html) | Perform a localized multitaper spectral analysis. |
| [SHMultiTaperCSE](pyshmultitapercse.html) | Perform a localized multitaper cross-spectral analysis. |
| [SHLocalizedAdmitCorr](pyshlocalizedadmitcorr.html) | Calculate the localized admittance and correlation spectra of two functions at a given location. |
| [SHReturnTapers](pyshreturntapers.html) | Calculate the eigenfunctions of the spherical-cap concentration problem. |
| [SHReturnTapersM](pyshreturntapersm.html) | Calculate the eigenfunctions of the spherical-cap concentration problem for a single angular order. |
| [ComputeDm](pycomputedm.html) | Compute the space-concentration kernel of a spherical cap. |
| [ComputeDG82](pycomputedg82.html) | Compute the tridiagonal matrix of *Gr&uuml;nbaum et al.* (1982) that commutes with the space-concentration kernel of a spherical cap. |
| [SHFindLWin](pyshfindlwin.html) | Determine the spherical-harmonic bandwidth that is necessary to achieve a certain concentration factor. |
| [SHBiasK](pyshbiask.html) | Calculate the multitaper (cross-)power spectrum expectation of a windowed function. |
| [SHMTCouplingMatrix](pyshmtcouplingmatrix.html) | Calculate the multitaper coupling matrix for a given set of localization windows. |
| [SHBiasAdmitCorr](pyshbiasadmitcorr.html) | Calculate the expected multitaper admittance and correlation spectra associated with the input global cross-power spectra of two functions. |
| [SHMTDebias](pyshmtdebias.html) | Invert for the global power spectrum given a localized multitaper spectrum estimate. |
| [SHMTVarOpt](pyshmtvaropt.html) | Calculate the theoretical minimum variance of a localized multitaper spectral estimate and the corresponding optimal weights to apply to each localized spectrum. |
| [SHMTVar](pyshmtvar.html) | Calculate the theoretical variance of a multitaper spectral estimate for a given input power spectrum. |
| [SHSjkPG](pyshsjkpg.html) | Calculate the expectation of the product of two functions, each multiplied by a different data taper, for a given spherical harmonic degree and two different angular orders. |
| [SHRotateTapers](pyshrotatetapers.html) | Rotate orthogonal spherical-cap Slepian functions centered at the North pole to a different location. |

## Localization windows (arbitrary domain)

| Function name | Description |
| ------------- | ----------- |
| [SHReturnTapersMap](pyshreturntapersmap.html) | Calculate the eigenfunctions of the concentration problem for an arbitrary concentration region. |
| [SHBiasKMask](pyshbiaskmask.html) | Calculate the multitaper (cross-)power spectrum expectation of a function localized by arbitrary windows derived from a mask. |
| [SHMultiTaperMaskSE](pyshmultitapermaskse.html) | Perform a localized multitaper spectral analysis using arbitrary windows. |
| [SHMultiTaperMaskCSE](pyshmultitapermaskcse.html) | Perform a localized multitaper cross-spectral analysis using arbitrary windows. |
| [ComputeDMap](pycomputedmap.html) | Compute the space-concentration kernel of a mask defined on the sphere. |
| [Curve2Mask](pycurve2mask.html) | Given a set of latitude and longitude coordinates representing a closed curve, output a gridded mask. |

## Localization bias (general)

| Function name | Description |
| ------------- | ----------- |
| [SHBias](pyshbias.html) | Calculate the (cross-)power spectrum expectation of a windowed function. |

## Slepian function expansions

| Routine name | Description |
| ------------ | ----------- |
| [SlepianCoeffs](pyslepiancoeffs.html) | Determine the expansion coefficients of a function for a given set of input Slepian functions. |
| [SlepianCoeffsToSH](pyslepiancoeffstosh.html) | Convert a function expressed in Slepian coefficients to spherical harmonic coefficients. |
| [SHSCouplingMatrix](pyshscouplingmatrix.html) | Compute the spherical harmonic coupling matrix for a given set of Slepian functions. |
| [SHSCouplingMatrixCap](pyshscouplingmatrixcap.html) | Compute the spherical harmonic coupling matrix for a given set of spherical-cap Slepian functions. |
| [SHSlepianVar](pyshslepianvar.html) | Calculate the theoretical variance of the power of a function expanded in spherical-cap Slepian functions. |

## Other routines

| Function name | Description |
| ------------- | ----------- |
| [SphericalCapCoef](pysphericalcapcoef.html) | Calculate the spherical harmonic coefficients of a spherical cap.|

## References

* Grunbaum, F. A., L. Longhi, and M. Perlstadt, Differential operators commuting with finite convolution integral operators: some non-abelian examples, SIAM J. Appl. Math., 42, 941-955, doi:[10.1137/0142067](https://doi.org/10.1137/0142067), 1982.
