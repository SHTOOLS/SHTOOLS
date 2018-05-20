---
title: "Localized spectral analysis"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: localized-spectral-analysis.html
summary: 
toc: true
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

## Multitaper spectral estimation (spherical cap domain)

| Routine name | Description |
| ------------ | ----------- |
| [SHMultiTaperSE](shmultitaperse.html) | Perform a localized multitaper spectral analysis. |
| [SHMultiTaperCSE](shmultitapercse.html) | Perform a localized multitaper cross-spectral analysis. |
| [SHLocalizedAdmitCorr](shlocalizedadmitcorr.html) | Calculate the localized admittance and correlation spectra of two functions at a given location. |
| [SHReturnTapers](shreturntapers.html) | Calculate the eigenfunctions of the spherical-cap concentration problem. |
| [SHReturnTapersM](shreturntapersm.html) | Calculate the eigenfunctions of the spherical-cap concentration problem for a single angular order. |
| [ComputeDm](computedm.html) | Compute the space-concentration kernel of a spherical cap. |
| [ComputeDG82](computedg82.html) | Compute the tridiagonal matrix of *Gr&uuml;nbaum et al.* (1982) that commutes with the space-concentration kernel of a spherical cap. |
| [SHFindLWin](shfindlwin.html) | Determine the spherical-harmonic bandwidth that is necessary to achieve a certain concentration factor. |
| [SHBiasK](shbiask.html) | Calculate the multitaper (cross-)power spectrum expectation of a windowed function. |
| [SHMTCouplingMatrix](shmtcouplingmatrix.html) | Calculate the multitaper coupling matrix for a given set of localization windows. |
| [SHBiasAdmitCorr](shbiasadmitcorr.html) | Calculate the expected multitaper admittance and correlation spectra associated with the input global cross-power spectra of two functions. |
| [SHMTDebias](shmtdebias.html) | Invert for the global power spectrum given a localized multitaper spectrum estimate. |
| [SHMTVarOpt](shmtvaropt.html) | Calculate the minimum variance and corresponding optimal weights of a localized multitaper spectral estimate. |
| [SHSjkPG](shsjkpg.html) | Calculate the expectation of the product of two functions, each multiplied by a different data taper, for a given spherical harmonic degree and two different angular orders. |

## Localization windows (arbitrary domain)

| Routine name | Description |
| ------------ | ----------- |
| [SHReturnTapersMap](shreturntapersmap.html) | Calculate the eigenfunctions of the concentration problem for an arbitrary concentration region. |
| [SHMultiTaperMaskSE](shmultitapermaskse.html) | Perform a localized multitaper spectral analysis using arbitrary windows. |
| [SHMultiTaperMaskCSE](shmultitapermaskcse.html) | Perform a localized multitaper cross-spectral analysis using arbitrary windows. |
| [SHBiasKMask](shbiaskmask.html) | Calculate the multitaper (cross-)power spectrum expectation of a function localized by arbitrary windows derived from a mask. |
| [ComputeDMap](computedmap.html) | Compute the space-concentration kernel of a mask defined on the sphere. |
| [Curve2Mask](curve2mask.html) | Given a set of latitude and longitude coordinates representing a closed curve, output a gridded mask. |

## Localization bias (general)

| Routine name | Description |
| ------------ | ----------- |
| [SHBias](shbias.html) | Calculate the (cross-)power spectrum expectation of a windowed function. |

## Other routines

| Routine name | Description |
| ------------ | ----------- |
| [SphericalCapCoef](sphericalcapcoef.html) | Calculate the spherical harmonic coefficients of a spherical cap. |

## References

* Grünbaum, F.A., L. Longhi, and M. Perlstadt, Differential operators commuting with finite convolution integral operators: some non-abelian examples, SIAM J. Appl. Math., 42, 941-955, doi:[10.1137/0142067](https://doi.org/10.1137/0142067), 1982.
