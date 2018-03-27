---
title: SHReturnTapersMap (Python)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: pyshreturntapersmap.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Calculate the eigenfunctions and eigenvalues of the space-concentration problem for an arbitrary region.

## Usage

`tapers`, `eigenvalues` = SHReturnTapersMap (`dh_mask`, `lmax`, [`n`, `ntapers`, `sampling`])

## Returns

`tapers` : float, dimension ((`lmax`+1)\*\*2, `ntapers`)
:   The spherical harmonic coefficients of the tapers, arranged in columns, from best to worst concentrated. The spherical harmonic coefficients in each column are indexed according to the scheme described in `YilmIndexVector`.

`eigenvalues` : float, dimension (`ntapers`)
:   The concentration factor for each localization window specified in the columns of `tapers`.

## Parameters

`dh_mask` : integer, dimension (`n`, `n`\*`sampling`)
:   A Driscoll and Healy (1994) sampled grid describing the concentration region R. All elements should either be 1 (for inside the concentration region) or 0 (for outside R).

`lmax` : integer
:   The spherical harmonic bandwidth of the localization windows.

`n` : optional, integer
:   The number of latitudinal samples in `dh_mask`. The effective spherical harmonic bandwidth of this grid is `L=n/2-1`.

`ntapers` : optional, integer, default = (`lmax`+1)\*\*2
:   The number of best concentrated eigenvalues and corresponding eigenfunctions to return in `tapers` and `eigenvalues`. The default value is to return all tapers.

`sampling` : optional, integer, default determined from `dh_mask`
:   For 1, `dh_mask` has `n x n` samples. For 2, `dh_mask` has `n x 2n` samples. The default it to determine `sampling` from the dimensions of `dh_mask`.

## Description

`SHReturnTapersMap` will calculate the eigenfunctions (i.e., localization windows) of the space-concentration problem for an arbitrary concentration region specified in `dh_mask` (see Simons et al. (2006) for further details). The input mask `dh_mask` must be sampled according to the Driscoll and Healy (1994) sampling theorem with `n` samples in latitude, and possess a value of 1 inside the concentration region, and 0 elsewhere. `dh_mask` can either possess `n` samples in longitude (`sampling=1`) or `2n` samples in longitude (`sampling=2`). Given the approximate way in which the elements of the space-concentration kernel are calculated (see `ComputeDMap` for details), `sampling=2` should be preferred. The effective spherical harmonic bandwidth (L=N/2-1) of the grid `dh_mask` determines the accuracy of the results, and experience shows that this should be about 4 times larger than `lmax`. 

The spherical harmonic coefficients of each window are given in the columns of `tapers`, and the corresponding concentration factors are given in `eigenvaules`. The spherical harmonic coefficients are ordered according to the scheme described in `YilmIndexVector`, which can be converted to matrix form using `SHVectorToCilm`, and the columns of `tapers` are ordered from best to worst concentrated. The localization windows are normalized such that they have unit power. If the optional parameter `ntapers` is specified, then only the `ntapers` largest eigenvalues and corresponding eigenfunctions will be calculated and returned.

## References

Driscoll, J.R. and D.M. Healy, Computing Fourier transforms and convolutions on the 2-sphere, Adv. Appl. Math., 15, 202-250, 1994.

Simons, F.J., F.A. Dahlen, and M.A. Wieczorek, Spatiospectral concentration on a sphere, SIAM Review, 48, 504-536, 2006.

## See also

[computedmap](pycomputedmap.html), [yilmindexvector](pyyilmindexvector.html), [shvectortocilm](pyshvectortocilm.html)
