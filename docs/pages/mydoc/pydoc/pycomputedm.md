---
title: ComputeDM (Python)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: pycomputedm.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Compute the space-concentration kernel of a spherical cap.

## Usage

`dm` = ComputeDM (`lmax`, `m`, `theta0`)

## Returns

`dm` : float, dimension (`lmax`+1, `lmax`+1)
:   The space-concentration kernel or angular order `m`.

## Parameters

`lmax` : integer
:   The spherical harmonic bandwidth of the windows.

`m` : integer
:   The angular order of the concentration problem.

`theta0` : float
:   The angular radius of the spherical cap in radians.

## Description

`ComputeDM` will calculate the space-concentration kernel of angular order `m` for the spherical-cap concentration problem. The eigenfunctions of this matrix correspond to a family of orthogonal windowing functions, and the eigenvalues correspond to the window's concentration factor (i.e., the power of the window within `theta0` divided by the total power of the function). It is assumed that the employed spherical harmonic functions are normalized to the same value for all degrees and angular orders, which is the case for both the geodesy 4-pi and orthonormalized harmonics. This kernel is symmetric and is computed exactly by Gauss-Legendre quadrature.

## References

Simons, F.J., F.A. Dahlen, and M.A. Wieczorek, Spatiospectral concentration on a sphere, SIAM Review, 48, 504-536, 2006.

## See also

[computedg82](pycomputedg82.html), [shreturntapers](pyshreturntapers.html), [shreturntapersm](pyshreturntapersm.html)
