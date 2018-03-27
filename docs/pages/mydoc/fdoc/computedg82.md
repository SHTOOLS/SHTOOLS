---
title: ComputeDG82 (Fortran)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: computedg82.html
summary:
tags: [fortran]
toc: false
editdoc: fdoc
---

Compute the tridiagonal matrix of Grunbaum et al. (1982) that commutes with the space-concentration kernel of a spherical cap.

## Usage

call ComputeDG82 (`dg82`, `lmax`, `m`, `theta0`, `exitstatus`)

## Parameters

`dg82` : output, real\*8, dimension (`lmax`-abs(`m`)+1, `lmax`-abs(`m`)+1)
:   The tridiagonal matrix of Grunbaum et al. (1982) that commutes with the space-concentration kernel of order M of a spherical cap.

`lmax` : input, integer
:   The spherical harmonic bandwidth of the windows.

`m` : input, integer
:   The angular order of the concentration problem.

`theta0` : input, real\*8
:   The angular radius of the spherical cap in radians.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

## Description

`ComputeDG82` will calculate the tridiagonal matrix of Grunbaum et al. (1982) that commutes with the space-concentration kernel of order `m` of a spherical cap. The eigenfunctions of this matrix correspond to a family of orthogonal windowing functions, and the eigenvalues correspond to the window's concentration factor (i.e., the power of the window within `theta0` divided by the total power of the function). It is assumed that the employed spherical harmonic functions are normalized to the same value for all degrees and angular orders, which is the case for both the geodesy 4-pi and orthonormalized harmonics. The returned matrix is symmetric, and the first element corresponds to (abs(`m`), abs(`m`)) as the values for elements less than this are identically zero.

## References

Grunbaum, F.A., L. Longhi, and M. Perlstadt, Differential operators commuting with finite convolution integral operators: some non-abelian examples, SIAM J. Appl. Math., 42, 941-955, 1982.

Simons, F.J., F.A. Dahlen, and M.A. Wieczorek, Spatiospectral concentration on a sphere, SIAM Review, 48, 504-536, 2006.

## See also

[computedm](computedm.html), [shreturntapers](shreturntapers.html), [shreturntapersm](shreturntapersm.html)
