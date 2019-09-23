---
title: SHSCouplingMatrixCap (Fortran)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: shscouplingmatrixcap.html
summary:
tags: [fortran]
toc: false
editdoc: fdoc
---

This routine returns the spherical harmonic coupling matrix for a given set of spherical-cap Slepian basis functions. This matrix relates the power spectrum expectation of the function expressed in a subset of the best-localized Slepian functions to the expectation of the global power spectrum.

## Usage

call SHSCouplingMatrixCap (`kij`, `galpha`, `galpha_order`, `lmax`, `nmax`, `exitstatus`)

## Parameters

`kij` : output, real(dp), dimension (`lmax`+1, `lmax`+1)
:   The coupling matrix that relates the power spectrum expectation of the function expressed in a subset of the best-localized spherical-cap Slepian functions to the expectation of the global power spectrum.

`galpha` : input, real(dp), dimension (`lmax`+1, `nmax`)
:   An array of spherical-cap Slepian functions arranged in columns from best to worst localized and obtained from a call to `SHReturnTapers`.

`galpha_order` : input, integer, dimension (`kmax`)
:   The angular orders of the spherical-cap Slepian functions in `galpha`.

`lmax` : input, integer
:   The spherical harmonic bandwidth of the Slepian functions.

`nmax` : input, integer
:   The number of Slepian functions used in reconstructing the function.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

## Description

`SHSCouplingMatrixCap` returns the spherical harmonic coupling matrix that relates the power spectrum expectation of the function expressed in a subset of the best-localized spherical-cap Slepian functions to the expectation of the global power spectrum (assumed to be stationary). The spherical-cap Slepian functions are determined by a call to `SHReturnTapers` and each row of `galpha` contains the (`lmax`+1) spherical harmonic coefficients for the single angular order as given in `galpha_order`.

The relationship between the global and localized power spectra is given by:

`< S_{\tilde{f}}(l) > = \sum_{l'=0}^lmax K_{ll'} S_{f}(l')`

where `S_{\tilde{f}}` is the expectation of the power spectrum at degree l of the function expressed in Slepian functions, `S_{f}(l')` is the expectation of the global power spectrum, and `< ... >` is the expectation operator. The coupling matrix is given explicitly by

`K_{ll'} = \frac{1}{2l'+1} Sum_{m=-mmax}^mmax ( Sum_{alpha=1}^nmax g_{l'm}(alpha) g_{lm}(alpha) )**2`

where mmax is min(l, l').

## See also

[shreturntapers](shreturntapers.html), [shscouplingmatrix](shscouplingmatrix.html), [shslepianvar](shslepianvar.html), [shmtcouplingmatrix](shmtcouplingmatrix.html)
