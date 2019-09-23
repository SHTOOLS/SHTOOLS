---
title: SHSCouplingMatrix (Fortran)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: shscouplingmatrix.html
summary:
tags: [fortran]
toc: false
editdoc: fdoc
---

This routine returns the spherical harmonic coupling matrix for a given set of Slepian basis functions. This matrix relates the power spectrum expectation of the function expressed in a subset of the best-localized Slepian functions to the expectation of the global power spectrum.

## Usage

call SHSCouplingMatrix (`kij`, `galpha`, `lmax`, `nmax`, `exitstatus`)

## Parameters

`kij` : output, real(dp), dimension (`lmax`+1, `lmax`+1)
:   The coupling matrix that relates the power spectrum expectation of the function expressed in a subset of the best-localized Slepian functions to the expectation of the global power spectrum.

`galpha` : input, real(dp), dimension ((`lmax`+1)**2, `nmax`)
:   An array of Slepian functions, arranged in columns from best to worst localized.

`lmax` : input, integer
:   The spherical harmonic bandwidth of the Slepian functions.

`nmax` : input, integer
:   The number of Slepian functions used in reconstructing the function.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

## Description

`SHSCouplingMatrix` returns the spherical harmonic coupling matrix that relates the power spectrum expectation of the function expressed in a subset of the best-localized Slepian functions to the expectation of the global power spectrum (assumed to be stationary). The Slepian functions are determined by a call to either (1) `SHReturnTapers` and then `SHRotateTapers`, or (2) `SHReturnTapersMap`. Each row of `galpha` contains the (`lmax`+1)**2 spherical harmonic coefficients of a Slepian function that can be unpacked using `SHVectorToCilm`. The Slepian functions must be normalized to have unit power (that is the sum of the coefficients squared is 1).

The relationship between the global and localized power spectra is given by:

`< S_{\tilde{f}}(l) > = \sum_{l'=0}^lmax K_{ll'} S_{f}(l')`

where `S_{\tilde{f}}` is the expectation of the power spectrum at degree l of the function expressed in Slepian functions, `S_{f}(l')` is the expectation of the global power spectrum, and `< ... >` is the expectation operator. The coupling matrix is given explicitly by

`K_{ll'} = \frac{1}{2l'+1} Sum_{m=-l}^l Sum_{m'=-l'}^l' ( Sum_{alpha=1}^nmax g_{l'm'}(alpha) g_{lm}(alpha) )**2`

## See also

[shreturntapers](shreturntapers.html), [shrotatetapers](shrotatetapers.html), [shreturntapersmap](shreturntapersmap.html), [shvectorrocilm](shvectortocilm.html), [shscouplingmatrixcap](shmtcouplingmatrixcap.html), [shmtcouplingmatrix](shmtcouplingmatrix.html)
