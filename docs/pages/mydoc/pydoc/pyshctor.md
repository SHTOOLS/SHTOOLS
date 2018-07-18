---
title: SHctor (Python)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: pyshctor.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Convert complex spherical harmonics to real form.

## Usage

`rcilm` = SHctor (`ccilm`, [`lmax`, `convention`, `switchcs`])

## Returns

`rcilm` : float, dimension (2, `lmax`+1, `lamx`+1)
:   The output real spherical harmonic coefficients. `rcilm[0,:,:]` and `rcilm[1,:,:]` correspond to the cosine and sine terms, respectively.

## Parameters

`ccilm` : float, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The input complex spherical harmonic coefficients. `ccilm[0,:,:]` and `ccilm[1,:,:]` correspond to the real and complex part of the coefficients, respectively. Only the positive angular orders are input; the negative orders are assumed to satisfy the relation `C_{l-m}=(-1)^m C_{lm}^*`.

`lmax` : optional, integer, default = `lmaxin`
:   The maximum degree of the output coefficients.

`convention` : optional, integer, default = 1
:   If 1 (default), the input and output coefficients will have the same normalization. If 2, orthonormalized coefficients will be converted to real geodesy 4-pi form.

`swtichcs` : optional, integer, default = 0
:   If 0 (default), the input and output coefficients will possess the same Condon-Shortley phase convention. If 1, the input coefficients will first be multiplied by (-1)^m.

## Description

`SHctor` will convert complex spherical harmonics of a real function to real form. The normalization of the input and output coefficients are by default the same, but if the optional argument `convention` is set to 2, this routine will convert from geodesy 4-pi normalized coefficients to orthonormalized coefficients. The Condon-Shortley phase convention between the input an output coefficients can be modified by the optional argument `switchcs`.

## See also

[shrtoc](pyshrtoc.html)
