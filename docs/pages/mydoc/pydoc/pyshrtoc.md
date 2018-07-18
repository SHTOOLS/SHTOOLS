---
title: SHrtoc (Python)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: pyshrtoc.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Convert real spherical harmonics to complex form.

## Usage

`ccilm` = SHrtoc (`rcilm`, [`lmax`, `convention`, `switchcs`])

## Returns

`ccilm` : float, dimension (2, `lmax`+1, `lmax`+1)
:   The output complex spherical harmonic coefficients. `ccilm[0,:,:]` and `ccilm[1,:,:]` correspond to the real and complex part of the coefficients, respectively. Only the positive angular orders are output; the negative orders can be calculated from the relation `C_{l-m}=(-1)^m C_{lm}^*`.

## Parameters

`rcilm` : float, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The input real spherical harmonic coefficients. `rcilm[0,:,:]` and `rcilm[1,:,:]` correspond to the cosine and sine terms, respectively.

`lmax` : optional, integer, default = `lmaxin`
:   The maximum degree of the output coefficients.

`convention` : optional, integer, default = 1
:   If 1 (default), the input and output coefficients will have the same normalization. If 2, real geodesy 4-pi coefficients will be converted to complex orthonormal form.

`swtichcs` : optional, integer default = 0
:   If 0 (default), the input and output coefficients will possess the same Condon-Shortley phase convention. If 1, the input coefficients will first be multiplied by (-1)^m.

## Description

`SHrtoc` will convert real spherical harmonics to complex form. The normalization of the input and output coefficients are by default the same, but if the optional argument `convention` is set to 2, this routine will convert from geodesy 4-pi normalized coefficients to orthonormalized coefficients. The Condon-Shortley phase convention between the input an output coefficients can be modified by the optional argument `switchcs`.

## See also

[shctor](pyshctor.html)
