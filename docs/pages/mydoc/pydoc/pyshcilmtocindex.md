---
title: SHCilmToCindex (Python)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: pyshcilmtocindex.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Convert a three-dimensional array of spherical harmonic coefficients to a two-dimensional indexed array.

## Usage

`cindex` = SHCilmToCindex (`cilm`, [`lmax`])

## Returns

`cindex` : float, dimension (2, (`lmax`+1)\*(`lmax`+2)/2)
:   The indexed output spherical harmonic coefficients.

## Parameters

`cilm` : float, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The input spherical harmonic coefficients. `cilm[0,:,:]` and `cilm[1,:,:]` correspond to either the real and imaginary components, or cosine and sine coefficients, respectively.

`lmax` : optional, integer, default = `lmaxin`
:   Maximum degree of input spherical harmonics to convert.

## Description

`SHCilmToCindex` will convert a three-dimensional array of spherical harmonic coefficients to a two-dimensional indexed array.  The degree `l` and order `m` corresponds to the index `l*(l+1)/2+m`.

## See also

[shcindextocilm](pyshcindextocilm.html), [shcilmtovector](pyshcilmtovector.html), [shvectortocilm](pyshvectortocilm.html)
