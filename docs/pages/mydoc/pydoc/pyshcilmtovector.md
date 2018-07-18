---
title: SHCilmToVector (Python)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: pyshcilmtovector.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Convert a three-dimensional array of real spherical harmonic coefficients to a 1-dimensional indexed vector.

## Usage

`vector` = SHCilmToVector (`cilm`, [`lmax`])

## Returns

`vector` : float, dimension ( (`lmax`+1)\*\*2 )
:   The indexed output real spherical harmonic coefficients.

## Parameters

`cilm` : float, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The input real spherical harmonic coefficients.

`lmax` : optional, integer, default = `lmaxin`
:   The maximum degree of the output coefficients to convert.

## Description

`SHCilmToVector` will convert a three-dimensional array of real spherical harmonic coefficients to a 1-dimensional indexed array.  The degree `l`, order `m`, and `i` (1 for cosine, 2 for sine) corresponds to the index `l**2+(i-1)*l+m`.

## See also

[shvectortocilm](pyshvectortocilm.html), [yilmindexvector](pyyilmindexvector.html), [shcindextocilm](pyshcindextocilm.html), [shcilmtocindex](pyshcilmtocindex.html)
